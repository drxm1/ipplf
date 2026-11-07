/** @file FieldRunOutput.h
 * @brief Shared VTK frames and energy CSV output for staggered R2 O1-VS runs.
 */
#pragma once
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "core/input/UnitSystem.h"
#include "core/output/FieldOutput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/ObservationSchedule.h"
#include "core/simulation/SimulationConfig.h"
#include "opt-in/diagnostics/AlignedFieldEnergy.h"
#include "opt-in/diagnostics/FieldEnergy.h"

namespace chdr::output {
    /// Resolved sampling choices for the two outputs supplied by FieldRunOutput.
    struct FieldRunSettings {
        simulation::ObservationSchedule fields_m, diagnostics_m;
        bool fieldsEnabled_m = true, diagnosticsEnabled_m = true;
    };

    /// Resolve input once so output dispatch and the saved record use identical settings.
    inline FieldRunSettings fieldRunSettings(const input::ChdrYamlRunInput& request, int every) {
        const simulation::ObservationSchedule defaults{.every_m = every, .includeFinal_m = true};
        return {input::observationSchedule(request, "fields", defaults),
                input::observationSchedule(request, "diagnostics", defaults),
                input::observationEnabled(request, "fields"),
                input::observationEnabled(request, "diagnostics")};
    }

    /** @brief Write one material map, field frames and two energy diagnostics for a run.
     * @details The experiment schedules fields and diagnostics independently through
     * writeFields() and writeDiagnostics().
     * E and raw B retain their separate sample times in VTK. The CSV records the discrete
     * quantities and the squared-field energy from alignedFieldEnergy. The vector and aligned-field
     * energies share the electric-field time; the scalar diagnostic has the integer time.
     * Every rank participates in construction, each scheduled observation and finish().
     * Mesh/layout owners must outlive this object. The application owns the run callback
     * and calls finish() explicitly after successful evolution.
     * @pre The solver uses O1-VS with staggered R2 reconstruction.
     */
    class FieldRunOutput {
    public:
        template <class Solver>
        FieldRunOutput(const std::string& directory,
                       const simulation::ChdrCartesianGridSettings& box, Solver& state,
                       const material::MaterialOwner& owner, const std::string& metadata,
                       input::UnitSystem units = {})
            : writer_m(directory, {box.cells_m, {0, 0, 0}, box.spacing_m},
                       state.A_n.getLayout().comm, "wave", true, fieldOutputScales(units))
            , units_m(std::move(units))
            , communicator_m(state.A_n.getLayout().comm) {
            raw_m.initialize(state.A_n.get_mesh(), state.A_n.getLayout());
            writer_m.writeMaterial(state.A_n, owner.epsilon());
            openEnergy(directory, metadata, state.A_n.getLayout().comm);
        }

        /// Record the existing energy definitions without also writing a field frame.
        template <class Solver>
        void writeDiagnostics(long long step, Solver& state, const material::MaterialOwner& owner) {
            checkStep(step);
            const auto energy    = diagnostics::potentialEnergy(state, owner);
            const double aligned = diagnostics::alignedFieldEnergy(state, owner);
            if (!std::isfinite(energy.vector_m.total()) || !std::isfinite(aligned))
                throw std::runtime_error("Nonfinite energy in material movie");
            writeEnergy(step, energy, aligned);
        }

        /// Write a raw-field frame without evaluating the optional energy diagnostics.
        template <class Solver>
        void writeFields(long long step, Solver& state) {
            checkStep(step);
            state.computeRawMagneticField(raw_m);
            const std::array snapshots{
                VectorFieldSnapshot{"E_raw", state.electricField(), Quantity::Electric,
                                    units_m.electricFactor()},
                VectorFieldSnapshot{"B_raw", raw_m, Quantity::RawMagnetic,
                                    units_m.magneticFactor()}};
            writer_m.writeFrame({static_cast<int>(step),
                                 state.labelTime(Quantity::VectorPotential, step), state.getDt()},
                                state.levelContract(), snapshots);
        }

        /// Close the diagnostic stream and complete the series before recording success.
        void finish() {
            int failed = 0;
            if (energy_m.is_open()) {
                energy_m.close();
                failed = !energy_m;
            }
            MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MAX, communicator_m);
            if (failed)
                throw std::runtime_error("Cannot complete the energy CSV");
            writer_m.writeSeries();
        }

        /** @brief Write whichever existing outputs select this completed solver state.
         * @details Called by the application's manager callback, including state zero.
         * This method owns no callback or lifecycle. Field and diagnostic schedules are
         * independent; their common default preserves the existing output cadence.
         */
        template <class Solver>
        void observe(long long step, Solver& state, const material::MaterialOwner& owner,
                     long long finalStep, const FieldRunSettings& settings) {
            if (settings.diagnosticsEnabled_m && settings.diagnostics_m.selects(step, finalStep))
                writeDiagnostics(step, state, owner);
            if (settings.fieldsEnabled_m && settings.fields_m.selects(step, finalStep))
                writeFields(step, state);
        }

    private:
        /// Keep the public observation labels within the integer range of the VTK writer.
        static void checkStep(long long step) {
            if (step < 0 || step > std::numeric_limits<int>::max())
                throw std::invalid_argument("The field writer's step index is outside its range");
        }

        /// Open the rank-zero CSV and propagate failure before any later collective observation.
        void openEnergy(const std::string& directory, const std::string& metadata, MPI_Comm comm) {
            int rank = 0;
            MPI_Comm_rank(comm, &rank);
            if (rank == 0) {
                energy_m.open(directory + "/energy.csv");
                energy_m << std::setprecision(17) << "# " << metadata
                         << " energy_units=normalized time_units="
                         << (units_m.isSI() ? "s" : "normalized") << '\n'
                         << "step,time,vacuum,dielectric,total,scalar,absolute,field_energy";
                if (units_m.isSI())
                    energy_m << ",field_energy_J";
                energy_m << '\n';
            }
            int failed = rank == 0 && !energy_m;
            MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MAX, comm);
            if (failed)
                throw std::runtime_error("Cannot open the energy CSV");
        }

        /// Write diagnostics on the rank that owns the stream; other ranks do nothing.
        void writeEnergy(long long step, const diagnostics::EnergyRecord& record, double aligned) {
            if (energy_m.is_open()) {
                energy_m << step << ',' << record.vectorTime_m * units_m.timeFactor() << ','
                         << record.vector_m.vacuum_m << ',' << record.vector_m.dielectric_m << ','
                         << record.vector_m.total() << ',' << record.scalar_m.total() << ','
                         << record.vector_m.absolute_m << ',' << aligned;
                if (units_m.isSI())
                    energy_m << ',' << aligned * units_m.energyFactor();
                energy_m << '\n';
            }
        }

        FieldOutput writer_m;
        input::UnitSystem units_m;
        MPI_Comm communicator_m;
        EMFieldsType raw_m;
        std::ofstream energy_m;
    };

    /// Record the resolved schedules and the intentionally mixed diagnostic/SI column units.
    inline void recordFieldRunOutput(RunRecordDescription& description,
                                     const input::ChdrYamlRunInput& request,
                                     const FieldRunSettings& settings) {
        recordObservation(description, "fields", settings.fields_m, settings.fieldsEnabled_m);
        recordObservation(description, "diagnostics", settings.diagnostics_m,
                          settings.diagnosticsEnabled_m);
        description.effective_m["output/energy_csv/energy_columns_units"] = "normalized";
        description.effective_m["output/energy_csv/time_units"] =
            request.units_m.isSI() ? "s" : "normalized";
        if (request.units_m.isSI())
            description.effective_m["output/energy_csv/physical_energy_column"] = "field_energy_J";
    }
}  // namespace chdr::output
