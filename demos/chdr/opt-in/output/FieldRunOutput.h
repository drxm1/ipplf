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

#include "core/output/FieldOutput.h"
#include "core/simulation/SimulationConfig.h"
#include "opt-in/diagnostics/AlignedFieldEnergy.h"
#include "opt-in/diagnostics/FieldEnergy.h"

namespace chdr::output {
    /// Record Cartesian component offsets in the writer's E-then-B component order.
    inline std::array<double, 18> fieldOffsets(const solver::LevelContract& contract) {
        std::array<double, 18> result{};
        for (unsigned c = 0; c < 6; ++c) {
            const auto offset = contract.positionOffset(
                c < 3 ? solver::Quantity::Electric : solver::Quantity::RawMagnetic, c % 3);
            for (unsigned d = 0; d < 3; ++d)
                result[3 * c + d] = offset[d];
        }
        return result;
    }

    /** @brief Write one material map, field frames and two energy diagnostics for a run.
     * @details The experiment chooses frame times by calling write(). E and raw B retain
     * their separate sample times in VTK. The CSV records the existing discrete quantities
     * and the squared-field energy from alignedFieldEnergy. The vector and aligned-field
     * energies share the electric-field time; the scalar diagnostic has the integer time.
     * Every rank participates in construction, write() and finish(). Mesh/layout owners
     * must outlive this object; captured observers must not outlive the output object.
     * @pre The solver uses O1-VS with staggered R2 reconstruction.
     */
    class FieldRunOutput {
    public:
        template <class Solver>
        FieldRunOutput(const std::string& directory, const simulation::RunBox& box, Solver& state,
                       const material::MaterialOwner& owner, const std::string& metadata)
            : writer_m(directory, {box.cells_m, {0, 0, 0}, box.spacing_m},
                       state.A_n.getLayout().comm, true) {
            raw_m.initialize(state.A_n.get_mesh(), state.A_n.getLayout());
            writer_m.writeMaterial(state.A_n, owner.epsilon(), fieldOffsets(state.levelContract()));
            openEnergy(directory, metadata, state.A_n.getLayout().comm);
        }

        /// Record diagnostics and raw field samples for one completed solver state.
        template <class Solver>
        void write(long long step, Solver& state, const material::MaterialOwner& owner) {
            if (step < 0 || step > std::numeric_limits<int>::max())
                throw std::invalid_argument("The field writer's step index is outside its range");
            const auto energy    = solver::potentialEnergy(state, owner);
            const double aligned = solver::alignedFieldEnergy(state, owner);
            if (!std::isfinite(energy.vector_m.total()) || !std::isfinite(aligned))
                throw std::runtime_error("Nonfinite energy in material movie");
            writeEnergy(step, energy, aligned);
            state.computeRawMagneticField(raw_m);
            writer_m.writeFrame(
                state.electricField(), raw_m,
                {static_cast<int>(step), state.labelTime(solver::Quantity::VectorPotential, step),
                 state.labelTime(solver::Quantity::Electric, step),
                 state.labelTime(solver::Quantity::RawMagnetic, step), state.getDt()});
        }

        /// Complete the VTK series after the last frame.
        void finish() { writer_m.writeSeries(); }

    private:
        /// Open the rank-zero CSV and propagate failure before any later collective observation.
        void openEnergy(const std::string& directory, const std::string& metadata, MPI_Comm comm) {
            int rank = 0;
            MPI_Comm_rank(comm, &rank);
            if (rank == 0) {
                energy_m.open(directory + "/energy.csv");
                energy_m << std::setprecision(17) << "# " << metadata << '\n'
                         << "step,time,vacuum,dielectric,total,scalar,absolute,field_energy\n";
            }
            int failed = rank == 0 && !energy_m;
            MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MAX, comm);
            if (failed)
                throw std::runtime_error("Cannot open the energy CSV");
        }

        /// Write diagnostics on the rank that owns the stream; other ranks do nothing.
        void writeEnergy(long long step, const solver::EnergyRecord& record, double aligned) {
            if (energy_m.is_open())
                energy_m << step << ',' << record.vectorTime_m << ',' << record.vector_m.vacuum_m
                         << ',' << record.vector_m.dielectric_m << ',' << record.vector_m.total()
                         << ',' << record.scalar_m.total() << ',' << record.vector_m.absolute_m
                         << ',' << aligned << '\n';
        }

        FieldOutput writer_m;
        solver::EMFieldsType raw_m;
        std::ofstream energy_m;
    };
}  // namespace chdr::output
