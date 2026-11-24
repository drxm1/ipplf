/** @file E16E17Output.h
 * @brief Measurement tables, native profiles and full-domain field snapshots for E16/E17.
 */
#pragma once

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "core/input/RunInput.h"
#include "core/output/FieldOutput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/SimulationConfig.h"
#include "core/solver/ChdrSolver.h"
#include "experiments/E16E17_normal_incidence_on_planar_dielectric_interface/E16E17Measurements.h"

namespace chdr::experiments::e16e17 {
    namespace outputDetail {
        /// Open an output table with round-trip double precision and a locale-independent format.
        inline void openCsv(std::ofstream& file, const std::string& directory, const char* name,
                            const char* header) {
            file.exceptions(std::ios::failbit | std::ios::badbit);
            file.imbue(std::locale::classic());
            file.open(std::filesystem::path(directory) / name);
            file << std::setprecision(std::numeric_limits<double>::max_digits10) << header << '\n';
        }

        /** @brief Gather the unique owned samples of one transverse line onto rank zero.
         * @details The lower of the two central cells selects the line on an even mesh.
         * A device kernel copies only the local line into a temporary device view.
         * Its host mirror is reduced across ranks. Each global index has one owner,
         * and the zero-initialized view supplies zero on all other ranks.
         */
        template <class Field>
        std::vector<double> gatherLine(Field& field, unsigned component) {
            const auto local  = field.getLayout().getLocalNDIndex();
            const auto global = field.getLayout().getDomain();
            const auto view   = field.getView();
            const int lineX   = global[0].first() + global[0].length() / 2 - 1;
            const int lineY   = global[1].first() + global[1].length() / 2 - 1;
            const int ghost = field.getNghost(), firstZ = local[2].first() - global[2].first();
            using ExecutionSpace = typename Field::view_type::execution_space;
            Kokkos::View<double*, typename Field::view_type::device_type> line(
                "E16/E17 profile samples", global[2].length());
            if (local[0].first() <= lineX && lineX <= local[0].last() && local[1].first() <= lineY
                && lineY <= local[1].last()) {
                const int i = lineX - local[0].first() + ghost;
                const int j = lineY - local[1].first() + ghost;
                Kokkos::parallel_for(
                    "E16/E17 extract profile",
                    Kokkos::RangePolicy<ExecutionSpace>(0, local[2].length()),
                    KOKKOS_LAMBDA(int k) { line(firstZ + k) = view(i, j, k + ghost)[component]; });
            }
            Kokkos::fence();
            const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, line);
            std::vector<double> values(global[2].length());
            MPI_Reduce(host.data(), values.data(), values.size(), MPI_DOUBLE, MPI_SUM, 0,
                       field.getLayout().comm);
            return values;
        }

        /// Identify the nonzero component of each profile channel in Cartesian numbering.
        inline unsigned component(Quantity quantity) {
            return quantity == Quantity::Magnetic ? 1 : 0;
        }

        /// Profile labels retain the E01 table format.
        inline const char* name(Quantity quantity) {
            if (quantity == Quantity::VectorPotential)
                return "Ax";
            return quantity == Quantity::Electric ? "Ex" : "By";
        }

        /** @brief Locate one profile sample using the solver's component placement.
         * @details Positive half-cell offsets in LevelContract specify the stored
         * staggered components [ryu2016potentialfdtd, Appendix A, Eqs. (19)--(24)].
         * The line uses the lower central transverse cell and the global longitudinal index.
         */
        template <class Field>
        auto position(Field& field, LevelContract contract, Quantity quantity, int k) {
            const auto h = field.get_mesh().getMeshSpacing(), origin = field.get_mesh().getOrigin();
            const auto global = field.getLayout().getDomain();
            const auto offset = contract.positionOffset(quantity, component(quantity));
            const int i       = global[0].first() + global[0].length() / 2 - 1;
            const int j       = global[1].first() + global[1].length() / 2 - 1;
            return std::array<double, 3>{origin[0] + (i + 0.5 + offset[0]) * h[0],
                                         origin[1] + (j + 0.5 + offset[1]) * h[1],
                                         origin[2] + (k + 0.5 + offset[2]) * h[2]};
        }
    }  // namespace outputDetail

    /** @brief Save the raw sums and energies required by E16.meas01--03.
     * @pre The output directory has been prepared collectively. All ranks construct
     * and call this writer in the same order. Only rank zero writes to disk.
     */
    class MeasurementCsv {
    public:
        MeasurementCsv(const std::string& directory, MPI_Comm communicator) {
            MPI_Comm_rank(communicator, &rank_m);
            if (rank_m == 0)
                outputDetail::openCsv(file_m, directory, "measurements.csv",
                                      "quantity,checkpoint,target_time,step,actual_time,value");
        }

        /// Store the actual native time separately from the requested observation target.
        void row(const char* quantity, const char* checkpoint, double target, long long step,
                 double time, double value) {
            if (rank_m != 0)
                return;
            if (!std::isfinite(value))
                throw std::runtime_error("E16/E17 measurement CSV contains a nonfinite value");
            file_m << quantity << ',' << checkpoint << ',' << target << ',' << step << ',' << time
                   << ',' << value << '\n';
            file_m.flush();
        }

        /// Close the table before publishing a completed run record.
        void finish() {
            if (rank_m == 0)
                file_m.close();
        }

    private:
        int rank_m = 0;
        std::ofstream file_m;
    };

    /** @brief Write Ax, Ex and R2 By along the native longitudinal line of E01.
     * @details The analytic reference uses each component's actual coordinate and time.
     * No field interpolation or reconstruction occurs here. Only the extracted line
     * is copied to the host. The full-domain arrays remain in their own memory space.
     * @pre The output directory has already been prepared collectively.
     */
    template <class Reference = E16E17AnalyticContinuumReferenceSampler>
    class ProfileCsv {
    public:
        ProfileCsv(const std::string& directory, const Reference& reference, LevelContract contract,
                   MPI_Comm communicator, const char* filename = "profiles.csv")
            : reference_m(reference)
            , contract_m(contract) {
            MPI_Comm_rank(communicator, &rank_m);
            if (rank_m == 0)
                outputDetail::openCsv(
                    file_m, directory, filename,
                    "component,checkpoint,target_time,step,actual_time,k,x,y,z,numerical,analytic");
        }

        /// Gather one component collectively and write the native values on rank zero.
        template <class Field>
        void write(Field& field, Quantity quantity, const char* checkpoint, double target,
                   long long step, double time) {
            const unsigned slot =
                outputDetail::component(quantity) + (quantity == Quantity::VectorPotential ? 1 : 0);
            const auto values = outputDetail::gatherLine(field, slot);
            if (rank_m == 0) {
                writeRows(field, quantity, checkpoint, target, step, time, values);
                file_m.flush();
            }
        }

        /// Close profiles before publishing the completed run record.
        void finish() {
            if (rank_m == 0)
                file_m.close();
        }

    private:
        /// Write one row per global longitudinal cell, preserving the periodic endpoint sample.
        template <class Field>
        void writeRows(Field& field, Quantity quantity, const char* checkpoint, double target,
                       long long step, double time, const std::vector<double>& values) {
            const int firstZ = field.getLayout().getDomain()[2].first();
            for (std::size_t index = 0; index < values.size(); ++index) {
                const int k        = firstZ + static_cast<int>(index);
                const auto point   = outputDetail::position(field, contract_m, quantity, k);
                const double exact = measurementDetail::referenceValue(
                    reference_m, quantity, point[0], point[1], point[2], time,
                    outputDetail::component(quantity));
                if (!std::isfinite(values[index]) || !std::isfinite(exact))
                    throw std::runtime_error("E16/E17 profile contains a nonfinite sample");
                file_m << outputDetail::name(quantity) << ',' << checkpoint << ',' << target << ','
                       << step << ',' << time << ',' << k << ',' << point[0] << ',' << point[1]
                       << ',' << point[2] << ',' << values[index] << ',' << exact << '\n';
            }
        }

        Reference reference_m;
        LevelContract contract_m;
        int rank_m = 0;
        std::ofstream file_m;
    };

    /** @brief Save an optional sequence of native Ex and R2 By profiles for movies.
     * @details Movie frames use the same sampling and CSV layout as the three prescribed
     * profiles. Each field's actual time is also its target time. The initial state and
     * final state are included, together with the configured cadence. The separate file
     * lets the analysis retain its three prescribed observations.
     */
    class MovieProfiles {
        using FieldSolver = solver::ChdrSolver<solver::StandardSolverPeriodic>;

    public:
        MovieProfiles(const input::ChdrYamlRunInput& request,
                      const E16E17AnalyticContinuumReferenceSampler& reference,
                      LevelContract contract, MPI_Comm communicator, int defaultEvery)
            : schedule_m(input::observationSchedule(
                  request, "movie_profiles", {.every_m = defaultEvery, .includeFinal_m = true}))
            , enabled_m(input::observationEnabled(
                  request, "movie_profiles", request.observations_m.contains("movie_profiles"))) {
            if (enabled_m)
                profiles_m =
                    std::make_unique<ProfileCsv<>>(request.outputDirectory_m, reference, contract,
                                                   communicator, "movie-profiles.csv");
        }

        /// Save the numerical and continuum fields at each selected native half-step.
        void observe(long long step, FieldSolver& state, int finalStep) {
            if (!enabled_m || !schedule_m.selects(step, finalStep))
                return;
            const double electricTime = state.labelTime(Quantity::Electric, step);
            const double magneticTime = state.labelTime(Quantity::Magnetic, step);
            profiles_m->write(state.electricField(), Quantity::Electric, "movie", electricTime,
                              step, electricTime);
            profiles_m->write(state.magneticField(), Quantity::Magnetic, "movie", magneticTime,
                              step, magneticTime);
        }

        /// Close the optional movie source before publishing the completed run record.
        void finish() {
            if (profiles_m)
                profiles_m->finish();
        }

        /// Record the effective movie cadence and its native field convention.
        void record(output::RunRecordDescription& description) const {
            if (!enabled_m)
                return;
            output::recordObservation(description, "movie_profiles", schedule_m, enabled_m);
            auto node                 = description.effective_m["movie_profiles"];
            node["file"]              = "movie-profiles.csv";
            node["components"]        = "Ex,By";
            node["sampling_metadata"] = "derived from solver LevelContract";
        }

    private:
        simulation::ObservationSchedule schedule_m;
        bool enabled_m;
        std::unique_ptr<ProfileCsv<>> profiles_m;
    };

    /** @brief Save E, R2 B and raw B at the three prescribed field observations.
     * @details FieldOutput records positions and times from LevelContract. Raw B is
     * reconstructed in independent storage, so magneticField() retains the R2 samples
     * used by the measurements. The caller provides the selected step indices.
     * @pre All ranks construct, observe and finish in the same order.
     */
    class FieldSnapshots {
        using FieldSolver = solver::ChdrSolver<solver::StandardSolverPeriodic>;

    public:
        FieldSnapshots(const input::ChdrYamlRunInput& request,
                       const simulation::ChdrCartesianGridSettings& grid, FieldSolver& state,
                       std::vector<int> prescribedSteps)
            : steps_m(std::move(prescribedSteps))
            , enabled_m(
                  input::observationEnabled(request, "fields", request.caseName_m != "E16-S")) {
            if (!enabled_m)
                return;
            writer_m = std::make_unique<output::FieldOutput>(
                std::filesystem::path(request.outputDirectory_m) / "fields",
                output::OutputGeometry{grid.cells_m, {0, 0, 0}, grid.spacing_m},
                state.A_n.getLayout().comm, "fields", true,
                output::fieldOutputScales(request.units_m));
            rawMagnetic_m.initialize(state.A_n.get_mesh(), state.A_n.getLayout());
        }

        /// Serialize the native arrays only at the supplied field checkpoints.
        void observe(long long step, FieldSolver& state) {
            if (!enabled_m || std::find(steps_m.begin(), steps_m.end(), step) == steps_m.end())
                return;
            state.computeRawMagneticField(rawMagnetic_m);
            const std::array fields{
                output::VectorFieldSnapshot{"E", state.electricField(), Quantity::Electric, 1.0},
                output::VectorFieldSnapshot{"B_R2", state.magneticField(), Quantity::Magnetic, 1.0},
                output::VectorFieldSnapshot{"B_raw", rawMagnetic_m, Quantity::RawMagnetic, 1.0}};
            const output::FrameRecord frame{static_cast<int>(step),
                                            state.labelTime(Quantity::VectorPotential, step),
                                            state.getDt()};
            writer_m->writeFrame(frame, state.levelContract(), fields);
        }

        /// Publish the completed ParaView series after the last frame.
        void finish() {
            if (writer_m)
                writer_m->writeSeries();
        }

        /// Keep the selected steps and array conventions in the run record.
        void record(output::RunRecordDescription& description) const {
            auto node       = description.effective_m["field_snapshots"];
            node["enabled"] = enabled_m;
            node["prescribed_steps"].set(steps_m);
            node["arrays"]            = "E,B_R2,B_raw";
            node["sampling_metadata"] = "derived from solver LevelContract";
            if (enabled_m)
                node["series"] = "fields/fields.pvti.series";
        }

    private:
        std::vector<int> steps_m;
        bool enabled_m;
        std::unique_ptr<output::FieldOutput> writer_m;
        EMFieldsType rawMagnetic_m;
    };
}  // namespace chdr::experiments::e16e17
