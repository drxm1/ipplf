/** @file E03.cpp
 * @brief Run the periodic spherical pulse and record E03.meas01–06.
 */
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <string_view>

#include "core/output/RunRecord.h"
#include "core/simulation/ChDRSimulationManager.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Measurements.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Plane.h"

namespace chdr::experiments::e03 {
    namespace {
        using Simulation = simulation::ChDRSimulationManager<>;
        using Solver     = Simulation::Solver;

        /// Validate the cubic resolution before converting the YAML number to an index.
        int cellsPerAxis(const input::ChdrYamlRunInput& request) {
            if (!request.parameters_m.has_path("cells_per_axis"))
                throw std::invalid_argument("E03 requires cells_per_axis");
            const auto node = request.parameters_m["cells_per_axis"];
            if (node.dtype().is_number() && node.number_of_elements() == 1) {
                const double cells = node.to_float64();
                if (std::find(E03Levels.begin(), E03Levels.end(), cells) != E03Levels.end())
                    return static_cast<int>(cells);
            }
            throw std::invalid_argument("E03 cells_per_axis must be 96, 144, 216 or 324");
        }

        /// Both cases use the same grid and histories at their respective native times.
        std::string polarization(const input::ChdrYamlRunInput& request) {
            if (request.parameters_m.has_path("polarization")) {
                const auto node = request.parameters_m["polarization"];
                if (node.dtype().is_string()) {
                    const auto name = node.as_string();
                    if (name == "TEr" || name == "TMr")
                        return name;
                }
            }
            throw std::invalid_argument("E03 polarization must be TEr or TMr");
        }

        /// Invert LevelContract's label map to select the nearest state to a target.
        int nearestStep(Solver& state, Quantity quantity, double target) {
            return static_cast<int>(
                std::floor((target - state.labelTime(quantity, 0)) / state.getDt() + 0.5));
        }

        /// The potential, scalar and reconstructed fields may select different states.
        struct Checkpoint {
            double target_m;
            int aStep_m, phiStep_m, ebStep_m;
        };

        std::array<Checkpoint, 4> e03Checkpoints(Solver& state, double v) {
            std::array<Checkpoint, 4> result{};
            for (unsigned i = 0; i < result.size(); ++i) {
                const double target = 0.5 * (i + 1) / v;
                result[i] = {target, nearestStep(state, Quantity::VectorPotential, target),
                             nearestStep(state, Quantity::Phi, target),
                             nearestStep(state, Quantity::Electric, target)};
            }
            return result;
        }

        /** @brief Record raw numerators and their scales once on rank zero.
         * @details All ranks measure the same quantities collectively. Ratios can be recovered
         * from the CSV without storing duplicates. Checkpoint times are in solver units.
         */
        class E03MeasurementsCsv {
        public:
            E03MeasurementsCsv(const input::ChdrYamlRunInput& request, MPI_Comm comm, int cells,
                               std::string_view polarization)
                : cells_m(cells)
                , polarization_m(polarization) {
                int rank;
                MPI_Comm_rank(comm, &rank);
                if (rank != 0)
                    return;
                file_m.exceptions(std::ios::failbit | std::ios::badbit);
                file_m.open(std::filesystem::path(request.outputDirectory_m) / "measurements.csv");
                file_m << std::setprecision(std::numeric_limits<double>::max_digits10)
                       << "measurement,quantity,cells_per_axis,polarization,target_time,step,"
                          "actual_time,value,reference_scale,samples\n";
            }

            void row(std::string_view measurement, std::string_view quantity, double target,
                     long long step, double time, double value, double scale,
                     std::uint64_t samples = 0) {
                if (file_m.is_open())
                    file_m << measurement << ',' << quantity << ',' << cells_m << ','
                           << polarization_m << ',' << target << ',' << step << ',' << time << ','
                           << value << ',' << scale << ',' << samples << '\n';
            }

            void finish() {
                if (file_m.is_open())
                    file_m.close();
            }

        private:
            int cells_m;
            std::string_view polarization_m;
            std::ofstream file_m;
        };

        /** @brief Preserve the fixed initial reference norms used by E03.meas06. */
        struct ConstraintScales {
            double qL2_m, qLinf_m, gL2_m, gLinf_m;
        };

        template <class Field>
        E03FieldNorms norms(Solver& state, Field& field, const E03Reference& ref, Quantity quantity,
                            long long step) {
            return measureE03Field(field, ref, state.levelContract(), quantity,
                                   state.labelTime(quantity, step));
        }

        /// Compare both loaded storage slots, including the TE scalar zeros, on short checks.
        void checkLoadedHistory(Solver& state, const E03Reference& ref) {
            for (const auto quantity : {Quantity::VectorPotential, Quantity::Phi}) {
                if (norms(state, state.A_n, ref, quantity, 0).errorLinf_m != 0
                    || norms(state, state.A_nm1, ref, quantity, -1).errorLinf_m != 0)
                    throw std::runtime_error("E03 loaded history differs from its native samples");
            }
        }

        /// Save both norms of one residual with E03's fixed full-domain reference scales.
        void writeResidual(E03MeasurementsCsv& csv, const E03ResidualNorms& residual,
                           std::string_view l2, std::string_view linf, double target,
                           long long step, double time, double scaleL2, double scaleLinf) {
            csv.row("E03.meas06", l2, target, step, time, residual.l2_m, scaleL2, residual.count_m);
            csv.row("E03.meas06", linf, target, step, time, residual.linf_m, scaleLinf,
                    residual.count_m);
        }

        /// Observe E03.meas01–03 and electric Gauss law at the shared R2 field time.
        void observeFields(Solver& state, const E03Reference& ref, E03MeasurementsCsv& csv,
                           const ConstraintScales& scales, long long step, double target) {
            const double time = state.labelTime(Quantity::Electric, step);
            const auto e      = norms(state, state.electricField(), ref, Quantity::Electric, step);
            const auto b      = norms(state, state.magneticField(), ref, Quantity::Magnetic, step);
            csv.row("E03.meas01", "E_L2", target, step, time, e.errorL2_m, e.referenceL2_m);
            csv.row("E03.meas02", "B_L2", target, step, time, b.errorL2_m, b.referenceL2_m);
            csv.row("E03.meas03", "E_Linf", target, step, time, e.errorLinf_m, e.referenceLinf_m);
            writeResidual(csv, measureE03Gauss(state, E03Epsilon), "Q_L2", "Q_Linf", target, step,
                          time, scales.qL2_m, scales.qLinf_m);
        }

        /// Observe E03.meas04 and Lorenz gauge at the vector-potential time.
        void observeVectorPotential(Solver& state, const E03Reference& ref, E03MeasurementsCsv& csv,
                                    const ConstraintScales& scales, long long step, double target) {
            const double time = state.labelTime(Quantity::VectorPotential, step);
            if (step > 0) {
                const auto a = norms(state, state.A_n, ref, Quantity::VectorPotential, step);
                csv.row("E03.meas04", "A_L2", target, step, time, a.errorL2_m, a.referenceL2_m);
            }
            writeResidual(csv, measureE03Gauge(state, E03Epsilon), "G_L2", "G_Linf", target, step,
                          time, scales.gL2_m, scales.gLinf_m);
        }

        /// Use the independent scalar clock of R2 for E03.meas05.
        void observeScalarPotential(Solver& state, const E03Reference& ref, E03MeasurementsCsv& csv,
                                    long long step, double target) {
            const auto phi = norms(state, state.A_n, ref, Quantity::Phi, step);
            csv.row("E03.meas05", "phi_L2", target, step, state.labelTime(Quantity::Phi, step),
                    phi.errorL2_m, phi.referenceL2_m);
        }

        /// The reference scales are fixed at the initial E and A times on each mesh.
        ConstraintScales constraintScales(Solver& state, const E03Reference& ref) {
            const auto e = norms(state, state.electricField(), ref, Quantity::Electric, 0);
            const auto a = norms(state, state.A_n, ref, Quantity::VectorPotential, 0);
            return {E03Epsilon * e.referenceL2_m, E03Epsilon * e.referenceLinf_m, a.referenceL2_m,
                    a.referenceLinf_m};
        }

        /// Publish provenance only when the requested schedule and output have finished.
        void recordE03(const input::ChdrYamlRunInput& request,
                       const simulation::ChdrCartesianGridSettings& grid, Solver& state, int steps,
                       int protocolSteps, bool plane) {
            auto record = output::describeRun("experiments-wave-pulses-in-homogeneous-dielectric",
                                              grid, state, steps, "E03 periodic spherical pulse");
            record.effective_m["protocol_steps"]                                = protocolSteps;
            record.effective_m["material/relative_permittivity"]                = E03Epsilon;
            record.effective_m["analytic_reference/centre_in_mesh_coordinates"] = E03BoxCentre;
            record.effective_m["analytic_reference/speed"] = 1 / std::sqrt(E03Epsilon);
            record.effective_m["measurements/file"]        = "measurements.csv";
            if (plane)
                record.effective_m["plane/file"] = "plane-z0.csv";
            output::writeRunRecord(request.outputDirectory_m, request, record,
                                   state.A_n.getLayout().comm);
        }

        /// Route each state to the observations whose native clocks select it.
        void observeE03(Solver& state, const E03Reference& ref, E03MeasurementsCsv& csv,
                        const ConstraintScales& scales, const std::array<Checkpoint, 4>& points,
                        bool tm, long long step) {
            if (step == 0) {
                observeFields(state, ref, csv, scales, 0, 0);
                observeVectorPotential(state, ref, csv, scales, 0, 0);
            }
            for (const auto& point : points) {
                if (step == point.ebStep_m)
                    observeFields(state, ref, csv, scales, step, point.target_m);
                if (step == point.aStep_m)
                    observeVectorPotential(state, ref, csv, scales, step, point.target_m);
                if (tm && step == point.phiStep_m)
                    observeScalarPotential(state, ref, csv, step, point.target_m);
            }
        }

        /// Drive the manager once, sampling only the state indices selected by each clock.
        void evolveE03(Simulation& run, const input::ChdrYamlRunInput& request,
                       const simulation::ChdrCartesianGridSettings& grid, const E03Reference& ref,
                       std::string_view polarizationName) {
            auto& state         = run.solver();
            const auto points   = e03Checkpoints(state, ref.v_m);
            const auto last     = points.back();
            const bool tm       = polarizationName == "TMr";
            const int fullSteps = std::max({last.aStep_m, last.ebStep_m, tm ? last.phiStep_m : 0});
            const int steps     = request.steps_m.value_or(fullSteps);
            if (steps < 0 || steps > fullSteps)
                throw std::invalid_argument("E03 run.steps must fit the prescribed interval");
            const bool plane = input::booleanParameter(request, "write_plane", false);
            const auto comm  = state.A_n.getLayout().comm;
            output::prepareRunDirectory(request.outputDirectory_m, comm);
            E03MeasurementsCsv csv(request, comm, grid.cells_m[0], polarizationName);
            const auto scales = constraintScales(state, ref);
            if (request.steps_m)
                checkLoadedHistory(state, ref);
            run.advance(steps, [&](long long step, Solver& solver, const auto&) {
                observeE03(solver, ref, csv, scales, points, tm, step);
                if (plane && (step == 0 || step == last.ebStep_m))
                    writeE03Plane(solver, ref, request.outputDirectory_m, step);
            });
            if (request.steps_m && !tm
                && norms(state, state.A_n, ref, Quantity::Phi, steps).errorLinf_m != 0)
                throw std::runtime_error("E03 TE scalar potential did not remain zero");
            csv.finish();
            recordE03(request, grid, state, steps, fullSteps, plane);
        }
    }  // namespace

    void runE03(const input::ChdrYamlRunInput& request) {
        input::validateParameters(request, {"cells_per_axis", "polarization", "write_plane"});
        input::validateObservationNames(request, {});
        if (request.units_m.system_m != "normalized" || request.duration_m
            || request.outputDirectory_m.empty())
            throw std::invalid_argument(
                "E03 requires normalized units, fixed duration and output.directory");
        const int n     = cellsPerAxis(request);
        const auto name = polarization(request);
        const double h  = 2 * E03BoxCentre / n;
        const simulation::ChdrCartesianGridSettings grid{
            {n, n, n}, {h, h, h}, 0, {true, true, true}, 2};
        const E03Reference ref{
            1 / std::sqrt(E03Epsilon),
            name == "TEr" ? reference::Polarization::TEr : reference::Polarization::TMr,
            E03BoxCentre};
        const auto dielectric =
            KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
            marked = false;
            return E03Epsilon;
        };
        Simulation run{grid, dielectric, {.timeStep_m = h / 2}};
        run.load(ref);
        evolveE03(run, request, grid, ref, name);
    }
}  // namespace chdr::experiments::e03
