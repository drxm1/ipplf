/** @file E04.cpp
 * @brief Observe spherical pulse exit, return and constraints in the all-Mur cube.
 */
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string_view>

#include "core/output/RunRecord.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03E04Setup.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Measurements.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04Control.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04Monitor.h"
#include "opt-in/diagnostics/AlignedFieldEnergy.h"

namespace chdr::experiments::e04 {
    namespace {
        using namespace e03;
        using Solver                    = SphericalSimulation::Solver;
        constexpr double w              = reference::w;
        constexpr double finalEta       = 16;
        constexpr double precursorStart = 4.5;
        constexpr std::array targets{2., 3., 4., 5., 6., 8., 10., 12., finalEta};

        /// Label one row at its actual native time, independently of the target eta.
        struct Observation {
            long long step_m;
            double time_m, eta_m, targetEta_m;
        };

        class E04Csv {
        public:
            E04Csv(const std::string& directory, MPI_Comm comm) {
                int rank;
                MPI_Comm_rank(comm, &rank);
                if (rank != 0)
                    return;
                const std::filesystem::path dir(directory);
                open(energy_m, dir / "energy.csv", "step,time,eta,remaining_fraction");
                open(monitor_m, dir / "monitor.csv",
                     "step,time,eta,xminus_in,xminus_out,xplus_in,xplus_out,"
                     "yminus_in,yminus_out,yplus_in,yplus_out,zminus_in,zminus_out,zplus_in,zplus_"
                     "out");
                open(measurements_m, dir / "measurements.csv",
                     "measurement,quantity,region,target_eta,step,actual_time,actual_eta,"
                     "value,reference_scale,samples,unavailable");
            }

            /// Save the single E04.meas01 energy trace, including the initial energy bias.
            void energy(const Observation& at, double fraction) {
                if (energy_m.is_open())
                    energy_m << at.step_m << ',' << at.time_m << ',' << at.eta_m << ',' << fraction
                             << '\n';
            }

            /// Save the E04.meas02 per-face, sign-split powers for exact-window integration.
            void monitor(const Observation& at, const std::array<double, 12>& power) {
                if (!monitor_m.is_open())
                    return;
                monitor_m << at.step_m << ',' << at.time_m << ',' << at.eta_m;
                for (const double facePower : power)
                    monitor_m << ',' << facePower;
                monitor_m << '\n';
            }

            void row(std::string_view meas, std::string_view quantity, std::string_view region,
                     const Observation& at, double value, double scale, std::uint64_t samples = 0,
                     std::uint64_t unavailable = 0) {
                if (measurements_m.is_open())
                    measurements_m << meas << ',' << quantity << ',' << region << ','
                                   << at.targetEta_m << ',' << at.step_m << ',' << at.time_m << ','
                                   << at.eta_m << ',' << value << ',' << scale << ',' << samples
                                   << ',' << unavailable << '\n';
            }

            void finish() {
                for (auto* stream : {&energy_m, &monitor_m, &measurements_m})
                    if (stream->is_open())
                        stream->close();
            }

        private:
            static void open(std::ofstream& stream, const std::filesystem::path& path,
                             const char* header) {
                stream.exceptions(std::ios::failbit | std::ios::badbit);
                stream.open(path);
                stream << std::setprecision(std::numeric_limits<double>::max_digits10) << header
                       << '\n';
            }

            std::ofstream energy_m, monitor_m, measurements_m;
        };

        /// E04.meas01 uses the same energy weights for the reference and the difference.
        double energyFromNorms(double electric, double magnetic) {
            // Unit permeability [ellingson2020ElectromagneticsVol2, Eqs. (3.20), (3.22)].
            return 0.5 * (E03Epsilon * electric * electric + magnetic * magnetic);
        }

        struct InitialNorms {
            E03FieldNorms e_m, b_m, a_m;

            double referenceEnergy() const {
                return energyFromNorms(e_m.referenceL2_m, b_m.referenceL2_m);
            }
        };

        template <class Field>
        E03FieldNorms norms(Solver& state, Field& field, const E03Reference& ref, Quantity quantity,
                            long long step) {
            return measureE03Field(field, ref, state.levelContract(), quantity,
                                   state.labelTime(quantity, step));
        }

        /// The constraint scales and energy denominator are fixed at initialization.
        InitialNorms initialNorms(Solver& state, const E03Reference& ref) {
            return {norms(state, state.electricField(), ref, Quantity::Electric, 0),
                    norms(state, state.magneticField(), ref, Quantity::Magnetic, 0),
                    norms(state, state.A_n, ref, Quantity::VectorPotential, 0)};
        }

        /// At startup and eta=2 compare with both the exact pulse and the periodic replay.
        void compareControl(Solver& state, const E04PeriodicControl& control, E04Csv& csv,
                            const Observation& at, const E03FieldNorms& e, const E03FieldNorms& b) {
            const auto difference = control.difference(state, at.step_m);
            csv.row("E04.control", "E_continuum_L2", "box", at, e.errorL2_m, e.referenceL2_m);
            csv.row("E04.control", "B_continuum_L2", "box", at, b.errorL2_m, b.referenceL2_m);
            csv.row("E04.control", "E_periodic_L2", "box", at, difference[0], e.referenceL2_m);
            csv.row("E04.control", "B_periodic_L2", "box", at, difference[1], b.referenceL2_m);
        }

        /// Observe the field difference before the outgoing continuum support leaves the cube.
        void observeFields(Solver& state, const E03Reference& ref,
                           const E04PeriodicControl& control, E04Csv& csv, const Observation& at,
                           double U0h) {
            const auto e = norms(state, state.electricField(), ref, Quantity::Electric, at.step_m);
            const auto b = norms(state, state.magneticField(), ref, Quantity::Magnetic, at.step_m);
            if (at.targetEta_m == 2) {
                compareControl(state, control, csv, at, e, b);
                return;
            }
            const double error = std::sqrt(energyFromNorms(e.errorL2_m, b.errorL2_m));
            csv.row("E04.meas01", "e_field", "box", at, error, std::sqrt(U0h));
        }

        /** @brief Save E04.meas03 at one native time with residual-specific availability.
         * @details Gauss needs physical E samples. Gauge can use the recorded quadratic
         * potential closure. Both use the initial full-box scales specified in E03.typ.
         */
        template <bool Gauss>
        void observeConstraint(Solver& state, E04Csv& csv, const Observation& at,
                               const InitialNorms& initial) {
            const double l2Scale =
                Gauss ? E03Epsilon * initial.e_m.referenceL2_m : initial.a_m.referenceL2_m;
            const double maxScale =
                Gauss ? E03Epsilon * initial.e_m.referenceLinf_m : initial.a_m.referenceLinf_m;
            for (auto region : {ConstraintRegion::Interior, ConstraintRegion::Boundary}) {
                const auto residual = Gauss ? measureE03Gauss(state, E03Epsilon, region)
                                            : measureE03Gauge(state, E03Epsilon, region);
                const char* name = region == ConstraintRegion::Interior ? "interior" : "boundary";
                csv.row("E04.meas03", Gauss ? "Q_L2" : "G_L2", name, at, residual.l2_m, l2Scale / w,
                        residual.count_m, residual.unavailable_m);
                csv.row("E04.meas03", Gauss ? "Q_Linf" : "G_Linf", name, at, residual.linf_m,
                        maxScale / w, residual.count_m, residual.unavailable_m);
            }
        }

        struct Checkpoint {
            double eta_m;
            int fieldStep_m, gaugeStep_m;
        };

        Observation observation(Solver& state, Quantity quantity, long long step, double v,
                                double targetEta) {
            const double time = state.labelTime(quantity, step);
            return {step, time, v * time / w, targetEta};
        }

        /// Select field and gauge observations with their own native clocks.
        std::array<Checkpoint, targets.size()> checkpoints(Solver& state, double v) {
            std::array<Checkpoint, targets.size()> result{};
            for (unsigned i = 0; i < targets.size(); ++i) {
                const double time = targets[i] * w / v;
                result[i] = {targets[i], nearestSphericalStep(state, Quantity::Electric, time),
                             nearestSphericalStep(state, Quantity::VectorPotential, time)};
            }
            return result;
        }

        /// Evaluate only the target observables selected by this state.
        void observeTargets(Solver& state, const E03Reference& ref,
                            const E04PeriodicControl& control, E04Csv& csv,
                            const InitialNorms& initial,
                            const std::array<Checkpoint, targets.size()>& points, long long step) {
            for (const auto& point : points) {
                if (step == point.fieldStep_m) {
                    const auto at =
                        observation(state, Quantity::Electric, step, ref.v_m, point.eta_m);
                    if (point.eta_m <= 8)
                        observeFields(state, ref, control, csv, at, initial.referenceEnergy());
                    if (point.eta_m != 5)
                        observeConstraint<true>(state, csv, at, initial);
                }
                if (step == point.gaugeStep_m && point.eta_m != 5) {
                    const auto at =
                        observation(state, Quantity::VectorPotential, step, ref.v_m, point.eta_m);
                    observeConstraint<false>(state, csv, at, initial);
                }
            }
        }

        /// Record the actual component locations of Mur's outermost owned potential layers.
        void recordMurLayers(conduit_cpp::Node node, Solver& state, double h) {
            for (unsigned c = 0; c < 4; ++c) {
                const auto offset = state.levelContract().positionOffset(
                    c == 0 ? Quantity::Phi : Quantity::VectorPotential, c == 0 ? 0 : c - 1);
                auto component                           = node.append();
                component["potential_storage_component"] = c;
                for (unsigned axis = 0; axis < 3; ++axis) {
                    auto planes     = component["positions"].append();
                    planes["axis"]  = axis;
                    planes["lower"] = -E03BoxCentre + (0.5 + offset[axis]) * h;
                    planes["upper"] = E03BoxCentre + (offset[axis] - 0.5) * h;
                }
            }
        }

        void recordE04(const input::ChdrYamlRunInput& request,
                       const simulation::ChdrCartesianGridSettings& grid, Solver& state, int steps,
                       int protocolSteps, double U0h, int controlStep, double sharpness) {
            auto record =
                output::describeRun("experiments-wave-pulses-in-homogeneous-dielectric", grid,
                                    state, steps, "E04 spherical pulse with all-face Mur");
            auto& node                                   = record.effective_m;
            node["protocol_steps"]                       = protocolSteps;
            node["material/relative_permittivity"]       = E03Epsilon;
            node["analytic_reference/profile_sharpness"] = sharpness;
            node["boundary/mur_mask"]                    = state.boundary().murMask();
            node["boundary/kernel_speed"]                = state.boundary().speed[0];
            node["boundary/potential_ghost_closure"]     = "Quadratic extrapolation";
            recordMurLayers(node["boundary/potential_layers_centred"], state, grid.spacing_m[0]);
            node["initial_reference_energy"]           = U0h;
            node["periodic_control/target_eta"]        = 2.0;
            node["periodic_control/target_field_step"] = controlStep;
            output::writeRunRecord(request.outputDirectory_m, request, record,
                                   state.A_n.getLayout().comm);
        }

        /// Collect the dense E04.meas01/02 histories at the common R2 half-step.
        void observeTraces(Solver& state, const material::MaterialOwner& owner, E04Csv& csv,
                           long long step, int firstMonitorStep, double v, double U0h) {
            const auto at         = observation(state, Quantity::Electric, step, v, 0);
            const double fraction = diagnostics::alignedFieldEnergy(state, owner) / U0h;
            csv.energy(at, fraction);
            if (!std::isfinite(fraction))
                throw std::runtime_error("Nonfinite E04 electromagnetic energy");
            if (step >= firstMonitorStep)
                csv.monitor(at, measureE04Monitor(state, v, U0h));
        }

        /// Bracket the exact endpoints without duplicating the R2 half-step offset.
        int bracketStep(Solver& state, double eta, double v, bool upper) {
            const double step =
                (eta * w / v - state.labelTime(Quantity::Electric, 0)) / state.getDt();
            return static_cast<int>(upper ? std::ceil(step) : std::floor(step));
        }

        void evolveE04(SphericalSimulation& run, const input::ChdrYamlRunInput& request,
                       const simulation::ChdrCartesianGridSettings& grid, const E03Reference& ref) {
            auto& state         = run.solver();
            const auto points   = checkpoints(state, ref.v_m);
            const int fullSteps = bracketStep(state, finalEta, ref.v_m, true);
            const int steps     = request.steps_m.value_or(fullSteps);
            if (steps < 0 || steps > fullSteps)
                throw std::invalid_argument("E04 run.steps must fit the prescribed interval");
            const int controlStep = points.front().fieldStep_m;
            const auto comm       = state.A_n.getLayout().comm;
            output::prepareRunDirectory(request.outputDirectory_m, comm);
            E04PeriodicControl control(grid, ref, std::min(steps, controlStep), controlStep);
            E04Csv csv(request.outputDirectory_m, comm);
            const auto initial         = initialNorms(state, ref);
            const double U0h           = initial.referenceEnergy();
            const int firstMonitorStep = bracketStep(state, precursorStart, ref.v_m, false);
            compareControl(state, control, csv,
                           observation(state, Quantity::Electric, 0, ref.v_m, 0), initial.e_m,
                           initial.b_m);
            observeConstraint<true>(state, csv,
                                    observation(state, Quantity::Electric, 0, ref.v_m, 0), initial);
            observeConstraint<false>(
                state, csv, observation(state, Quantity::VectorPotential, 0, ref.v_m, 0), initial);
            run.advance(steps, [&](long long step, Solver& solver, const auto& owner) {
                observeTraces(solver, owner, csv, step, firstMonitorStep, ref.v_m, U0h);
                observeTargets(solver, ref, control, csv, initial, points, step);
            });
            csv.finish();
            recordE04(request, grid, state, steps, fullSteps, U0h, controlStep, ref.sharpness_m);
        }
    }  // namespace

    void runE04(const input::ChdrYamlRunInput& request) {
        input::validateParameters(request, {"cells_per_axis", "profile_sharpness"});
        e03::validateSphericalRun(request);
        const auto grid = e03::sphericalGrid(e03::sphericalCells(request), false);
        const auto ref  = e03::sphericalReference(e03::sphericalSharpness(request));
        auto run        = e03::makeSphericalSimulation(grid, ref);
        run.load(ref);
        evolveE04(run, request, grid, ref);
    }
}  // namespace chdr::experiments::e04
