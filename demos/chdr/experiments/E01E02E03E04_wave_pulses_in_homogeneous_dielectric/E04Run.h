/** @file E04Run.h
 * @brief Observe spherical pulse exit, return and constraints in the all-Mur cube.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <string_view>

#include "core/input/RunInput.h"
#include "core/output/RunRecord.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03E04FieldSnapshots.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03E04Setup.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Measurements.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04Control.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04Monitor.h"
#include "opt-in/diagnostics/AlignedFieldEnergy.h"

namespace chdr::experiments::e04 {
    namespace e04RunDetail {
        using namespace e03;
        using Solver                           = SphericalSimulation::Solver;
        inline constexpr double w              = reference::w;
        inline constexpr double finalEta       = 16;
        inline constexpr double precursorStart = 4.5;
        inline constexpr std::array targets{2., 3., 4., 5., 6., 8., 10., 12., finalEta};

        /// Label one row at its actual native time, independently of the target eta.
        struct Observation {
            long long step_m;
            double time_m, eta_m, targetEta_m;
        };

        /// The three E04 CSV files on rank zero. The constructor first prepares the directory.
        class E04Csv {
        public:
            E04Csv(const std::string& directory, MPI_Comm comm) {
                output::prepareRunDirectory(directory, comm);
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
        inline double energyFromNorms(double electric, double magnetic) {
            // Unit permeability [ellingson2020ElectromagneticsVol2, Eqs. (3.20), (3.22)].
            return 0.5 * (E03Epsilon * electric * electric + magnetic * magnetic);
        }

        struct InitialNorms {
            E03FieldNorms e_m, b_m, a_m;

            double referenceEnergy() const {
                return energyFromNorms(e_m.referenceL2_m, b_m.referenceL2_m);
            }
        };

        template <class Field, DipolePulseReference Reference>
        E03FieldNorms norms(Solver& state, Field& field, const Reference& ref, Quantity quantity,
                            long long step) {
            return measureE03Field(field, ref, state.levelContract(), quantity,
                                   state.labelTime(quantity, step));
        }

        /// The constraint scales and energy denominator are fixed at initialization.
        template <DipolePulseReference Reference>
        InitialNorms initialNorms(Solver& state, const Reference& ref) {
            return {norms(state, state.electricField(), ref, Quantity::Electric, 0),
                    norms(state, state.magneticField(), ref, Quantity::Magnetic, 0),
                    norms(state, state.A_n, ref, Quantity::VectorPotential, 0)};
        }

        /// At startup and eta=2 compare with both the exact pulse and the periodic replay.
        inline void compareControl(Solver& state, const E04PeriodicControl& control, E04Csv& csv,
                                   const Observation& at, const E03FieldNorms& e,
                                   const E03FieldNorms& b) {
            const auto difference = control.difference(state, at.step_m);
            csv.row("E04.control", "E_continuum_L2", "box", at, e.errorL2_m, e.referenceL2_m);
            csv.row("E04.control", "B_continuum_L2", "box", at, b.errorL2_m, b.referenceL2_m);
            csv.row("E04.control", "E_periodic_L2", "box", at, difference[0], e.referenceL2_m);
            csv.row("E04.control", "B_periodic_L2", "box", at, difference[1], b.referenceL2_m);
        }

        /// Observe the field difference before the outgoing continuum support leaves the cube.
        template <DipolePulseReference Reference>
        void observeFields(Solver& state, const Reference& ref, const E04PeriodicControl& control,
                           E04Csv& csv, const Observation& at, double U0h) {
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

        inline Observation observation(Solver& state, Quantity quantity, long long step, double v,
                                       double targetEta) {
            const double time = state.labelTime(quantity, step);
            return {step, time, v * time / w, targetEta};
        }

        /// Select field and gauge observations with their own native clocks.
        inline std::array<Checkpoint, targets.size()> checkpoints(Solver& state, double v) {
            std::array<Checkpoint, targets.size()> result{};
            for (unsigned i = 0; i < targets.size(); ++i) {
                const double time = targets[i] * w / v;
                result[i] = {targets[i], nearestSphericalStep(state, Quantity::Electric, time),
                             nearestSphericalStep(state, Quantity::VectorPotential, time)};
            }
            return result;
        }

        /// Evaluate only the target observables selected by this state.
        template <DipolePulseReference Reference>
        void observeTargets(Solver& state, const Reference& ref, const E04PeriodicControl& control,
                            E04Csv& csv, const InitialNorms& initial,
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
        inline void recordMurLayers(conduit_cpp::Node node, Solver& state, double h) {
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

        /// The values that recordE04 writes besides the reference sampler's own entries.
        struct E04RecordValues {
            int steps_m, protocolSteps_m;
            double U0h_m;
            int controlStep_m;
            bool planeFrames_m;
        };

        template <DipolePulseReference Reference>
        void recordE04(const input::ChdrYamlRunInput& request,
                       const simulation::ChdrCartesianGridSettings& grid, Solver& state,
                       const E04RecordValues& values, const Reference& ref,
                       const SphericalExperiment& experiment,
                       const SphericalFieldSnapshots& snapshots) {
            auto record = output::describeRun(
                experiment.application_m, grid, state, values.steps_m,
                std::string(experiment.name_m) + " spherical pulse with all-face Mur");
            snapshots.record(record);
            auto& node                             = record.effective_m;
            node["protocol_steps"]                 = values.protocolSteps_m;
            node["material/relative_permittivity"] = E03Epsilon;
            ref.record(node["analytic_reference"]);
            node["boundary/mur_mask"]                = state.boundary().murMask();
            node["boundary/kernel_speed"]            = state.boundary().speed[0];
            node["boundary/potential_ghost_closure"] = "Quadratic extrapolation";
            recordMurLayers(node["boundary/potential_layers_centred"], state, grid.spacing_m[0]);
            node["initial_reference_energy"]           = values.U0h_m;
            node["periodic_control/target_eta"]        = 2.0;
            node["periodic_control/target_field_step"] = values.controlStep_m;
            if (values.planeFrames_m)
                node["plane/file"] = "plane-z0.csv";
            output::writeRunRecord(request.outputDirectory_m, request, record,
                                   state.A_n.getLayout().comm);
        }

        /// Collect the dense E04.meas01/02 histories at the common R2 half-step.
        inline void observeTraces(Solver& state, const material::MaterialOwner& owner, E04Csv& csv,
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
        inline int bracketStep(Solver& state, double eta, double v, bool upper) {
            const double step =
                (eta * w / v - state.labelTime(Quantity::Electric, 0)) / state.getDt();
            return static_cast<int>(upper ? std::ceil(step) : std::floor(step));
        }

        /// Save the startup control, field and constraint rows before the first update.
        template <DipolePulseReference Reference>
        void observeStartup(Solver& state, const Reference& ref, const E04PeriodicControl& control,
                            E04Csv& csv, const InitialNorms& initial) {
            compareControl(state, control, csv,
                           observation(state, Quantity::Electric, 0, ref.v_m, 0), initial.e_m,
                           initial.b_m);
            observeConstraint<true>(state, csv,
                                    observation(state, Quantity::Electric, 0, ref.v_m, 0), initial);
            observeConstraint<false>(
                state, csv, observation(state, Quantity::VectorPotential, 0, ref.v_m, 0), initial);
        }

        template <DipolePulseReference Reference>
        void evolveE04(SphericalSimulation& run, const input::ChdrYamlRunInput& request,
                       const simulation::ChdrCartesianGridSettings& grid, const Reference& ref,
                       const SphericalExperiment& experiment) {
            auto& state           = run.solver();
            const auto points     = checkpoints(state, ref.v_m);
            const int fullSteps   = bracketStep(state, finalEta, ref.v_m, true);
            const int steps       = sphericalRunSteps(request, fullSteps, "E04");
            const int controlStep = points.front().fieldStep_m;
            E04Csv csv(request.outputDirectory_m, state.A_n.getLayout().comm);
            E04PeriodicControl control(grid, ref, std::min(steps, controlStep), controlStep);
            const auto initial         = initialNorms(state, ref);
            const double U0h           = initial.referenceEnergy();
            const int firstMonitorStep = bracketStep(state, precursorStart, ref.v_m, false);
            observeStartup(state, ref, control, csv, initial);
            SphericalFieldSnapshots snapshots(request, grid, state);
            const SphericalPlaneFrames planes(request);
            run.advance(steps, [&](long long step, Solver& solver, const auto& owner) {
                observeTraces(solver, owner, csv, step, firstMonitorStep, ref.v_m, U0h);
                observeTargets(solver, ref, control, csv, initial, points, step);
                snapshots.observe(step, solver, steps);
                planes.observe(step, solver, ref, steps);
            });
            csv.finish();
            snapshots.finish();
            recordE04(request, grid, state, {steps, fullSteps, U0h, controlStep, planes.writes()},
                      ref, experiment, snapshots);
        }
    }  // namespace e04RunDetail

    /** @brief Compare periodic propagation and matched-speed Mur exit as specified in E04.typ.
     * @param ref Continuum sampler of the spherical pulse, for example e03::E03Reference.
     * @param experiment Names for the run record.
     * @pre The caller has validated the experiment's parameter names.
     */
    template <DipolePulseReference Reference>
    void runE04(const input::ChdrYamlRunInput& request, const Reference& ref,
                const e03::SphericalExperiment& experiment) {
        e03::validateSphericalRun(request);
        const auto grid = e03::sphericalGrid(e03::sphericalCells(request), false);
        auto run        = e03::makeSphericalSimulation(grid, ref);
        run.load(ref);
        e04RunDetail::evolveE04(run, request, grid, ref, experiment);
    }
}  // namespace chdr::experiments::e04
