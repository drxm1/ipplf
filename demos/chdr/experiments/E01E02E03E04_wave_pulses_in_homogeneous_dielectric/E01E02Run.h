/** @file E01E02Run.h
 * @brief Shared E01/E02 observation planning, evolution and output.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/input/RunInput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/ChDRSimulationManager.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02FieldSnapshots.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Measurements.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Profiles.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Protocol.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E02Diagnostics.h"

namespace chdr::experiments::e01e02 {
    using PulseSimulation = simulation::ChDRSimulationManager<PulseSolver>;

    /// Require an exact positive integer before narrowing a YAML scalar.
    inline int pulseInteger(const input::ChdrYamlRunInput& request, const char* key) {
        if (!request.parameters_m.has_path(key))
            throw std::invalid_argument(std::string("Missing pulse parameter: ") + key);
        const auto node = request.parameters_m[key];
        if (!node.dtype().is_number() || node.number_of_elements() != 1)
            throw std::invalid_argument(std::string(key) + " must be an integer");
        const double value = node.to_float64();
        if (!std::isfinite(value) || value < 2 || value > std::numeric_limits<int>::max()
            || std::floor(value) != value)
            throw std::invalid_argument(std::string(key) + " must be an integer >= 2");
        return static_cast<int>(value);
    }

    /// Select only resolutions belonging to the named experiment's ladder.
    inline int pulseCellsPerPulse(const input::ChdrYamlRunInput& request,
                                  const PulseProtocol& protocol) {
        const int value = pulseInteger(request, "cells_per_pulse");
        if (std::find(protocol.resolutions_m.begin(), protocol.resolutions_m.end(), value)
            == protocol.resolutions_m.end())
            throw std::invalid_argument(std::string(protocol.name_m)
                                        + " cells_per_pulse is outside its refinement ladder");
        return value;
    }

    /// Preserve E01's selectable transverse count and E02's fixed eight-by-eight mesh.
    inline int pulseTransverseCells(const input::ChdrYamlRunInput& request,
                                    const PulseProtocol& protocol) {
        const int value = pulseInteger(request, "transverse_cells_per_axis");
        if (protocol.transverseCells_m && value != protocol.transverseCells_m)
            throw std::invalid_argument(std::string(protocol.name_m)
                                        + " requires transverse_cells_per_axis = "
                                        + std::to_string(protocol.transverseCells_m));
        return value;
    }

    /// Validate common input choices before allocating fields or opening output.
    inline void validatePulseInput(const input::ChdrYamlRunInput& request,
                                   const PulseProtocol& protocol) {
        if (request.caseName_m != protocol.name_m || request.units_m.system_m != "normalized")
            throw std::invalid_argument(std::string(protocol.name_m)
                                        + " requires its named case and normalized units");
        input::validateParameters(request, {"cells_per_pulse", "transverse_cells_per_axis"});
        input::validateObservationNames(request, {"fields"});
        pulseCellsPerPulse(request, protocol);
        pulseTransverseCells(request, protocol);
        if (request.outputDirectory_m.empty()
            && (!protocol.periodicZ_m
                || (request.observations_m.contains("fields")
                    && input::observationEnabled(request, "fields"))))
            throw std::invalid_argument("Pulse field/probe output requires output.directory");
        if (request.duration_m || (request.steps_m && *request.steps_m < 0))
            throw std::invalid_argument("Pulse runs use fixed targets and nonnegative run.steps");
    }

    namespace pulseRunDetail {
        /// One prescribed target has separate nearest integer-potential and half-step-field states.
        struct Checkpoint {
            const char* name_m;
            double travelWidths_m, targetTime_m;
            int potentialStep_m, fieldStep_m;
            bool continuumErrors_m, fieldSnapshot_m;
        };

        /// The full protocol schedule stays fixed when run.steps requests an earlier stop.
        struct RunPlan {
            std::vector<Checkpoint> checkpoints_m;
            std::vector<int> snapshotSteps_m;
            int fullSteps_m, steps_m, defaultSnapshotEvery_m, lastErrorStep_m = 0;
        };

        /** @brief Select the nearest available state using the actual sampling contract.
         * @details The fixed-step labels are supplied by LevelContract. Subtracting the
         * zero-state label before rounding inverts that map without duplicating its offsets.
         */
        inline int nearestObservationStep(double target, Quantity quantity,
                                          const LevelContract& contract, double dt) {
            return static_cast<int>(
                std::floor((target - contract.labelTime(quantity, 0, dt)) / dt + 0.5));
        }

        /// Convert an experiment target using the reference speed and the actual level contract.
        inline Checkpoint observationCheckpoint(const PulseTarget& target,
                                                const E01Reference& reference,
                                                const LevelContract& contract, double dt) {
            const double time = target.travelWidths_m * e01Reference::w / reference.v_m;
            return {target.name_m,
                    target.travelWidths_m,
                    time,
                    nearestObservationStep(time, Quantity::VectorPotential, contract, dt),
                    nearestObservationStep(time, Quantity::Electric, contract, dt),
                    target.continuumErrors_m,
                    target.fieldSnapshot_m};
        }

        /// Resolve the protocol independently of an optional short-run limit.
        inline void appendCheckpoint(RunPlan& plan, const Checkpoint& point) {
            plan.checkpoints_m.push_back(point);
            if (point.continuumErrors_m)
                plan.lastErrorStep_m = point.fieldStep_m;
            if (point.fieldSnapshot_m)
                plan.snapshotSteps_m.push_back(point.fieldStep_m);
        }

        /// Resolve the protocol independently of an optional short-run limit.
        inline RunPlan runPlan(const input::ChdrYamlRunInput& request,
                               const E01Reference& reference, const LevelContract& contract,
                               double dt, const PulseProtocol& protocol) {
            RunPlan plan{};
            for (const auto& target : protocol.targets_m)
                appendCheckpoint(plan, observationCheckpoint(target, reference, contract, dt));
            const auto& last = plan.checkpoints_m.back();
            plan.fullSteps_m = std::max(last.potentialStep_m, last.fieldStep_m);
            plan.defaultSnapshotEvery_m =
                protocol.periodicZ_m
                    ? nearestObservationStep(last.targetTime_m / 4, Quantity::VectorPotential,
                                             contract, dt)
                    : plan.fullSteps_m + 1;
            plan.steps_m = request.steps_m.value_or(plan.fullSteps_m);
            if (plan.steps_m < 0 || plan.steps_m > plan.fullSteps_m)
                throw std::invalid_argument("Pulse run.steps exceeds the prescribed run interval");
            return plan;
        }

        /** @brief Write rank-zero error rows while every rank participates in measurements.
         * @details A requested output directory is prepared collectively. Filesystem failures
         * throw to the application's MPI-aware failure handler; earlier CSV rows are retained.
         */
        class ErrorCsv {
        public:
            ErrorCsv(const std::string& directory,
                     const simulation::ChdrCartesianGridSettings& grid, double dt,
                     MPI_Comm communicator, const char* caseName)
                : caseName_m(caseName)
                , grid_m(grid)
                , dt_m(dt) {
                MPI_Comm_rank(communicator, &rank_m);
                output::prepareRunDirectory(directory, communicator);
                if (rank_m == 0 && !directory.empty()) {
                    file_m.exceptions(std::ios::failbit | std::ios::badbit);
                    file_m.imbue(std::locale::classic());
                    file_m.open(std::filesystem::path(directory) / "errors.csv");
                    stream_m = &file_m;
                }
                if (rank_m == 0)
                    *stream_m << std::setprecision(std::numeric_limits<double>::max_digits10)
                              << "measurement,checkpoint,mesh_id,nx,ny,nz,hx,hy,hz,dt,"
                                 "target_time,step,actual_time,relative_error\n";
            }

            /// A CSV row records the continuum comparison time independently of its target.
            void row(const char* measurement, const char* checkpoint, double target, long long step,
                     double actual, double error) {
                if (rank_m != 0)
                    return;
                const auto& n = grid_m.cells_m;
                const auto& h = grid_m.spacing_m;
                *stream_m << measurement << ',' << checkpoint << ',' << caseName_m << '_' << n[0]
                          << 'x' << n[1] << 'x' << n[2] << ',' << n[0] << ',' << n[1] << ',' << n[2]
                          << ',' << h[0] << ',' << h[1] << ',' << h[2] << ',' << dt_m << ','
                          << target << ',' << step << ',' << actual << ',' << error << std::endl;
                if (!*stream_m)
                    throw std::runtime_error("Cannot write pulse error CSV");
            }

            /// Finalize the error file before publishing a completed run record.
            void finish() {
                if (rank_m != 0)
                    return;
                stream_m->flush();
                if (!*stream_m)
                    throw std::runtime_error("Cannot flush pulse error CSV");
                if (file_m.is_open())
                    file_m.close();
            }

        private:
            const char* caseName_m;
            const simulation::ChdrCartesianGridSettings& grid_m;
            double dt_m;
            int rank_m = 0;
            std::ofstream file_m;
            std::ostream* stream_m = &std::cout;
        };

        /// Measure E01.meas04 at the integer time of the current vector-potential slot.
        inline void observePotential(PulseSolver& state, const E01Reference& reference,
                                     const char* checkpoint, double target, long long step,
                                     ErrorCsv& csv, PulseProfileCsv& profiles,
                                     bool continuumErrors = true) {
            const double time = state.labelTime(Quantity::VectorPotential, step);
            if (continuumErrors) {
                const auto error = measurePulseFieldError(
                    state.A_n, reference, state.levelContract(), Quantity::VectorPotential, time);
                csv.row("A_L2", checkpoint, target, step, time, error.relativeL2_m);
            }
            profiles.write(state.A_n, Quantity::VectorPotential, checkpoint, target, step, time);
        }

        /// Measure E01.meas01-03 using each reconstructed field's actual R2 time label.
        inline void observeFields(PulseSolver& state, const E01Reference& reference,
                                  const char* checkpoint, double target, long long step,
                                  ErrorCsv& csv, PulseProfileCsv& profiles,
                                  bool continuumErrors = true) {
            const double electricTime = state.labelTime(Quantity::Electric, step);
            const double magneticTime = state.labelTime(Quantity::Magnetic, step);
            if (continuumErrors) {
                const auto electric =
                    measurePulseFieldError(state.electricField(), reference, state.levelContract(),
                                           Quantity::Electric, electricTime);
                const auto magnetic =
                    measurePulseFieldError(state.magneticField(), reference, state.levelContract(),
                                           Quantity::Magnetic, magneticTime);
                csv.row("E_L2", checkpoint, target, step, electricTime, electric.relativeL2_m);
                csv.row("E_Linf", checkpoint, target, step, electricTime, electric.relativeLinf_m);
                csv.row("B_L2", checkpoint, target, step, magneticTime, magnetic.relativeL2_m);
            }
            profiles.write(state.electricField(), Quantity::Electric, checkpoint, target, step,
                           electricTime);
            profiles.write(state.magneticField(), Quantity::Magnetic, checkpoint, target, step,
                           magneticTime);
        }

        /// Track which quantity groups already have a scheduled row at the current state.
        struct ObservedGroups {
            bool potential_m = false, fields_m = false;
        };

        /// All ranks enter the same measurement branches, including their global reductions.
        inline ObservedGroups observeScheduled(PulseSolver& state, const E01Reference& reference,
                                               const RunPlan& plan, long long step, ErrorCsv& csv,
                                               PulseProfileCsv& profiles) {
            ObservedGroups observed;
            for (const auto& checkpoint : plan.checkpoints_m) {
                if (step == checkpoint.potentialStep_m) {
                    observePotential(state, reference, checkpoint.name_m, checkpoint.targetTime_m,
                                     step, csv, profiles, checkpoint.continuumErrors_m);
                    observed.potential_m = true;
                }
                if (step == checkpoint.fieldStep_m) {
                    observeFields(state, reference, checkpoint.name_m, checkpoint.targetTime_m,
                                  step, csv, profiles, checkpoint.continuumErrors_m);
                    observed.fields_m = true;
                }
            }
            return observed;
        }

        /** @brief Keep initialization and short-run endpoints separate from protocol targets.
         * @details Initial fields have time -dt/2. A partial_final diagnostic uses the final
         * stored state, with integer potential time as its target and actual per-field times
         * for comparisons. Scheduled quantity rows at that state are not duplicated.
         */
        inline void observe(PulseSolver& state, const E01Reference& reference, const RunPlan& plan,
                            long long step, ErrorCsv& csv, PulseProfileCsv& profiles) {
            if (step == 0) {
                observePotential(state, reference, "initial", 0, step, csv, profiles);
                observeFields(state, reference, "initial", 0, step, csv, profiles);
                return;
            }
            const auto observed = observeScheduled(state, reference, plan, step, csv, profiles);
            if (step != plan.steps_m || plan.steps_m == plan.fullSteps_m)
                return;
            const double target        = state.labelTime(Quantity::VectorPotential, step);
            const bool continuumErrors = step <= plan.lastErrorStep_m;
            if (!observed.potential_m)
                observePotential(state, reference, "partial_final", target, step, csv, profiles,
                                 continuumErrors);
            if (!observed.fields_m)
                observeFields(state, reference, "partial_final", target, step, csv, profiles,
                              continuumErrors);
        }

        /// Record the analytic input values without inferring them from numerical field samples.
        inline void recordReference(output::RunRecordDescription& description,
                                    const E01Reference& reference) {
            auto node               = description.effective_m["analytic_reference"];
            node["name"]            = reference.periodicZ_m ? "E01 periodic compact continuum pulse"
                                                            : "E02 unwrapped compact continuum pulse";
            node["periodic_z"]      = reference.periodicZ_m;
            node["domain_length_z"] = reference.Lz_m;
            if (reference.periodicZ_m)
                node["period"] = reference.Lz_m;
            node["initial_centre_z"]    = reference.z0_m;
            node["wave_speed"]          = reference.v_m;
            node["potential_amplitude"] = e01Reference::A0;
            node["support_half_width"]  = e01Reference::w;
        }

        /// Record each target's name and distance alongside its selected native state indices.
        inline void recordTargets(output::RunRecordDescription& description, const RunPlan& plan,
                                  const char* key, bool errorsOnly) {
            auto node = description.effective_m[key];
            std::vector<double> targets, distances;
            std::vector<int> potentialSteps, fieldSteps;
            for (const auto& checkpoint : plan.checkpoints_m) {
                if (errorsOnly && !checkpoint.continuumErrors_m)
                    continue;
                node["protocol_checkpoint_names"].append().set(std::string(checkpoint.name_m));
                targets.push_back(checkpoint.targetTime_m);
                distances.push_back(checkpoint.travelWidths_m);
                potentialSteps.push_back(checkpoint.potentialStep_m);
                fieldSteps.push_back(checkpoint.fieldStep_m);
            }
            node["protocol_travel_distances_in_w"].set(distances);
            node["protocol_target_times"].set(targets);
            node["protocol_potential_steps"].set(potentialSteps);
            node["protocol_field_steps"].set(fieldSteps);
        }

        /// Keep the prescribed schedule visible even when a partial run stops before its targets.
        inline void recordSchedule(output::RunRecordDescription& description, const RunPlan& plan) {
            recordTargets(description, plan, "error_observations", true);
            recordTargets(description, plan, "profile_observations", false);
            description.effective_m["protocol_steps"] = plan.fullSteps_m;
            auto node                   = description.effective_m["error_observations"];
            node["initial_step"]        = 0;
            node["partial_final"]       = plan.steps_m < plan.fullSteps_m;
            node["file"]                = "errors.csv";
            node["norms"]               = "A_L2,E_L2,E_Linf,B_L2; full-domain relative errors";
            node["magnetic_convention"] = "R2 averaged magnetic field";
        }

        /// Export native comparison times so analysis need not duplicate LevelContract rules.
        inline void recordSampleTimes(output::RunRecordDescription& description,
                                      const PulseSolver& state, const RunPlan& plan,
                                      const char* key, bool errorsOnly) {
            const std::array quantities{Quantity::VectorPotential, Quantity::Electric,
                                        Quantity::Magnetic};
            const std::array names{"A", "E", "B"};
            for (unsigned c = 0; c < quantities.size(); ++c) {
                std::vector<double> times{state.labelTime(quantities[c], 0)};
                for (const auto& checkpoint : plan.checkpoints_m) {
                    if (errorsOnly && !checkpoint.continuumErrors_m)
                        continue;
                    const int step = quantities[c] == Quantity::VectorPotential
                                         ? checkpoint.potentialStep_m
                                         : checkpoint.fieldStep_m;
                    times.push_back(state.labelTime(quantities[c], step));
                }
                description.effective_m[std::string(key) + "/sample_times/" + names[c]].set(times);
            }
        }

        /// Identify the fixed native line and its connection to the error observations.
        inline void recordProfiles(output::RunRecordDescription& description,
                                   const simulation::ChdrCartesianGridSettings& grid) {
            auto profiles                   = description.effective_m["profiles"];
            profiles["file"]                = "profiles.csv";
            profiles["global_cell_i"]       = pulseProfileDetail::lineCell(grid.cells_m[0]);
            profiles["global_cell_j"]       = pulseProfileDetail::lineCell(grid.cells_m[1]);
            profiles["components"]          = "Ax,Ex,By";
            profiles["magnetic_convention"] = "R2 averaged magnetic field";
            profiles["schedule"] = "profile_observations supplies native targets and times";
        }

        /// Record the actual face mask and kernel speeds resolved by the constructed solver.
        inline void recordBoundary(output::RunRecordDescription& description,
                                   const PulseSolver& state) {
            const auto& boundary = state.boundary();
            auto node            = description.effective_m["boundary"];
            node["face_order"]   = "-x,+x,-y,+y,-z,+z";
            node["mur_mask"]     = boundary.murMask();
            node["kernel_speeds"].set(
                std::vector<double>(boundary.speed.begin(), boundary.speed.end()));
            node["speed_units"] = "physical solver units for O1VS; periodic-face entries unused";
        }

        /// Publish the common provenance record only after the selected run and CSVs finish.
        inline void recordRun(PulseSimulation& run, const input::ChdrYamlRunInput& request,
                              const simulation::ChdrCartesianGridSettings& grid,
                              const E01Reference& reference, const RunPlan& plan,
                              const PulseFieldSnapshots& snapshots, const PulseProtocol& protocol,
                              const std::optional<E02Diagnostics>& diagnostics) {
            if (request.outputDirectory_m.empty())
                return;
            auto description = output::describeRun(
                "experiments-wave-pulses-in-homogeneous-dielectric", grid, run.solver(),
                plan.steps_m, std::string(protocol.name_m) + " homogeneous pulse refinement");
            description.effective_m["material/relative_permittivity"] =
                run.material().homogeneousValue();
            recordReference(description, reference);
            recordBoundary(description, run.solver());
            recordSchedule(description, plan);
            recordSampleTimes(description, run.solver(), plan, "error_observations", true);
            recordSampleTimes(description, run.solver(), plan, "profile_observations", false);
            recordProfiles(description, grid);
            snapshots.record(description);
            if (diagnostics)
                diagnostics->record(description);
            output::writeRunRecord(request.outputDirectory_m, request, description,
                                   run.solver().A_n.getLayout().comm);
        }
    }  // namespace pulseRunDetail

    /** @brief Advance one pulse experiment with its prescribed observations.
     * @details The manager observes state zero and every completed solve. Expensive reductions
     * occur only at the independently selected potential/field checkpoints and, for a short
     * run, its additional endpoint. E01 and E02 select their own schedules and additional
     * diagnostics.
     * @pre validatePulseInput has succeeded and the analytic histories are already loaded.
     */
    inline void evolvePulse(PulseSimulation& run, const input::ChdrYamlRunInput& request,
                            const simulation::ChdrCartesianGridSettings& grid,
                            const E01Reference& reference, const PulseProtocol& protocol) {
        const auto plan = pulseRunDetail::runPlan(request, reference, run.solver().levelContract(),
                                                  run.solver().getDt(), protocol);
        pulseRunDetail::ErrorCsv csv(request.outputDirectory_m, grid, run.solver().getDt(),
                                     run.solver().A_n.getLayout().comm, protocol.name_m);
        PulseProfileCsv profiles(request.outputDirectory_m, reference, run.solver().levelContract(),
                                 run.solver().A_n.getLayout().comm);
        PulseFieldSnapshots snapshots(request, grid, run.solver(), plan.defaultSnapshotEvery_m,
                                      plan.snapshotSteps_m);
        std::optional<E02Diagnostics> diagnostics;
        if (!protocol.periodicZ_m)
            diagnostics.emplace(request.outputDirectory_m, run.solver(), grid);
        run.advance(plan.steps_m, [&](long long step, PulseSolver& state, const auto&) {
            pulseRunDetail::observe(state, reference, plan, step, csv, profiles);
            snapshots.observe(step, state, plan.steps_m);
            if (diagnostics)
                diagnostics->observe(step, state, plan.steps_m, plan.fullSteps_m,
                                     plan.checkpoints_m.back().targetTime_m);
        });
        csv.finish();
        profiles.finish();
        snapshots.finish();
        if (diagnostics)
            diagnostics->finish();
        pulseRunDetail::recordRun(run, request, grid, reference, plan, snapshots, protocol,
                                  diagnostics);
    }

    /// Build either physical box with the same transverse placement and z-only decomposition.
    inline simulation::ChdrCartesianGridSettings pulseGrid(const input::ChdrYamlRunInput& request,
                                                           const PulseProtocol& protocol) {
        const int longitudinal = pulseCellsPerPulse(request, protocol);
        const int transverse   = pulseTransverseCells(request, protocol);
        return {{transverse, transverse, protocol.longitudinalCellsPerPulse_m * longitudinal},
                {1.0 / transverse, 1.0 / transverse, 2.0 / longitudinal},
                0,
                {true, true, protocol.periodicZ_m},
                2};
    }

    /** @brief Select the common discretization and E02's explicit physical Mur speed.
     * @details O1VS uses the physical timestep. Its dielectric wave speed is
     * \f$1/\sqrt{\epsilon_r}\f$ [chew2014generalizedgauge, Eq. (13)], so the
     * unscaled kernel default of one is not the E02 boundary speed.
     */
    inline simulation::ChdrSolverSettings pulseSolverSettings(double dt, double speed,
                                                              const PulseProtocol& protocol) {
        return {.reconstruction_m = EBReconstructionConvention::R2,
                .spatial_m        = SpatialDiscretization::Staggered,
                .formulation_m =
                    solver::FormulationParameters{
                        .formulation_m = solver::MaterialFormulation::O1VS, .alpha_m = 1.0},
                .murSpeed_m = protocol.periodicZ_m ? solver::MurSpeed::matched()
                                                   : solver::MurSpeed::physical(speed),
                .timeStep_m = dt};
    }

    /// Construct, initialize and evolve one named experiment without duplicating solver setup.
    inline void runPulse(const input::ChdrYamlRunInput& request, const PulseProtocol& protocol) {
        validatePulseInput(request, protocol);
        const auto grid = pulseGrid(request, protocol);
        const E01Reference reference{grid.lengths()[2], 2.0, 1.0 / std::sqrt(2.13),
                                     protocol.periodicZ_m};
        const auto homogeneousDielectric =
            KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
            marked = false;
            return 2.13;
        };
        PulseSimulation run{grid, homogeneousDielectric,
                            pulseSolverSettings(grid.spacing_m[2] / 2, reference.v_m, protocol)};
        run.load(reference);
        evolvePulse(run, request, grid, reference, protocol);
    }
}  // namespace chdr::experiments::e01e02
