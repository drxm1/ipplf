/** @file E01Run.h
 * @brief E01 observation times, evolution, global errors, profiles and field snapshots.
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
#include <stdexcept>
#include <string>
#include <vector>

#include "core/input/RunInput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/ChDRSimulationManager.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01FieldSnapshots.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01Measurements.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01Profiles.h"

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric {
    using E01Simulation = simulation::ChDRSimulationManager<E01Solver>;

    /** @brief Select the number of longitudinal cells across the full E01 pulse width.
     * @details E01.typ specifies 32, 64, 128 and 256, with 512 as an optional extension.
     * The parameter is required. No rounding or narrowing of other values is used.
     */
    inline int e01CellsPerPulse(const input::ChdrYamlRunInput& request) {
        if (!request.parameters_m.has_path("cells_per_pulse"))
            throw std::invalid_argument("E01 requires parameters.cells_per_pulse");
        const auto node = request.parameters_m["cells_per_pulse"];
        if (!node.dtype().is_number() || node.number_of_elements() != 1)
            throw std::invalid_argument("E01 parameters.cells_per_pulse must be a numeric scalar");
        const double value = node.to_float64();
        for (int allowed : std::array{32, 64, 128, 256, 512})
            if (value == allowed)
                return allowed;
        throw std::invalid_argument(
            "E01 parameters.cells_per_pulse must be 32, 64, 128, 256 or 512");
    }

    /// Reject unsupported E01 case, unit, parameter and observation options.
    inline void validateE01Input(const input::ChdrYamlRunInput& request) {
        if (request.caseName_m != "E01")
            throw std::invalid_argument("This executable currently supports case E01 only");
        if (request.units_m.system_m != "normalized")
            throw std::invalid_argument("E01 currently requires normalized units");
        input::validateParameters(request, {"cells_per_pulse"});
        input::validateObservationNames(request, {"fields"});
        if (request.outputDirectory_m.empty() && request.observations_m.contains("fields")
            && input::observationEnabled(request, "fields"))
            throw std::invalid_argument("E01 field snapshots require an output directory");
        if (request.duration_m)
            throw std::invalid_argument("E01 fixes its target time; use run.steps to stop early");
        if (request.steps_m && *request.steps_m < 0)
            throw std::invalid_argument("E01 run.steps must be nonnegative");
    }

    namespace e01_run_detail {
        /// One prescribed target has separate nearest integer-potential and half-step-field states.
        struct Checkpoint {
            const char* name_m;
            double targetTime_m;
            int potentialStep_m, fieldStep_m;
        };

        /// The full protocol schedule stays fixed when run.steps requests an earlier stop.
        struct RunPlan {
            std::array<Checkpoint, 4> checkpoints_m;
            int fullSteps_m, steps_m;
        };

        /** @brief Resolve the four targets and the last required solver state.
         * @details E01.typ prescribes travel distance 40 and targets at its quarter intervals.
         * Its nearest-time rule uses \f$n_A=\lfloor\tau/h_t+1/2\rfloor\f$ and
         * \f$n_{EB}=\lfloor\tau/h_t+1\rfloor\f$. These follow from the integer potential
         * and R2 half-step labels [fallahi2020mithra20fullwavesimulation, Eqs. (3.50)-(3.55)].
         */
        inline RunPlan runPlan(const input::ChdrYamlRunInput& request,
                               const E01Reference& reference, double dt) {
            RunPlan plan{};
            constexpr std::array names{"quarter_1", "quarter_2", "quarter_3", "final"};
            const double finalTarget = 40.0 / reference.v_m;
            for (unsigned q = 0; q < plan.checkpoints_m.size(); ++q) {
                const double target   = (q + 1) * finalTarget / 4.0;
                plan.checkpoints_m[q] = {names[q], target, int(std::floor(target / dt + 0.5)),
                                         int(std::floor(target / dt + 1.0))};
            }
            const auto& last = plan.checkpoints_m.back();
            plan.fullSteps_m = std::max(last.potentialStep_m, last.fieldStep_m);
            plan.steps_m     = request.steps_m.value_or(plan.fullSteps_m);
            if (plan.steps_m < 0 || plan.steps_m > plan.fullSteps_m)
                throw std::invalid_argument("E01 run.steps exceeds the prescribed run interval");
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
                     MPI_Comm communicator)
                : grid_m(grid)
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
                *stream_m << measurement << ',' << checkpoint << ",E01_" << n[0] << 'x' << n[1]
                          << 'x' << n[2] << ',' << n[0] << ',' << n[1] << ',' << n[2] << ',' << h[0]
                          << ',' << h[1] << ',' << h[2] << ',' << dt_m << ',' << target << ','
                          << step << ',' << actual << ',' << error << std::endl;
                if (!*stream_m)
                    throw std::runtime_error("Cannot write E01 error CSV");
            }

            /// Finalize the error file before publishing a completed run record.
            void finish() {
                if (rank_m != 0)
                    return;
                stream_m->flush();
                if (!*stream_m)
                    throw std::runtime_error("Cannot flush E01 error CSV");
                if (file_m.is_open())
                    file_m.close();
            }

        private:
            const simulation::ChdrCartesianGridSettings& grid_m;
            double dt_m;
            int rank_m = 0;
            std::ofstream file_m;
            std::ostream* stream_m = &std::cout;
        };

        /// Measure E01.meas04 at the integer time of the current vector-potential slot.
        inline void observePotential(E01Solver& state, const E01Reference& reference,
                                     const char* checkpoint, double target, long long step,
                                     ErrorCsv& csv, E01ProfileCsv& profiles) {
            const double time = state.labelTime(Quantity::VectorPotential, step);
            const auto error  = measureE01FieldError(state.A_n, reference, state.levelContract(),
                                                     Quantity::VectorPotential, time);
            csv.row("A_L2", checkpoint, target, step, time, error.relativeL2_m);
            profiles.write(state.A_n, Quantity::VectorPotential, checkpoint, target, step, time);
        }

        /// Measure E01.meas01-03 using each reconstructed field's actual R2 time label.
        inline void observeFields(E01Solver& state, const E01Reference& reference,
                                  const char* checkpoint, double target, long long step,
                                  ErrorCsv& csv, E01ProfileCsv& profiles) {
            const double electricTime = state.labelTime(Quantity::Electric, step);
            const double magneticTime = state.labelTime(Quantity::Magnetic, step);
            const auto electric =
                measureE01FieldError(state.electricField(), reference, state.levelContract(),
                                     Quantity::Electric, electricTime);
            const auto magnetic =
                measureE01FieldError(state.magneticField(), reference, state.levelContract(),
                                     Quantity::Magnetic, magneticTime);
            csv.row("E_L2", checkpoint, target, step, electricTime, electric.relativeL2_m);
            csv.row("E_Linf", checkpoint, target, step, electricTime, electric.relativeLinf_m);
            csv.row("B_L2", checkpoint, target, step, magneticTime, magnetic.relativeL2_m);
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
        inline ObservedGroups observeScheduled(E01Solver& state, const E01Reference& reference,
                                               const RunPlan& plan, long long step, ErrorCsv& csv,
                                               E01ProfileCsv& profiles) {
            ObservedGroups observed;
            for (const auto& checkpoint : plan.checkpoints_m) {
                if (step == checkpoint.potentialStep_m) {
                    observePotential(state, reference, checkpoint.name_m, checkpoint.targetTime_m,
                                     step, csv, profiles);
                    observed.potential_m = true;
                }
                if (step == checkpoint.fieldStep_m) {
                    observeFields(state, reference, checkpoint.name_m, checkpoint.targetTime_m,
                                  step, csv, profiles);
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
        inline void observe(E01Solver& state, const E01Reference& reference, const RunPlan& plan,
                            long long step, ErrorCsv& csv, E01ProfileCsv& profiles) {
            if (step == 0) {
                observePotential(state, reference, "initial", 0, step, csv, profiles);
                observeFields(state, reference, "initial", 0, step, csv, profiles);
                return;
            }
            const auto observed = observeScheduled(state, reference, plan, step, csv, profiles);
            if (step != plan.steps_m || plan.steps_m == plan.fullSteps_m)
                return;
            const double target = state.labelTime(Quantity::VectorPotential, step);
            if (!observed.potential_m)
                observePotential(state, reference, "partial_final", target, step, csv, profiles);
            if (!observed.fields_m)
                observeFields(state, reference, "partial_final", target, step, csv, profiles);
        }

        /// Record the analytic input values without inferring them from numerical field samples.
        inline void recordReference(output::RunRecordDescription& description,
                                    const E01Reference& reference) {
            auto node                   = description.effective_m["analytic_reference"];
            node["name"]                = "E01 periodic compact continuum pulse";
            node["period"]              = reference.Lz_m;
            node["initial_centre_z"]    = reference.z0_m;
            node["wave_speed"]          = reference.v_m;
            node["potential_amplitude"] = E01_analytic_continuum_reference::A0;
            node["support_half_width"]  = E01_analytic_continuum_reference::w;
        }

        /// Keep the prescribed schedule visible even when a partial run stops before its targets.
        inline void recordSchedule(output::RunRecordDescription& description, const RunPlan& plan) {
            std::vector<double> targets;
            std::vector<int> potentialSteps, fieldSteps;
            for (const auto& checkpoint : plan.checkpoints_m) {
                targets.push_back(checkpoint.targetTime_m);
                potentialSteps.push_back(checkpoint.potentialStep_m);
                fieldSteps.push_back(checkpoint.fieldStep_m);
            }
            auto node = description.effective_m["error_observations"];
            node["protocol_target_times"].set(targets);
            node["protocol_potential_steps"].set(potentialSteps);
            node["protocol_field_steps"].set(fieldSteps);
            node["initial_step"]        = 0;
            node["partial_final"]       = plan.steps_m > 0 && plan.steps_m < plan.fullSteps_m;
            node["file"]                = "errors.csv";
            node["norms"]               = "A_L2,E_L2,E_Linf,B_L2; full-domain relative errors";
            node["magnetic_convention"] = "R2 averaged magnetic field";
        }

        /// Identify the fixed native line and its connection to the error observations.
        inline void recordProfiles(output::RunRecordDescription& description) {
            auto profiles                   = description.effective_m["profiles"];
            profiles["file"]                = "profiles.csv";
            profiles["global_cell_i"]       = e01ProfileDetail::LineCellX;
            profiles["global_cell_j"]       = e01ProfileDetail::LineCellY;
            profiles["components"]          = "Ax,Ex,By";
            profiles["magnetic_convention"] = "R2 averaged magnetic field";
            profiles["schedule"] = "same observations as the corresponding global errors";
        }

        /// Publish the common provenance record only after the selected run and CSVs finish.
        inline void recordRun(E01Simulation& run, const input::ChdrYamlRunInput& request,
                              const simulation::ChdrCartesianGridSettings& grid,
                              const E01Reference& reference, const RunPlan& plan,
                              const E01FieldSnapshots& snapshots) {
            if (request.outputDirectory_m.empty())
                return;
            auto description = output::describeRun(
                "experiments-wave-pulses-in-homogeneous-dielectric", grid, run.solver(),
                plan.steps_m,
                "E01 longitudinal refinement; errors, profiles and optional field snapshots");
            description.effective_m["material/relative_permittivity"] =
                run.material().homogeneousValue();
            recordReference(description, reference);
            recordSchedule(description, plan);
            recordProfiles(description);
            snapshots.record(description);
            output::writeRunRecord(request.outputDirectory_m, request, description,
                                   run.solver().A_n.getLayout().comm);
        }
    }  // namespace e01_run_detail

    /** @brief Advance E01 and save errors, longitudinal profiles and scheduled field frames.
     * @details The manager observes state zero and every completed solve. Expensive reductions
     * occur only at the independently selected potential/field checkpoints and, for a short
     * run, its additional endpoint. The output and schedule implement E01.typ, meas01-04.
     * @pre validateE01Input has succeeded and the analytic histories are already loaded.
     */
    inline void evolveE01(E01Simulation& run, const input::ChdrYamlRunInput& request,
                          const simulation::ChdrCartesianGridSettings& grid,
                          const E01Reference& reference) {
        const auto plan = e01_run_detail::runPlan(request, reference, run.solver().getDt());
        e01_run_detail::ErrorCsv csv(request.outputDirectory_m, grid, run.solver().getDt(),
                                     run.solver().A_n.getLayout().comm);
        E01ProfileCsv profiles(request.outputDirectory_m, reference, run.solver().levelContract(),
                               run.solver().A_n.getLayout().comm);
        E01FieldSnapshots snapshots(request, grid, run.solver(),
                                    plan.checkpoints_m.front().potentialStep_m);
        run.advance(plan.steps_m, [&](long long step, E01Solver& state, const auto&) {
            e01_run_detail::observe(state, reference, plan, step, csv, profiles);
            snapshots.observe(step, state, plan.steps_m);
        });
        csv.finish();
        profiles.finish();
        snapshots.finish();
        e01_run_detail::recordRun(run, request, grid, reference, plan, snapshots);
    }
}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric
