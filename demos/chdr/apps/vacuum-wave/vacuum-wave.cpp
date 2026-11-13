#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include "VacuumChecks.h"
#include "VacuumDiagnostics.h"
#include "VacuumEnergy.h"
#include "VacuumSetup.h"
#include "core/input/Application.h"
#include "core/input/RunInput.h"
#include "core/output/FieldOutput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/ObservationSchedule.h"

namespace {
    void runConfigured(const chdr::input::ChdrYamlRunInput& input);
}  // namespace

/** @brief Run the demo and finalize IPPL after simulation objects are destroyed.
 * @param argc Number of command-line arguments, including the program name.
 * @param argv Command-line strings passed to IPPL initialization.
 * @return Status from the shared application lifecycle and this case's checks.
 */
int main(int argc, char* argv[]) {
    return chdr::input::runApplication(argc, argv, runConfigured);
}

namespace {
    using namespace chdr::apps::vacuumWave;

    /// Select an existing benchmark grid; its fixed acceptance checks require these cases.
    int selectedResolution(const chdr::input::ChdrYamlRunInput& input) {
        chdr::input::validateParameters(input, {});
        chdr::input::validateObservationNames(input, {"diagnostics", "fields"});
        if (input.outputDirectory_m.empty() && input.observations_m.contains("fields")
            && chdr::input::observationEnabled(input, "fields", true))
            throw std::invalid_argument("Enabled fields observation requires an output directory");
        if (input.caseName_m == "coarse")
            return 64;
        if (input.caseName_m == "fine")
            return 128;
        throw std::invalid_argument("Vacuum-wave case must be coarse or fine");
    }

    /// @brief Copy the actual field geometry into the writer's unit-neutral metadata.
    chdr::output::OutputGeometry outputGeometry(Fields& fields) {
        chdr::output::OutputGeometry geometry;
        for (unsigned d = 0; d < Dim; ++d) {
            geometry.cells_m[d]   = fields.getFL().getDomain()[d].length();
            geometry.origin_m[d]  = fields.getMesh().getOrigin()[d];
            geometry.spacing_m[d] = fields.getMesh().getMeshSpacing()[d];
        }
        return geometry;
    }

    /// @brief Create optional collective output after all simulation storage exists.
    std::unique_ptr<chdr::output::FieldOutput> makeOutput(
        Fields& fields, const chdr::input::ChdrYamlRunInput& input) {
        if (input.outputDirectory_m.empty())
            return nullptr;
        if (!chdr::input::observationEnabled(input, "fields", true)) {
            chdr::output::prepareRunDirectory(input.outputDirectory_m,
                                              fields.getFL().comm.getCommunicator());
            return nullptr;
        }
        return std::make_unique<chdr::output::FieldOutput>(
            input.outputDirectory_m, outputGeometry(fields), fields.getFL().comm.getCommunicator(),
            "wave", false, chdr::output::fieldOutputScales(input.units_m));
    }

    /// @brief Bound a requested partial run by the one-period benchmark duration.
    int numberOfSteps(const chdr::input::ChdrYamlRunInput& input, double dt) {
        const int periodSteps = static_cast<int>(std::llround(1.0 / dt));
        const int steps       = chdr::input::resolvedSteps(input, dt, periodSteps);
        if (steps > periodSteps)
            throw std::invalid_argument("Requested run exceeds the one-period benchmark");
        return steps;
    }

    /// @brief Select physical-error checkpoints at quarter periods and the final state.
    bool isCheckpoint(int step, int lastStep, double dt) {
        const int quarterSteps = static_cast<int>(std::llround(0.25 / dt));
        return step % quarterSteps == 0 || step == lastStep;
    }

    /// @brief Append the energy columns and end the row; W_E and W_B have different times.
    void printEnergy(const chdr::apps::vacuumWave::EnergySample& energy,
                     const chdr::input::UnitSystem& units) {
        const double factor = units.energyFactor();
        std::cout << ',' << energy.electric_m * factor << ',' << energy.magnetic_m * factor << ','
                  << energy.total_m * factor << ',' << energy.invariant_m << std::endl;
    }

    /// Convert timestamps only; the recorded errors and normalized mode coefficient stay unitless.
    void printStateAndErrors(const Diagnostics& d, const chdr::input::UnitSystem& units) {
        const double time = units.timeFactor();
        std::cout << std::setprecision(17) << d.step_m << ',' << d.dt_m * time << ','
                  << d.timeA_m * time << ',' << d.timeE_m * time << ',' << d.timeB_m * time << ','
                  << d.relativeE_m << ',' << d.relativeB_m << ',' << d.discreteError_m << ','
                  << d.zeroError_m << ',' << d.sourceError_m << ',' << d.phaseReal_m << ','
                  << d.phaseImag_m << ',';
    }

    /// Label SI publication without changing the normalized CSV header or values.
    void printCsvHeader(const chdr::input::UnitSystem& units) {
        if (units.isSI())
            std::cout << std::setprecision(17)
                      << "# time_units=s phase_speed_units=m/s "
                         "physical_energy_columns=energyE,energyB,energy physical_energy_units=J "
                         "energy_invariant_units=normalized length_m="
                      << units.length_m << " charge_density_C_m3=" << units.chargeDensityScale_m
                      << '\n';
        std::cout << "step,dt,timeA,timeE,timeB,relativeE,relativeB,discreteError,zeroError,"
                     "sourceError,phaseReal,phaseImag,phaseAmplitude,phaseIncrement,phaseSpeed,"
                     "count,energyE,energyB,energy,energyInvariant\n";
    }

    /// @brief Emit owned-cell diagnostics; field times and undefined initial speed are explicit.
    void printDiagnostics(const chdr::apps::vacuumWave::Diagnostics& d,
                          const chdr::apps::vacuumWave::PhaseSample& phase,
                          const chdr::apps::vacuumWave::EnergySample& energy,
                          const chdr::input::UnitSystem& units) {
        if (ippl::Comm->rank() != 0)
            return;
        if (d.step_m == 0)
            printCsvHeader(units);
        printStateAndErrors(d, units);
        std::cout << phase.amplitude_m << ',' << phase.increment_m << ',';
        if (phase.hasSpeed_m)
            std::cout << phase.speed_m * units.velocityFactor();
        std::cout << ',' << d.count_m;
        printEnergy(energy, units);
        if (!std::cout)
            throw std::runtime_error("Failed to write diagnostics");
    }

    /// Publish contract-labelled fields; this benchmark's scalar potential is identically zero.
    void writeCheckpoint(chdr::output::FieldOutput& output, VacuumSimulation& run,
                         const Diagnostics& d, const chdr::input::UnitSystem& units) {
        const chdr::output::FrameRecord frame{d.step_m, d.timeA_m, d.dt_m};
        const std::array snapshots{
            chdr::output::VectorFieldSnapshot{"E_raw", run.fields().getE(),
                                              chdr::Quantity::Electric, units.electricFactor()},
            chdr::output::VectorFieldSnapshot{"B_raw", run.fields().getB(),
                                              chdr::Quantity::RawMagnetic, units.magneticFactor()}};
        output.writeFrame(frame, run.solver().levelContract(), snapshots, /*phiIsZero=*/true);
        output.writeSeries();
    }

    /** @brief Require every state for phase tracking and the existing benchmark checks.
     * @details Output cadence is independent. Skipping diagnostic states would change the
     * phase-increment checks and cease to run the existing benchmark protocol.
     */
    void validateDiagnostics(const chdr::input::ChdrYamlRunInput& input, int lastStep) {
        const auto schedule = chdr::input::observationSchedule(input, "diagnostics", {});
        if (!chdr::input::observationEnabled(input, "diagnostics", true) || schedule.every_m != 1
            || schedule.first_m != 0 || (schedule.last_m >= 0 && schedule.last_m < lastStep))
            throw std::invalid_argument("Vacuum-wave diagnostics must observe every run state");
    }

    /// Measure, print and check the existing diagnostics without changing their arithmetic.
    Diagnostics observeState(int step, VacuumSimulation& run, PhaseTracker& tracker, int lastStep,
                             const chdr::input::ChdrYamlRunInput& input) {
        const auto diagnostic = measure(run.solver(), run.fields(), step);
        const auto energy     = measureEnergy(run.solver(), run.fields());
        const auto phase      = tracker.update(diagnostic);
        printDiagnostics(diagnostic, phase, energy, input.units_m);
        validateStep(diagnostic, phase, isCheckpoint(step, lastStep, run.solver().getDt()),
                     input.caseName_m == "fine");
        return diagnostic;
    }

    /// Resolve the default quarter-period cadence and any requested field-output overrides.
    chdr::simulation::ObservationSchedule fieldSchedule(const chdr::input::ChdrYamlRunInput& input,
                                                        double dt) {
        const chdr::simulation::ObservationSchedule fallback{
            .every_m = static_cast<int>(std::llround(0.25 / dt)), .includeFinal_m = true};
        return chdr::input::observationSchedule(input, "fields", fallback);
    }

    /** @brief Observe every manager state and write the selected field frames.
     * @details Diagnostics run before field output and retain the existing benchmark checks.
     * Each completed frame publishes its series; no separate finalizer is required.
     */
    void evolve(VacuumSimulation& run, const chdr::input::ChdrYamlRunInput& input, int lastStep,
                chdr::simulation::ObservationSchedule schedule) {
        auto output = makeOutput(run.fields(), input);
        PhaseTracker tracker;
        run.advance(lastStep, [&](long long step, auto&, const auto&) {
            const auto diagnostic =
                observeState(static_cast<int>(step), run, tracker, lastStep, input);
            if (output && schedule.selects(step, lastStep))
                writeCheckpoint(*output, run, diagnostic, input.units_m);
        });
    }

    /// Record the CSV's units, including the deliberately unscaled discrete invariant.
    void recordCsvUnits(chdr::output::RunRecordDescription& description,
                        const chdr::input::UnitSystem& units) {
        auto csv                       = description.effective_m["output/stdout_csv"];
        csv["time_units"]              = units.isSI() ? "s" : "normalized";
        csv["phase_speed_units"]       = units.isSI() ? "m/s" : "normalized";
        csv["physical_energy_columns"] = "energyE,energyB,energy";
        csv["physical_energy_units"]   = units.isSI() ? "J" : "normalized";
        csv["energy_invariant_units"]  = "normalized";
        if (units.isSI()) {
            csv["length_m"]            = units.length_m;
            csv["charge_density_C_m3"] = units.chargeDensityScale_m;
        }
    }

    /// Record the actual schedules and distinguish normalized mode coefficients from fields.
    void recordRun(VacuumSimulation& run, const chdr::simulation::ChdrCartesianGridSettings& box,
                   const chdr::input::ChdrYamlRunInput& input, int steps,
                   chdr::simulation::ObservationSchedule schedule) {
        auto description = chdr::output::describeRun(
            "chdr-vacuum-wave", box, run.solver(), steps,
            "CSV errors and modal phase coefficients are dimensionless normalized quantities; "
            "raw E/B energies use distinct times; energyInvariant remains normalized.");
        recordCsvUnits(description, input.units_m);
        chdr::output::recordObservation(description, "diagnostics",
                                        chdr::input::observationSchedule(input, "diagnostics", {}));
        chdr::output::recordObservation(
            description, "fields", schedule,
            !input.outputDirectory_m.empty()
                && chdr::input::observationEnabled(input, "fields", true));
        chdr::output::writeRunRecord(input.outputDirectory_m, input, description,
                                     run.fields().getFL().comm.getCommunicator());
    }

    /// @brief Own one complete benchmark run; field lifetime encloses the solver's lifetime.
    void runConfigured(const chdr::input::ChdrYamlRunInput& input) {
        const auto box = vacuumBox(selectedResolution(input));
        VacuumSimulation run(box, VacuumMaterial{}, vacuumSettings(box));
        run.load(VacuumPotential{});
        const int lastStep  = numberOfSteps(input, run.solver().getDt());
        const auto schedule = fieldSchedule(input, run.solver().getDt());
        validateDiagnostics(input, lastStep);
        evolve(run, input, lastStep, schedule);
        recordRun(run, box, input, lastStep, schedule);
    }
}  // namespace
