#include <cmath>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <locale>
#include <memory>
#include <stdexcept>
#include <string>

#include "VacuumChecks.h"
#include "VacuumDiagnostics.h"
#include "VacuumFailure.h"
#include "VacuumOutput.h"
#include "VacuumSetup.h"

namespace {
    int runVacuumWave(int argc, char* argv[]);
}  // namespace

/** @brief Run the demo and finalize IPPL after simulation objects are destroyed.
 * @param argc Number of command-line arguments, including the program name.
 * @param argv Command-line strings passed to IPPL initialization.
 * @return Status returned by runVacuumWave().
 */
int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    const int result = runVacuumWave(argc, argv);
    ippl::finalize();
    return result;
}

namespace {
    using namespace chdr::vacuum;

    /// Bounded choices for the agreed coarse/fine benchmark; no adjustable timestep.
    struct Options {
        int cellsZ_m = 64;
        int steps_m  = -1;  ///< -1 selects one full period; 0/1 support initialization checks.
        std::filesystem::path output_m;  ///< Empty disables VTK output; CSV goes to stdout.
    };

    /// @brief Read a required next CLI value or throw before field construction.
    std::string argumentValue(int& index, int argc, char* argv[]) {
        if (index + 1 >= argc)
            throw std::invalid_argument("Missing value after " + std::string(argv[index]));
        return argv[++index];
    }

    /// @brief Parse a nonnegative step count with no ignored suffix.
    int parseSteps(const std::string& text) {
        std::size_t consumed = 0;
        const int steps      = std::stoi(text, &consumed);
        if (consumed != text.size() || steps < 0)
            throw std::invalid_argument("--steps needs a nonnegative integer");
        return steps;
    }

    /// @brief Consume one demo argument; skip IPPL options already handled by initialize.
    void parseArgument(Options& options, int& index, int argc, char* argv[]) {
        const std::string arg(argv[index]);
        if (arg == "--fine")
            options.cellsZ_m = 128;
        else if (arg == "--steps")
            options.steps_m = parseSteps(argumentValue(index, argc, argv));
        else if (arg == "--output")
            options.output_m = argumentValue(index, argc, argv);
        else if (arg == "--info" || arg == "-i" || arg == "--overallocate" || arg == "-b"
                 || arg == "--timer-fences")
            argumentValue(index, argc, argv);
        else if (arg == "--debug" || arg == "-g")
            return;  // IPPL has already handled the debugger pause.
        else
            throw std::invalid_argument("Unknown argument: " + arg);
    }

    /// @brief Parse identical per-rank options; output must be new/empty when supplied.
    Options parseOptions(int argc, char* argv[]) {
        Options options;
        for (int index = 1; index < argc; ++index)
            parseArgument(options, index, argc, argv);
        return options;
    }

    /// @brief Copy the actual field geometry into the writer's unit-neutral metadata.
    chdr::vacuum::OutputGeometry outputGeometry(Fields& fields) {
        chdr::vacuum::OutputGeometry geometry;
        for (unsigned d = 0; d < Dim; ++d) {
            geometry.cells_m[d]   = fields.getFL().getDomain()[d].length();
            geometry.origin_m[d]  = fields.getMesh().getOrigin()[d];
            geometry.spacing_m[d] = fields.getMesh().getMeshSpacing()[d];
        }
        return geometry;
    }

    /// @brief Create optional collective output after all simulation storage exists.
    std::unique_ptr<chdr::vacuum::VacuumOutput> makeOutput(Fields& fields, const Options& options) {
        if (options.output_m.empty())
            return nullptr;
        return std::make_unique<chdr::vacuum::VacuumOutput>(
            options.output_m, outputGeometry(fields), fields.getFL().comm.getCommunicator());
    }

    /// @brief Bound a requested partial run by the one-period benchmark duration.
    int numberOfSteps(const Options& options, double dt) {
        const int periodSteps = static_cast<int>(std::llround(1.0 / dt));
        if (options.steps_m > periodSteps)
            throw std::invalid_argument("--steps exceeds the one-period benchmark");
        return options.steps_m < 0 ? periodSteps : options.steps_m;
    }

    /// @brief Select physical-error checkpoints at quarter periods and the final state.
    bool isCheckpoint(int step, int lastStep, double dt) {
        const int quarterSteps = static_cast<int>(std::llround(0.25 / dt));
        return step % quarterSteps == 0 || step == lastStep;
    }

    /// @brief Emit owned-cell diagnostics; field times and undefined initial speed are explicit.
    void printDiagnostics(const chdr::vacuum::Diagnostics& d,
                          const chdr::vacuum::PhaseSample& phase) {
        if (ippl::Comm->rank() != 0)
            return;
        if (d.step_m == 0)
            std::cout << "step,dt,timeA,timeE,timeB,relativeE,relativeB,discreteError,zeroError,"
                         "sourceError,phaseReal,phaseImag,phaseAmplitude,phaseIncrement,phaseSpeed,"
                         "count\n";
        std::cout << std::setprecision(17) << d.step_m << ',' << d.dt_m << ',' << d.timeA_m << ','
                  << d.timeE_m << ',' << d.timeB_m << ',' << d.relativeE_m << ',' << d.relativeB_m
                  << ',' << d.discreteError_m << ',' << d.zeroError_m << ',' << d.sourceError_m
                  << ',' << d.phaseReal_m << ',' << d.phaseImag_m << ',' << phase.amplitude_m << ','
                  << phase.increment_m << ',';
        if (phase.hasSpeed_m)
            std::cout << phase.speed_m;
        std::cout << ',' << d.count_m << std::endl;
        if (!std::cout)
            throw std::runtime_error("Failed to write diagnostics");
    }

    /// @brief Publish measured array times and update the series after a complete frame.
    void writeCheckpoint(chdr::vacuum::VacuumOutput& output, Fields& fields,
                         const chdr::vacuum::Diagnostics& d) {
        const chdr::vacuum::FrameRecord frame{d.step_m, d.timeA_m, d.timeE_m, d.timeB_m, d.dt_m};
        output.writeFrame(fields.getE(), fields.getB(), frame);
        output.writeSeries();
    }

    /// @brief Measure/check every state; optionally publish owned fields at checkpoints.
    void evolve(VacuumSolver& solver, Fields& fields, const Options& options) {
        const int lastStep = numberOfSteps(options, solver.getDt());
        auto output        = makeOutput(fields, options);
        chdr::vacuum::PhaseTracker phaseTracker;
        // Observe states 0 through lastStep; perform exactly lastStep advances.
        for (int step = 0; step <= lastStep; ++step) {
            const auto diagnostic = chdr::vacuum::measure(solver, fields, step);
            const auto phase      = phaseTracker.update(diagnostic);
            const bool checkpoint = isCheckpoint(step, lastStep, solver.getDt());
            printDiagnostics(diagnostic, phase);
            validateStep(diagnostic, phase, checkpoint, options.cellsZ_m == 128);
            if (output && checkpoint)
                writeCheckpoint(*output, fields, diagnostic);
            if (step < lastStep)
                solver.solve();
        }
    }

    /// @brief Own one complete benchmark run; field lifetime encloses the solver's lifetime.
    void runConfigured(const Options& options) {
        auto fields = makePeriodicFieldContainer(options.cellsZ_m);
        checkLocalDomain(fields.getFL());
        initializeFieldStorage(fields);
        VacuumSolver solver(fields.getJ(), fields.getE(), fields.getB());
        initializePotentials(solver);
        solver.evaluate_EB();
        evolve(solver, fields, options);
    }

    /// @brief Return success or abort MPI on a rank-local exception.
    int runVacuumWave(int argc, char* argv[]) {
        try {
            std::cout.imbue(std::locale::classic());
            runConfigured(parseOptions(argc, argv));
            return 0;
        } catch (const IpplException& error) {
            return abortFailure(error.where(), error.what());
        } catch (const std::exception& error) {
            return abortFailure("runVacuumWave", error.what());
        }
    }
}  // namespace
