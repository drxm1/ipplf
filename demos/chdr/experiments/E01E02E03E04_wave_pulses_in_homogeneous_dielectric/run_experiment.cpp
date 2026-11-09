#include <optional>

#include "core/input/Application.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02AnalysisProcess.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E02.h"

namespace chdr::experiments::e01e02 {
    /// Periodic propagation uses the shared pulse implementation with E01's explicit protocol.
    void runE01(const input::ChdrYamlRunInput& request) {
        runPulse(request, E01);
    }

    /// Mur exit uses the same implementation with E02's protocol and additional diagnostics.
    void runE02(const input::ChdrYamlRunInput& request) {
        runPulse(request, E02);
    }
}  // namespace chdr::experiments::e01e02

namespace {
    namespace experiment = chdr::experiments::e01e02;

    // Only host data crosses the IPPL lifetime boundary; Python starts after MPI finalization.
    std::optional<experiment::PulseAnalysisCommand> pendingAnalysis;

    /// Keep experiment and study selection visible at the application's entry point.
    void runConfigured(const chdr::input::ChdrYamlRunInput& input) {
        if (input.caseName_m == "E01-study")
            pendingAnalysis = experiment::runPulseStudy(input, experiment::E01, experiment::runE01);
        else if (input.caseName_m == "E02-study")
            pendingAnalysis = experiment::runPulseStudy(input, experiment::E02, experiment::runE02);
        else if (input.caseName_m == "E01")
            experiment::runE01(input);
        else if (input.caseName_m == "E02")
            experiment::runE02(input);
        else
            throw std::invalid_argument("Supported cases are E01, E02, E01-study and E02-study");
    }
}  // namespace

int main(int argc, char** argv) {
    const int status = chdr::input::runApplication(argc, argv, runConfigured);
    if (status != 0 || !pendingAnalysis)
        return status;
    return experiment::runPulseAnalysis(*pendingAnalysis);
}
