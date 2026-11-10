#include "core/input/Application.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Study.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E02.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04.h"

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

    /// Keep experiment and study selection visible at the application's entry point.
    void runConfigured(const chdr::input::ChdrYamlRunInput& input) {
        if (input.caseName_m == "E01-study")
            experiment::runPulseStudy(input, experiment::E01, experiment::runE01);
        else if (input.caseName_m == "E02-study")
            experiment::runPulseStudy(input, experiment::E02, experiment::runE02);
        else if (input.caseName_m == "E01")
            experiment::runE01(input);
        else if (input.caseName_m == "E02")
            experiment::runE02(input);
        else if (input.caseName_m == "E03")
            chdr::experiments::e03::runE03(input);
        else if (input.caseName_m == "E04")
            chdr::experiments::e04::runE04(input);
        else
            throw std::invalid_argument(
                "Supported cases are E01 through E04, E01-study and E02-study");
    }
}  // namespace

int main(int argc, char** argv) {
    return chdr::input::runApplication(argc, argv, runConfigured);
}
