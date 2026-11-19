#include "core/input/Application.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Study.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E02.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Run.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04Run.h"

namespace chdr::experiments::e01e02 {
    /// Periodic propagation uses the shared pulse implementation with E01's explicit protocol.
    void runE01(const input::ChdrYamlRunInput& request) {
        runPulse<E01Reference>(request, E01);
    }

    /// Mur exit uses the same implementation with E02's protocol and additional diagnostics.
    void runE02(const input::ChdrYamlRunInput& request) {
        runPulse<E01Reference>(request, E02);
    }
}  // namespace chdr::experiments::e01e02

namespace {
    namespace experiment = chdr::experiments::e01e02;
    namespace spherical  = chdr::experiments::e03;

    /// E03 and E04 use the compact dipole pulse with the sharpness selected in the input.
    spherical::E03Reference compactDipolePulse(const chdr::input::ChdrYamlRunInput& input) {
        return spherical::sphericalReference(spherical::sphericalSharpness(input));
    }

    /// Periodic spherical pulse of E03.typ.
    void runE03(const chdr::input::ChdrYamlRunInput& input) {
        chdr::input::validateParameters(input, {"cells_per_axis", "write_plane",
                                                "profile_sharpness", "single_precision_fields"});
        spherical::runE03(input, compactDipolePulse(input),
                          {"E03", experiment::WavePulsesApplication, spherical::E03Epsilon,
                           chdr::experiments::O1VSConfiguration});
    }

    /// The same pulse leaving the all-Mur cube of E04.typ.
    void runE04(const chdr::input::ChdrYamlRunInput& input) {
        chdr::input::validateParameters(
            input, {"cells_per_axis", "profile_sharpness", "single_precision_fields"});
        chdr::experiments::e04::runE04(
            input, compactDipolePulse(input),
            {"E04", experiment::WavePulsesApplication, spherical::E03Epsilon,
             chdr::experiments::O1VSConfiguration});
    }

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
            runE03(input);
        else if (input.caseName_m == "E04")
            runE04(input);
        else
            throw std::invalid_argument(
                "Supported cases are E01 through E04, E01-study and E02-study");
    }
}  // namespace

int main(int argc, char** argv) {
    return chdr::input::runApplication(argc, argv, runConfigured);
}
