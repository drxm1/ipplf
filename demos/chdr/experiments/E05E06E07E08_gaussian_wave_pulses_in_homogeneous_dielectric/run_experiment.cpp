#include "core/input/Application.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Study.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Run.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04Run.h"
#include "experiments/E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/E05.h"
#include "experiments/E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/E06.h"
#include "experiments/E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/GaussianProfiles.h"

namespace chdr::experiments::e05e08 {
    namespace gaussian = E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric;
    using gaussian::GaussianDipolePulse;
    using gaussian::GaussianPlanePulse;

    /// E01 with the Gaussian profile and the ladder 32 to 1024.
    void runE05(const input::ChdrYamlRunInput& request) {
        e01e02::runPulse<GaussianPlanePulse>(request, E05);
    }

    /// E02 with the Gaussian profile, without the comparison with the periodic run.
    void runE06(const input::ChdrYamlRunInput& request) {
        e01e02::runPulse<GaussianPlanePulse>(request, E06);
    }

    /// E03 with the Gaussian profile. E07 has no `profile_sharpness` parameter.
    void runE07(const input::ChdrYamlRunInput& request) {
        input::validateParameters(request,
                                  {"cells_per_axis", "write_plane", "single_precision_fields"});
        e03::runE03(request, GaussianDipolePulse{1 / std::sqrt(e03::E03Epsilon), e03::E03BoxCentre},
                    {"E07", GaussianPulsesApplication, e03::E03Epsilon, O1VSConfiguration});
    }

    /// E04 with the pulse of E07.
    void runE08(const input::ChdrYamlRunInput& request) {
        input::validateParameters(request, {"cells_per_axis", "single_precision_fields"});
        e04::runE04(request, GaussianDipolePulse{1 / std::sqrt(e03::E03Epsilon), e03::E03BoxCentre},
                    {"E08", GaussianPulsesApplication, e03::E03Epsilon, O1VSConfiguration});
    }
}  // namespace chdr::experiments::e05e08

namespace {
    namespace experiment = chdr::experiments::e05e08;
    namespace pulse      = chdr::experiments::e01e02;

    /// Keep experiment and study selection visible at the application's entry point.
    void runConfigured(const chdr::input::ChdrYamlRunInput& input) {
        if (input.caseName_m == "E05-study")
            pulse::runPulseStudy(input, experiment::E05, experiment::runE05);
        else if (input.caseName_m == "E06-study")
            pulse::runPulseStudy(input, experiment::E06, experiment::runE06);
        else if (input.caseName_m == "E05")
            experiment::runE05(input);
        else if (input.caseName_m == "E06")
            experiment::runE06(input);
        else if (input.caseName_m == "E07")
            experiment::runE07(input);
        else if (input.caseName_m == "E08")
            experiment::runE08(input);
        else
            throw std::invalid_argument(
                "Supported cases are E05 through E08, E05-study and E06-study");
    }
}  // namespace

int main(int argc, char** argv) {
    return chdr::input::runApplication(argc, argv, runConfigured);
}
