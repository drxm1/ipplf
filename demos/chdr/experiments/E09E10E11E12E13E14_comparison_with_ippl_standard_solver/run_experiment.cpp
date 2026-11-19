#include <cmath>

#include "core/input/Application.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Study.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Run.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04Run.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/SolverConfiguration.h"
#include "experiments/E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/GaussianProfiles.h"
#include "experiments/E09E10E11E12E13E14_comparison_with_ippl_standard_solver/E09.h"

namespace chdr::experiments::e09e14 {
    namespace gaussian = E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric;

    /// E05 in vacuum with the solver configuration of the input.
    void runE09(const input::ChdrYamlRunInput& request) {
        e01e02::runPulse<gaussian::GaussianPlanePulse>(request, E09);
    }

    /// E06 in vacuum with the solver configuration of the input.
    void runE10(const input::ChdrYamlRunInput& request) {
        e01e02::runPulse<gaussian::GaussianPlanePulse>(request, E10);
    }

    /** @brief The Gaussian dipole pulse of E07 in a medium of relative permittivity epsilon.
     * @details E11 and E12 use the vacuum, E13 and E14 the dielectric of E03 to E08.
     */
    e03::SphericalExperiment dipoleExperiment(const input::ChdrYamlRunInput& request,
                                              const char* name, double epsilon) {
        return {name, SolverComparisonApplication, epsilon, solverConfiguration(request)};
    }

    /// The Gaussian dipole pulse of E07 at the speed \f$1/\sqrt{\varepsilon_r}\f$ of the medium.
    gaussian::GaussianDipolePulse dipolePulse(double epsilon) {
        return {1 / std::sqrt(epsilon), e03::E03BoxCentre};
    }

    /// E07 (periodic cube) in the medium of epsilon; E11 in vacuum, E13 in the dielectric.
    void runPeriodicDipole(const input::ChdrYamlRunInput& request, const char* name,
                           double epsilon) {
        input::validateParameters(request, {"cells_per_axis", "write_plane", "write_native_lines",
                                            "single_precision_fields", "solver_configuration"});
        e03::runE03(request, dipolePulse(epsilon), dipoleExperiment(request, name, epsilon));
    }

    /// E08 (all-face Mur cube) in the medium of epsilon; E12 in vacuum, E14 in the dielectric.
    void runMurDipole(const input::ChdrYamlRunInput& request, const char* name, double epsilon) {
        input::validateParameters(request, {"cells_per_axis", "single_precision_fields",
                                            "solver_configuration", "write_energy_parts"});
        e04::runE04(request, dipolePulse(epsilon), dipoleExperiment(request, name, epsilon));
    }
}  // namespace chdr::experiments::e09e14

namespace {
    namespace experiment = chdr::experiments::e09e14;
    namespace pulse      = chdr::experiments::e01e02;
    using chdr::experiments::e03::E03Epsilon;

    /// Keep experiment and study selection visible at the application's entry point.
    void runConfigured(const chdr::input::ChdrYamlRunInput& input) {
        const auto& name = input.caseName_m;
        if (name == "E09-study")
            pulse::runPulseStudy(input, experiment::E09, experiment::runE09);
        else if (name == "E10-study")
            pulse::runPulseStudy(input, experiment::E10, experiment::runE10);
        else if (name == "E09")
            experiment::runE09(input);
        else if (name == "E10")
            experiment::runE10(input);
        else if (name == "E11" || name == "E13")
            experiment::runPeriodicDipole(input, name.c_str(), name == "E11" ? 1.0 : E03Epsilon);
        else if (name == "E12" || name == "E14")
            experiment::runMurDipole(input, name.c_str(), name == "E12" ? 1.0 : E03Epsilon);
        else
            throw std::invalid_argument(
                "Supported cases are E09 through E14, E09-study and E10-study");
    }
}  // namespace

int main(int argc, char** argv) {
    return chdr::input::runApplication(argc, argv, runConfigured);
}
