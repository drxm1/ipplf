#include <cmath>

#include "core/input/Application.h"
#include "core/simulation/ChDRSimulationManager.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01AnalyticContinuumReference.h"

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric {

    void runE01([[maybe_unused]] const chdr::input::ChdrYamlRunInput& input) {
        // Grid settings for the coarsest E01 mesh.
        const simulation::ChdrCartesianGridSettings gridSettings{
            {16, 16, 704}, {1.0 / 16, 1.0 / 16, 1.0 / 16}, 0, {true, true, true}, 2};

        simulation::ChdrSolverSettings solverSettings;
        solverSettings.timeStep_m = gridSettings.spacing_m[2] / 2;  //< h_t = h_z / 2

        // Construct a permittivity sampler for the homogeneous dielectric
        const auto homogeneousDielectric =
            KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
            marked = false;  //< Setting marked must be done for material statistics
            return 2.13;
        };

        // Construct the ChDRSimulationManager for the E01 simulation from the settings and
        // permittivity sampler
        simulation::ChDRSimulationManager<> simulationManager{gridSettings, homogeneousDielectric,
                                                              solverSettings};

        // Construct the analytic reference sampler
        const E01_analytic_continuum_reference::E01AnalyticContinuumReferenceSampler reference{
            gridSettings.lengths()[2], 2., 1. / std::sqrt(2.13)};

        // Load the analytic continuum wave / pulse into the staggered storage
        simulationManager.load(reference);
    }
}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric

int main(int argc, char** argv) {
    return chdr::input::runApplication(
        argc, argv, chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::runE01);
}
