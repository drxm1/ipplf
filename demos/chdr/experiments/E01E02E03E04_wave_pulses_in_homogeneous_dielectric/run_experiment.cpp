#include "core/input/Application.h"
#include "core/simulation/ChDRSimulationManager.h"

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric {

    void runE01([[maybe_unused]] const chdr::input::ChdrYamlRunInput& input) {
        // Grid settings for the coarsest E01 mesh.
        const simulation::ChdrCartesianGridSettings gridSettings{
            {16, 16, 704}, {1.0 / 16, 1.0 / 16, 1.0 / 16}, 0, {true, true, true}, 2};

        simulation::ChdrSolverSettings solverSettings;
        solverSettings.timeStep_m = gridSettings.spacing_m[2] / 2;

        // Construct a material map for the homogeneous dielectric
        const auto homogeneousDielectric =
            KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
            marked = false;
            return 2.13;
        };
        simulation::ChDRSimulationManager<> simulationManager{gridSettings, homogeneousDielectric,
                                                              solverSettings};
    }
}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric

int main(int argc, char** argv) {
    return chdr::input::runApplication(
        argc, argv, chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::runE01);
}
