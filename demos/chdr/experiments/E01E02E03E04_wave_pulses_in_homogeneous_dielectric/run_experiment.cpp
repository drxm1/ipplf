#include <cmath>

#include "core/input/Application.h"
#include "core/simulation/ChDRSimulationManager.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01AnalyticContinuumReference.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01Run.h"

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric {

    void runE01(const chdr::input::ChdrYamlRunInput& input) {
        validateE01Input(input);
        const int cellsPerPulse = e01CellsPerPulse(input);
        // Keep the physical box and transverse mesh fixed; refine only z as specified in E01.typ.
        const simulation::ChdrCartesianGridSettings gridSettings{
            {16, 16, 22 * cellsPerPulse},
            {1.0 / 16, 1.0 / 16, 2.0 / cellsPerPulse},
            0,
            {true, true, true},
            2};

        // Solver settings selected explicitly as described in the experiment E01 specs
        const simulation::ChdrSolverSettings solverSettings{
            .reconstruction_m = EBReconstructionConvention::R2,
            .spatial_m        = SpatialDiscretization::Staggered,
            .formulation_m =
                solver::FormulationParameters{.formulation_m = solver::MaterialFormulation::O1VS,
                                              .alpha_m       = 1.},
            .timeStep_m = gridSettings.spacing_m[2] / 2.};

        // Construct a permittivity sampler for the homogeneous dielectric
        const auto homogeneousDielectric =
            KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
            marked = false;  //< Setting marked must be done for material statistics
            return 2.13;
        };

        // Construct the analytic reference sampler
        const E01_analytic_continuum_reference::E01AnalyticContinuumReferenceSampler reference{
            gridSettings.lengths()[2], 2., 1. / std::sqrt(2.13)};

        // Construct the ChDRSimulationManager for the E01 simulation from the settings and
        // permittivity sampler
        simulation::ChDRSimulationManager<solver::ChdrSolver<solver::StandardSolverPeriodic>>
            simulationManager{gridSettings, homogeneousDielectric, solverSettings};

        // Load the analytic continuum wave / pulse into the staggered storage
        simulationManager.load(reference);
        evolveE01(simulationManager, input, gridSettings, reference);
    }
}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric

int main(int argc, char** argv) {
    return chdr::input::runApplication(
        argc, argv, chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::runE01);
}
