/** @file SimulationConfig.h
 * @brief Experiment-supplied mesh, partition, boundary and solver settings.
 */
#pragma once

#include <array>
#include <optional>

#include "ChdrSolver.h"

namespace chdr::simulation {
    /** @brief Define the Cartesian box and its MPI partition independently of the material map.
     * @details Arrays use x, y, z order. A zero seam selects automatic decomposition;
     * a positive seam is a global cell boundary along splitAxis_m and requires two ranks.
     * Periodicity is passed to the shared layout used by fields and material.
     */
    struct RunBox {
        std::array<int, 3> cells_m;
        std::array<double, 3> spacing_m;
        int seam_m = 0;
        std::array<bool, 3> periodic_m{true, true, true};
        unsigned splitAxis_m = 0;

        /// Return the physical extent from the global cell counts and spacings.
        std::array<double, 3> lengths() const {
            return {cells_m[0] * spacing_m[0], cells_m[1] * spacing_m[1],
                    cells_m[2] * spacing_m[2]};
        }
    };

    /** @brief Choose the solver's runtime settings for one experiment.
     * @details ChdrSolver validates supported combinations. The potential stencil is
     * selected separately through ChDRSimulationManager's template argument. An empty
     * formulation selects the owner-free vacuum update and requires unit permittivity.
     * safetyFactor_m is passed to ChdrSolver::setSafetyFactor before histories are loaded.
     */
    struct SolverSettings {
        double safetyFactor_m = 0.9;
        solver::EBReconstructionConvention reconstruction_m =
            solver::EBReconstructionConvention::R2;
        solver::SpatialDiscretization spatial_m = solver::SpatialDiscretization::Staggered;
        std::optional<solver::FormulationParameters> formulation_m =
            solver::FormulationParameters{solver::MaterialFormulation::O1VS, 1};
        solver::MurSpeed murSpeed_m = solver::MurSpeed::matched();
    };
}  // namespace chdr::simulation
