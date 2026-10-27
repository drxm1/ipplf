/** @file InterfaceRun.h
 * @brief Compatibility names for interface experiments using the shared field-simulation manager.
 */
#pragma once

#include "core/material/PlanarSlab.h"
#include "core/simulation/ChDRSimulationManager.h"

namespace chdr::apps::planarInterface {
    using RunBox          = simulation::RunBox;
    using SolverSettings  = simulation::SolverSettings;
    using InterfaceRun    = simulation::ChDRSimulationManager<>;
    using InterfaceSolver = InterfaceRun::Solver;
}  // namespace chdr::apps::planarInterface
