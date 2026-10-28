/** @file InterfaceRun.h
 * @brief Compatibility names for interface experiments using the shared field-simulation manager.
 */
#pragma once

#include "core/simulation/ChDRSimulationManager.h"

namespace chdr::apps::planarInterface {
    // TODO: This file doesn't really have a point anymore,
    // we could just use the actual names without these
    // using declarations.
    using RunBox          = simulation::RunBox;
    using SolverSettings  = simulation::SolverSettings;
    using InterfaceRun    = simulation::ChDRSimulationManager<>;
    using InterfaceSolver = InterfaceRun::Solver;
}  // namespace chdr::apps::planarInterface
