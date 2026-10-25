/** @file InterfaceRun.h
 * @brief Compatibility names for interface experiments using the shared field-simulation manager.
 */
#pragma once

#include "../simulation/ChDRSimulationManager.h"
#include "PlanarSlab.h"
#include "opt-in/diagnostics/FieldEnergy.h"

namespace chdr::interface {
    using RunBox          = simulation::RunBox;
    using SolverSettings  = simulation::SolverSettings;
    using InterfaceRun    = simulation::ChDRSimulationManager<>;
    using InterfaceSolver = InterfaceRun::Solver;
}  // namespace chdr::interface
