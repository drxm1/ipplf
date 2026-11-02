/** @file VacuumSetup.h
 * @brief Fixed vacuum-wave configuration for the shared simulation manager.
 */
#pragma once

#include "Ippl.h"

#include "VacuumWave.h"
#include "core/base/FieldContainer.h"
#include "core/simulation/ChDRSimulationManager.h"

/// Vacuum-wave setup; callers initialize IPPL before constructing the manager.
namespace chdr::apps::vacuumWave {
    constexpr unsigned Dim = 3;  ///< Cartesian coordinates, ordered x, y, z.
    /// The shared field container of `core/base`, under the example's historical name.
    using Fields = ChDRFieldContainer;
    /// Ownership check and storage allocation come from `core/base`.
    using chdr::checkLocalDomain;
    using chdr::initializeFieldStorage;
    using VacuumSimulation = simulation::ChDRSimulationManager<>;
    using VacuumSolver     = VacuumSimulation::Solver;

    /** @brief Supply uniform unit permittivity without selecting a material update.
     * @details The manager owns a material map for every run; the vacuum settings select
     * its owner-free Standard update. No cell belongs to a separate material region.
     */
    struct VacuumMaterial {
        KOKKOS_INLINE_FUNCTION double operator()(material::CellIndex, material::Position,
                                                 bool& marked) const {
            marked = false;
            return 1.0;
        }
    };

    /** @brief Adapt the existing analytic wave to the manager's component sampler.
     * @details The shared history loader supplies the current and previous sample times.
     * referencePotential retains the wave of Fallahi, MITHRA 2.0, Eqs. (3.6)-(3.9).
     */
    struct VacuumPotential {
        KOKKOS_INLINE_FUNCTION double operator()(double, double, double z, double time,
                                                 unsigned component) const {
            return referencePotential(z, time)[component];
        }
    };

    simulation::RunBox vacuumBox(int cellsZ = 64);
    simulation::SolverSettings vacuumSettings(const simulation::RunBox& box);
}  // namespace chdr::apps::vacuumWave
