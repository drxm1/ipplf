/** @file VacuumSetup.h
 * @brief Periodic grid and potential initialization of the vacuum-wave example.
 */
#pragma once

#include "Ippl.h"

#include "FELFieldContainer.hpp"
#include "MaxwellSolvers/StandardFDTDSolver.h"
#include "core/base/FieldContainer.h"

/// Vacuum-wave setup; callers initialize IPPL and retain the owning fields.
namespace chdr::apps::vacuumWave {
    constexpr unsigned Dim = 3;  ///< Cartesian coordinates, ordered x, y, z.
    /// The shared field container of `core/base`, under the example's historical name.
    using Fields = ChDRFieldContainer;
    /// Ownership check and storage allocation come from `core/base`.
    using chdr::checkLocalDomain;
    using chdr::initializeFieldStorage;
    using VacuumSolver =
        ippl::StandardFDTDSolver<VField_t<double, Dim>, SourceField_t<double, Dim>, ippl::periodic>;

    Fields makePeriodicFieldContainer(int cellsZ = 64);
    void initializePotentials(VacuumSolver& solver);
}  // namespace chdr::apps::vacuumWave
