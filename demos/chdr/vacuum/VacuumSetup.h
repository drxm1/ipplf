/** @file VacuumSetup.h
 * @brief Shared periodic-grid, field-allocation and potential-initialization interface.
 */
#pragma once

#include "Ippl.h"

#include "FELFieldContainer.hpp"
#include "MaxwellSolvers/StandardFDTDSolver.h"

/// Shared demo/test setup; callers initialize IPPL and retain the owning fields.
namespace chdr::vacuum {
    constexpr unsigned Dim = 3;  ///< Cartesian coordinates, ordered x, y, z.
    using Fields           = FELFieldContainer<double, Dim>;
    using VacuumSolver =
        ippl::StandardFDTDSolver<VField_t<double, Dim>, SourceField_t<double, Dim>, ippl::periodic>;

    Fields makePeriodicFieldContainer(int cellsZ = 64);
    void checkLocalDomain(const ippl::FieldLayout<Dim>& layout);
    void initializeFieldStorage(Fields& fields);
    /** @brief Synchronize allocated one-cell potential halos and periodic faces.
     * @pre All layout ranks participate; local configuration is checked in debug builds.
     */
    void preparePotentialHalo(SourceField_t<double, Dim>& potential);
    void initializePotentials(VacuumSolver& solver);
}  // namespace chdr::vacuum
