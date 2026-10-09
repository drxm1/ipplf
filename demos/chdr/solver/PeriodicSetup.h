/** @file PeriodicSetup.h
 * @brief Runtime periodic boxes with a single selectable decomposition axis.
 */
#pragma once

#include <array>

#include "ChdrSolver.h"

namespace chdr::solver {
    /** @brief Construct mesh/layout ownership; allocate with initializeFieldStorage.
     * @pre All ranks must call it with identical arguments; lengths are finite and positive.
     * @param cells     Global cell counts in x, y and z; each must be at least two.
     * @param lengths   Physical box lengths in x, y and z; finite and positive.
     *                  Size of the whole container, in each axis.
     * @param splitAxis Axis eligible for MPI decomposition: 0=x, 1=y, 2=z.
     *                  The other two axes remain undistributed. Defaults to x.
     * @return Prvalue container; retain its address once its fields are allocated.
     * @details Construction defines mesh and MPI ownership only. Allocation is a separate
     * operation because the field objects retain pointers to this owner's mesh and layout.
     * Do not return or move an already allocated container through a setup helper.
     */
    ChDRFieldContainer makePeriodicFieldContainer(std::array<int, 3> cells,
                                                  std::array<double, 3> lengths,
                                                  unsigned splitAxis = 0);

    /// Reusing from the vacuum wave example
    using chdr::vacuum::checkLocalDomain;

    /// Reusing from the vacuum wave example
    using chdr::vacuum::initializeFieldStorage;

    /// Reusing from the vacuum wave example
    using chdr::vacuum::preparePotentialHalo;

}  // namespace chdr::solver
