/** @file testing/Util.h
 * @brief Utilities shared by the ChDR test executables.
 */
#pragma once

#include "Ippl.h"

#include <stdexcept>

#include "core/base/FieldContainer.h"

namespace chdr::testing::util {
    /** @brief Finish pending writes and synchronize one potential's periodic halos.
     * @param potential Allocated field with one-cell halos and periodic layout/faces.
     * @pre All layout ranks call in matching order; IPPL/Kokkos are initialized.
     * @details Debug builds check the local configuration; the fence completes kernel work.
     * Production code prepares both histories of a solver with
     * FDTDSolverBase::prepareHistories(); tests use this helper for one field at a time.
     */
    inline void preparePotentialHalo(AJFourFieldsType& potential) {
        Kokkos::fence();
#ifndef NDEBUG
        if (!potential.getView().data())
            throw std::logic_error("Potential storage must be allocated before halo preparation");
        if (potential.getNghost() != 1 || !potential.getLayout().isAllPeriodic_m)
            throw std::logic_error("Expected one-cell halos on a periodic layout");
        for (unsigned face = 0; face < 2 * SpatialDim; ++face) {
            const auto& bc = potential.getFieldBC()[face];
            if (!bc || bc->getBCType() != ippl::PERIODIC_FACE || bc->getFace() != face)
                throw std::logic_error("Expected a matching periodic condition on every face");
        }
#endif
        potential.fillHalo();
        potential.getFieldBC().apply(potential);
    }
}  // namespace chdr::testing::util
