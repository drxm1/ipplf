/** @file MonitorTestSupport.h
 * @brief Shared final-transition checks for evolved periodic test cases.
 *
 * TODO: REVIEW.
 */
#pragma once
#include <gtest/gtest.h>

#include "Monitors.h"
#include "SpatialDimension.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "MonitorTestSupport.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    /** @brief Check a prescribed continuity residual at every owned cell.
     * @details The source fixture constructs charge with independent neighbour
     * differences. This subtracts its prescribed value before normalization,
     * so a local downward error cannot hide behind the global maximum or norm.
     * The continuity expression and time weights are documented in Monitors::gauge.
     * @param currentSource Source consumed by the transition under test: charge at
     * component 0 and current density at components 1 through SpatialDim.
     * @param previousSource Source consumed one timestep earlier, with the same layout.
     * @param dt Positive timestep separating the two source histories.
     * @param weight Current-source interpolation weight from LevelContract.
     * @param expected Prescribed continuity residual to subtract at every owned cell.
     * @param spatialDiscretization Source-component positions and difference directions.
     * @pre Current and previous source halos are valid for the selected differences.
     */
    inline ResidualMetrics expectedContinuityMismatch(
        const AJFourFieldsType& currentSource, const AJFourFieldsType& previousSource, double dt,
        double weight, double expected,
        SpatialDiscretization spatialDiscretization = SpatialDiscretization::Collocated) {
        const auto now = currentSource.getView(), old = previousSource.getView();
        const auto h        = currentSource.get_mesh().getMeshSpacing();
        const int direction = spatialDiscretization == SpatialDiscretization::Staggered ? -1 : 0;
        return reduceResidual(
            currentSource, KOKKOS_LAMBDA(int i, int j, int k) {
                return (monitorDetail::Expression(now(i, j, k)[0])
                        - monitorDetail::Expression(old(i, j, k)[0]))
                           * (1.0 / dt)
                       + monitorDetail::divergence(now, 1, i, j, k, h, direction) * weight
                       + monitorDetail::divergence(old, 1, i, j, k, h, direction) * (1.0 - weight)
                       - monitorDetail::Expression(expected);
            });
    }

}  // namespace chdr::solver::test
