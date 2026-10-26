/** @file PeriodicMode.h
 * @brief Shared host validation for analytic periodic mode generators.
 */
#pragma once
#include "Ippl.h"

#include <array>
#include <cmath>
#include <stdexcept>

#include "core/base/SpatialDimension.h"

namespace chdr::solver {
    /** @brief Reject a wave whose phase fails to close on a declared periodic box.
     * @details Each axis requires an integer number of cycles. Checking this before
     * initialization prevents a seam discontinuity from being mistaken for solver
     * error. The tolerance tests host-generated mode metadata, not a field residual.
     */
    inline void requirePeriodicMode(ippl::Vector<double, SpatialDim> waveNumber,
                                    std::array<double, SpatialDim> lengths) {
        for (unsigned d = 0; d < SpatialDim; ++d) {
            const double mode = waveNumber[d] * lengths[d] / (2 * Kokkos::numbers::pi_v<double>);
            if (!std::isfinite(mode) || std::abs(mode - std::round(mode)) > 1e-12)
                throw std::invalid_argument("Wave does not fit the periodic box");
        }
    }
}  // namespace chdr::solver
