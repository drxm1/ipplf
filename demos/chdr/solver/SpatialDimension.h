/** @file SpatialDimension.h
 * @brief Share the spatial dimension without depending on solver definitions.
 */
#pragma once

namespace chdr::solver {
    /// Dimension of a cauchy slice.
    /// Must match the amount of components that \f$\vec{E}, \vec{B}, \vec{A}, \vec{J}\f$ have.
    constexpr unsigned SpatialDim = 3;
}  // namespace chdr::solver
