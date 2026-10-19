/** @file RectangleXZ.h
 * @brief Rectangular material region in the x-z plane, uniform along y.
 */
#pragma once
#include <Kokkos_Core.hpp>

namespace chdr::material {
    /** @brief Select a rectangle by its closed x and z intervals.
     * @details The material map samples this predicate at cell centres. Bounds on cell
     * faces therefore select whole material cells; ties are included, as in PlanarSlab.
     * The y coordinate is unrestricted so the same map serves a thin periodic extrusion.
     */
    struct RectangleXZ {
        double lowerX_m, upperX_m, lowerZ_m, upperZ_m;
        KOKKOS_INLINE_FUNCTION bool contains(double x, double, double z) const {
            return lowerX_m <= x && x <= upperX_m && lowerZ_m <= z && z <= upperZ_m;
        }
    };
}  // namespace chdr::material
