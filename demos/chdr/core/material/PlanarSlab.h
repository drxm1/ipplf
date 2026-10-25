/** @file PlanarSlab.h
 * @brief Closed planar slab sampled at material cell centres.
 */
#pragma once
#include <Kokkos_Core.hpp>

namespace chdr::material {
    /** @brief Select \f$\mathrm{lower}\le x_{\mathrm{axis}}\le\mathrm{upper}\f$, with ties inside.
     * @details axis_m is 0, 1 or 2. A half-space in a periodic box is a slab
     * extending to a box end; its second interface lies at the periodic wrap.
     * Requested faces are retained here; realized faces follow cell sampling.
     */
    struct PlanarSlab {
        unsigned axis_m;
        double lower_m, upper_m;
        KOKKOS_INLINE_FUNCTION bool contains(double x, double y, double z) const {
            const double position = axis_m == 0 ? x : axis_m == 1 ? y : z;
            return lower_m <= position && position <= upper_m;
        }
    };
}  // namespace chdr::material
