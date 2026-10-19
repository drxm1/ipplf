/** @file PacketStart.h
 * @brief Localized transverse packet for qualitative periodic interface demonstrations.
 */
#pragma once
#include "PlanarSlab.h"

namespace chdr::interface {
    /** @brief Sample a Gaussian carrier, or its weighted discrete curl, at two time levels.
     * @details TE uses \f$A_y\f$ in the x-z plane. TM uses the adjoint forward curl of a scalar
     * edge field, divided by the harmonic face permittivity. The backward divergence of this
     * weighted vector cancels identically (commuting differences), so \f$\Phi=0\f$ has zero initial
     * gauge residual [chew2014generalizedgauge, section 2, Eqs. (15)-(18); ryu2016potentialfdtd,
     * Appendix A]. The translating Gaussian is an approximate narrow-band start, suitable for field
     * pictures; PlanePulse supplies exact lattice modes.
     * @pre The demonstration uses unit grid spacing and periodic x and z extents
     * both equal to length_m, with coordinates starting at zero. The fixed half-cell
     * offsets and omitted spacing divisors in the TM differences assume this spacing.
     * @pre slab_m.axis_m is zero: the material varies only along x. slab_m and
     * epsilon_m must match the solver's two-material map, with permittivity one
     * outside the slab. The sampler's material lookup wraps x periodically.
     * @pre directionX_m and directionY_m form a unit direction vector; packet
     * widths and length_m are positive. The translating start assumes vacuum
     * propagation speed, so the incident packet should be localized in vacuum.
     * @details These conditions are supplied by the demonstration driver; this
     * sampler does not validate them. It is independent of y and describes a
     * packet extruded along that periodic axis.
     */
    struct PacketStart {
        material::PlanarSlab slab_m;
        double epsilon_m = 2.13, centreX_m = 371, centreZ_m = 256, length_m = 512;
        double directionX_m = -1, directionZ_m = 0, width_m = 29, transverseWidth_m = 45;
        double waveNumber_m = 0.215258381712;
        bool tm_m           = false;
        KOKKOS_INLINE_FUNCTION double psi(double x, double z, double time) const {
            double dx = x - centreX_m, dz = z - centreZ_m;
            dx -= Kokkos::floor(dx / length_m + 0.5) * length_m;
            dz -= Kokkos::floor(dz / length_m + 0.5) * length_m;
            const double along  = directionX_m * dx + directionZ_m * dz - time;
            const double across = -directionZ_m * dx + directionX_m * dz;
            return Kokkos::exp(-0.5
                               * (along * along / (width_m * width_m)
                                  + across * across / (transverseWidth_m * transverseWidth_m)))
                   * Kokkos::cos(waveNumber_m * along);
        }
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double t,
                                                 unsigned c) const {
            if (!tm_m)
                return c == 2 ? psi(x, z, t) : 0;
            if (c != 1 && c != 3)
                return 0;
            double left = x - (c == 1 ? 0.5 : 0), right = x + (c == 1 ? 0.5 : 0);
            left -= Kokkos::floor(left / length_m) * length_m;
            right -= Kokkos::floor(right / length_m) * length_m;
            const double a    = slab_m.contains(left, y, z) ? epsilon_m : 1;
            const double b    = slab_m.contains(right, y, z) ? epsilon_m : 1;
            const double curl = c == 1 ? psi(x, z - 0.5, t) - psi(x, z + 0.5, t)
                                       : psi(x + 0.5, z, t) - psi(x - 0.5, z, t);
            return curl * (1 / a + 1 / b) / 2;
        }
    };
}  // namespace chdr::interface
