/** @file AnnulusStart.h
 * @brief Localized, approximately outgoing TE annulus in the x-z plane.
 */
#pragma once
#include <Kokkos_Core.hpp>

namespace chdr::apps::interfaceMovies {
    /** @brief Sample an annulus in the y component of the vector potential.
     * @details With radial distance \f$r\f$ and \f$s=r-r_0-t\f$, the profile is
     * \f$A_y=a(s/\lambda)\exp[-(\pi s/\lambda)^2]\sqrt{r_0/\max(r,h)}\f$.
     * Its radial profile is an antiderivative of a Ricker pulse
     * [schneider2010UnderstandingFdtd, section 5.2.3, Eq. (5.20)]. The radial factor
     * is a large-radius outgoing approximation: away from the centre, substituting
     * \f$r^{-1/2}g(r-t)\f$ into the radial vacuum wave equation leaves a residual
     * \f$-A_y/(4r^2)\f$ (deduction from [chew2014generalizedgauge, Eqs. (16)-(18)]).
     * This is initial data, not an exact outgoing solution of the discrete scheme.
     * Small Gaussian tails outside \f$|s|\le1.5\lambda\f$ are set to zero so both
     * histories can be placed strictly inside vacuum. guard_m keeps the centre finite.
     * Only component 2 is nonzero; component 0 is the scalar potential.
     * @pre Positive radius, wavelength and guard; both histories' support lies in
     * vacuum and clear of open faces. Coordinates and time use normalized vacuum speed one.
     */
    struct AnnulusStart {
        double centreX_m = 150, centreZ_m = 300;
        double radius_m = 72, wavelength_m = 36, guard_m = 1, amplitude_m = 1;

        /// Radius enclosing every nonzero sample at the supplied potential time.
        double outerRadius(double time) const { return radius_m + time + 1.5 * wavelength_m; }

        KOKKOS_INLINE_FUNCTION double operator()(double x, double, double z, double time,
                                                 unsigned component) const {
            if (component != 2)
                return 0;
            const double dx = x - centreX_m, dz = z - centreZ_m;
            const double r = Kokkos::sqrt(dx * dx + dz * dz), s = r - radius_m - time;
            if (Kokkos::abs(s) > 1.5 * wavelength_m)
                return 0;
            const double u = Kokkos::numbers::pi * s / wavelength_m;
            return amplitude_m * (s / wavelength_m) * Kokkos::exp(-u * u)
                   * Kokkos::sqrt(radius_m / Kokkos::max(r, guard_m));
        }
    };
}  // namespace chdr::apps::interfaceMovies
