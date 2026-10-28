/** @file TMAnnulusStart.h
 * @brief Source-free TM initial histories from the discrete curl of a radial profile.
 */
#pragma once
#include "AnnulusStart.h"

namespace chdr::apps::interfaceMovies {
    /** @brief Initialize in-plane vector potential and out-of-plane magnetic field.
     * @details Sample the radial profile on y-directed edges and form
     * \f$A_x=-D_z^-\psi\f$, \f$A_z=D_x^-\psi\f$ at their natural face positions.
     * Commuting differences cancel the backward divergence of A; because the entire
     * support is in vacuum, its permittivity-weighted divergence also vanishes initially
     * [ryu2016potentialfdtd, Appendix A; chew2014generalizedgauge, Eqs. (15)-(18)].
     * Both histories use this construction, with scalar potential zero. The profile
     * remains an approximate outgoing start (AnnulusStart). This initial cancellation
     * does not assert constraint preservation by later componentwise Mur updates.
     * The extra spatial difference gives TM fields one higher radial derivative than
     * TE. At large radius their peak wavelength is approximately
     * \f$\lambda/\sqrt{3/2}\f$, where \f$\lambda\f$ is the profile wavelength parameter.
     * This is a deduction from the extra wavenumber factor multiplying the Ricker spectrum
     * [schneider2010UnderstandingFdtd, section 5.2.3, Eq. (5.20)]. The field waveform
     * therefore differs from TE even when the same profile parameters are supplied.
     * @pre Positive spacings; the profile support enlarged by half the larger spacing
     * lies in vacuum and outside the first three owned layers at each open face,
     * at both history times. This keeps quadratic open-ghost extrapolation zero initially.
     */
    struct TMAnnulusStart {
        AnnulusStart profile_m;
        double hx_m = 1, hz_m = 1;

        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double time,
                                                 unsigned component) const {
            if (component == 1)
                return (profile_m(x, y, z - hz_m / 2, time, 2)
                        - profile_m(x, y, z + hz_m / 2, time, 2))
                       / hz_m;
            if (component == 3)
                return (profile_m(x + hx_m / 2, y, z, time, 2)
                        - profile_m(x - hx_m / 2, y, z, time, 2))
                       / hx_m;
            return 0;
        }
    };
}  // namespace chdr::apps::interfaceMovies
