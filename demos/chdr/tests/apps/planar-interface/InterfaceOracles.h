/** @file InterfaceOracles.h
 * @brief Expected values of the planar-interface tests: lattice dispersion and normal root of the
 * staggered Standard scheme in a uniform medium, exact lattice reflection and transmission of a
 * planar dielectric interface on a cell boundary (TE and TM), the continuum Fresnel coefficients in
 * the same bases. PlanePresets.h documents the three recorded run configurations and their pulse,
 * probe and gate conventions; these formulas remain independent of those configurations.
 * @details Units \f$c = \varepsilon_0 = \mu_0 = 1\f$, \f$\mu_r = 1\f$, \f$\varepsilon \ge 1\f$,
 * \f$n = \sqrt\varepsilon\f$. The interface normal is x. "Medium 1" (\f$\varepsilon_1\f$) holds the
 * incident wave and "medium 2" the transmitted one, on whichever side of the interface they lie:
 * the lattice is mirror symmetric in x, so r and t do not depend on the side. A mode is
 * \f$e^{\mathrm i(k_xx + k_yy + k_zz - \omega t)}\f$; the Fourier probes use the kernel
 * \f$e^{+\mathrm i\omega t}\f$, which extracts this amplitude. Symbols of the Standard scheme:
 * \f$K_u = (2/h_u)\sin(k_uh_u/2)\f$, \f$W = (2/\Delta t)\sin(\omega\Delta t/2)\f$,
 * \f$K_t^2 = K_y^2 + K_z^2\f$ [fallahi2020mithra20fullwavesimulation, section 3.1.3, Eq. (3.17),
 * p. 14]. TE: \f$\mathbf E\f$ along the discrete transverse direction
 * \f$\hat{\mathbf s} = \hat{\mathbf x}\times\mathbf K_t/K_t\f$ (\f$\hat{\mathbf z}\f$ for
 * \f$k_z = 0\f$); TM: \f$\mathbf H\f$ along it. r and t are ratios of \f$E_s\f$ (TE) or \f$H_s\f$
 * (TM) at the material interface, which is their phase origin. A common sign change of the
 * amplitude basis leaves these ratios unchanged; the x-z application uses fixed +y.
 * These host functions are written
 * independently of the solver; the tests take their expected values from here.
 */
#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <numbers>

#include "apps/planar-interface/PlaneRow.h"

namespace chdr::apps::planarInterface {
    using Complex = std::complex<double>;
    using Vec3    = std::array<double, 3>;
    enum class LatticeRoot {
        Propagating,
        Evanescent,
        Sawtooth
    };

    /// \f$K = (2/h)\sin(kh/2)\f$: the symbol of one compact second difference is \f$-K^2\f$.
    inline Complex symbolK(Complex k, double h) {
        return 2.0 * std::sin(k * (h / 2)) / h;
    }
    inline double symbolK(double k, double h) {
        return 2 * std::sin(k * h / 2) / h;
    }
    /// \f$W = (2/\Delta t)\sin(\omega\Delta t/2)\f$: the leapfrog second difference gives
    /// \f$-W^2\f$.
    inline double symbolW(double omega, double dt) {
        return 2 * std::sin(omega * dt / 2) / dt;
    }

    /// Normal wavenumber of one lattice mode: real positive, imaginary positive, or sawtooth.
    struct LatticeNormal {
        Complex kx;
        LatticeRoot root;
    };

    /** @brief Lattice normal wavenumber of a uniform medium at \f$(\omega, k_y, k_z)\f$.
     * @details For a transverse mode (zero discrete divergence) the generalized-gauge potential
     * update read staggered [chew2014generalizedgauge, section 2, Eqs. (15)-(18), p. 3;
     * ryu2016potentialfdtd, Appendix A, Eqs. (19)-(24)] has the symbol of the Standard vacuum
     * update with the mass \f$\varepsilon\f$ (Deduction: its divergence correction vanishes):
     * \f$\varepsilon W^2 = K_x^2 + K_y^2 + K_z^2\f$ [fallahi2020mithra20fullwavesimulation,
     * Eq. (3.17)], so \f$\sin^2(k_xh_x/2) = s := h_x^2(\varepsilon W^2 - K_t^2)/4\f$. Branches:
     * \f$0 < s < 1\f$ propagating, \f$k_x = (2/h_x)\arcsin\sqrt s > 0\f$, the root with positive
     * group velocity \f$\sin(k_xh_x)/(h_x\varepsilon W\cos(\omega\Delta t/2))\f$
     * [trefethen1982GroupVelocity, section 1, Eqs. (1.4)-(1.5), p. 115]; \f$s < 0\f$ evanescent,
     * \f$k_x = (2\mathrm i/h_x)\,\mathrm{asinh}\sqrt{-s}\f$, decaying toward +x; \f$s > 1\f$
     * sawtooth, \f$k_x = \pi/h_x + (2\mathrm i/h_x)\,\mathrm{acosh}\sqrt s\f$. This is the root
     * continued from \f$\mathrm{Im}\,\omega > 0\f$ [gustafsson1972StabilityTheoryII, section 5,
     * Lemma 5.1, p. 659]; the caller gives \f$k_x\f$ the sign of the travel or decay direction,
     * away from the material interface for reflected and transmitted waves.
     */
    inline LatticeNormal latticeNormal(double omega, double ky, double kz, double eps, Vec3 h,
                                       double dt) {
        const double W = symbolW(omega, dt), Ky = symbolK(ky, h[1]), Kz = symbolK(kz, h[2]);
        const double s = (eps * W * W - Ky * Ky - Kz * Kz) * h[0] * h[0] / 4;
        if (s < 0)
            return {Complex(0, 2 * std::asinh(std::sqrt(-s)) / h[0]), LatticeRoot::Evanescent};
        if (s > 1)
            return {Complex(std::numbers::pi / h[0], 2 * std::acosh(std::sqrt(s)) / h[0]),
                    LatticeRoot::Sawtooth};
        return {Complex(2 * std::asin(std::sqrt(s)) / h[0], 0), LatticeRoot::Propagating};
    }

    /// Interface-referred r, t and the normal lattice roots in the incident and transmitted media.
    struct LatticeCoefficients {
        Complex r, t, k1, k2;
    };

    /** @brief Exact lattice reflection and transmission of a material interface on a cell boundary.
     * @details With zero discrete gauge residual the staggered potentials give fields that obey
     * Yee's scheme [yee1966MaxwellIsotropicMedia, Eqs. (5)-(6), p. 303] with the harmonic
     * permittivity \f$\varepsilon_x = 2\varepsilon_1\varepsilon_2/(\varepsilon_1 +
     * \varepsilon_2)\f$ on the coincident cell face [eymard2000FiniteVolume, section 7.1,
     * Example 7.1, p. 22]. Per mode, TE and TM decouple into three-point recurrences along x
     * (Deduction from the discrete Faraday and Ampere laws): TE, \f$E_s\f$ on cell centres with the
     * cell permittivity,
     * \f$E_{j+1} + E_{j-1} - 2E_j + h_x^2(\varepsilon_jW^2 - K_t^2)E_j = 0\f$; TM, \f$H_s\f$ on
     * x-faces, \f$(H_{m+1} - H_m)/\varepsilon_{m+1/2} - (H_m - H_{m-1})/\varepsilon_{m-1/2}
     * + h_x^2(W^2 - K_t^2/\varepsilon_{x,m})H_m = 0\f$, \f$\varepsilon_{x,m}\f$ the cell value
     * except \f$\varepsilon_x\f$ on the interface cell face. Away from the interface each
     * recurrence is homogeneous, so every sample is a sum of its region's two modes, and both
     * regional forms hold on the two cells (TE) or the one cell face (TM) adjoining the interface.
     * Eliminating the overlap gives, with \f$h = h_x\f$ (Deduction; for TM the mass term at that
     * cell face equals
     * \f$(1 - \cos k_1h)/\varepsilon_1 + (1 - \cos k_2h)/\varepsilon_2\f$ by the two dispersion
     * relations):
     * \f[ r_{\mathrm{TE}} = \frac{\sin((k_1 - k_2)h/2)}{\sin((k_1 + k_2)h/2)},\quad
     *     t_{\mathrm{TE}} = \frac{\sin(k_1h)}{\sin((k_1 + k_2)h/2)};\qquad
     *     r_{\mathrm{TM}} = \frac{\varepsilon_2\sin k_1h - \varepsilon_1\sin k_2h}
     *                            {\varepsilon_2\sin k_1h + \varepsilon_1\sin k_2h},\quad
     *     t_{\mathrm{TM}} = 1 + r_{\mathrm{TM}}. \f]
     * TM has one stored value at the interface cell face, so \f$t = 1 + r\f$ exactly; TE has not.
     * For \f$h \to 0\f$
     * both become Fresnel's forms in the wavenumbers; an evanescent side has complex \f$k\f$ and
     * \f$\lvert r\rvert = 1\f$. Serves the type-A coefficient checks at the band centre.
     */
    inline LatticeCoefficients latticeCoefficients(Polarization p, double eps1, double eps2, Vec3 h,
                                                   double dt, double ky, double kz, double omega) {
        const Complex k1 = latticeNormal(omega, ky, kz, eps1, h, dt).kx;
        const Complex k2 = latticeNormal(omega, ky, kz, eps2, h, dt).kx;
        if (p == Polarization::TE) {
            const Complex sum = std::sin((k1 + k2) * (h[0] / 2));
            return {std::sin((k1 - k2) * (h[0] / 2)) / sum, std::sin(k1 * h[0]) / sum, k1, k2};
        }
        const Complex s1 = eps2 * std::sin(k1 * h[0]), s2 = eps1 * std::sin(k2 * h[0]);
        const Complex r = (s1 - s2) / (s1 + s2);
        return {r, 1.0 + r, k1, k2};
    }

    /** @brief Continuum Fresnel coefficients in the same bases, with the normal wavenumbers.
     * @details \f$k_{jx} = \sqrt{\varepsilon_j\omega^2 - k_t^2}\f$, branch \f$\mathrm{Im}\,k \ge
     * 0\f$ (decay away from the material interface under total reflection; the conjugate of the
     * \f$e^{j\omega t}\f$ branch [ellingson2020ElectromagneticsVol2, section 5.12, Eq. (5.288), p.
     * 91]). TE:
     * \f$\Gamma_{\mathrm{TE}} = (\eta_2\cos\psi^i - \eta_1\cos\psi^t)/(\eta_2\cos\psi^i +
     * \eta_1\cos\psi^t)\f$ [ibid., section 5.6, Eq. (5.147), p. 75] with \f$\eta_j = \eta_0/n_j\f$,
     * \f$\cos\psi = k_x/(n\omega)\f$ gives \f$r = (k_{1x} - k_{2x})/(k_{1x} + k_{2x})\f$,
     * \f$t = 1 + r\f$ [ibid., Eq. (5.142)]. TM: \f$\Gamma_{\mathrm{TM}} = (-\eta_1\cos\psi^i +
     * \eta_2\cos\psi^t)/(\eta_1\cos\psi^i + \eta_2\cos\psi^t)\f$ [ibid., section 5.7, Eq. (5.191),
     * p. 79] refers the reflected \f$\mathbf H\f$ to the opposite direction [ibid., Eq. (5.156)];
     * in the \f$H_s\f$ basis \f$r = -\Gamma_{\mathrm{TM}} = (\varepsilon_2k_{1x} -
     * \varepsilon_1k_{2x})
     * /(\varepsilon_2k_{1x} + \varepsilon_1k_{2x})\f$ and \f$t = 1 - \Gamma_{\mathrm{TM}} = 1 +
     * r\f$ [ibid., Eq. (5.187), p. 78] (Deduction). Serves the continuum comparisons.
     */
    inline LatticeCoefficients fresnelCoefficients(Polarization p, double eps1, double eps2,
                                                   double ky, double kz, double omega) {
        const double kt2 = ky * ky + kz * kz;
        const Complex k1 = std::sqrt(Complex(eps1 * omega * omega - kt2));
        const Complex k2 = std::sqrt(Complex(eps2 * omega * omega - kt2));
        const Complex a  = p == Polarization::TE ? k1 : eps2 * k1;
        const Complex b  = p == Polarization::TE ? k2 : eps1 * k2;
        const Complex r  = (a - b) / (a + b);
        return {r, 1.0 + r, k1, k2};
    }

}  // namespace chdr::apps::planarInterface
