/** @file InterfaceOracles.h
 * @brief Expected values of the planar-interface tests: lattice dispersion and normal root of the
 * staggered Standard scheme in a uniform medium, exact lattice reflection and transmission of a
 * planar dielectric face on a cell boundary (TE and TM), the continuum Fresnel coefficients in the
 * same bases, and the frozen parameters of the plane-pulse runs.
 * @details Units \f$c = \varepsilon_0 = \mu_0 = 1\f$, \f$\mu_r = 1\f$, \f$\varepsilon \ge 1\f$,
 * \f$n = \sqrt\varepsilon\f$. The face normal is x. "Medium 1" (\f$\varepsilon_1\f$) holds the
 * incident wave and "medium 2" the transmitted one, on whichever side of the face they lie: the
 * lattice is mirror symmetric in x, so r and t do not depend on the side. A mode is
 * \f$e^{\mathrm i(k_xx + k_yy + k_zz - \omega t)}\f$; the Fourier probes use the kernel
 * \f$e^{+\mathrm i\omega t}\f$, which extracts this amplitude. Symbols of the Standard scheme:
 * \f$K_u = (2/h_u)\sin(k_uh_u/2)\f$, \f$W = (2/\Delta t)\sin(\omega\Delta t/2)\f$,
 * \f$K_t^2 = K_y^2 + K_z^2\f$ [fallahi2020mithra20fullwavesimulation, section 3.1.3, Eq. (3.17),
 * p. 14]. TE: \f$\mathbf E\f$ along the discrete transverse direction
 * \f$\hat{\mathbf s} = \hat{\mathbf x}\times\mathbf K_t/K_t\f$ (\f$\hat{\mathbf z}\f$ for
 * \f$k_z = 0\f$); TM: \f$\mathbf H\f$ along it. r and t are ratios of \f$E_s\f$ (TE) or \f$H_s\f$
 * (TM) at the face, so the phase is referred to the face. These host functions are written
 * independently of the solver; the tests take their expected values from here.
 */
#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <numbers>

namespace chdr::interface {
    using Complex = std::complex<double>;
    using Vec3    = std::array<double, 3>;
    enum class Polarization {
        TE,
        TM
    };
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
     * away from the face for reflected and transmitted waves.
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

    /// r, t at the face and the lattice normal wavenumbers of the incident and the other medium.
    struct LatticeCoefficients {
        Complex r, t, k1, k2;
    };

    /** @brief Exact lattice reflection and transmission of a planar face on a cell boundary.
     * @details With zero discrete gauge residual the staggered potentials give fields that obey
     * Yee's scheme [yee1966MaxwellIsotropicMedia, Eqs. (5)-(6), p. 303] with the harmonic
     * permittivity \f$\varepsilon_x = 2\varepsilon_1\varepsilon_2/(\varepsilon_1 +
     * \varepsilon_2)\f$ on the material face [eymard2000FiniteVolume, section 7.1, Example 7.1, p.
     * 22]. Per mode, TE and TM decouple into three-point recurrences along x (Deduction from the
     * discrete Faraday and Ampere laws): TE, \f$E_s\f$ on cell centres with the cell permittivity,
     * \f$E_{j+1} + E_{j-1} - 2E_j + h_x^2(\varepsilon_jW^2 - K_t^2)E_j = 0\f$; TM, \f$H_s\f$ on
     * x-faces, \f$(H_{m+1} - H_m)/\varepsilon_{m+1/2} - (H_m - H_{m-1})/\varepsilon_{m-1/2}
     * + h_x^2(W^2 - K_t^2/\varepsilon_{x,m})H_m = 0\f$, \f$\varepsilon_{x,m}\f$ the cell value
     * except \f$\varepsilon_x\f$ on the material face. Beyond the face each recurrence is
     * homogeneous, so every sample is a sum of its region's two modes, and both regional forms
     * hold on the two cells (TE) or the one face (TM) adjoining the material face. Eliminating
     * the overlap gives, with \f$h = h_x\f$ (Deduction; for TM the face mass term equals
     * \f$(1 - \cos k_1h)/\varepsilon_1 + (1 - \cos k_2h)/\varepsilon_2\f$ by the two dispersion
     * relations):
     * \f[ r_{\mathrm{TE}} = \frac{\sin((k_1 - k_2)h/2)}{\sin((k_1 + k_2)h/2)},\quad
     *     t_{\mathrm{TE}} = \frac{\sin(k_1h)}{\sin((k_1 + k_2)h/2)};\qquad
     *     r_{\mathrm{TM}} = \frac{\varepsilon_2\sin k_1h - \varepsilon_1\sin k_2h}
     *                            {\varepsilon_2\sin k_1h + \varepsilon_1\sin k_2h},\quad
     *     t_{\mathrm{TM}} = 1 + r_{\mathrm{TM}}. \f]
     * TM has one stored face value, so \f$t = 1 + r\f$ exactly; TE has not. For \f$h \to 0\f$
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
     * 0\f$ (decay away from the face under total reflection; the conjugate of the \f$e^{j\omega
     * t}\f$ branch [ellingson2020ElectromagneticsVol2, section 5.12, Eq. (5.288), p. 91]). TE:
     * \f$\Gamma_{\mathrm{TE}} = (\eta_2\cos\psi^i - \eta_1\cos\psi^t)/(\eta_2\cos\psi^i +
     * \eta_1\cos\psi^t)\f$ [ibid., section 5.6, Eq. (5.147), p. 75] with \f$\eta_j = \eta_0/n_j\f$,
     * \f$\cos\psi = k_x/(n\omega)\f$ gives \f$r = (k_{1x} - k_{2x})/(k_{1x} + k_{2x})\f$,
     * \f$t = 1 + r\f$ [ibid., Eq. (5.142)]. TM: \f$\Gamma_{\mathrm{TM}} = (-\eta_1\cos\psi^i +
     * \eta_2\cos\psi^t)/(\eta_1\cos\psi^i + \eta_2\cos\psi^t)\f$ [ibid., section 5.7, Eq. (5.191),
     * p. 79] refers the reflected \f$\mathbf H\f$ to the opposite direction [ibid., Eq. (5.156)];
     * in the \f$H_s\f$ basis \f$r = -\Gamma_{\mathrm{TM}} = (\varepsilon_2k_{1x} -
     * \varepsilon_1k_{2x})
     * /(\varepsilon_2k_{1x} + \varepsilon_1k_{2x})\f$ and \f$t = 1 - \Gamma_{\mathrm{TM}} = 1 +
     * r\f$ [ibid., Eq. (5.187), p. 78] (Deduction). Serves the type-B checks.
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

    /** @brief Frozen parameters of one plane-pulse run in a periodic box; \f$h_x = 1\f$ is the unit
     * of length.
     * @details The dielectric fills the cells \f$[0, D)\f$, the face is at \f$x_I = D\f$. The start
     * is the exact discrete sum of the box modes \f$k_m = 2\pi m/N_x\f$ with one sign of \f$k_x\f$
     * and the weight \f$e^{-(\lvert k_m\rvert - k_{x0})^2/(2\sigma_k^2)}\f$ for
     * \f$\lvert\lvert k_m\rvert - k_{x0}\rvert < 9\sigma_k\f$ (boxModes of them), each at its
     * lattice frequency in the incident medium; \f$\sigma_k = (1 - f)(k_{x0} -
     * k_{\mathrm{cut}})/q\f$ with \f$k_{\mathrm{cut}}\f$ the incident-medium \f$k_x\f$ at the
     * frequency where the vacuum side reaches its fixed-\f$k_y\f$ cutoff (0 for incidence from
     * vacuum). The envelope half-width is \f$w = q/\sigma_k\f$; the pulse centre lies probeCells +
     * w beyond the face on the incident side. The probe samples sit probeCells cells beyond the
     * face on each side (TE:
     * \f$E_z\f$ at the cell centre, physical distance probeCells + 1/2, index \f$x - 1/2\f$; TM:
     * averaged \f$B_z\f$ on the x-face, index \f$x - 1\f$). A gate holds the samples 0 to gate-1
     * (label \f$(k - 1/2)\Delta t\f$); steps is the number of solves. Box lengths, gates and steps
     * come from timing every contaminating path (second-face reflection, periodic wraps, the
     * periodic image of the closed-form incident) per spectral component with the lattice group
     * velocity; each such path must arrive after the end of the gate it would disturb.
     */
    struct PlaneRun {
        Polarization polarization;
        double epsilon;
        bool incidentFromVacuum;
        int cellsX, cellsY, cellsZ, dielectricCells, transversePeriods, probeCells;
        Vec3 h;
        double courant, q, f;
        int gateIncident, gateOther, steps, boxModes;
        double dt = 0, ky = 0, omega0 = 0, cellsPerWavelength = 0, kx0 = 0, omegaCut = 0, kCut = 0,
               sigmaK = 0, envelopeHalfWidth = 0, pulseCentre = 0, incidentProbe = 0,
               otherProbe      = 0;
        int incidentProbeIndex = 0, otherProbeIndex = 0;
    };

    /// Derived values of a run from its frozen integers and ratios (omega0 and ky set before).
    inline PlaneRun completeRun(PlaneRun r) {
        const double epsIn = r.incidentFromVacuum ? 1 : r.epsilon, D = r.dielectricCells;
        const double side = r.incidentFromVacuum ? 1 : -1;
        r.dt              = r.courant
               / std::sqrt(1 / (r.h[0] * r.h[0]) + 1 / (r.h[1] * r.h[1]) + 1 / (r.h[2] * r.h[2]));
        r.cellsPerWavelength = 2 * std::numbers::pi / (std::sqrt(r.epsilon) * r.omega0 * r.h[0]);
        r.kx0                = latticeNormal(r.omega0, r.ky, 0, epsIn, r.h, r.dt).kx.real();
        r.omegaCut = r.ky == 0 ? 0 : 2 / r.dt * std::asin(r.dt * symbolK(r.ky, r.h[1]) / 2);
        r.kCut     = r.incidentFromVacuum
                         ? 0
                         : latticeNormal(r.omegaCut, r.ky, 0, epsIn, r.h, r.dt).kx.real();
        r.sigmaK   = (1 - r.f) * (r.kx0 - r.kCut) / r.q;
        r.envelopeHalfWidth  = r.q / r.sigmaK;
        r.pulseCentre        = D + side * (r.probeCells + r.envelopeHalfWidth);
        const double offset  = r.polarization == Polarization::TE ? 0.5 : 1.0;
        r.incidentProbe      = D + side * (r.probeCells + 1 - offset);
        r.otherProbe         = D - side * (r.probeCells + 1 - offset);
        r.incidentProbeIndex = static_cast<int>(std::lround(r.incidentProbe - offset));
        r.otherProbeIndex    = static_cast<int>(std::lround(r.otherProbe - offset));
        return r;
    }

    /// TE at normal incidence from vacuum onto \f$\varepsilon_r = 2.13\f$, 20 cells per dielectric
    /// wavelength (\f$\omega_0 = 2\pi/(20n)\f$), cubic cells, \f$S = 0.9\f$, \f$q = 8\f$, \f$f =
    /// 0\f$.
    inline PlaneRun normalIncidenceRun() {
        PlaneRun r{Polarization::TE, 2.13, true, 812, 2,    2,    211,  0, 3,
                   {1, 1, 1},        0.9,  8,    0,   1175, 1178, 1180, 59};
        r.ky     = 0;
        r.omega0 = 2 * std::numbers::pi / (20 * std::sqrt(2.13));
        return completeRun(r);
    }
    /// TM from vacuum at the Brewster angle \f$\arctan n\f$ at the band centre
    /// (\f$\omega_0 = k_y\sqrt{1 + \varepsilon}/n\f$), one period over 18 cells of \f$h_y = 2\f$,
    /// cells \f$(1, 2, 2)\f$, \f$S = 0.9\f$, \f$q = 7\f$, \f$f = 0.5\f$.
    inline PlaneRun brewsterIncidenceRun() {
        PlaneRun r{Polarization::TM, 2.13, true, 3581, 18,   2,    1484, 1, 3,
                   {1, 2, 2},        0.9,  7,    0.5,  6759, 6755, 6761, 88};
        r.ky     = 2 * std::numbers::pi / (18 * 2);
        r.omega0 = r.ky * std::sqrt(1 + 2.13) / std::sqrt(2.13);
        return completeRun(r);
    }
    /// TM from the dielectric \f$\varepsilon_r = 3.4\f$ at 18 degrees at the band centre
    /// (\f$\omega_0 = k_y/(n\sin18^\circ)\f$), one period over 52 cells of \f$h_y = 1.25\f$, cells
    /// \f$(1, 1.25, 0.85)\f$, \f$S = 0.9\f$, \f$q = 7\f$, \f$f = 0.3\f$.
    inline PlaneRun obliqueDielectricRun() {
        PlaneRun r{Polarization::TM, 3.4, false, 1828, 52,   2,    949,  1, 3,
                   {1, 1.25, 0.85},  0.9, 7,     0.3,  3787, 3785, 3789, 78};
        r.ky     = 2 * std::numbers::pi / (52 * 1.25);
        r.omega0 = r.ky / (std::sqrt(3.4) * std::sin(std::numbers::pi / 10));
        return completeRun(r);
    }
}  // namespace chdr::interface
