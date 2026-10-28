/** @file ReflectionOracles.h
 * @brief Expected values of the reflection tests: the discrete reflection coefficient of Mur's
 * second-order face on the leapfrog lattice, its continuum value and observed order, the field
 * rows, the edge and corner plane-wave families and the reflected gauge residual.
 * @details One potential component, interior speed \f$v_u\f$, Mur speed \f$v_m\f$, physical step
 * \f$\Delta t\f$ (the kernel of a medium runs at \f$v_u\Delta t\f$ with speed 1; only
 * \f$v\Delta t\f$ enters). At a face with normal axis \f$\nu\f$, layer \f$j = 0\f$ is the outermost
 * owned layer (written by Mur) and j grows into the box. A lattice mode is
 * \f$\kappa^j e^{\mathrm i(\mathbf k_t\cdot\mathbf x_t - \omega k\Delta t)}\f$ with
 * \f$\kappa = e^{\mathrm i\theta}\f$, \f$\theta = k_\nu h_\nu\f$; the incident root is
 * \f$-\theta\f$ of the reflected root \f$\theta\f$ (the interior symbols are even in
 * \f$\theta\f$). Symbols: \f$W = (2/\Delta t)\sin(\omega\Delta t/2)\f$,
 * \f$C_t = \cos(\omega\Delta t/2)\f$, \f$K_d^2 = (4/h_d^2)\sin^2(k_dh_d/2)\f$.
 */
#pragma once
#include <array>
#include <cmath>
#include <complex>
#include <numbers>

namespace chdr::solver {
    using OracleComplex = std::complex<double>;
    enum class OracleStencil {
        Standard,
        NonStandard
    };
    enum class LatticeModeClass {
        Propagating,
        Evanescent,
        GridBranch
    };

    /// Reflection data at one face: classification, roots \f$\theta\f$, \f$R_{1/2}\f$, \f$R_0\f$.
    struct DiscreteReflection {
        LatticeModeClass modeClass;
        OracleComplex reflectedRoot, incidentRoot, centre, layerZero;
    };

    inline double squaredSymbol(double k, double h) {
        const double s = 2 * std::sin(k * h / 2) / h;
        return s * s;
    }

    /** @brief Discrete reflection coefficient of Mur's face for one tangential mode.
     * @details Interior: substituting the mode in the leapfrog update gives
     * \f$W^2/v_u^2 = D'K_\nu^2 + \mathrm{rest}\f$; Standard \f$D' = 1\f$, rest the side
     * \f$K^2\f$; NonStandard with \f$\mathcal A = \frac14(1 + 0.02/((h_z/h_x)^2 + (h_z/h_y)^2))\f$:
     * face normal to z, \f$D' = 1 - \mathcal Ah_z^2(K_x^2 + K_y^2)\f$, rest \f$K_x^2 + K_y^2\f$;
     * otherwise \f$D' = 1 - \mathcal Ah_z^2K_z^2\f$, rest \f$D'K_s^2 + K_z^2\f$
     * [fallahi2020mithra20fullwavesimulation, section 3.1.3, Eq. (3.24), p. 15] (Deduction).
     * Hence \f$q = \sin^2(\theta/2) = (W^2/v_u^2 - \mathrm{rest})h_\nu^2/(4D')\f$. The reflected
     * root carries energy into the box, group velocity
     * \f$v_u^2D'\Delta t\sin\theta/(h_\nu\sin\omega\Delta t) > 0\f$, or decays into it,
     * \f$|\kappa| < 1\f$ [trefethen1982GroupVelocity, section 5, p. 131]:
     * \f$0 < q < 1\f$ propagating, \f$\theta = \mathrm{sgn}(D')2\arcsin\sqrt q\f$; \f$q < 0\f$
     * evanescent, \f$\theta = 2\mathrm i\,\mathrm{asinh}\sqrt{-q}\f$; \f$q > 1\f$ grid branch,
     * \f$\theta = \pi + 2\mathrm i\,\mathrm{acosh}\sqrt q\f$. Face: Mur's centred row
     * [mur1981AbsorbingBoundaryFdtd, section IV, Eq. (16), p. 380] applied to one mode is
     * \f$P(\kappa) = (e^{\mathrm i\theta/2}/v_m)(X + Y)\f$ with
     * \f$X = v_mWK_\nu C_t\f$, \f$Y = (W^2 - v_m^2K_f^2/2)C_\nu\f$,
     * \f$K_\nu = (2/h_\nu)\sin(\theta/2)\f$, \f$C_\nu = \cos(\theta/2)\f$ and \f$K_f^2\f$ the
     * five-point side symbol of both stencils (Deduction). The reflected amplitude at layer 0 is
     * \f$R_0 = -P(1/\kappa)/P(\kappa)\f$
     * [higdon1986AbsorbingDifferenceMultidim, section 3, Eq. (3.6), p. 440], so
     * \f$R_{1/2} = R_0\kappa = (X - Y)/(X + Y)\f$ at the Mur centre plane. sin and cos of
     * \f$\theta/2\f$ are evaluated from q without trigonometry.
     * @pre \f$0 < \omega\Delta t < \pi\f$, \f$q \ne 0, 1\f$, \f$D' \ne 0\f$.
     */
    inline DiscreteReflection discreteReflection(std::array<double, 3> h, double dt, double vu,
                                                 double vm, OracleStencil stencil, unsigned axis,
                                                 std::array<double, 3> kt, double omega) {
        const unsigned s1 = (axis + 1) % 3, s2 = (axis + 2) % 3;
        const double W = 2 / dt * std::sin(omega * dt / 2), ct = std::cos(omega * dt / 2);
        const double k1 = squaredSymbol(kt[s1], h[s1]), k2 = squaredSymbol(kt[s2], h[s2]);
        double dp = 1, rest = k1 + k2;
        if (stencil == OracleStencil::NonStandard) {
            const double a =
                0.25 * (1 + 0.02 / ((h[2] / h[0]) * (h[2] / h[0]) + (h[2] / h[1]) * (h[2] / h[1])));
            // On a face normal to x or y, the side axes are z and one transverse axis.
            const double kz = s1 == 2 ? k1 : k2, kOther = s1 == 2 ? k2 : k1;
            dp   = axis == 2 ? 1 - a * h[2] * h[2] * (k1 + k2) : 1 - a * h[2] * h[2] * kz;
            rest = axis == 2 ? k1 + k2 : dp * kOther + kz;
        }
        const double q  = (W * W / (vu * vu) - rest) * h[axis] * h[axis] / (4 * dp);
        const double rq = std::sqrt(std::fabs(q)), sign = dp > 0 ? 1.0 : -1.0;
        DiscreteReflection r;
        OracleComplex sn, cn;  // sin and cos of theta / 2, from q
        if (q < 0) {
            r.modeClass     = LatticeModeClass::Evanescent;
            sn              = {0, rq};
            cn              = std::sqrt(1 - q);
            r.reflectedRoot = {0, 2 * std::asinh(rq)};
        } else if (q > 1) {
            r.modeClass     = LatticeModeClass::GridBranch;
            sn              = rq;
            cn              = {0, -std::sqrt(q - 1)};
            r.reflectedRoot = {std::numbers::pi, 2 * std::acosh(rq)};
        } else {
            r.modeClass     = LatticeModeClass::Propagating;
            sn              = sign * rq;
            cn              = std::sqrt(1 - q);
            r.reflectedRoot = sign * 2 * std::asin(rq);
        }
        const OracleComplex x        = vm * W * (2 / h[axis]) * sn * ct;
        const OracleComplex y        = (W * W - vm * vm * (k1 + k2) / 2) * cn;
        const OracleComplex rootHalf = cn - OracleComplex(0, 1) * sn;  // exp(-i theta / 2)
        r.incidentRoot = -r.reflectedRoot, r.centre = (x - y) / (x + y);
        r.layerZero = r.centre * rootHalf * rootHalf;
        return r;
    }

    /** @brief Continuum coefficient \f$R_c\f$ of the condition at speed \f$v_m\f$ for interior
     * speed \f$v_u\f$, at \f$s = v_u|\mathbf k_t|/\omega\f$ (\f$s = \sin\theta\f$ below 1).
     * @details With \f$r = v_m/v_u\f$ and \f$c = \cos\theta = \sqrt{1 - s^2}\f$, or
     * \f$c = \mathrm i\sqrt{s^2 - 1}\f$ for the evanescent wave that decays into the box,
     * \f$R_c = (rc - 1 + r^2s^2/2)/(rc + 1 - r^2s^2/2)\f$ (Deduction: the condition applied to
     * the incident and reflected plane waves). At \f$r = 1\f$ this is
     * \f$-((1 - c)/(1 + c))^2\f$, Engquist and Majda's reflection with its sign
     * [engquist1977AbsorbingBoundaryWaves, section 1.B, p. 633]. Evaluated with
     * \f$1 - rc = (1 - r^2 + r^2s^2)/(1 + rc)\f$ and \f$1 - s^2 = -(s - 1)(s + 1)\f$ (s - 1 from
     * one fused multiply-add) to avoid cancellation near normal and grazing incidence.
     */
    inline OracleComplex continuumReflection(double vu, double vm, double ktAbs, double omega) {
        const double d = std::fma(vu, ktAbs, -omega) / omega;  // s - 1 with one rounding
        const double s = vu * ktAbs / omega, r = vm / vu, rs = r * r * s * s, c2 = -d * (2 + d);
        const OracleComplex c =
            c2 > 0 ? OracleComplex(std::sqrt(c2), 0) : OracleComplex(0, std::sqrt(-c2));
        const OracleComplex one = 1.0 + r * c, oneMinus = (1 - r * r + rs) / one;
        return (2 * (r * r - 1) - rs * oneMinus) / (2.0 * one) / (one - rs / 2);
    }

    /// Observed order \f$\log(e_{\mathrm{coarse}}/e_{\mathrm{fine}})/\log\rho\f$, ratio \f$\rho\f$.
    inline double observedOrder(double errorCoarse, double errorFine, double refinementRatio) {
        return std::log(errorCoarse / errorFine) / std::log(refinementRatio);
    }

    /** @brief Reflected-to-incident ratio of TE \f$E_y\f$ or TM0 averaged \f$B_y\f$ at sample j.
     * @details Face normal to x, \f$\mathbf k_t\f$ in the x-z plane; TE: only \f$A_y\f$; TM0:
     * \f$A_x = 0\f$. Each potential array reflects with \f$R_0\f$, so array value j of the
     * reflected wave is \f$R_0\kappa^{2j}\f$ times that of the incident one. \f$E_y\f$ uses the
     * time difference of \f$A_y\f$ at the same position: ratio \f$R_0\kappa^{2j}\f$. \f$B_y =
     * -\partial_xA_z\f$ changes sign with \f$k_x\f$: centred (R2) \f$-R_0\kappa^{2j}\f$, forward
     * (R2s; the edge of array value j lies half a cell inward of it at the lower face and outward
     * at the upper face) \f$-R_0\kappa^{2j + 1}\f$ and \f$-R_0\kappa^{2j - 1}\f$. The time label is
     * common to both waves. Deduction from the R2 reconstructions with
     * \f$\mathbf B = \nabla\times\mathbf A\f$ [fallahi2020mithra20fullwavesimulation, Eq. (3.8)].
     */
    inline OracleComplex fieldRowRatio(const DiscreteReflection& r, bool magneticTM0,
                                       bool staggered, bool lowerFace, int j) {
        double power = 2.0 * j;
        if (magneticTM0 && staggered) {
            power += lowerFace ? 1 : -1;
        }
        const OracleComplex ratio =
            r.layerZero * std::exp(OracleComplex(0, 1) * power * r.reflectedRoot);
        return magneticTM0 ? -ratio : ratio;
    }

    /// Standard lattice frequency of \f$|k_d|\f$: \f$\omega = (2/\Delta t)\arcsin(v\Delta t
    /// K/2)\f$.
    inline double latticeFrequency(std::array<double, 3> h, double dt, double v,
                                   std::array<double, 3> kAbs) {
        double k2 = 0;
        for (unsigned d = 0; d < 3; ++d) {
            k2 += squaredSymbol(kAbs[d], h[d]);
        }
        return 2 / dt * std::asin(v * dt * std::sqrt(k2) / 2);
    }

    /** @brief Amplitude of family wave \f$m\f$ at an edge (two open axes) or corner (three).
     * @details Bit d of m set: the wave moves into the box along open axis d (reflected there),
     * root \f$+|k_d|h_d\f$ in the inward layer index; clear: outgoing, \f$-|k_d|h_d\f$. Periodic
     * axes keep \f$+|k_d|\f$ for every wave. Amplitude \f$\prod_{d \in m}R_0^{(d)}\f$, with
     * \f$R_0^{(d)}\f$ the Standard face coefficient of axis d for the other two components and
     * the lattice frequency; the outgoing wave has amplitude 1. The sum solves the interior and
     * each selected face row exactly, so the edge and corner rows too
     * [bamberger1987SecondOrderCorner, section II, p. 05] (Deduction: linear superposition).
     * @pre Bits of m only on open axes; matched speed v.
     */
    inline OracleComplex familyAmplitude(std::array<double, 3> h, double dt, double v,
                                         std::array<double, 3> kAbs, std::array<bool, 3> open,
                                         unsigned m) {
        const double omega      = latticeFrequency(h, dt, v, kAbs);
        OracleComplex amplitude = 1;
        for (unsigned d = 0; d < 3; ++d) {
            if (open[d] && ((m >> d) & 1u)) {
                amplitude *=
                    discreteReflection(h, dt, v, v, OracleStencil::Standard, d, kAbs, omega)
                        .layerZero;
            }
        }
        return amplitude;
    }

    /** @brief Gauge residual \f$G^{\mathrm r}\f$ of the reflected wave, referred to x = 0.
     * @details \f$G = \varepsilon\delta_t\phi + \mathrm{div}_h\mathbf A\f$ at cell centres; face
     * normal to x at x = 0, centres at \f$(j + 1/2)h\f$. \f$G^{\mathrm i}\f$ and
     * \f$A_x^{\mathrm i}\f$ are the incident amplitudes referred to x = 0. Componentwise Mur gives
     * every centred component \f$\rho = R_0e^{-\mathrm i\theta}\f$; the staggered \f$A_x\f$ (at
     * \f$(j + 1)h\f$) \f$\rho e^{-\mathrm i\theta}\f$; the x-difference symbol is odd in
     * \f$k_x\f$. Hence \f$G^{\mathrm r} = \rho(G^{\mathrm i} + 2\mathrm ig_xA_x^{\mathrm
     * i}\beta)\f$,
     * \f$g_x = \sin\theta/h\f$, \f$\beta = 1\f$ collocated, \f$e^{-\mathrm i\theta/2}\f$ staggered
     * (Deduction; continuum form \f$2\mathrm ik_xRA_x\f$ with the Lorenz gauge of
     * [chew2014generalizedgauge, section 2, Eq. (11), p. 3]).
     */
    inline OracleComplex reflectedGaugeResidual(const DiscreteReflection& r, double hx,
                                                OracleComplex gaugeIncident,
                                                OracleComplex axIncident, bool staggered) {
        const OracleComplex i(0, 1), theta = r.reflectedRoot;
        const OracleComplex rho  = r.layerZero * std::exp(-i * theta);
        const OracleComplex beta = staggered ? std::exp(-i * theta / 2.0) : OracleComplex(1);
        return rho * (gaugeIncident + 2.0 * i * std::sin(theta) / hx * axIncident * beta);
    }
}  // namespace chdr::solver
