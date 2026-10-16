/** @file OuterLayerOracles.h
 * @brief Expected values of the outer-layer tests: the closed ghost value beyond an open face and
 * the analytic fields of the polynomial potentials at their sample positions and times.
 * @details Normalized units \f$c = \varepsilon_0 = \mu_0 = 1\f$. Component 0 is \f$\phi\f$,
 * components 1-3 are \f$A_x, A_y, A_z\f$; \f$\mathbf E = -\partial_t\mathbf A - \nabla\phi\f$,
 * \f$\mathbf B = \nabla\times\mathbf A\f$ [fallahi2020mithra20fullwavesimulation, section 3.1.1,
 * Eqs. (3.6)-(3.8), p. 12]. Owned cell \f$i_d = 0 \ldots n_d - 1\f$ has its centre at
 * \f$(i_d + 1/2)h_d\f$ (global origin 0). Readings: R1 is the inherited reconstruction; R2 and its
 * staggered form R2s sample \f$\mathbf E\f$ and the averaged \f$\mathbf B\f$ at
 * \f$t_{k-1/2}\f$ and the raw \f$\mathbf B\f$ at \f$t_k\f$. In R2s, \f$A_d\f$ and \f$E_d\f$ lie on
 * the positive \f$d\f$ face and \f$B_d\f$ on the edge displaced half a cell along the other two
 * axes [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(24)].
 */
#pragma once
#include <Kokkos_Core.hpp>

#include <array>
#include <cmath>

namespace open_boundary_oracle {
    /** @brief Ghost value beyond an open face from the three outermost owned samples.
     * @details \f$g = 3u_0 - 3u_1 + u_2\f$ is the unique value with a vanishing third difference
     * \f$g - 3u_0 + 3u_1 - u_2 = 0\f$, so it equals \f$p(x_0 - h)\f$ for every polynomial \f$p\f$
     * of degree at most two in the normal coordinate (Deduction: the third difference of a
     * quadratic vanishes). The centred difference at layer 0 then reads
     * \f$(u_1 - g)/(2h) = (-3u_0 + 4u_1 - u_2)/(2h)\f$, the second-order one-sided first
     * derivative [fornberg1988FiniteDifferenceFormulas, Table 3, p. 704].
     */
    KOKKOS_INLINE_FUNCTION double closedGhostValue(double u0, double u1, double u2) {
        return 3 * u0 - 3 * u1 + u2;
    }

    /// Per component: \f$P_x\f$, \f$P_y\f$, \f$P_z\f$ (ascending powers) and the time slope a.
    using PolynomialCoefficients = std::array<std::array<double, 10>, 4>;

    /** @brief Forty frozen coefficients of the polynomial potentials
     * \f$\psi_q = P_x^q(x)P_y^q(y)P_z^q(z)(1 + a_qt)\f$, degree two in space, one in time.
     * @details Standard-normal draws (NumPy default_rng seed 7, in the order given), stored to
     * round-trip precision. With periodicY, \f$P_y^q = 1\f$, so the potentials are periodic in y.
     */
    inline PolynomialCoefficients polynomialFixture(bool periodicY) {
        PolynomialCoefficients c = {{
            {0.0012301533574825742, 0.2987455375084699, -0.2741378553622176, -0.8905918387572742,
             -0.45467078517172255, -0.9916465549964624, 0.060143602597438485, 1.3402152455545335,
             -0.49220651855132963, -0.6204748998199404},
            {0.4898420501851982, 0.35688700816006075, 0.10541424899789856, -0.9304680447082047,
             -0.02925182246327349, 0.6953031944582878, -1.344214547285082, -0.45761576104021817,
             -1.901222739800844, -1.289537739784976},
            {-1.8417350377917323, -0.23509113107468127, -1.2674464814437032, 0.2712643588217015,
             0.15675108662422516, -0.18693094462995438, -2.516759710820513, -0.5386928958466366,
             -0.048500945401071985, 0.11330898600330756},
            {-1.5301357655053935, -0.47775327603393064, -0.9785190780566395, -0.8088372394255993,
             1.0608986233860787, -0.8075346753318965, -0.0325217049455206, 0.8843898673831739,
             -0.583600432743302, -0.11170194958415963},
        }};
        for (auto& q : c) {
            if (periodicY) {
                q[3] = 1, q[4] = 0, q[5] = 0;
            }
        }
        return c;
    }

    /** @brief \f$\partial^m\psi_q/\partial x_{\mathrm{axis}}^m\f$ at x and t, m = 1 when axis is
     * 0-2 and m = 0 when axis is -1; timeDerivative gives \f$\partial_t\f$ instead of the value.
     */
    inline double polynomialPotential(const PolynomialCoefficients& c, unsigned q,
                                      std::array<double, 3> x, double t, int axis,
                                      bool timeDerivative) {
        double value = timeDerivative ? c[q][9] : 1 + c[q][9] * t;
        for (int d = 0; d < 3; ++d) {
            const double* p = &c[q][3 * d];
            value *= d == axis ? p[1] + 2 * p[2] * x[d] : p[0] + x[d] * (p[1] + x[d] * p[2]);
        }
        return value;
    }

    /** @brief Analytic \f$E_d = -\partial_tA_d - \partial_d\phi\f$ or
     * \f$B_d = \partial_aA_b - \partial_bA_a\f$, \f$(a, b) = (d + 1, d + 2) \bmod 3\f$, at x, t.
     * @details Time label t: R1 uses 0 for E and B (its E combines \f$\partial_t\mathbf A\f$ at
     * \f$t_{-1/2}\f$ with \f$\phi\f$ at \f$t_0\f$; exact because \f$\partial_t\mathbf A\f$ does not
     * depend on t); R2 and R2s use \f$-\Delta t/2\f$ for E and averaged B (the mean of two linear
     * samples is the midpoint value) and 0 for raw B. Histories: \f$\mathbf A\f$ at
     * \f$0, -\Delta t\f$; \f$\phi\f$ at \f$0, -\Delta t\f$ (R1) or \f$\Delta t/2, -\Delta t/2\f$
     * (R2, R2s).
     */
    inline double polynomialField(const PolynomialCoefficients& c, bool magnetic, unsigned d,
                                  std::array<double, 3> x, double t) {
        if (!magnetic) {
            return -polynomialPotential(c, d + 1, x, t, -1, true)
                   - polynomialPotential(c, 0, x, t, int(d), false);
        }
        const unsigned a = (d + 1) % 3, b = (d + 2) % 3;
        return polynomialPotential(c, b + 1, x, t, int(a), false)
               - polynomialPotential(c, a + 1, x, t, int(b), false);
    }

    /// Sample position of field component d of owned cell i (offsets of the file description).
    inline std::array<double, 3> fieldPosition(std::array<int, 3> i, std::array<double, 3> h,
                                               bool staggered, bool magnetic, unsigned d) {
        std::array<double, 3> x;
        for (unsigned e = 0; e < 3; ++e) {
            const bool shifted = staggered && (magnetic ? e != d : e == d);
            x[e]               = (i[e] + (shifted ? 1.0 : 0.5)) * h[e];
        }
        return x;
    }

    /// Positive normalization: maximum of \f$|E_d|\f$ or \f$|B_d|\f$ over owned cells and d.
    inline double polynomialFieldMaximum(const PolynomialCoefficients& c, std::array<int, 3> n,
                                         std::array<double, 3> h, bool staggered, bool magnetic,
                                         double t) {
        double maximum = 0;
        for (int i = 0; i < n[0]; ++i) {
            for (int j = 0; j < n[1]; ++j) {
                for (int k = 0; k < n[2]; ++k) {
                    for (unsigned d = 0; d < 3; ++d) {
                        const auto x = fieldPosition({i, j, k}, h, staggered, magnetic, d);
                        maximum =
                            std::fmax(maximum, std::fabs(polynomialField(c, magnetic, d, x, t)));
                    }
                }
            }
        }
        return maximum;
    }
}  // namespace open_boundary_oracle
