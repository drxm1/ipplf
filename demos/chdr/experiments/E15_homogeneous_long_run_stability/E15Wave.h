/** @file E15Wave.h
 * @brief Diagonal periodic mode, discrete initial histories and continuum reference of E15.
 */
#pragma once
#include <Kokkos_Core.hpp>

#include <cmath>

namespace chdr::experiments::e15 {
    constexpr double Epsilon = 2.13;
    constexpr double TwoPi   = 2 * Kokkos::numbers::pi_v<double>;

    /** @brief Wave of E15.typ, "Analytic continuum reference" and "Expected behaviour".
     * @details Chew (2014), Eqs. (5)-(6), gives the fields of the potentials.
     * The leapfrog dispersion relation follows by the harmonic substitution of
     * Schneider (2010), section 7.4, on the three equal mesh axes.
     */
    struct Wave {
        int cells_m;
        double dt_m, omega_m, discreteOmega_m, kappaE_m, kappaB_m;

        explicit Wave(int cells)
            : cells_m(cells)
            , dt_m(0.5 / cells) {
            omega_m         = TwoPi * std::sqrt(3.0 / Epsilon);
            const double x  = 0.5 * TwoPi / cells;
            const double y  = std::asin(0.5 * std::sqrt(3.0 / Epsilon) * std::sin(x));
            discreteOmega_m = 2 * y / dt_m;
            kappaE_m        = std::sin(x) / x;
            kappaB_m        = kappaE_m * std::cos(y);
        }

        /// Amplitudes of the positive-frequency complex fields, excluding the factor i.
        KOKKOS_INLINE_FUNCTION double amplitude(bool electric, unsigned d) const {
            if (electric)
                return omega_m * (d == 0 ? 1 : d == 1 ? -1 : 0) / Kokkos::sqrt(2.0);
            return TwoPi * (d == 2 ? -2 : 1) / Kokkos::sqrt(2.0);
        }

        /// Core loadHistory supplies each slot's physical position and time from LevelContract.
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double time,
                                                 unsigned component) const {
            const double wave = Kokkos::cos(TwoPi * (x + y + z) - discreteOmega_m * time);
            if (component == 0)
                return wave / Kokkos::sqrt(Epsilon);
            return wave * (amplitude(true, component - 1) / omega_m + 1 / Kokkos::sqrt(3.0));
        }
    };
}  // namespace chdr::experiments::e15
