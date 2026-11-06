#pragma once

#include <Kokkos_Macros.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric {

    namespace E01_analytic_continuum_reference {
        constexpr double A0 = 1.0;
        constexpr double w  = 1.0;

        /** @brief pulse profile defined for the E01 reference potential.
         * @param s \f$:= z - z_0 - v t\f$, signed displacement from pulse center.
         * @returns \f$p(s)\f$, zero for \f$|s| \ge w\f$.
         */
        KOKKOS_INLINE_FUNCTION double p(const double s) {
            return Kokkos::abs(s) >= w ? 0. : A0 * Kokkos::exp(1 - 1 / (1 - Kokkos::pow(s / w, 2)));
        }

        /** @brief Derivative of the pulse profile defined for the E01 reference potential.
         * @param s := z - z_0 - v t
         * @returns \f$p'(s)\f$, the derivative of p wrt s
         */
        KOKKOS_INLINE_FUNCTION double p_prime(const double s) {
            return Kokkos::abs(s) >= w
                       ? 0.
                       : -2 * s * p(s) / (Kokkos::pow(w * (1 - Kokkos::pow(s / w, 2)), 2));
        }

        KOKKOS_INLINE_FUNCTION double r(const double s, const double L_z) {
            return s - L_z * Kokkos::floor(s / L_z + 0.5);
        }

        /// @brief Periodic extension of \f$p(s)\f$
        KOKKOS_INLINE_FUNCTION double a(const double s, const double L_z) {
            return p(r(s, L_z));
        }

        /// @brief Derivative of \f$a(s)\f$ wrt \f$s\f$
        KOKKOS_INLINE_FUNCTION double a_prime(const double s, const double L_z) {
            return p_prime(r(s, L_z));
        }

    }  // namespace E01_analytic_continuum_reference

}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric
