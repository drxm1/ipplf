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

        /** @brief Use this struct to sample the continuum analytic four-potential A reference of
         * experiment E01. You can also sample the corresponding E and B fields by explicitly
         * calling their functions.
         */
        struct E01AnalyticContinuumReferenceSampler {
            /// The length of the computational mesh domain
            const double Lz_m;

            /// The center z-coordinate of the initialized vector potential pulse
            const double z0_m;

            /// The wave speed inside the homogeneous medium, \f$v = 1/\sqrt \epsilon_r\f$
            const double v_m;

            /// @brief sample the four analytic continuum reference four potential
            KOKKOS_INLINE_FUNCTION double operator()([[maybe_unused]] const double _x,
                                                     [[maybe_unused]] const double _y,
                                                     const double z, const double t,
                                                     const unsigned component) const {
                return component == 1 ? /* A_x = */ Ax(z, t) : /* phi, Ay, Az = */ 0.;
            }

            KOKKOS_INLINE_FUNCTION double referenceFourPotentialA(const double x, const double y,
                                                                  const double z, const double t,
                                                                  const unsigned component) const {
                return this->operator()(x, y, z, t, component);
            }

            KOKKOS_INLINE_FUNCTION double referenceE([[maybe_unused]] const double _x,
                                                     [[maybe_unused]] const double _y,
                                                     const double z, const double t,
                                                     const unsigned component) const {
                return component == 0 /* x */ ? /* E_x = */ v_m * a_prime(s(z, t))
                                              : /* E_y, E_z = */ 0.;
            }

            KOKKOS_INLINE_FUNCTION double referenceB([[maybe_unused]] const double _x,
                                                     [[maybe_unused]] const double _y,
                                                     const double z, const double t,
                                                     const unsigned component) const {
                return component == 1 ? a_prime(s(z, t)) : 0.;
            }

            KOKKOS_INLINE_FUNCTION double s(const double z, const double t) const {
                return z - z0_m - v_m * t;
            }

            KOKKOS_INLINE_FUNCTION double r(const double s) const {
                return s - Lz_m * Kokkos::floor(s / Lz_m + 0.5);
            }

            /// @brief Periodic extension of \f$p(s)\f$
            KOKKOS_INLINE_FUNCTION double a(const double s) const { return p(r(s)); }

            /// @brief Derivative of \f$a(s)\f$ wrt \f$s\f$
            KOKKOS_INLINE_FUNCTION double a_prime(const double s) const { return p_prime(r(s)); }

            KOKKOS_INLINE_FUNCTION double Ax(const double z, const double t) const {
                return a(s(z, t));
            }
        };

    }  // namespace E01_analytic_continuum_reference

}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric
