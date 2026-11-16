#pragma once

#include <Kokkos_Macros.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <catalyst_conduit.hpp>
#include <concepts>
#include <string>

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

        /** @brief Potential profile \f$p\f$ of one copy of the plane pulse.
         * @details It supplies \f$p(s)\f$ and \f$p'(s)\f$ at the displacement \f$s\f$ from the
         * pulse center, a static word `name` for the run record ("compact", "Gaussian") and
         * `record(node)`, which writes its parameters to the run record.
         */
        template <class T>
        concept PlanePulseProfile = requires(const T& profile, double s, conduit_cpp::Node& node) {
            { profile.value(s) } -> std::same_as<double>;
            { profile.derivative(s) } -> std::same_as<double>;
            { T::name } -> std::convertible_to<const char*>;
            profile.record(node);
        };

        /** @brief The compact potential profile of E01.typ as a profile of the plane-pulse sampler.
         * @details A profile supplies \f$p(s)\f$ and \f$p'(s)\f$ of one pulse copy and the
         * run-record entries of its parameters. The plane-pulse sampler adds the travelling
         * coordinate and the periodic extension, so another profile (for example the Gaussian
         * of E05.typ) reuses the whole E01/E02 run path without changing it.
         */
        struct E01CompactProfile {
            /// Word in the run record's reference name, as in "E01 periodic compact ...".
            static constexpr const char* name = "compact";

            KOKKOS_INLINE_FUNCTION double value(double s) const { return p(s); }

            KOKKOS_INLINE_FUNCTION double derivative(double s) const { return p_prime(s); }

            /** @brief Record the amplitude and the support half-width of E01.typ.
             * @details The node is taken by reference: a copied conduit node is a separate
             * tree, and entries written to it would not reach the run record.
             */
            void record(conduit_cpp::Node& node) const {
                node["potential_amplitude"] = A0;
                node["support_half_width"]  = w;
            }
        };

        /** @brief Use this struct to sample the continuum analytic four-potential A reference of
         * experiment E01. You can also sample the corresponding E and B fields by explicitly
         * calling their functions.
         * @tparam Profile Potential profile of one pulse copy, for example E01CompactProfile.
         */
        template <PlanePulseProfile Profile = E01CompactProfile>
        struct E01AnalyticContinuumReferenceSampler {
            /// The length of the computational mesh domain
            const double Lz_m;

            /// The center z-coordinate of the initialized vector potential pulse
            const double z0_m;

            /// The wave speed inside the homogeneous medium, \f$v = 1/\sqrt \epsilon_r\f$
            const double v_m;

            /// Disable periodic wrapping when reusing this reference for E02.
            const bool periodicZ_m = true;

            /// Potential profile \f$p\f$ of one pulse copy.
            const Profile profile_m{};

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
                return periodicZ_m ? s - Lz_m * Kokkos::floor(s / Lz_m + 0.5) : s;
            }

            /// @brief Periodic extension of \f$p(s)\f$
            KOKKOS_INLINE_FUNCTION double a(const double s) const { return profile_m.value(r(s)); }

            /// @brief Derivative of \f$a(s)\f$ wrt \f$s\f$
            KOKKOS_INLINE_FUNCTION double a_prime(const double s) const {
                return profile_m.derivative(r(s));
            }

            KOKKOS_INLINE_FUNCTION double Ax(const double z, const double t) const {
                return a(s(z, t));
            }

            /** @brief Write the reference's name, geometry, speed and profile to the run record.
             * @param node The record's `analytic_reference` node. Pass the node that
             * `operator[]` returns: a copy of a conduit node is a separate tree.
             * @param experiment The experiment name, for example "E01".
             */
            void record(conduit_cpp::Node node, const std::string& experiment) const {
                node["name"] = experiment + (periodicZ_m ? " periodic " : " unwrapped ")
                               + Profile::name + " continuum pulse";
                node["periodic_z"]      = periodicZ_m;
                node["domain_length_z"] = Lz_m;
                if (periodicZ_m)
                    node["period"] = Lz_m;
                node["initial_centre_z"] = z0_m;
                node["wave_speed"]       = v_m;
                profile_m.record(node);
            }
        };

    }  // namespace E01_analytic_continuum_reference

}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric
