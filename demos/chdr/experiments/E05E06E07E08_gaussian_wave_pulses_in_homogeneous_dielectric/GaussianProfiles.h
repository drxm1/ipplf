/** @file GaussianProfiles.h
 * @brief The Gaussian potential profile of E05.typ for the plane and the spherical pulse.
 */
#pragma once

#include <Kokkos_Macros.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01AnalyticContinuumReference.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03AnalyticContinuumReference.h"

namespace chdr::experiments::E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric {
    namespace e01Reference =
        E01E02E03E04_wave_pulses_in_homogeneous_dielectric::E01_analytic_continuum_reference;
    namespace e03Reference =
        E01E02E03E04_wave_pulses_in_homogeneous_dielectric::E03_analytic_continuum_reference;

    /// The amplitude \f$A_0\f$ and the length unit \f$w\f$ of E01 to E04 stay unchanged.
    inline constexpr double A0 = e01Reference::A0;
    inline constexpr double w  = e01Reference::w;

    /** @brief Standard deviation \f$\sigma=w/4\f$ of the Gaussian.
     * @details Near its center the sharpened profile of E03.typ with \f$\beta=8\f$ is close to
     * \f$A_0\exp(-8s^2/w^2)\f$, which is this Gaussian. Where \f$|s|\ge2w=8\sigma\f$,
     * \f$p\f$, \f$\sigma|p'|\f$ and \f$\sigma^2|p''|\f$ are below \f$10^{-12}\f$ of their
     * maxima, so the geometries of E01 to E04 stay unchanged.
     */
    inline constexpr double Sigma = w / 4;

    /// \f$p(s)=A_0\exp(-s^2/(2\sigma^2))\f$ of E05.typ.
    KOKKOS_INLINE_FUNCTION double gaussian(double s) {
        return A0 * Kokkos::exp(-s * s / (2 * Sigma * Sigma));
    }

    /** @brief The Gaussian as a profile of the plane-pulse sampler of E01 and E02.
     * @details The sampler evaluates it at the displacement from the nearest periodic copy.
     * Every other copy is at least \f$L_z/2=22w=88\sigma\f$ away, where
     * \f$p<\exp(-3872)\f$ is below the smallest positive double, so the nearest copy alone
     * is exact in double precision (E05.typ).
     */
    struct GaussianPotentialProfile {
        /// Word in the run record's reference name, as in "E05 periodic Gaussian ...".
        static constexpr const char* name = "Gaussian";

        KOKKOS_INLINE_FUNCTION double value(double s) const { return gaussian(s); }

        /// \f$p'(s)=-s\,p(s)/\sigma^2\f$.
        KOKKOS_INLINE_FUNCTION double derivative(double s) const {
            return -s / (Sigma * Sigma) * gaussian(s);
        }

        /// Record the amplitude and the standard deviation (by reference; copies are detached).
        void record(conduit_cpp::Node& node) const {
            node["potential_amplitude"] = A0;
            node["standard_deviation"]  = Sigma;
        }
    };

    /** @brief The Gaussian as the radial profile of the dipole-pulse sampler of E03 and E04.
     * @details The generator is \f$f=w^2p\f$ with \f$f'=-s f/\sigma^2\f$ and
     * \f$f''=(s^2/\sigma^4-1/\sigma^2)f\f$. It has no support boundary, so `contains` is
     * always true. E07.typ bounds the inward-traveling term that this generator omits near
     * the center.
     */
    struct GaussianRadialProfile {
        KOKKOS_INLINE_FUNCTION bool contains(double) const { return true; }

        KOKKOS_INLINE_FUNCTION e03Reference::RadialGenerator generator(double s) const {
            const double variance = Sigma * Sigma;
            const double f        = w * w * gaussian(s);
            return {f, -s / variance * f, (s * s / (variance * variance) - 1 / variance) * f};
        }

        /// Record the profile name and its standard deviation (by reference).
        void record(conduit_cpp::Node& node) const {
            node["profile"]            = "Gaussian";
            node["standard_deviation"] = Sigma;
        }
    };

    static_assert(e01Reference::PlanePulseProfile<GaussianPotentialProfile>);
    static_assert(e03Reference::RadialPulseProfile<GaussianRadialProfile>);
}  // namespace chdr::experiments::E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric
