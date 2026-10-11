/** @file LorenzWave.h
 * @brief Source-free and driven continuum potential, source and field samples for Waves.h.
 *
 * An anylatic test generator that supplies potential samples and source sample.
 * Also provides corresponding E and B reference fields.
 *
 * TODO: REVIEW.
 */
#pragma once
#include "PeriodicMode.h"
#include "SpatialDimension.h"

static_assert(chdr::solver::SpatialDim == 3,
              "LorenzWave.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    /** @brief Real part of a complex Lorenz-gauge plane wave.
     * @details \f$\widehat\phi=(k\cdot\widehat A)/\Omega\f$ and
     * \f$s=(|k|^2-\Omega^2)\psi\f$. Deduction from
     * fallahi2020mithra20fullwavesimulation, Eqs. (3.6)-(3.9).
     */
    struct LorenzWave {
        ippl::Vector<double, SpatialDim> waveNumber_m, realA_m, imagA_m;
        double omega_m, sourceFactor_m;
        /// @brief Obtain one part of the scalar amplitude from the continuum Lorenz relation.
        KOKKOS_INLINE_FUNCTION double phiAmplitude(bool imaginary) const {
            const auto a = imaginary ? imagA_m : realA_m;
            return (waveNumber_m[0] * a[0] + waveNumber_m[1] * a[1] + waveNumber_m[2] * a[2])
                   / omega_m;
        }
        KOKKOS_INLINE_FUNCTION double phase(double x, double y, double z, double t) const {
            return waveNumber_m[0] * x + waveNumber_m[1] * y + waveNumber_m[2] * z - omega_m * t;
        }
        /** @brief Return the real plane-wave component without changing its supplied label.
         * @details The mode convention is \f$\exp(i(k\cdot x-\Omega t))\f$. Component zero
         * is phi; the other components are A. Device callers provide each component's
         * own spatial/time location through the common history/source sampling helpers.
         */
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double t,
                                                 unsigned c) const {
            const double real = c == 0 ? phiAmplitude(false) : realA_m[c - 1];
            const double imag = c == 0 ? phiAmplitude(true) : imagA_m[c - 1];
            return real * Kokkos::cos(phase(x, y, z, t)) - imag * Kokkos::sin(phase(x, y, z, t));
        }
        /// @brief Apply the continuum wave-operator multiplier, retaining IPPL's positive source
        /// sign.
        KOKKOS_INLINE_FUNCTION double source(double x, double y, double z, double t,
                                             unsigned c) const {
            return sourceFactor_m * (*this)(x, y, z, t, c);
        }
        /// @brief Copy a device-callable source functor; no reference to the host wave object
        /// escapes.
        auto sources() const {
            const auto wave = *this;
            return KOKKOS_LAMBDA(double x, double y, double z, double t, unsigned c)->double {
                return wave.source(x, y, z, t, c);
            };
        }
        /** @brief Continuum electric reference at the requested component location.
         * @details \f$\widehat E=i\Omega\widehat A-i k\widehat\phi\f$;
         * source: [fallahi2020mithra20fullwavesimulation, Eq. (3.9)].
         * The real/imaginary arithmetic below
         * evaluates that amplitude directly, without differencing sampled potentials.
         */
        KOKKOS_INLINE_FUNCTION double electric(double x, double y, double z, double t,
                                               unsigned c) const {
            const double real = omega_m * realA_m[c] - waveNumber_m[c] * phiAmplitude(false);
            const double imag = omega_m * imagA_m[c] - waveNumber_m[c] * phiAmplitude(true);
            return -imag * Kokkos::cos(phase(x, y, z, t)) - real * Kokkos::sin(phase(x, y, z, t));
        }
        /** @brief Continuum magnetic reference from the analytic curl amplitude.
         * @details \f$\widehat B=i k\times\widehat A\f$, from
         * [fallahi2020mithra20fullwavesimulation, Eq. (3.8)].
         * Cyclic component indices preserve the right-handed curl.
         */
        KOKKOS_INLINE_FUNCTION double magnetic(double x, double y, double z, double t,
                                               unsigned c) const {
            const unsigned a = (c + 1) % SpatialDim, b = (c + 2) % SpatialDim;
            const double real = waveNumber_m[a] * realA_m[b] - waveNumber_m[b] * realA_m[a];
            const double imag = waveNumber_m[a] * imagA_m[b] - waveNumber_m[b] * imagA_m[a];
            return -imag * Kokkos::cos(phase(x, y, z, t)) - real * Kokkos::sin(phase(x, y, z, t));
        }
        void requirePeriodic(std::array<double, SpatialDim> lengths) const {
            requirePeriodicMode(waveNumber_m, lengths);
        }
    };

    /** @brief Construct a sourced polarization with caller-selectable wave number and frequency.
     * @details The source multiplier is computed from those actual parameters.
     * Constructing the functor does not run the holdout; its only use is sampling.
     */
    inline LorenzWave makeSourcedWave(
        ippl::Vector<double, SpatialDim> waveNumber = 2 * Kokkos::numbers::pi_v<double>
                                                      * ippl::Vector<double, SpatialDim>(1, 2, 2),
        double omega = 4 * Kokkos::numbers::pi_v<double>) {
        const double norm2 = waveNumber[0] * waveNumber[0] + waveNumber[1] * waveNumber[1]
                             + waveNumber[2] * waveNumber[2];
        return {waveNumber, {0.6, 0.9, -0.3}, {0.2, -0.4, 0.5}, omega, norm2 - omega * omega};
    }

    /** @brief Construct a transverse polarization plus a unit gauge component.
     * @details \f$\widehat A=e_1+0.5i e_2+ik/|k|\f$, with a right-handed transverse
     * basis. Source: chew2014generalizedgauge, section 7, Eqs. (76)-(80); specific
     * polarization is a manufactured test input, not a claim from that source.
     */
    inline LorenzWave makeMixedWave() {
        const auto k =
            2 * Kokkos::numbers::pi_v<double> * ippl::Vector<double, SpatialDim>(1, 2, 2);
        const double norm       = std::sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);
        const double transverse = std::sqrt(k[1] * k[1] + k[2] * k[2]);
        const ippl::Vector<double, SpatialDim> e1(0, k[2] / transverse, -k[1] / transverse);
        const ippl::Vector<double, SpatialDim> e2((k[1] * e1[2] - k[2] * e1[1]) / norm,
                                                  -k[0] * e1[2] / norm, k[0] * e1[1] / norm);
        const ippl::Vector<double, SpatialDim> imaginary = 0.5 * e2 + k / norm;
        return {k, e1, imaginary, norm, 0.0};
    }
}  // namespace chdr::solver
