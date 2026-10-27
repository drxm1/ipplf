/** @file LorenzWave.h
 * @brief Source-free and driven continuum potential, source and field samples for Waves.h.
 *
 * An anylatic test generator that supplies potential samples and source sample.
 * Also provides corresponding E and B reference fields.
 */
#pragma once
#include "MediumOracles.h"
#include "PeriodicMode.h"
#include "core/base/SpatialDimension.h"

static_assert(chdr::SpatialDim == 3,
              "LorenzWave.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    /** @brief Analytic potentials, physical sources and fields for one continuum plane wave.
     * @details Defaults describe vacuum; makeMediumWave supplies the uniform relative
     * permittivity and caches the scalar amplitude and separate charge/current factors.
     * In vacuum, \f$\widehat\phi=(k\cdot\widehat A)/\Omega\f$ and
     * \f$s=(|k|^2-\Omega^2)\psi\f$ (deduction from
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.6)-(3.9)]).
     * The medium amplitudes instead obey [chew2014generalizedgauge, Eqs. (11)-(13)].
     */
    struct LorenzWave {
        ippl::Vector<double, SpatialDim> waveNumber_m, realA_m, imagA_m;
        double omega_m, sourceFactor_m;
        /// Uniform relative permittivity and host-computed medium scalar/source data.
        double epsilon_m = 1, realPhi_m = 0, imagPhi_m = 0, chargeFactor_m = 0;
        /** @brief Obtain one part of the scalar amplitude from the continuum Lorenz relation.
         * @details Medium values come from the independent host oracle. The vacuum branch
         * retains the existing arithmetic and also supports directly constructed waves.
         */
        KOKKOS_INLINE_FUNCTION double phiAmplitude(bool imaginary) const {
            if (epsilon_m != 1)
                return imaginary ? imagPhi_m : realPhi_m;
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
        /** @brief Sample physical charge (component zero) or current (components one to three).
         * @details The medium charge and current use different multipliers, as derived in
         * drivenSourceFactor from [chew2014generalizedgauge, Eqs. (12)-(13)]. These values
         * are input sources; the material update applies its own dielectric scaling.
         */
        KOKKOS_INLINE_FUNCTION double source(double x, double y, double z, double t,
                                             unsigned c) const {
            return (epsilon_m != 1 && c == 0 ? chargeFactor_m : sourceFactor_m)
                   * (*this)(x, y, z, t, c);
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

    /** @brief Prepare device samples from the independent medium-wave oracles.
     * @details The host supplies all medium-dependent amplitudes and source factors.
     * Device sampling retains the field definitions of LorenzWave; the gauge and
     * forcing follow [chew2014generalizedgauge, Eqs. (11)-(13)].
     * A free wave uses the medium frequency and gamma sets its longitudinal amplitude.
     * A driven wave uses a fixed polarization and two thirds of that frequency, so its
     * physical source is nonzero; gamma is unused in that branch. This detuning is a
     * test input, not an approximation to the source-free dispersion relation.
     */
    inline LorenzWave makeMediumWave(ippl::Vector<double, SpatialDim> k, double eps,
                                     double gamma = 1, bool sourced = false) {
        const double omega = mediumFrequency(k, eps) * (sourced ? 2.0 / 3 : 1);
        LorenzWave wave;
        if (sourced) {
            wave = makeSourcedWave(k, omega);
        } else {
            const auto amplitude = mediumWaveAmplitude(k, gamma);
            wave                 = LorenzWave{k, amplitude.real_m, amplitude.imag_m, omega, 0};
        }
        wave.epsilon_m      = eps;
        wave.realPhi_m      = mediumScalarAmplitude(k, wave.realA_m, eps, omega);
        wave.imagPhi_m      = mediumScalarAmplitude(k, wave.imagA_m, eps, omega);
        wave.sourceFactor_m = sourced ? drivenSourceFactor(1, k, eps, omega) : 0;
        wave.chargeFactor_m = sourced ? drivenSourceFactor(0, k, eps, omega) : 0;
        return wave;
    }
}  // namespace chdr::solver
