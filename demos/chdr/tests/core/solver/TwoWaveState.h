/** @file TwoWaveState.h
 * @brief Exact lattice histories and a row-scaled two-root reflection fit.
 */
#pragma once
#include <vector>

#include "LevelContract.h"
#include "OpenBoxCases.h"

namespace chdr::solver::test::reflection {
    using Complex = OracleComplex;

    /** @brief Device-copyable real or imaginary part of an exact two-wave state.
     * @details The root and reflection coefficient come only from ReflectionOracles.
     * Layer zero is the first owned sample. TE uses only \f$A_y\f$; TM0 uses
     * \f$A_z=1\f$ and \f$\phi=S_z/(\varepsilon W)\f$, the discrete incident polarization from the
     * Lorenz gauge [chew2014generalizedgauge, section 2, equation (11)] specialized to the symbols
     * of LevelContract and ReflectionOracles (deduction). Retaining the scalar potential makes
     * the incident loaded state gauge-consistent, including when only its magnetic field is read.
     * Normal coordinates of these three
     * components coincide for both readings.
     */
    struct TwoWaveState {
        double hx_m, ky_m, kz_m, omega_m, thetaReal_m, thetaImag_m, rr_m, ri_m, phi_m;
        bool magnetic_m;
        int quadrature_m;
        TwoWaveState(const Case& row, SpatialDiscretization reading, bool tm, int part)
            : hx_m(row.h_m[0])
            , ky_m(row.tangent()[1])
            , kz_m(row.tangent()[2])
            , omega_m(row.omega_m)
            , magnetic_m(tm)
            , quadrature_m(part) {
            const auto r    = row.oracle();
            thetaReal_m     = r.reflectedRoot.real();
            thetaImag_m     = r.reflectedRoot.imag();
            rr_m            = r.layerZero.real();
            ri_m            = r.layerZero.imag();
            const double sz = reading == SpatialDiscretization::Staggered
                                  ? 2 * std::sin(kz_m * row.h_m[2] / 2) / row.h_m[2]
                                  : std::sin(kz_m * row.h_m[2]) / row.h_m[2];
            phi_m = sz / (row.epsilon_m * (2 / row.dt_m * std::sin(omega_m * row.dt_m / 2)));
        }

        /** @brief Sample the selected quadrature at the component's physical label. */
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double t,
                                                 unsigned q) const {
            const double amplitude =
                q == (magnetic_m ? 3u : 2u) ? 1 : (magnetic_m && q == 0 ? phi_m : 0);
            const double j = x / hx_m - .5, phase = ky_m * y + kz_m * z - omega_m * t;
            const double a = phase - thetaReal_m * j, b = phase + thetaReal_m * j;
            const double inc = Kokkos::exp(thetaImag_m * j), ref = Kokkos::exp(-thetaImag_m * j);
            return amplitude
                   * (quadrature_m == 0
                          ? inc * Kokkos::cos(a)
                                + ref * (rr_m * Kokkos::cos(b) - ri_m * Kokkos::sin(b))
                          : inc * Kokkos::sin(a)
                                + ref * (rr_m * Kokkos::sin(b) + ri_m * Kokkos::cos(b)));
        }
    };

    /** @brief Fit both exact lattice roots and return reflected/incident amplitude.
     * @details Rows are divided by their larger root magnitude. The two-column
     * least-squares model is the wave decomposition of Jones and Parrott (1988),
     * "Evaluation of a multi-point method for determining acoustic impedance",
     * Analysis, equations (1)-(3). Accumulation and the 2x2 Hermitian solve use
     * long double; the requested coefficient refers to the first input layer.
     */
    inline Complex fit(const std::vector<Complex>& values, Complex theta) {
        using Wide     = std::complex<long double>;
        long double aa = 0, bb = 0;
        Wide ab{}, ay{}, by{};
        for (std::size_t j = 0; j < values.size(); ++j) {
            const Wide phase =
                Wide(0, 1) * Wide(theta.real(), theta.imag()) * static_cast<long double>(j);
            Wide a = std::exp(-phase), b = std::exp(phase);
            const long double scale = std::max(std::abs(a), std::abs(b));
            a /= scale;
            b /= scale;
            const Wide y = Wide(values[j].real(), values[j].imag()) / scale;
            aa += std::norm(a);
            bb += std::norm(b);
            ab += std::conj(a) * b;
            ay += std::conj(a) * y;
            by += std::conj(b) * y;
        }
        const auto ratio = (aa * by - std::conj(ab) * ay) / (bb * ay - ab * by);
        return {double(ratio.real()), double(ratio.imag())};
    }
}  // namespace chdr::solver::test::reflection
