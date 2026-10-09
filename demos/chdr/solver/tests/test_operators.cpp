/** @file test_operators.cpp
 * @brief Independent Fourier oracles for the spatial stencil symbols.
 * @details Wave symbols follow [fallahi2020mithra20fullwavesimulation,
 * Eqs. (3.17), (3.24)]; the centered symbol follows
 * [trefethen1996FiniteDifferenceSpectral, Eqs. (5.1.8)-(5.1.9), p. 195].
 *
 * TODO: REVIEW.
 */
#include <gtest/gtest.h>
#include <limits>

#include "Operators.h"
#include "PeriodicSetup.h"

namespace {
    using namespace chdr::solver;

    /// Put cosine and sine in adjacent components so each derivative has an independent oracle.
    void loadFourierMode(AJFourFieldsType& potential, ippl::Vector<double, 3> waveNumber) {
        const auto h     = potential.get_mesh().getMeshSpacing();
        const auto local = potential.getLayout().getLocalNDIndex();
        const ippl::Vector<int, 3> first(local[0].first(), local[1].first(), local[2].first());
        const auto view = potential.getView();
        Kokkos::parallel_for(
            potential.getFieldRangePolicy(), KOKKOS_LAMBDA(int i, int j, int k) {
                const double phase = waveNumber[0] * (i + first[0] - 0.5) * h[0]
                                     + waveNumber[1] * (j + first[1] - 0.5) * h[1]
                                     + waveNumber[2] * (k + first[2] - 0.5) * h[2];
                view(i, j, k)[0] = Kokkos::cos(phase);
                view(i, j, k)[1] = Kokkos::sin(phase);
            });
        preparePotentialHalo(potential);
    }

    /// One owned-cell Fourier comparison; every nonfinite error is a failure, not a skipped
    /// maximum.
    template <Stencil Scheme>
    struct FourierError {
        AJFourFieldsType::view_type view_m;
        ippl::Vector<double, 3> h_m, g_m;
        double symbol_m, infinity_m;
        KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& maximum) const {
            const auto f = component(view_m, 0);
            const double lapError =
                Kokkos::abs(laplacian<Scheme>(f, i, j, k, h_m) - symbol_m * view_m(i, j, k)[0])
                / Kokkos::abs(symbol_m);
            maximum = Kokkos::max(maximum, Kokkos::isfinite(lapError) ? lapError : infinity_m);
            for (int d = 0; d < 3; ++d) {
                const double derivativeError = Kokkos::abs(centeredDifference(f, i, j, k, d, h_m[d])
                                                           + g_m[d] * view_m(i, j, k)[1])
                                               / Kokkos::abs(g_m[d]);
                maximum = Kokkos::max(
                    maximum, Kokkos::isfinite(derivativeError) ? derivativeError : infinity_m);
            }
        }
    };

    /** @brief Independently evaluate the trigonometric dispersion polynomial.
     * @details Direct substitution in [fallahi2020mithra20fullwavesimulation,
     * Eqs. (3.17), (3.24)-(3.25)] avoids production symbol helpers. Unequal transverse
     * spacings and distinct phase increments make axis/coefficient swaps visible.
     */
    double independentSymbol(ippl::Vector<double, 3> waveNumber, ippl::Vector<double, 3> h,
                             Stencil stencil) {
        double transverse = 0.0;
        for (unsigned d = 0; d < 2; ++d)
            transverse += 4.0 * std::pow(std::sin(0.5 * waveNumber[d] * h[d]), 2) / (h[d] * h[d]);
        const double z = std::pow(std::sin(0.5 * waveNumber[2] * h[2]), 2);
        const double alpha =
            stencil == Stencil::Standard
                ? 0.0
                : 1.0 + 0.02 / (std::pow(h[2] / h[0], 2) + std::pow(h[2] / h[1], 2));
        return -(transverse * (1.0 - alpha * z) + 4.0 * z / (h[2] * h[2]));
    }

    template <Stencil Scheme>
    double localFourierError(AJFourFieldsType& potential, ippl::Vector<double, 3> waveNumber) {
        const auto h        = potential.get_mesh().getMeshSpacing();
        const double symbol = independentSymbol(waveNumber, h, Scheme);
        const ippl::Vector<double, 3> g(std::sin(waveNumber[0] * h[0]) / h[0],
                                        std::sin(waveNumber[1] * h[1]) / h[1],
                                        std::sin(waveNumber[2] * h[2]) / h[2]);
        double error = 0;
        Kokkos::parallel_reduce(potential.getFieldRangePolicy(),
                                FourierError<Scheme>{potential.getView(), h, g, symbol,
                                                     std::numeric_limits<double>::infinity()},
                                Kokkos::Max<double>(error));
        return error;
    }

    template <Stencil Scheme>
    void checkFourierStencil() {
        auto fields = makePeriodicFieldContainer({12, 16, 24}, {1, 1, 1});
        initializeFieldStorage(fields);
        StandardSolverPeriodic storage(fields.getJ(), fields.getE(), fields.getB());
        const ippl::Vector<double, 3> waveNumber =
            2 * Kokkos::numbers::pi_v<double> * ippl::Vector<double, 3>(1, 2, -5);
        loadFourierMode(storage.A_n, waveNumber);
        double error = localFourierError<Scheme>(storage.A_n, waveNumber);
        MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, fields.getFL().comm);
        EXPECT_LE(error, 1e-12);
    }

    TEST(Operators, MatchFourierSymbols) {
        checkFourierStencil<Stencil::Standard>();
        checkFourierStencil<Stencil::NonStandard>();
    }

}  // namespace
