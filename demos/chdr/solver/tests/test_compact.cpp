/** @file test_compact.cpp
 * @brief Natural-position staggered derivative symbols and component-offset metadata.
 * @details Signed half-cell symbols follow by Fourier substitution into
 * [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(20)] (deduction).
 *
 * TODO: REVIEW.
 */
#include <Kokkos_MathematicalConstants.hpp>
#include <gtest/gtest.h>
#include <limits>

#include "PeriodicSetup.h"
#include "StaggeredOperators.h"

namespace {
    using namespace chdr::solver;
    struct SampledMode {
        ippl::Vector<double, 3> waveNumber_m, spacing_m;
        KOKKOS_INLINE_FUNCTION double operator()(int i, int j, int k) const {
            return Kokkos::sin(waveNumber_m[0] * (i + 0.5) * spacing_m[0]
                               + waveNumber_m[1] * (j + 0.5) * spacing_m[1]
                               + waveNumber_m[2] * (k + 0.5) * spacing_m[2]);
        }
    };

    /// Compare both signed staggered derivative symbols at their respective half-cell positions.
    struct CompactSymbolError {
        SampledMode sample_m;
        ippl::Vector<double, 3> symbol_m;
        static constexpr double INFINITY_VALUE = std::numeric_limits<double>::infinity();
        KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& maximum) const {
            const auto waveNumber = sample_m.waveNumber_m, spacing = sample_m.spacing_m;
            const double phase = waveNumber[0] * (i + 0.5) * spacing[0]
                                 + waveNumber[1] * (j + 0.5) * spacing[1]
                                 + waveNumber[2] * (k + 0.5) * spacing[2];
            for (unsigned d = 0; d < 3; ++d) {
                const double forward = staggeredDifference(sample_m, i, j, k, d, spacing[d], true);
                const double backward =
                    staggeredDifference(sample_m, i, j, k, d, spacing[d], false);
                const double expectedForward =
                    symbol_m[d] * Kokkos::cos(phase + waveNumber[d] * spacing[d] / 2.0);
                const double expectedBackward =
                    symbol_m[d] * Kokkos::cos(phase - waveNumber[d] * spacing[d] / 2.0);
                const double forwardError =
                    Kokkos::abs(forward - expectedForward) / Kokkos::abs(symbol_m[d]);
                const double backwardError =
                    Kokkos::abs(backward - expectedBackward) / Kokkos::abs(symbol_m[d]);
                maximum = Kokkos::max(
                    maximum, Kokkos::isfinite(forwardError) ? forwardError : INFINITY_VALUE);
                maximum = Kokkos::max(
                    maximum, Kokkos::isfinite(backwardError) ? backwardError : INFINITY_VALUE);
            }
        }
    };

    TEST(CompactOperators, NaturalPositionSingleModeSymbols) {
        auto fields = makePeriodicFieldContainer({12, 10, 8}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        const auto spacing = fields.getMesh().getMeshSpacing();
        const auto waveNumber =
            2.0 * Kokkos::numbers::pi_v<double> * ippl::Vector<double, 3>(1, 2, -1);
        const SampledMode sample{waveNumber, spacing};
        ippl::Vector<double, 3> symbol;
        for (unsigned d = 0; d < 3; ++d)
            symbol[d] = 2.0 * std::sin(waveNumber[d] * spacing[d] / 2.0) / spacing[d];
        double error = 0.0;
        Kokkos::parallel_reduce("Compact single-mode symbols", fields.getJ().getFieldRangePolicy(),
                                CompactSymbolError{sample, symbol}, Kokkos::Max<double>(error));
        MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, fields.getFL().comm);
        EXPECT_LE(error, 1e-12);
    }

}  // namespace
