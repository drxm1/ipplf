/** @file test_comoving.cpp
 * @brief Compare reconstructed co-moving fields with independent analytic predictions.
 * @details Output rows tagged T7a describe one Fourier mode; T7b describes the
 * neutralized Gaussian train. P1 columns contain exact-symbol predictions and
 * P2 columns contain their leading-order expansions, defined by comovingSymbols().
 *
 * TODO: REVIEW.
 */
#include <algorithm>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>

#include "ComovingPotentials.h"
#include "MonitorTestSupport.h"
#include "RankRecords.h"
#include "SpatialDimension.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_comoving.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;
    using namespace chdr::solver::test;

    struct ErrorNorms {
        double l2_m, maximum_m;
    };
    struct ComovingComparison {
        ErrorNorms measured_m, exact_m, leading_m;
    };

    /** @brief Accumulate owned-cell errors without allowing a NaN to hide in a maximum.
     * @details The reference is the analytic E component at the reconstruction's
     * comparison time and position. Long-double local positive sums are converted
     * to plain doubles for MPI; halos never contribute to any metric.
     */
    struct MeasuredNorms {
        long double errorSum_m = 0.0L, exactSum_m = 0.0L;
        double maxima_m[2]{};

        void add(double actual, double reference) {
            const double error = actual - reference;
            if (!std::isfinite(actual)) {
                errorSum_m  = std::numeric_limits<long double>::infinity();
                maxima_m[0] = std::numeric_limits<double>::infinity();
            } else
                errorSum_m += static_cast<long double>(error) * error;
            exactSum_m += static_cast<long double>(reference) * reference;
            maxima_m[0] = std::max(maxima_m[0], std::abs(error));
            maxima_m[1] = std::max(maxima_m[1], std::abs(reference));
        }
        ErrorNorms global(MPI_Comm comm) {
            double sums[2]{double(errorSum_m), double(exactSum_m)};
            MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(MPI_IN_PLACE, maxima_m, 2, MPI_DOUBLE, MPI_MAX, comm);
            return {std::sqrt(sums[0] / sums[1]), maxima_m[0] / maxima_m[1]};
        }
    };

    /// Measure the actual solver kernel on owned cells; the reference arrays are independent.
    ErrorNorms measuredError(EMFieldsType& field, unsigned component,
                             const ComovingPotentials& wave, const std::vector<double>& exact) {
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, field.getView());
        const auto local = field.getLayout().getLocalNDIndex();
        const int ghost  = field.getNghost();
        MeasuredNorms norms;
        for (int i = local[0].first(); i <= local[0].last(); ++i)
            for (int j = local[1].first(); j <= local[1].last(); ++j)
                for (int k = local[2].first(); k <= local[2].last(); ++k)
                    norms.add(host(i - local[0].first() + ghost, j - local[1].first() + ghost,
                                   k - local[2].first() + ghost)[component],
                              exact[wave.index(i, j, k)]);
        return norms.global(field.getLayout().comm);
    }

    /** @brief Use Fourier-coefficient norms for relative L2 and sampled series for maxima.
     * @details Orthogonality of distinct integer box modes cancels the common
     * transform normalization in the ratio of squared coefficient norms. This is
     * the discrete Parseval identity, deduced from the DFT conventions in
     * [heinzel2002SpectrumSpectralDensity, Sec. 5, Eqs. (3)-(5)]. The maximum ratio
     * instead uses the independently evaluated error and reference series at the
     * actual sample positions. comovingSymbols() defines both prediction arrays;
     * no measured field enters this calculation.
     */
    ErrorNorms predictedError(const std::vector<Complex>& error,
                              const std::vector<Complex>& reference,
                              const std::vector<double>& errorGrid,
                              const std::vector<double>& referenceGrid) {
        long double errorNorm = 0.0L, referenceNorm = 0.0L;
        double errorMaximum = 0.0, referenceMaximum = 0.0;
        for (std::size_t i = 0; i < error.size(); ++i) {
            errorNorm += std::norm(error[i]);
            referenceNorm += std::norm(reference[i]);
            errorMaximum     = std::max(errorMaximum, std::abs(errorGrid[i]));
            referenceMaximum = std::max(referenceMaximum, std::abs(referenceGrid[i]));
        }
        return {std::sqrt(double(errorNorm / referenceNorm)), errorMaximum / referenceMaximum};
    }

    /// Preserve global-component traversal within each owned block for the bitwise R2 check.
    std::vector<double> ownedElectric(EMFieldsType& field) {
        const auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, field.getView());
        const auto local = field.getLayout().getLocalNDIndex();
        const int ghost  = field.getNghost();
        std::vector<double> values;
        for (int i = 0; i < int(local[0].length()); ++i)
            for (int j = 0; j < int(local[1].length()); ++j)
                for (int k = 0; k < int(local[2].length()); ++k)
                    for (unsigned c = 0; c < SpatialDim; ++c)
                        values.push_back(host(i + ghost, j + ghost, k + ghost)[c]);
        return values;
    }

}  // namespace
