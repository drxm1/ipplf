/** @file RankWitness.h
 * @brief Compare globally ordered serial and distributed samples with explicit scales.
 */
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>  // IWYU pragma: keep
#include <cstring>
#include <gtest/gtest.h>
#include <iostream>
#include <limits>
#include <mpi.h>

namespace chdr::solver::test {
    /// Maximum relative difference and binary equality flag, transported together across ranks.
    using RankAgreement = std::array<double, 2>;

    /** @brief Compare real or complex binary64 samples without discarding signed zeros.
     * @details A positive supplied scale is the serial run maximum or an uncancelled
     * reduction-input scale. Zero selects the maximum magnitude in the reference array.
     * Arrays contain owned samples only, in the same global order. Nonfinite differences
     * fail the comparison; binary equality is recorded separately from numerical tolerance.
     * Only the rank holding the serial reference calls this helper.
     */
    template <class Range>
    RankAgreement rankDifference(const Range& observed, const Range& reference, double scale = 0) {
        EXPECT_EQ(observed.size(), reference.size());
        if (observed.size() != reference.size())
            return {std::numeric_limits<double>::infinity(), 0};
        double error = 0, maximum = 0;
        for (std::size_t i = 0; i < reference.size(); ++i) {
            const double difference = std::abs(observed[i] - reference[i]);
            error                   = std::max(error, std::isfinite(difference)
                                                          ? difference
                                                          : std::numeric_limits<double>::infinity());
            maximum                 = std::max(maximum, std::abs(reference[i]));
        }
        if (scale == 0)
            scale = maximum;
        EXPECT_GT(scale, 0);
        EXPECT_TRUE(std::isfinite(scale));
        return {error / scale, std::memcmp(observed.data(), reference.data(),
                                           reference.size() * sizeof(typename Range::value_type))
                                       == 0
                                   ? 1.
                                   : 0.};
    }

    /** @brief Broadcast rank zero's comparison and enforce its scale-aware tolerance.
     * @details Cell comparisons use 1e-13, above binary64 round-off at the serial amplitude
     * scale; Fourier reductions pass their separately justified tolerance and input scale.
     * Every rank participates, including ranks that do not retain a serial reference.
     */
    inline void checkRankAgreement(RankAgreement report, MPI_Comm comm, const char* channel,
                                   double tolerance = 1e-13) {
        MPI_Bcast(report.data(), 2, MPI_DOUBLE, 0, comm);
        EXPECT_LE(report[0], tolerance);
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0)
            std::cout << channel << " rank_error=" << report[0] << " bitwise=" << report[1] << '\n';
    }
}  // namespace chdr::solver::test
