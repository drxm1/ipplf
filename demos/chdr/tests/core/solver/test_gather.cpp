/** @file test_gather.cpp
 * @brief Global index order, exact ownership and NaN halo exclusion in rank dumps.
 */
#include <gtest/gtest.h>

#include "RankRecords.h"
#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"

static_assert(chdr::SpatialDim == 3,
              "test_gather.cpp requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    /// Fill owned values with their expected flattened index; leave every halo as NaN.
    void fillGatherPattern(AJFourFieldsType& field) {
        field            = std::numeric_limits<double>::quiet_NaN();
        const auto view  = field.getView();
        const auto first = field.getLayout().getLocalNDIndex().first();
        Kokkos::parallel_for(
            "gather ownership pattern", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                const int cell =
                    ((i + first[0] - 1) * 10 + j + first[1] - 1) * 12 + k + first[2] - 1;
                for (unsigned c = 0; c < 4; ++c)
                    view(i, j, k)[c] = 4 * cell + c;
            });
    }

    /** @brief Compute this rank's verdict without reporting a test failure yet.
     * @details Only rank zero owns the gathered array; all other arrays must be empty.
     * Return an integer flag so every rank can participate in the subsequent verdict
     * reduction before any GoogleTest assertion is evaluated.
     */
    int localGatherMismatch(const std::vector<double>& values, int rank) {
        if (rank != 0)
            return !values.empty();
        int mismatch = values.size() != 8u * 10 * 12 * 4;
        for (std::size_t n = 0; n < values.size(); ++n)
            mismatch |= values[n] != double(n);
        return mismatch;
    }

    void checkGatherSplit(unsigned split) {
        auto fields = makeBoxFieldContainer({8, 10, 12}, {1, 1, 1}, split);
        initializeFieldStorage(fields);
        auto& field = fields.getJ();
        fillGatherPattern(field);
        const auto values = gatherOwned(field);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int mismatch = localGatherMismatch(values, rank);
        // Every rank finishes this split before any test failure is reported.
        MPI_Allreduce(MPI_IN_PLACE, &mismatch, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        EXPECT_EQ(mismatch, 0);
    }

    TEST(RankGather, OwnedIndexOrderAndGhostExclusion) {
        for (unsigned split = 0; split < SpatialDim; ++split)
            checkGatherSplit(split);
    }
    TEST(RankGather, RepeatedNaNIndexStillFailsOwnership) {
        const rankDetail::PackedValues repeated{{0, 0},
                                                {std::numeric_limits<double>::quiet_NaN(), 1.0}};
        EXPECT_THROW(rankDetail::orderValues(repeated, 2), std::runtime_error);
    }
}  // namespace chdr::solver::test
