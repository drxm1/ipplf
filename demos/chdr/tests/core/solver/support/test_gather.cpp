/** @file test_gather.cpp
 * @brief Global index order, exact ownership and NaN halo exclusion in rank records.
 * @details gatherOwned (RankRecords.h) collects the values of every rank record that
 * compare_ranks.py reads. OwnedIndexOrderAndGhostExclusion fills an 8 by 10 by 12
 * four-component field with its flattened global index on owned cells and NaN in the halo,
 * and splits it along x, y and z in turn. Rank 0 must receive exactly the values 0 to 3839
 * in x, y, z, component order, and every other rank an empty array.
 * RepeatedNaNIndexStillFailsOwnership requires that orderValues throws std::runtime_error
 * for the index list {0, 0}; it makes no MPI call. Both checks are exact. CTest runs the
 * suite RankGather as chdr.solver.gather.serial and chdr.solver.gather.mpi2.
 */
#include <gtest/gtest.h>

#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "tests/core/solver/support/RankRecords.h"

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

    /// For each split axis, rank 0 receives the owned values 0 to 3839 in global order and no
    /// halo value; other ranks receive nothing.
    TEST(RankGather, OwnedIndexOrderAndGhostExclusion) {
        for (unsigned split = 0; split < SpatialDim; ++split)
            checkGatherSplit(split);
    }
    /// A repeated global index throws std::runtime_error, also when its value is NaN.
    TEST(RankGather, RepeatedNaNIndexStillFailsOwnership) {
        const rankDetail::PackedValues repeated{{0, 0},
                                                {std::numeric_limits<double>::quiet_NaN(), 1.0}};
        EXPECT_THROW(rankDetail::orderValues(repeated, 2), std::runtime_error);
    }
}  // namespace chdr::solver::test
