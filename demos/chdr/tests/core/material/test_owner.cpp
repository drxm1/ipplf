/** @file test_owner.cpp
 * @brief Check rasterization, exact global statistics and unrestricted material storage.
 */
#include "OwnerTestSupport.h"

namespace chdr::material::test {
    /// Reproduce the specified occupied-ID fingerprint in an independent host scan.
    std::uint64_t referenceHash(std::uint64_t value) {
        value += 0x9e3779b97f4a7c15ULL;
        value ^= value >> 30;
        value *= 0xbf58476d1ce4e5b9ULL;
        value ^= value >> 27;
        value *= 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }

    /// A coordinate-dependent slab exercises the centre expression and geometry adapter.
    struct Slab {
        KOKKOS_INLINE_FUNCTION bool contains(double, double y, double) const {
            return y >= 1.5 && y < 3.0;
        }
    };

    /** @brief Compare the map and its statistics against a serial global host scan.
     * @details Each one-rank and two-rank case compares the same complete reference,
     * so integer equality across decompositions follows without gathering a field.
     * The allocation comparison includes ghosts and uses the separate sourceCell()
     * oracle. The statistics scan visits global owned cells exactly once, excluding
     * ghosts. Both maps supplied here contain only finite values.
     */
    template <class Map>
    void checkMap(Map map, unsigned axis) {
        auto grid = makeGrid(axis, false);
        const MaterialOwner owner(grid.getMesh(), grid.getFL(), map);
        const auto sample = [&](CellIndex global, bool& marked) {
            Position centre;
            for (unsigned d = 0; d < 3; ++d)
                centre[d] = grid.getMesh().getOrigin()[d]
                            + (global[d] + 0.5) * grid.getMesh().getMeshSpacing()[d];
            return map(global, centre, marked);
        };
        EXPECT_EQ(mismatches(owner.epsilon(), grid.getFL(),
                             [&](CellIndex global) {
                                 bool marked = false;
                                 return sample(sourceCell(global, Cells, {false, false, false}),
                                               marked);
                             }),
                  0U);
        MapStatistics expected{std::numeric_limits<double>::max(),
                               std::numeric_limits<double>::lowest()};
        for (int i = 0; i < Cells[0]; ++i)
            for (int j = 0; j < Cells[1]; ++j)
                for (int k = 0; k < Cells[2]; ++k) {
                    bool marked        = false;
                    const double value = sample(CellIndex(i, j, k), marked);
                    const auto hash    = referenceHash((i * 14 + j) * 16 + k);
                    expected.marked_m += marked;
                    if (marked)
                        expected.markedChecksum_m += hash;
                    expected.mapChecksum_m +=
                        referenceHash(hash ^ Kokkos::bit_cast<std::uint64_t>(value));
                    expected.minimum_m = std::min(expected.minimum_m, value);
                    expected.maximum_m = std::max(expected.maximum_m, value);
                }
        const auto& actual = owner.statistics();
        EXPECT_EQ(actual.marked_m, expected.marked_m);
        EXPECT_EQ(actual.markedChecksum_m, expected.markedChecksum_m);
        EXPECT_EQ(actual.mapChecksum_m, expected.mapChecksum_m);
        EXPECT_EQ(actual.minimum_m, expected.minimum_m);
        EXPECT_EQ(actual.maximum_m, expected.maximum_m);
        EXPECT_GT(actual.marked_m, 0U);
        EXPECT_LT(actual.marked_m, 12U * 14U * 16U);
        EXPECT_TRUE(actual.allFinite_m);
        EXPECT_FALSE(actual.homogeneous());
    }

    TEST(MaterialOwner, RasterizationAndRankIndependentStatistics) {
        for (unsigned axis = 0; axis < 3; ++axis) {
            checkMap(IndexMap{1}, axis);
            checkMap(TwoMaterialMap<Slab>{Slab{}, 2.13, 1.0}, axis);
        }
    }

    // Storage accepts values independently of a solver's material restrictions:
    // 0.5 is valid stored data, while NaN is retained but cannot be homogeneous.
    TEST(MaterialOwner, HomogeneousAndNonfiniteValues) {
        auto grid = makeGrid(0, false);
        for (double value : {0.5, 1.0, 2.13, std::numeric_limits<double>::quiet_NaN()}) {
            const auto map = KOKKOS_LAMBDA(CellIndex, Position, bool& marked) {
                marked = false;
                return value;
            };
            const MaterialOwner owner(grid.getMesh(), grid.getFL(), map);
            EXPECT_EQ(owner.statistics().allFinite_m, std::isfinite(value));
            if (std::isfinite(value))
                EXPECT_EQ(owner.homogeneousValue(), value);
            else
                EXPECT_THROW(owner.homogeneousValue(), std::invalid_argument);
        }
    }
}  // namespace chdr::material::test
