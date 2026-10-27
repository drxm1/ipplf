/** @file test_statistics.cpp
 * @brief Independent small-array checks of distributed field statistics.
 */
#include "SerialReference.h"
#include "core/base/SpatialDimension.h"

static_assert(
    chdr::SpatialDim == 3,
    "test_statistics.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr;
    using namespace chdr::solver;
    using namespace chdr::solver::test;
    using StatisticPairs = Kokkos::View<double* [2], Kokkos::LayoutLeft>;

    /** @brief A short positive sum for which a sequential traversal loses small terms.
     * @details The exact first result is representable. Five leaves exercise an
     * odd tail at more than one level; the second channel checks every leaf is
     * retained. EXPECT_EQ deliberately requires exact equality, not an ULP window.
     */
    TEST(StableStatistics, OddTreeKeepsSmallTermsAndTheTail) {
        StatisticPairs pairs("odd statistic fixture", 5);
        auto host          = Kokkos::create_mirror_view(pairs);
        const double large = std::ldexp(1.0, 53);
        for (unsigned i = 0; i < 5; ++i) {
            host(i, 0) = i == 0 ? large : 1.0;
            host(i, 1) = i + 1;
        }
        Kokkos::deep_copy(pairs, host);
        const auto sums = statisticsDetail::pairwiseSums(pairs);
        EXPECT_EQ(sums[0], large + 4.0);
        EXPECT_EQ(sums[1], 15.0);
    }

    TEST(StableStatistics, EmptyAndSingleLeaf) {
        EXPECT_EQ(statisticsDetail::pairwiseSums(StatisticPairs("empty statistic fixture", 0)),
                  (std::array<double, 2>{0.0, 0.0}));
        StatisticPairs pairs("single statistic fixture", 1);
        auto host  = Kokkos::create_mirror_view(pairs);
        host(0, 0) = 3.0;
        host(0, 1) = 0.25;
        Kokkos::deep_copy(pairs, host);
        EXPECT_EQ(statisticsDetail::pairwiseSums(pairs), (std::array<double, 2>{3.0, 0.25}));
    }

    /// Plant signed integer multiples of each global cell number and leave every ghost NaN.
    template <class FieldType>
    void plantIntegerValues(FieldType& field) {
        field             = std::numeric_limits<double>::quiet_NaN();
        const auto view   = field.getView();
        const auto first  = field.getLayout().getLocalNDIndex().first();
        const auto global = field.getLayout().getDomain();
        const int ny = global[1].length(), nz = global[2].length(), ghost = field.getNghost();
        constexpr unsigned COMPONENTS = FieldType::value_type::dim;
        Kokkos::parallel_for(
            "plant exact statistic fixture", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                const int cell = ((i - ghost + first[0]) * ny + j - ghost + first[1]) * nz + k
                                 - ghost + first[2];
                for (unsigned c = 0; c < COMPONENTS; ++c) {
                    const double value = double(cell + 1) * (c + 1);
                    view(i, j, k)[c]   = c % 2 ? -value : value;
                }
            });
    }

    /// Closed integer totals check both consumers of the shared helper independently of its tree.
    template <class FieldType>
    void checkIntegerValues(FieldType& field, int cells) {
        constexpr unsigned COMPONENTS = FieldType::value_type::dim;
        const double c = COMPONENTS, n = cells;
        const std::array<double, 3> exact{
            c * n, c * (c + 1) / 2 * n * (n + 1) / 2,
            c * (c + 1) * (2 * c + 1) / 6 * n * (n + 1) * (2 * n + 1) / 6};
        plantIntegerValues(field);
        EXPECT_EQ(ownedStatistics(field).values(), exact);
        EXPECT_EQ(rankDetail::statistics(field), exact);
        EXPECT_EQ(serialSnapshot(field).reductions_m, exact);
    }

    TEST(StableStatistics, OddOwnedBlocksExcludeNaNGhosts) {
        const int ranks = ippl::Comm->size();
        auto fields     = makeBoxFieldContainer({3, 3, 3 * ranks}, {1.0, 1.0, 1.0}, 2);
        initializeFieldStorage(fields);
        checkIntegerValues(fields.getJ(), 27 * ranks);
        checkIntegerValues(fields.getE(), 27 * ranks);
    }

    TEST(StableStatistics, OwnedNaNCannotDisappearFromMaximum) {
        auto fields = makeBoxFieldContainer({4, 4, 4}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        fields.getJ()     = std::numeric_limits<double>::quiet_NaN();
        const auto result = ownedStatistics(fields.getJ());
        EXPECT_TRUE(std::isinf(result.maxAbs_m));
        EXPECT_FALSE(std::isfinite(result.sumAbs_m));
        EXPECT_FALSE(std::isfinite(result.sumSquares_m));
    }

}  // namespace
