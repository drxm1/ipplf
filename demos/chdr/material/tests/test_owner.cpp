/** @file test_owner.cpp
 * @brief Check rasterization, exact global statistics and unrestricted material storage.
 */
#include "OwnerTestSupport.h"

namespace chdr::material::test {
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
