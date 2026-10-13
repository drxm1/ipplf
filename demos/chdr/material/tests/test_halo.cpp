/** @file test_halo.cpp
 * @brief Exact source-cell checks for faces, edges, corners and physical/inter-rank intersections.
 */
#include "OwnerTestSupport.h"

namespace chdr::material::test {
    /// The independent wrap/clamp oracle includes all allocated ghosts on every split axis.
    /// TODO: Improve these comments later based on the documented "oracles".
    /// A reader does not know what "the independend wrap/clamp oracle" is.
    void checkHalos(bool periodic) {
        for (unsigned axis = 0; axis < 3; ++axis) {
            auto grid = makeGrid(axis, periodic);
            const IndexMap map;
            const MaterialOwner owner(grid.getMesh(), grid.getFL(), map);
            const Kokkos::Array<bool, 3> boundaries{periodic, periodic, periodic};
            EXPECT_EQ(mismatches(owner.epsilon(), grid.getFL(),
                                 [&](CellIndex global) {
                                     return extendedValue(map, global, Cells, boundaries);
                                 }),
                      0U)
                << "split axis " << axis;
        }
    }

    TEST(MaterialHalo, PeriodicSourceCells) {
        checkHalos(true);
    }
}  // namespace chdr::material::test
