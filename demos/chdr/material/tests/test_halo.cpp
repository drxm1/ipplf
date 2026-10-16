/** @file test_halo.cpp
 * @brief Exact source-cell checks for faces, edges, corners and physical/inter-rank intersections.
 */
#include "OwnerTestSupport.h"

namespace chdr::material::test {
    /** @brief Require every owned and ghost cell to hold its source cell's value.
     * @details IndexMap assigns a distinct value to each domain cell, so copying
     * from a wrong neighbour changes the result. The expected value comes from
     * extendedValue() in OwnerOracles.h, independently of MaterialOwner's halo
     * implementation: it wraps indices at periodic boundaries and clamps them
     * to the nearest domain cell at physical boundaries. Inter-rank ghosts use
     * the global index of their owning cell, which remains inside the domain.
     *
     * Repeat with each axis selected for MPI decomposition. mismatches() visits
     * the complete allocation, including faces, edges and corners, and sums
     * bitwise differences across ranks; the assertion requires exactly zero.
     * On multiple ranks this also checks physical/inter-rank intersections.
     * @param periodic True for periodic global boundaries on all axes, false
     * for physical global boundaries on all axes; internal rank seams remain.
     */
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
    TEST(MaterialHalo, PhysicalSourceCells) {
        checkHalos(false);
    }

    /** @brief Check mixed wrap/copy edges with the independent source-cell oracle.
     * Every mixed topology and decomposition axis is included. Distinct integer
     * samples expose an accidental clamp at a periodic face or wrap at an open face.
     */
    TEST(MaterialHalo, MixedSourceCells) {
        for (unsigned mask = 1; mask < 7; ++mask) {
            const std::array<bool, 3> axes{bool(mask & 1), bool(mask & 2), bool(mask & 4)};
            const Kokkos::Array<bool, 3> expected{axes[0], axes[1], axes[2]};
            for (unsigned axis = 0; axis < 3; ++axis) {
                auto grid = makeGrid(axis, axes);
                const IndexMap map;
                const MaterialOwner owner(grid.getMesh(), grid.getFL(), map);
                EXPECT_EQ(mismatches(owner.epsilon(), grid.getFL(),
                                     [&](CellIndex global) {
                                         return extendedValue(map, global, Cells, expected);
                                     }),
                          0U);
            }
        }
    }
}  // namespace chdr::material::test
