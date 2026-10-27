/** @file test_derived.cpp
 * @brief Compare every cached coefficient cell with its independently authored closed form.
 */
#include <gtest/gtest.h>

#include "OwnerOracles.h"
#include "OwnerTestSupport.h"
#include "core/material/O1VSCoefficients.h"

namespace chdr::material::test {
    /** @brief Compare a face field over its complete allocation with the independent oracle.
     * @details faceValue() resolves the two source cells separately before applying
     * its coefficient formula. It therefore distinguishes the lower boundary pair
     * from a copied first-interior face. mismatches() checks bits, including ghosts.
     */
    template <unsigned Axis>
    void checkFace(const MaterialOwner& owner, IndexMap map, bool periodic) {
        const auto face = owner.derive(HarmonicFace<Axis>{});
        EXPECT_EQ(mismatches(face, owner.layout(),
                             [&](CellIndex global) {
                                 return faceValue(map, global, Axis, Cells,
                                                  {periodic, periodic, periodic});
                             }),
                  0U)
            << "face axis " << Axis;
    }

    TEST(MaterialDerived, CompleteExtendedMap) {
        // Vary the decomposition axis independently of the three face directions;
        // physical runs also cover edges shared by a rank boundary and a global face.
        for (bool periodic : {false, true})
            for (unsigned axis = 0; axis < 3; ++axis) {
                auto grid = makeGrid(axis, periodic);
                const IndexMap map{2};
                const MaterialOwner owner(grid.getMesh(), grid.getFL(), map);
                const auto q = owner.derive(InverseSquare{});
                EXPECT_EQ(mismatches(q, owner.layout(),
                                     [&](CellIndex global) {
                                         return inverseSquare(extendedValue(
                                             map, global, Cells, {periodic, periodic, periodic}));
                                     }),
                          0U);
                checkFace<0>(owner, map, periodic);
                checkFace<1>(owner, map, periodic);
                checkFace<2>(owner, map, periodic);
                // Repeated requests must reuse storage. Five allocations account for
                // epsilon, the inverse square and one face coefficient per direction.
                EXPECT_EQ(owner.derive(InverseSquare{}).data(), q.data());
                EXPECT_EQ(owner.allocatedBytes(), 5 * q.span() * sizeof(double));
            }
    }

    /// The same key with another type must not silently return a different rule's field.
    struct NameCollision : InverseSquare {};

    TEST(MaterialDerived, CacheIdentityAndHomogeneousFaces) {
        auto grid      = makeGrid(0, false);
        const auto map = KOKKOS_LAMBDA(CellIndex, Position, bool& marked) {
            marked = false;
            return 2.13;
        };
        const MaterialOwner owner(grid.getMesh(), grid.getFL(), map);
        EXPECT_THROW(owner.derived<InverseSquare>(), std::out_of_range);
        owner.derive(InverseSquare{});
        EXPECT_THROW(owner.derive(NameCollision{}), std::invalid_argument);
        // Require the original homogeneous value in ghosts as well as owned cells.
        EXPECT_EQ(mismatches(owner.derive(HarmonicFace<0>{}), owner.layout(),
                             [](CellIndex) {
                                 return 2.13;
                             }),
                  0U);
    }
}  // namespace chdr::material::test
