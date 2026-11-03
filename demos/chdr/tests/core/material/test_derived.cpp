/** @file test_derived.cpp
 * @brief Compare every cached coefficient cell with its independently authored closed form.
 * @details IndexMap{2} gives values in \f$[1, 16)\f$ on the 12 by 14 by 16 grid of
 * OwnerTestSupport.h, periodic or open, with each split axis. Every allocated cell of
 * \f$Q = \varepsilon^{-2}\f$ and of the three harmonic face fields must be bitwise equal to the
 * closed form on the extended map (OwnerOracles.h). In the lower ghost layer of an open axis
 * the face value is \f$H(\varepsilon_0, \varepsilon_0) = \varepsilon_0\f$, with
 * \f$\varepsilon_0\f$ the value of cell 0 on that axis. A repeated derive() must return the
 * same storage, and the owner must hold five fields. derived() before derive() must throw
 * std::out_of_range, and a rule with a colliding name must throw std::invalid_argument. A
 * uniform 2.13 must give 2.13 in every face cell, ghosts included.
 *
 * Reference: exact lattice value; self-consistency (cache, refusal). CTest entries:
 * chdr.material.Derived.serial and chdr.material.Derived.mpi2 (filter MaterialDerived.*).
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

    /// Periodic and open grids, 3 split axes: Q and the three face fields are bitwise equal to
    /// their closed forms in every allocated cell, a repeated derive() reuses the storage, and
    /// the owner holds five fields.
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

    /// Uniform 2.13 on an open grid split along x: derived() before derive() throws, a name
    /// collision throws, and every allocated cell of the x face field equals 2.13.
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
