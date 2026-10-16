#include "Ippl.h"

#include "MaxwellSolvers/FDTDBoundary.h"

#include <limits>

#include "BoundaryOracles.h"
#include "gtest/gtest.h"

namespace {
    /// Split along z; the small default box has two normal samples per rank on two ranks.
    ippl::FieldLayout<3> makeLayout(std::array<bool, 3> periodic,
                                    std::array<int, 3> size = {6, 5, 4}, int first = 0) {
        ippl::NDIndex<3> domain;
        for (unsigned d = 0; d < 3; ++d)
            domain[d] = ippl::Index(first, first + size[d] - 1);
        return {MPI_COMM_WORLD, domain, {false, false, true}, periodic};
    }

    /** @brief Check paired face kinds and interior ownership against the independent cell rule.
     * @details Visit owned global indices only. A cell belongs to the interior exactly when
     * the oracle assigns it no Mur face, including the periodic faces of mixed boxes.
     */
    void checkDescription(std::array<bool, 3> periodic, std::array<int, 3> size) {
        auto layout         = makeLayout(periodic, size);
        const auto boundary = ippl::FDTDBoundary::fromLayout(layout);
        const ippl::InteriorBounds interior(layout.getDomain(), boundary);
        for (unsigned face = 0; face < 6; ++face) {
            EXPECT_EQ(boundary.kind[face],
                      periodic[face / 2] ? ippl::FaceKind::periodic : ippl::FaceKind::mur);
            EXPECT_EQ((boundary.murMask() >> face) & 1u, unsigned(!periodic[face / 2]));
            EXPECT_EQ(boundary.speed[face], 1.);
        }
        for (int i = 1; i <= size[0]; ++i)
            for (int j = 1; j <= size[1]; ++j)
                for (int k = 1; k <= size[2]; ++k)
                    EXPECT_EQ(
                        interior.contains(i, j, k),
                        open_boundary_oracle::murFaceMask({i, j, k}, {size[0], size[1], size[2]},
                                                          {periodic[0], periodic[1], periodic[2]})
                            == 0);
    }
}  // namespace

TEST(FDTDBoundary, DescriptionAndInteriorOwnership) {
    for (unsigned bits = 0; bits < 8; ++bits)
        checkDescription({bool(bits & 1), bool(bits & 2), bool(bits & 4)}, {6, 5, 4});
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
