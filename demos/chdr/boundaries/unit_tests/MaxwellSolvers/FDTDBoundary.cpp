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

/** Opposite faces of a single open axis may have distinct valid speeds. Periodic speeds have
 * no meaning. Intersecting open faces require the matched kernel speed, as validateBoundary
 * documents from the edge/corner compatibility equations.
 */
TEST(FDTDBoundary, ValidateSpeedsAndIntersections) {
    auto layout   = makeLayout({false, true, true});
    auto boundary = ippl::FDTDBoundary::fromLayout(layout, .8);
    EXPECT_EQ(boundary.speed[0], .8);
    EXPECT_EQ(boundary.speed[1], .8);
    boundary.speed = {.8,
                      1.2,
                      0.,
                      -1.,
                      std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::quiet_NaN()};
    EXPECT_NO_THROW(ippl::validateBoundary(boundary, layout, 1, 1, 1));
    for (double speed : {0., -1., std::sqrt(2.), 2., std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
        boundary.speed[0] = speed;
        EXPECT_THROW(ippl::validateBoundary(boundary, layout, 1, 1, 1), IpplException);
    }
    for (const auto periodic : {std::array{false, false, true}, std::array{false, false, false}}) {
        auto intersecting = makeLayout(periodic);
        auto matched      = ippl::FDTDBoundary::fromLayout(intersecting);
        EXPECT_NO_THROW(ippl::validateBoundary(matched, intersecting, 1, 1, 1));
        matched.speed[1] = .9;
        EXPECT_THROW(ippl::validateBoundary(matched, intersecting, 1, 1, 1), IpplException);
    }
}

/// Invalid descriptions and storage widths must be rejected on every rank before exchange.
TEST(FDTDBoundary, RejectFaceMismatchAndGhostWidths) {
    auto layout         = makeLayout({false, true, true});
    const auto boundary = ippl::FDTDBoundary::fromLayout(layout);
    for (unsigned face = 0; face < 6; ++face) {
        auto mismatched       = boundary;
        mismatched.kind[face] = face < 2 ? ippl::FaceKind::periodic : ippl::FaceKind::mur;
        EXPECT_THROW(ippl::validateBoundary(mismatched, layout, 1, 1, 1), IpplException);
    }
    for (int width : {0, 2}) {
        EXPECT_THROW(ippl::validateBoundary(boundary, layout, width, 1, 1), IpplException);
        EXPECT_THROW(ippl::validateBoundary(boundary, layout, 1, width, 1), IpplException);
        EXPECT_THROW(ippl::validateBoundary(boundary, layout, 1, 1, width), IpplException);
    }
}

/** Three global samples support extrapolation; two local samples support Mur's next-level
 * reads. The three-cell split gives one rank only one cell, and both ranks must reject it.
 */
TEST(FDTDBoundary, RejectUnsupportedDomains) {
    auto tooShort     = makeLayout({false, true, true}, {2, 5, 4});
    auto thinPeriodic = makeLayout({true, true, true}, {1, 5, 4});
    auto shifted      = makeLayout({false, false, false}, {6, 5, 4}, 1);
    for (const auto* layout : {&tooShort, &thinPeriodic, &shifted})
        EXPECT_THROW(
            ippl::validateBoundary(ippl::FDTDBoundary::fromLayout(*layout), *layout, 1, 1, 1),
            IpplException);
    if (ippl::Comm->size() == 2) {
        auto split = makeLayout({false, false, false}, {6, 5, 3});
        EXPECT_THROW(ippl::validateBoundary(ippl::FDTDBoundary::fromLayout(split), split, 1, 1, 1),
                     IpplException);
    }
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
