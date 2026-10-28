#include "Ippl.h"

#include <bit>
#include <cstdint>
#include <limits>

#include "BoundaryOracles.h"
#include "MaxwellSolvers/AbsorbingBC.h"
#include "gtest/gtest.h"

namespace {
    const ippl::Vector<double, 3> spacing{1.25, .8, 1.};
    constexpr double dt = .3;

    /** @brief Replay the original unit-speed constructor arithmetic under the compiler under test.
     * @details These expressions are copied from AbsorbingBC.h at IPPL 5ce8e37ed.
     * Keep their grouping and order: this compatibility reference checks bits,
     * while BoundaryOracles supplies the independent mathematical reference.
     */
    template <unsigned Normal, unsigned Side1, unsigned Side2>
    std::array<double, 5> shippedFace(ippl::Vector<double, 3> hr, double dt) {
        constexpr double c           = 1;
        constexpr unsigned main_axis = Normal, side_axes[2] = {Side1, Side2};
        std::array<double, 5> weights;
        double d = 1.0 / (2.0 * dt * hr[main_axis]) + 1.0 / (2.0 * c * dt * dt);

        weights[0] = (1.0 / (2.0 * dt * hr[main_axis]) - 1.0 / (2.0 * c * dt * dt)) / d;
        weights[1] = (-1.0 / (2.0 * dt * hr[main_axis]) - 1.0 / (2.0 * c * dt * dt)) / d;
        weights[2] = (1.0 / (c * dt * dt) - c / (2.0 * hr[side_axes[0]] * hr[side_axes[0]])
                      - c / (2.0 * hr[side_axes[1]] * hr[side_axes[1]]))
                     / d;
        weights[3] = (c / (4.0 * hr[side_axes[0]] * hr[side_axes[0]])) / d;
        weights[4] = (c / (4.0 * hr[side_axes[1]] * hr[side_axes[1]])) / d;
        return weights;
    }

    /// Preserve the shipped edge arithmetic, including its asymmetric normal terms.
    template <unsigned Edge, unsigned Normal1, unsigned Normal2>
    std::array<double, 5> shippedEdge(ippl::Vector<double, 3> hr, double dt) {
        constexpr double c0_         = 1;
        constexpr unsigned edge_axis = Edge, normal_axis1 = Normal1, normal_axis2 = Normal2;
        std::array<double, 5> weights;
        double d = (1.0 / hr[normal_axis1] + 1.0 / hr[normal_axis2]) / (4.0 * dt)
                   + 3.0 / (8.0 * c0_ * dt * dt);

        weights[0] = (-(1.0 / hr[normal_axis2] - 1.0 / hr[normal_axis1]) / (4.0 * dt)
                      - 3.0 / (8.0 * c0_ * dt * dt))
                     / d;
        weights[1] = ((1.0 / hr[normal_axis2] - 1.0 / hr[normal_axis1]) / (4.0 * dt)
                      - 3.0 / (8.0 * c0_ * dt * dt))
                     / d;
        weights[2] = ((1.0 / hr[normal_axis2] + 1.0 / hr[normal_axis1]) / (4.0 * dt)
                      - 3.0 / (8.0 * c0_ * dt * dt))
                     / d;
        weights[3] =
            (3.0 / (4.0 * c0_ * dt * dt) - c0_ / (4.0 * hr[edge_axis] * hr[edge_axis])) / d;
        weights[4] = c0_ / (8.0 * hr[edge_axis] * hr[edge_axis]) / d;
        return weights;
    }

    /// Preserve all seventeen shipped corner expressions before normalizing any coefficient.
    std::array<double, 17> shippedCorner(ippl::Vector<double, 3> hr, double dt) {
        constexpr double c0_ = 1;
        std::array<double, 17> weights;
        weights[0] =
            (-1.0 / hr[0] - 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[1] =
            (1.0 / hr[0] - 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[2] =
            (-1.0 / hr[0] + 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[3] =
            (-1.0 / hr[0] - 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[4] =
            (1.0 / hr[0] + 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[5] =
            (1.0 / hr[0] - 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[6] =
            (-1.0 / hr[0] + 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[7] =
            (1.0 / hr[0] + 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[8] =
            -(-1.0 / hr[0] - 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[9] =
            -(1.0 / hr[0] - 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[10] =
            -(-1.0 / hr[0] + 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[11] =
            -(-1.0 / hr[0] - 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[12] =
            -(1.0 / hr[0] + 1.0 / hr[1] - 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[13] =
            -(1.0 / hr[0] - 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[14] =
            -(-1.0 / hr[0] + 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[15] =
            -(1.0 / hr[0] + 1.0 / hr[1] + 1.0 / hr[2]) / (8.0 * dt) - 1.0 / (4.0 * c0_ * dt * dt);
        weights[16] = 1.0 / (2.0 * c0_ * dt * dt);
        return weights;
    }

    /** @brief Check every weight against its independent oracle and the unit-speed baseline.
     * @details The tolerance is four machine epsilons times the largest absolute oracle weight
     * of that struct; small or vanishing individual coefficients do not define the scale.
     */
    template <std::size_t N, typename Actual>
    void checkWeights(const Actual& actual, const std::array<double, N>& expected,
                      const std::array<double, N>& shipped, bool unitSpeed) {
        double scale = 0.;
        for (double value : expected)
            scale = std::max(scale, std::abs(value));
        for (unsigned q = 0; q < N; ++q) {
            EXPECT_NEAR(actual[q], expected[q], 4 * std::numeric_limits<double>::epsilon() * scale);
            if (unitSpeed) {
                EXPECT_EQ(std::bit_cast<std::uint64_t>(actual[q]),
                          std::bit_cast<std::uint64_t>(shipped[q]));
            }
        }
    }

    /// Mur's face weights: unequal tangential spacings, all three normals and both side signs.
    template <unsigned Normal, unsigned Side1, unsigned Side2>
    void checkFace() {
        using Face = second_order_abc_face<double, Normal, Side1, Side2>;
        for (double speed : {.6, 1., 1.2})
            for (int sign : {-1, 1}) {
                const Face face(spacing, dt, sign, speed);
                checkWeights(face.Cweights,
                             open_boundary_oracle::murFaceWeights(speed, dt, spacing[Normal],
                                                                  spacing[Side1], spacing[Side2]),
                             shippedFace<Normal, Side1, Side2>(spacing, dt), speed == 1.);
            }
        const Face defaults(spacing, dt, 1);
        checkWeights(defaults.Cweights,
                     open_boundary_oracle::murFaceWeights(1., dt, spacing[Normal], spacing[Side1],
                                                          spacing[Side2]),
                     shippedFace<Normal, Side1, Side2>(spacing, dt), true);
    }

    /// Low-level edge formulas at arbitrary speed; solver validation restricts intersecting faces.
    template <unsigned Edge, unsigned Normal1, unsigned Normal2>
    void checkEdge() {
        using EdgeType = second_order_abc_edge<double, Edge, Normal1, Normal2, true, true>;
        for (double speed : {.6, 1., 1.2}) {
            const EdgeType edge(spacing, dt, speed);
            checkWeights(edge.Eweights,
                         open_boundary_oracle::murEdgeWeights(speed, dt, spacing[Normal1],
                                                              spacing[Normal2], spacing[Edge]),
                         shippedEdge<Edge, Normal1, Normal2>(spacing, dt), speed == 1.);
        }
        const EdgeType defaults(spacing, dt);
        checkWeights(defaults.Eweights,
                     open_boundary_oracle::murEdgeWeights(1., dt, spacing[Normal1],
                                                          spacing[Normal2], spacing[Edge]),
                     shippedEdge<Edge, Normal1, Normal2>(spacing, dt), true);
    }
}  // namespace

TEST(MurSpeed, FaceAndEdgeWeights) {
    checkFace<0, 1, 2>();
    checkFace<1, 0, 2>();
    checkFace<2, 0, 1>();
    checkEdge<0, 1, 2>();
    checkEdge<1, 2, 0>();
    checkEdge<2, 0, 1>();
}

/// The corner oracle supplies both raw and normalized coefficients in the constructor's order.
TEST(MurSpeed, CornerWeightsAndRatios) {
    using Corner       = second_order_abc_corner<double, true, true, true>;
    const auto shipped = shippedCorner(spacing, dt);
    auto shippedRatios = shipped;
    for (auto& value : shippedRatios)
        value /= shipped[0];
    for (double speed : {.6, 1., 1.2}) {
        std::array<double, 17> raw, normalized, ratios;
        open_boundary_oracle::murCornerWeights(speed, dt, {spacing[0], spacing[1], spacing[2]}, raw,
                                               normalized);
        const Corner corner(spacing, dt, speed);
        for (unsigned q = 0; q < 17; ++q)
            ratios[q] = corner.Cweights[q] / corner.Cweights[0];
        checkWeights(corner.Cweights, raw, shipped, speed == 1.);
        checkWeights(ratios, normalized, shippedRatios, speed == 1.);
        if (speed == 1.)
            checkWeights(Corner(spacing, dt).Cweights, raw, shipped, true);
    }
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
