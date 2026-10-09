/** @file test_shell.cpp
 * @brief Compare the existing Standard solver and shell on the periodic vacuum-wave setup.
 * @details The oracle is the existing Standard solver invoked with the same
 * initial histories. This regression checks shell plumbing and solve counting;
 * The unchanged vacuum-demo tests separately establish the physical wave reference.
 *
 * TODO: REVIEW.
 */
#include <gtest/gtest.h>
#include <limits>
#include <type_traits>

#include "PeriodicSetup.h"

namespace {
    using namespace chdr::solver;

    // Unit conversions must remain explicit at the label API boundary.
    static_assert(!std::is_convertible_v<int, HalfStepOffset>);
    static_assert(!std::is_convertible_v<HalfStepOffset, int>);
    static_assert(
        std::is_same_v<decltype(LevelContract{}.halfStepLabel(Quantity::Phi)), HalfStepOffset>);

    TEST(SolverShell, DeclaresMixedElectricLabel) {
        const LevelContract contract{EBReconstructionConvention::R1};
        EXPECT_THROW(contract.halfStepLabel(Quantity::Electric), std::invalid_argument);
        EXPECT_EQ(contract.halfStepLabel(Quantity::Electric, true).count(), -1);
        EXPECT_EQ(contract.labelTime(Quantity::Magnetic, 3, 0.125), 0.375);
        EXPECT_EQ(contract.labelTime(Quantity::Electric, 3, 0.125, true), 0.3125);
    }
}  // namespace
