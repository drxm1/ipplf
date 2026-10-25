/** @file test_open_medium.cpp
 * @brief Compare matched homogeneous absorption with the internal vacuum step.
 */
#include <cstring>
#include <gtest/gtest.h>

#include "BoxSetup.h"
#include "CoreTestSupport.h"
#include "OpenRunSupport.h"
#include "RandomPotentialHistoryFixture.h"

namespace chdr::solver::test {
    /** @brief Compare a medium step with the same vacuum kernel at its internal timestep.
     * Both source-free trajectories load the same independent random histories.
     * Identical histories and closures isolate the speed policy;
     * setting the kernel speed to the inverse refractive index is a control.
     * Potentials must agree bitwise, since both calls execute the same recurrence
     * [chew2014generalizedgauge, equations (12)-(13); HomogeneousUpdate].
     */
    TEST(OpenShell, MatchedMediumEqualsInternalVacuumStep) {
        for (double epsilon : {2.13, 3.7}) {
            auto mediumFields =
                makeBoxFieldContainer({16, 16, 16}, {16, 16, 16}, 2, {false, false, false});
            auto baseFields =
                makeBoxFieldContainer({16, 16, 16}, {16, 16, 16}, 2, {false, false, false});
            auto wrongFields =
                makeBoxFieldContainer({16, 16, 16}, {16, 16, 16}, 2, {false, false, false});
            for (auto* fields : {&mediumFields, &baseFields, &wrongFields})
                initializeFieldStorage(*fields);
            auto medium = makeOpenTestSolver(mediumFields.getJ(), mediumFields.getE(),
                                             mediumFields.getB(), epsilon);
            OpenStandardSolver base(baseFields.getJ(), baseFields.getE(), baseFields.getB(),
                                    EBReconstructionConvention::R1);
            OpenStandardSolver wrong(wrongFields.getJ(), wrongFields.getE(), wrongFields.getB(),
                                     EBReconstructionConvention::R1);
            medium->setTimeStepUnchecked(.5);
            base.setTimeStepUnchecked(.5 / std::sqrt(epsilon));
            wrong.setTimeStepUnchecked(base.getDt());
            wrong.setKernelMurSpeedForTest(1 / std::sqrt(epsilon));
            for (auto* solver : {medium.get(), &base, &wrong}) {
                loadRandomPotentialHistoryFixture(*solver, "medium", {16, 16, 16});
            }
            bool wrongDiffers = false;
            for (int step = 0; step < 256; ++step) {
                medium->step();
                base.step();
                wrong.step();
                const auto observed  = gatherOwned(medium->A_np1),
                           reference = gatherOwned(base.A_np1), control = gatherOwned(wrong.A_np1);
                if (ippl::Comm->rank() == 0) {
                    EXPECT_EQ(std::memcmp(observed.data(), reference.data(),
                                          observed.size() * sizeof(double)),
                              0);
                    wrongDiffers |= std::memcmp(observed.data(), control.data(),
                                                observed.size() * sizeof(double))
                                    != 0;
                }
                medium->timeShift();
                base.timeShift();
                wrong.timeShift();
            }
            if (ippl::Comm->rank() == 0) {
                EXPECT_TRUE(wrongDiffers);
            }
        }
    }
    /** @brief Detect whether step completion already published and closed every ghost.
     * Repeating the production halo/closure operations must leave a completed field
     * bitwise unchanged. The shipped-order control deliberately fails this check on
     * open boxes; the independently loaded histories make the omitted values nonzero.
     */
    bool controlGhostsComplete(OpenStandardSolver& solver) {
        const auto before = coreTest::copyToHost(solver.A_np1);
        solver.A_np1.fillHalo();
        ippl::fillOpenFaceGhosts(solver.A_np1, solver.boundary().murMask());
        const auto after = coreTest::copyToHost(solver.A_np1);
        return std::memcmp(before.data(), after.data(),
                           before.size() * sizeof(AJFourFieldsType::value_type))
               == 0;
    }

    /** @brief Check immutable configuration and both completion orders with a speed control.
     * Include a non-unit planar speed to detect accidental descriptor normalization.
     * The control speed exceeds the production guard deliberately; this test exercises
     * the existing negative-control bypass, not an admissible propagation configuration.
     */
    TEST(OpenShell, SpeedControlPreservesStateAndClosure) {
        for (auto axes : {std::array{true, true, true}, std::array{false, true, true},
                          std::array{false, false, false}}) {
            for (bool shipped : {false, true}) {
                auto fields = makeBoxFieldContainer({16, 16, 16}, {16, 16, 16}, 2, axes);
                initializeFieldStorage(fields);
                OpenStandardSolver solver(fields.getJ(), fields.getE(), fields.getB(),
                                          EBReconstructionConvention::R1,
                                          SpatialDiscretization::Collocated,
                                          axes[1] ? MurSpeed::physical(.8) : MurSpeed::matched());
                solver.setTimeStepUnchecked(.25);
                loadRandomPotentialHistoryFixture(solver, "medium", {16, 16, 16});
                const auto boundary = solver.boundary();
                const auto step     = solver.getDt();
                solver.useShippedOrder(shipped);
                solver.setKernelMurSpeedForTest(1.6);
                solver.step();
                EXPECT_EQ(solver.getDt(), step);
                EXPECT_EQ(solver.boundary().kind, boundary.kind);
                EXPECT_EQ(solver.boundary().speed, boundary.speed);
                EXPECT_EQ(controlGhostsComplete(solver), !shipped || boundary.murMask() == 0);
                solver.initialize();
                EXPECT_EQ(solver.getDt(), step);
                EXPECT_EQ(solver.boundary().kind, boundary.kind);
                EXPECT_EQ(solver.boundary().speed, boundary.speed);
                EXPECT_THROW(solver.setKernelMurSpeedForTest(0), IpplException);
                EXPECT_THROW(
                    solver.setKernelMurSpeedForTest(std::numeric_limits<double>::infinity()),
                    IpplException);
            }
        }
    }
}  // namespace chdr::solver::test
