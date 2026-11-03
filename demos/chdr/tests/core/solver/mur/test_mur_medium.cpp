/** @file test_mur_medium.cpp
 * @brief Compare matched homogeneous absorption with the internal vacuum step.
 * @details Box: 16 x 16 x 16 cells, unit spacing, split along z on two ranks, collocated R1.
 * Initial potentials come from the NumPy fixture medium (seed 271828), see
 * RandomPotentialHistoryFixture.h.
 * - MatchedMediumEqualsInternalVacuumStep: with Mur on all six faces, the Homogeneous update
 *   at \f$\Delta t = 0.5\f$ and \f$\varepsilon_r\f$ = 2.13 or 3.7 must give A_np1 bit for bit
 *   equal to the vacuum kernel at \f$0.5/\sqrt{\varepsilon_r}\f$, after each of 256 steps.
 *   A vacuum run with the kernel Mur speed \f$1/\sqrt{\varepsilon_r}\f$ must differ at least
 *   once. Reference: the vacuum kernel (bitwise regression).
 * - SpeedControlPreservesStateAndClosure: in three periodicity patterns and both finish
 *   orders, the test-only kernel Mur speed 1.6 leaves the step, the face kinds and the speeds
 *   unchanged, and only the old finish order leaves open-face ghosts incomplete.
 *   Reference: self-consistency.
 *
 * CTest: chdr.solver.mur.medium.serial and chdr.solver.mur.medium.mpi2 (filter MurMedium.*).
 */
#include <cstring>
#include <gtest/gtest.h>

#include "OpenRunSupport.h"
#include "core/solver/BoxSetup.h"
#include "tests/core/solver/support/CoreTestSupport.h"
#include "tests/core/solver/support/RandomPotentialHistoryFixture.h"
#include "tests/core/solver/support/RankRecords.h"

namespace chdr::solver::test {
    /** @brief Compare a medium step with the same vacuum kernel at its internal timestep.
     * All three source-free trajectories (medium, vacuum reference, vacuum control) load the
     * same independent random histories and run 256 steps with Mur on all faces.
     * Identical histories and closures isolate the speed policy;
     * setting the kernel speed to the inverse refractive index is a control.
     * Potentials must agree bitwise, since both calls execute the same recurrence
     * [chew2014generalizedgauge, equations (12)-(13); HomogeneousUpdate].
     */
    TEST(MurMedium, MatchedMediumEqualsInternalVacuumStep) {
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
     * Layouts: all periodic and x open (physical speed 0.8), all open (matched speed);
     * step 0.25, control speed 1.6. After step() and after initialize(), the step, the face
     * kinds and the speeds are unchanged. Ghosts are complete except after the old finish
     * order on a layout with a Mur face. Control speeds 0 and infinity throw IpplException.
     */
    TEST(MurMedium, SpeedControlPreservesStateAndClosure) {
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
