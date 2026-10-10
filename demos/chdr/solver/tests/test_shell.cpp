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
#include "SpatialDimension.h"

static_assert(chdr::solver::SpatialDim == 3,
              "test_shell.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;

    // Unit conversions must remain explicit at the label API boundary.
    static_assert(!std::is_convertible_v<int, HalfStepOffset>);
    static_assert(!std::is_convertible_v<HalfStepOffset, int>);
    static_assert(
        std::is_same_v<decltype(LevelContract{}.halfStepLabel(Quantity::Phi)), HalfStepOffset>);

    /** @brief Compare owned E or B values, including explicit nonfinite rejection.
     * @details The reduction is over the same field layout on both sides. A global
     * maximum precedes assertions, so every rank observes the same outcome.
     */
    double difference(EMFieldsType& first, EMFieldsType& second) {
        const auto a = first.getView(), b = second.getView();
        const double infinity = std::numeric_limits<double>::infinity();
        double error          = 0.0;
        Kokkos::parallel_reduce(
            "Shell field difference", first.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& maximum) {
                for (unsigned c = 0; c < SpatialDim; ++c) {
                    const double value = Kokkos::abs(a(i, j, k)[c] - b(i, j, k)[c]);
                    maximum = Kokkos::max(maximum, Kokkos::isfinite(value) ? value : infinity);
                }
            },
            Kokkos::Max<double>(error));
        MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, first.getLayout().comm);
        return error;
    }

    /// @brief Exercise the same base solve path for 128 steps, then check explicit reset counting.
    void compareBaseEvolution(StandardSolverPeriodic& base,
                              ChdrSolver<StandardSolverPeriodic>& shell,
                              ChDRFieldContainer& baseFields, ChDRFieldContainer& shellFields) {
        EXPECT_EQ(shell.getDt(), base.getDt());
        for (int step = 1; step <= 128; ++step) {
            base.solve();
            shell.solve();
            EXPECT_EQ(shell.completedSolves(), step);
            EXPECT_EQ(difference(baseFields.getE(), shellFields.getE()), 0.0);
            EXPECT_EQ(difference(baseFields.getB(), shellFields.getB()), 0.0);
        }
        shell.initialize();
        EXPECT_EQ(shell.completedSolves(), 0);
    }

    TEST(SolverShell, MatchesBaseOnM001) {
        auto baseFields  = chdr::vacuum::makePeriodicFieldContainer();
        auto shellFields = chdr::vacuum::makePeriodicFieldContainer();
        chdr::vacuum::initializeFieldStorage(baseFields);
        chdr::vacuum::initializeFieldStorage(shellFields);
        StandardSolverPeriodic base(baseFields.getJ(), baseFields.getE(), baseFields.getB());
        ChdrSolver<StandardSolverPeriodic> shell(shellFields.getJ(), shellFields.getE(),
                                                 shellFields.getB(),
                                                 EBReconstructionConvention::R1);
        chdr::vacuum::initializePotentials(base);
        chdr::vacuum::initializePotentials(shell);
        compareBaseEvolution(base, shell, baseFields, shellFields);
    }

    TEST(SolverShell, DeclaresMixedElectricLabel) {
        const LevelContract contract{EBReconstructionConvention::R1};
        EXPECT_THROW(contract.halfStepLabel(Quantity::Electric), std::invalid_argument);
        EXPECT_EQ(contract.halfStepLabel(Quantity::Electric, true).count(), -1);
        EXPECT_EQ(contract.labelTime(Quantity::Magnetic, 3, 0.125), 0.375);
        EXPECT_EQ(contract.labelTime(Quantity::Electric, 3, 0.125, true), 0.3125);
    }

    TEST(SolverShell, RefusesUnsupportedHaloWidth) {
        auto fields = chdr::vacuum::makePeriodicFieldContainer();
        chdr::vacuum::initializeFieldStorage(fields);
        EMFieldsType wrong;
        wrong.initialize(fields.getMesh(), fields.getFL(), 2);
        EXPECT_THROW((ChdrSolver<StandardSolverPeriodic>(fields.getJ(), wrong, fields.getB(),
                                                         EBReconstructionConvention::R1)),
                     IpplException);
    }

    TEST(SolverShell, RequiresExplicitChoiceAndFixedBindings) {
        static_assert(!std::is_constructible_v<ChdrSolver<StandardSolverPeriodic>,
                                               AJFourFieldsType&, EMFieldsType&, EMFieldsType&>);
        EXPECT_EQ(LevelContract{}.ebReconstructionConvention_m, EBReconstructionConvention::R2);
        auto fields = makePeriodicFieldContainer({8, 8, 8}, {1, 1, 1});
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2);
        EXPECT_THROW(solver.setSources(fields.getJ()), IpplException);
        EXPECT_THROW(solver.setEMFields(fields.getE(), fields.getB()), IpplException);
        EXPECT_THROW(solver.computeRawMagneticField(fields.getE()), IpplException);
        EXPECT_THROW(solver.computeRawMagneticField(fields.getB()), IpplException);
    }

    TEST(SolverShell, RejectsDistinctLayoutOwners) {
        auto first  = makePeriodicFieldContainer({8, 8, 8}, {1, 1, 1});
        auto second = makePeriodicFieldContainer({8, 8, 8}, {1, 1, 1});
        initializeFieldStorage(first);
        initializeFieldStorage(second);
        EXPECT_THROW((ChdrSolver<StandardSolverPeriodic>(first.getJ(), second.getE(), first.getB(),
                                                         EBReconstructionConvention::R2)),
                     IpplException);
    }

    TEST(SolverShell, InitialFieldsHonorSelectedReconstruction) {
        auto fields = makePeriodicFieldContainer({8, 8, 8}, {1, 1, 1});
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2);
        chdr::vacuum::initializePotentials(solver);
        solver.evaluateInitialFields();
        EMFieldsType expected;
        expected.initialize(fields.getMesh(), fields.getFL());
        Kokkos::deep_copy(expected.getView(), fields.getB().getView());
        solver.evaluate_EB();
        EXPECT_EQ(difference(expected, fields.getB()), 0.0);
        EXPECT_EQ(solver.completedSolves(), 0);
    }

    TEST(SolverShell, CollectivelyRefusesInsufficientLocalOwnership) {
        const int ranks = ippl::Comm->size();
        if (ranks < 2)
            GTEST_SKIP() << "Requires a split with unequal local extents";
        auto fields = makePeriodicFieldContainer({2 * ranks - 1, 8, 8}, {1, 1, 1});
        initializeFieldStorage(fields);
        EXPECT_THROW(
            (ChdrSolver<StandardSolverPeriodic>(fields.getJ(), fields.getE(), fields.getB(),
                                                EBReconstructionConvention::R2)),
            std::runtime_error);
    }
}  // namespace
