/** @file test_solver_contract.cpp
 * @brief The construction contract of ChdrSolver, and its evolution against IPPL's Standard
 * solver on the periodic vacuum-wave setup.
 * @details The oracle is IPPL's Standard solver invoked with the same initial histories. This
 * regression checks the field bindings of ChdrSolver and its solve counting; the vacuum-wave
 * tests separately establish the physical wave reference.
 */
#include <Kokkos_MathematicalConstants.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <type_traits>

#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"

static_assert(
    chdr::SpatialDim == 3,
    "test_solver_contract.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr;
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
            "ChdrSolver field difference", first.getFieldRangePolicy(),
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

    /// The periodic box of the vacuum plane wave: 8 by 8 by 64 cells on [0, 1/8]^2 x [0, 1].
    ChDRFieldContainer makeVacuumWaveBox() {
        return makeBoxFieldContainer({8, 8, 64}, {0.125, 0.125, 1.0}, 2);
    }

    /** @brief Load the +z vacuum plane wave into both histories of a configured solver.
     * @details The four-potential is \f$A_x=\sin(2\pi(z-t))/(2\pi)\f$ with the other
     * components zero, a solution of the vacuum wave equation for the potentials
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.6)-(3.9)] in units \f$c=1\f$ with
     * wavelength one. Both histories are sampled at collocated cell centres, at times 0
     * and \f$-\Delta t\f$.
     */
    void initializeVacuumHistories(StandardSolverPeriodic& solver) {
        const auto sample = KOKKOS_LAMBDA(double, double, double z, double t, unsigned component) {
            constexpr double WAVE_NUMBER = 2.0 * Kokkos::numbers::pi_v<double>;
            return component == 1 ? Kokkos::sin(WAVE_NUMBER * (z - t)) / WAVE_NUMBER : 0.0;
        };
        const double dt = solver.getDt();
        writeSamples(solver.A_n, sample, {0, 0, 0, 0});
        writeSamples(solver.A_nm1, sample, {-dt, -dt, -dt, -dt});
        solver.prepareHistories();
    }

    /// @brief Exercise the same base solve path for 128 steps, then check explicit reset counting.
    void compareBaseEvolution(StandardSolverPeriodic& base,
                              ChdrSolver<StandardSolverPeriodic>& chdrSolver,
                              ChDRFieldContainer& baseFields, ChDRFieldContainer& chdrFields) {
        EXPECT_EQ(chdrSolver.getDt(), base.getDt());
        for (int step = 1; step <= 128; ++step) {
            base.solve();
            chdrSolver.solve();
            EXPECT_EQ(chdrSolver.completedSolves(), step);
            EXPECT_EQ(difference(baseFields.getE(), chdrFields.getE()), 0.0);
            EXPECT_EQ(difference(baseFields.getB(), chdrFields.getB()), 0.0);
        }
        chdrSolver.initialize();
        EXPECT_EQ(chdrSolver.completedSolves(), 0);
    }

    TEST(ChdrSolverContract, MatchesBaseSolverOnTheVacuumWave) {
        auto baseFields = makeVacuumWaveBox();
        auto chdrFields = makeVacuumWaveBox();
        initializeFieldStorage(baseFields);
        initializeFieldStorage(chdrFields);
        StandardSolverPeriodic base(baseFields.getJ(), baseFields.getE(), baseFields.getB(),
                                    ippl::FDTDBoundary::fromLayout(baseFields.getFL()));
        ChdrSolver<StandardSolverPeriodic> chdrSolver(chdrFields.getJ(), chdrFields.getE(),
                                                      chdrFields.getB(),
                                                      EBReconstructionConvention::R1);
        initializeVacuumHistories(base);
        initializeVacuumHistories(chdrSolver);
        compareBaseEvolution(base, chdrSolver, baseFields, chdrFields);
    }

    TEST(ChdrSolverContract, DeclaresMixedElectricLabel) {
        const LevelContract contract{EBReconstructionConvention::R1};
        EXPECT_THROW(contract.halfStepLabel(Quantity::Electric), std::invalid_argument);
        EXPECT_EQ(contract.halfStepLabel(Quantity::Electric, true).count(), -1);
        EXPECT_EQ(contract.labelTime(Quantity::Magnetic, 3, 0.125), 0.375);
        EXPECT_EQ(contract.labelTime(Quantity::Electric, 3, 0.125, true), 0.3125);
    }

    TEST(ChdrSolverContract, RefusesUnsupportedHaloWidth) {
        auto fields = makeVacuumWaveBox();
        initializeFieldStorage(fields);
        EMFieldsType wrong;
        wrong.initialize(fields.getMesh(), fields.getFL(), 2);
        EXPECT_THROW((ChdrSolver<StandardSolverPeriodic>(fields.getJ(), wrong, fields.getB(),
                                                         EBReconstructionConvention::R1)),
                     IpplException);
    }

    TEST(ChdrSolverContract, RequiresExplicitChoiceAndFixedBindings) {
        static_assert(!std::is_constructible_v<ChdrSolver<StandardSolverPeriodic>,
                                               AJFourFieldsType&, EMFieldsType&, EMFieldsType&>);
        EXPECT_EQ(LevelContract{}.ebReconstructionConvention_m, EBReconstructionConvention::R2);
        auto fields = makeBoxFieldContainer({8, 8, 8}, {1, 1, 1});
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2);
        EXPECT_THROW(solver.setSources(fields.getJ()), IpplException);
        EXPECT_THROW(solver.setEMFields(fields.getE(), fields.getB()), IpplException);
        EXPECT_THROW(solver.computeRawMagneticField(fields.getE()), IpplException);
        EXPECT_THROW(solver.computeRawMagneticField(fields.getB()), IpplException);
    }

    TEST(ChdrSolverContract, RejectsDistinctLayoutOwners) {
        auto first  = makeBoxFieldContainer({8, 8, 8}, {1, 1, 1});
        auto second = makeBoxFieldContainer({8, 8, 8}, {1, 1, 1});
        initializeFieldStorage(first);
        initializeFieldStorage(second);
        EXPECT_THROW((ChdrSolver<StandardSolverPeriodic>(first.getJ(), second.getE(), first.getB(),
                                                         EBReconstructionConvention::R2)),
                     IpplException);
    }

    TEST(ChdrSolverContract, InitialFieldsHonorSelectedReconstruction) {
        auto fields = makeBoxFieldContainer({8, 8, 8}, {1, 1, 1});
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2);
        initializeVacuumHistories(solver);
        solver.evaluateInitialFields();
        EMFieldsType expected;
        expected.initialize(fields.getMesh(), fields.getFL());
        Kokkos::deep_copy(expected.getView(), fields.getB().getView());
        solver.evaluate_EB();
        EXPECT_EQ(difference(expected, fields.getB()), 0.0);
        EXPECT_EQ(solver.completedSolves(), 0);
    }

    TEST(ChdrSolverContract, CollectivelyRefusesInsufficientLocalOwnership) {
        const int ranks = ippl::Comm->size();
        if (ranks < 2)
            GTEST_SKIP() << "Requires a split with unequal local extents";
        auto fields = makeBoxFieldContainer({2 * ranks - 1, 8, 8}, {1, 1, 1});
        initializeFieldStorage(fields);
        EXPECT_THROW(
            (ChdrSolver<StandardSolverPeriodic>(fields.getJ(), fields.getE(), fields.getB(),
                                                EBReconstructionConvention::R2)),
            std::runtime_error);
    }
}  // namespace
