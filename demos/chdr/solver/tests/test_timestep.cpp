/** @file test_timestep.cpp
 * @brief Check time-step selection and stable or unstable evolution.
 * @details Bounds are derived in Operators.h from fallahi2020mithra20fullwavesimulation, Eqs.
 * (3.17), (3.19)-(3.20), (3.24). Tests compare the same initial data below and above the spectral
 * stability bound.
 *
 * TODO: REVIEW.
 */
#include <gtest/gtest.h>
#include <iostream>
#include <limits>

#include "CoreTestSupport.h"
#include "SpatialDimension.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_timestep.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;
    using namespace chdr::solver::coreTest;

    struct ZeroSamples {
        KOKKOS_INLINE_FUNCTION double operator()(double, double, double, double, unsigned) const {
            return 0.0;
        }
    };

    /// @brief Compare accepted factors with an independently evaluated host bound.
    template <class SolverType>
    void checkAllowedFactors(SolverType& solver, double bound) {
        for (double safety : {0.2, 0.9, 0.98}) {
            solver.setSafetyFactor(safety);
            EXPECT_LE(std::abs(solver.getDt() / (safety * bound) - 1), 1e-15);
            EXPECT_LE(std::abs(solver.safetyFactor() / safety - 1), 1e-15);
        }
    }

    /// @brief Exercise shared invalid inputs and the two distinct equality rules.
    template <class PotentialSolverBase>
    void checkRejectedFactors(ChdrSolver<PotentialSolverBase>& solver, double bound) {
        for (double invalid : {0.0, -1.0, 1.02, std::numeric_limits<double>::quiet_NaN()})
            EXPECT_THROW(solver.setSafetyFactor(invalid), IpplException);
        if constexpr (std::is_same_v<PotentialSolverBase, StandardSolverPeriodic>) {
            EXPECT_THROW(solver.setSafetyFactor(1.0), IpplException);
        } else {
            EXPECT_THROW(solver.setSafetyFactor(1.0 + std::ldexp(1.0, -50)), IpplException);
            solver.setSafetyFactor(1.0);
            EXPECT_EQ(solver.getDt(), bound);
        }
    }

    /// @brief A reset retains the explicit step; evolved histories refuse any later change.
    template <class SolverType>
    void checkResetAndLock(SolverType& solver, double bound) {
        for (double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()})
            EXPECT_THROW(solver.setTimeStepUnchecked(invalid), IpplException);
        solver.setSafetyFactor(0.98);
        solver.initialize();
        EXPECT_EQ(solver.getDt(), 0.98 * bound);
        solver.solve();
        EXPECT_THROW(solver.setSafetyFactor(0.9), IpplException);
        EXPECT_THROW(solver.setTimeStepUnchecked(bound), IpplException);
    }

    /// @brief Keep field ownership in this scope while the solver holds mesh/layout references.
    template <class PotentialSolverBase>
    void checkSetters(std::array<int, SpatialDim> cells, std::array<double, SpatialDim> lengths) {
        auto fields = makePeriodicFieldContainer(cells, lengths);
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               EBReconstructionConvention::R1);
        const auto h = fields.getMesh().getMeshSpacing();
        const double bound =
            std::is_same_v<PotentialSolverBase, StandardSolverPeriodic>
                ? 1 / std::sqrt(1 / (h[0] * h[0]) + 1 / (h[1] * h[1]) + 1 / (h[2] * h[2]))
                : h[2];
        checkAllowedFactors(solver, bound);
        checkRejectedFactors(solver, bound);
        checkResetAndLock(solver, bound);
    }

    TEST(TimeStep, ValidatesRangesAndInitialization) {
        checkSetters<StandardSolverPeriodic>({16, 16, 16}, {1, 1, 1});
        checkSetters<StandardSolverPeriodic>({20, 24, 16}, {1, 0.96, 1});
        checkSetters<NonStandardSolverPeriodic>({16, 16, 24}, {1, 1, 1});
    }
    TEST(TimeStep, LoadedAndReconstructedHistoriesLockStep) {
        auto fields = makePeriodicFieldContainer({8, 8, 8}, {1, 1, 1});
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2);
        const auto revision = solver.historyRevision();
        loadHistory(solver, ZeroSamples{});
        EXPECT_GT(solver.historyRevision(), revision);
        EXPECT_THROW(solver.setSafetyFactor(0.9), IpplException);
        EXPECT_THROW(solver.setTimeStepUnchecked(solver.getDt()), IpplException);
        solver.initialize();
        EXPECT_NO_THROW(solver.setSafetyFactor(0.9));
        solver.evaluateInitialFields();
        EXPECT_THROW(solver.setSafetyFactor(0.8), IpplException);
        solver.initialize();
        EXPECT_NO_THROW(solver.setSafetyFactor(0.8));
    }
}  // namespace
