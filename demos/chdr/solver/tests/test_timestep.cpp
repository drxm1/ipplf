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
#include "MonitorTestSupport.h"
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

    /** @brief Fill the scalar component with the original deterministic global-index hash.
     * @details Each rank samples the same integer rule at a given global cell.
     * The two history slots have distinct seeds; vector components remain zero.
     * No halo is read here, and the caller synchronizes after removing the mean.
     */
    void fillRandomPotential(AJFourFieldsType& field, unsigned seed) {
        const auto view   = field.getView();
        const auto domain = field.getLayout().getLocalNDIndex();
        const ippl::Vector<int, SpatialDim> first(domain[0].first(), domain[1].first(),
                                                  domain[2].first());
        Kokkos::parallel_for(
            field.getFieldRangePolicy(), KOKKOS_LAMBDA(int i, int j, int k) {
                unsigned hash = static_cast<unsigned>((i + first[0]) * 104729 + (j + first[1]) * 257
                                                      + k + first[2])
                                + seed;
                hash ^= hash >> 16;
                hash *= 0x7feb352dU;
                hash ^= hash >> 15;
                view(i, j, k)    = 0.0;
                view(i, j, k)[0] = static_cast<double>(hash % 10007) / 10007.0 - 0.5;
            });
    }

    /** @brief Remove the globally constant mode from one history before evolution.
     * @details The sum uses owned cells and the field communicator. Subtraction
     * retains the original reduction and arithmetic order; it is followed by a
     * true-copy halo synchronization. A constant spatial mode has zero Laplacian,
     * so unequal history means would produce linear time growth (deduction from
     * fallahi2020mithra20fullwavesimulation, Eqs. (3.11)-(3.15)) unrelated to the stability-bound
     * comparison.
     */
    void removeGlobalMean(AJFourFieldsType& field) {
        const auto view = field.getView();
        double sum      = 0;
        Kokkos::parallel_reduce(
            field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& result) { result += view(i, j, k)[0]; },
            sum);
        MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, field.getLayout().comm);
        const auto global = field.getLayout().getDomain();
        const double mean = sum / (global[0].length() * global[1].length() * global[2].length());
        Kokkos::parallel_for(
            field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) { view(i, j, k)[0] -= mean; });
        preparePotentialHalo(field);
    }

    /// @brief Initialize both histories independently, retaining the original seed order.
    template <class SolverType>
    void randomSlots(SolverType& solver) {
        for (auto* field : {&solver.A_n, &solver.A_nm1}) {
            fillRandomPotential(*field, field == &solver.A_n ? 4321 : 7890);
            removeGlobalMean(*field);
        }
    }

    /// @brief Store maxima over the frozen first and last observation windows.
    struct Envelope {
        double first_m = 0, last_m = 0;
        void observe(int step, bool bounded, double value) {
            if (step < (bounded ? 1000 : 1))
                first_m = std::max(first_m, value);
            if (step >= (bounded ? 9000 : 1))
                last_m = std::max(last_m, value);
        }
    };

    /// @brief Advance unchanged updates and observe the scalar maximum after each solve.
    template <class SolverType>
    void observeEnvelope(SolverType& solver, bool bounded, Envelope& envelope) {
        const int steps = bounded ? 10000 : 60;
        Monitors monitors(solver.sourceField());
        for (int step = 0; step < steps; ++step) {
            chdr::solver::test::advanceWithFinalMonitor(solver, monitors, step + 1 == steps);
            const auto view    = solver.A_n.getView();
            const double local = localMaximum(
                solver.A_n,
                KOKKOS_LAMBDA(int i, int j, int k) { return Kokkos::abs(view(i, j, k)[0]); });
            const double maximum = globalMaximum(solver.A_n, local);
            ASSERT_TRUE(std::isfinite(maximum));
            envelope.observe(step, bounded, maximum);
        }
    }

    /// @brief Apply the same frozen growth/upper-envelope comparisons after global observations.
    void expectEnvelope(const Envelope& envelope, bool bounded) {
        // "A8" labels the observed stable or unstable potential-amplitude envelope.
        if (ippl::Comm->rank() == 0)
            std::cout << "A8 bounded=" << bounded << " ratio=" << envelope.last_m / envelope.first_m
                      << '\n';
        if (bounded) {
            EXPECT_LE(envelope.last_m, 2 * envelope.first_m);
        } else {
            EXPECT_GE(envelope.last_m, 1e3 * envelope.first_m);
        }
    }

    /// @brief Build a stable owner, set the deliberately checked or unchecked experiment, and run.
    template <class PotentialSolverBase>
    void stability(std::array<int, SpatialDim> cells, std::array<double, SpatialDim> lengths,
                   bool bounded) {
        auto fields = makePeriodicFieldContainer(cells, lengths);
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               EBReconstructionConvention::R1);
        solver.setTimeStepUnchecked((bounded ? 0.98 : 1.02) * solver.maxStableTimeStep());
        randomSlots(solver);
        Envelope envelope;
        observeEnvelope(solver, bounded, envelope);
        expectEnvelope(envelope, bounded);
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
    TEST(TimeStep, GrowsAboveBound) {
        stability<StandardSolverPeriodic>({16, 16, 16}, {1, 1, 1}, false);
        stability<StandardSolverPeriodic>({20, 24, 16}, {1, 0.96, 1}, false);
        stability<NonStandardSolverPeriodic>({16, 16, 24}, {1, 1, 1}, false);
    }
    TEST(TimeStepLong, StandardBounded) {
        stability<StandardSolverPeriodic>({16, 16, 16}, {1, 1, 1}, true);
    }
    TEST(TimeStepLong, AnisotropicBounded) {
        stability<StandardSolverPeriodic>({20, 24, 16}, {1, 0.96, 1}, true);
    }
    TEST(TimeStepLong, NonStandardBounded) {
        stability<NonStandardSolverPeriodic>({16, 16, 24}, {1, 1, 1}, true);
    }
}  // namespace
