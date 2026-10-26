/** @file test_monitors.cpp
 * @brief Magnetic-divergence, Faraday and continuity checks with Nyquist content.
 * @details The identity and time-weight pairings are discrete deductions from
 * [christlieb2024GaugeConserving, Lemma 3.3, Theorem 3.5] and
 * [clemens2001FiniteIntegration, Eqs. (18)-(24)]. Fixtures use global integer
 * indices so decomposition cannot change the input bit patterns.
 *
 * TODO: REVIEW.
 */
#include <cstdint>
#include <gtest/gtest.h>
#include <iostream>

#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "opt-in/diagnostics/Monitors.h"
#include "tests/core/solver/MonitorTestSupport.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_monitors.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;

    KOKKOS_INLINE_FUNCTION double noiseValue(unsigned x, unsigned y, unsigned z, unsigned c,
                                             unsigned seed) {
        std::uint32_t bits = seed + 73856093u * x + 19349663u * y + 83492791u * z + 2654435761u * c;
        bits ^= bits >> 16;
        bits *= 0x7feb352du;
        bits ^= bits >> 15;
        const double random  = double(bits & 0xffffu) / 65535.0 - 0.5;
        const double checker = ((x + y + z) & 1u) ? -0.31 : 0.31;
        const double plane   = (x & 1u) ? -0.23 : 0.23;
        return random + checker + plane;
    }

    /// Preserve random, checkerboard and single-Nyquist-plane inputs across rank counts.
    void fillNoise(AJFourFieldsType& field, unsigned seed) {
        const auto view  = field.getView();
        const auto local = field.getLayout().getLocalNDIndex();
        const int firstX = local[0].first() - field.getNghost();
        const int firstY = local[1].first() - field.getNghost();
        const int firstZ = local[2].first() - field.getNghost();
        Kokkos::parallel_for(
            "Monitor random and Nyquist data", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                for (unsigned c = 0; c < 4; ++c)
                    view(i, j, k)[c] = noiseValue(i + firstX, j + firstY, k + firstZ, c, seed);
            });
        Kokkos::fence();
        field.fillHalo();
    }

    double checkIdentityReport(const MonitorReport& result, unsigned step) {
        EXPECT_LE(result.divMagnetic_m.cell_m, 1e-12);
        EXPECT_LE(result.divRawMagnetic_m.cell_m, 1e-12);
        double maximum = std::max(result.divMagnetic_m.cell_m, result.divRawMagnetic_m.cell_m);
        for (const auto& component : result.faraday_m) {
            EXPECT_LE(component.cell_m, 1e-12);
            maximum = std::max(maximum, component.cell_m);
        }
        EXPECT_EQ(result.continuityAvailable_m, step != 0);
        EXPECT_TRUE(std::isfinite(result.gauge_m.l2_m));
        return maximum;
    }

    /// Perturb actual output E to prove the monitor is sensitive to reconstruction errors.
    template <class SolverType>
    void checkElectricPerturbation(SolverType& solver, Monitors& monitors) {
        const auto electric = solver.electricField().getView();
        Kokkos::parallel_for(
            "Plant output error", solver.electricField().getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                electric(i, j, k)[1] += (i & 1) ? 0.0 : double(k);
            });
        const auto altered = monitors.evaluate(solver);
        EXPECT_GT(std::max(altered.faraday_m[0].cell_m, altered.faraday_m[2].cell_m), 1e-6);
    }

    template <class PotentialSolverBase>
    void checkReconstruction(std::array<int, SpatialDim> cells,
                             EBReconstructionConvention ebReconstructionConvention) {
        auto fields = makeBoxFieldContainer(cells, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               ebReconstructionConvention);
        fillNoise(solver.A_n, 19);
        fillNoise(solver.A_nm1, 41);
        Monitors monitors(fields.getJ());
        double maximum = 0.0;
        for (unsigned step = 0; step < 6; ++step) {
            fillNoise(fields.getJ(), 71 + step);
            maximum = std::max(maximum, checkIdentityReport(monitors.advance(solver), step));
        }
        // The "A2/A3" tag identifies magnetic-divergence and Faraday residual checks.
        if (ippl::Comm->rank() == 0)
            std::cout << "A2/A3 stencil=" << int(solver.stencil())
                      << " reconstruction=" << int(ebReconstructionConvention)
                      << " maximum=" << maximum << " threshold=1e-12\n";
        checkElectricPerturbation(solver, monitors);
    }

    template <class PotentialSolverBase>
    void checkBasicIdentities(std::array<int, SpatialDim> cells) {
        for (auto r : {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
                       EBReconstructionConvention::R2Prime})
            checkReconstruction<PotentialSolverBase>(cells, r);
    }
    TEST(ChdrMonitors, RandomNyquistStandard) {
        checkBasicIdentities<StandardSolverPeriodic>({12, 10, 8});
    }
    TEST(ChdrMonitors, RandomNyquistNonStandard) {
        checkBasicIdentities<NonStandardSolverPeriodic>({8, 8, 12});
    }

    /** @brief Independent source oracle: choose charge so the declared continuity is planted.
     * @details Explicit neighboring samples, rather than production operator calls,
     * construct the source. The device kernel is kept together so its three-axis
     * sum and arithmetic grouping remain visible. The resulting continuity
     * expression equals the specified planted scalar.
     */
    void setConservingCharge(AJFourFieldsType& source, AJFourFieldsType& previous, double dt,
                             double weight, double planted) {
        const auto now = source.getView(), old = previous.getView();
        const auto spacing = source.get_mesh().getMeshSpacing();
        Kokkos::parallel_for(
            "Construct conserving charge", source.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                double div = 0.0;
                for (unsigned d = 0; d < SpatialDim; ++d) {
                    const int di = d == 0, dj = d == 1, dk = d == 2;
                    const double currentDifference =
                        now(i + di, j + dj, k + dk)[d + 1] - now(i - di, j - dj, k - dk)[d + 1];
                    const double previousDifference =
                        old(i + di, j + dj, k + dk)[d + 1] - old(i - di, j - dj, k - dk)[d + 1];
                    div += (weight * currentDifference + (1.0 - weight) * previousDifference)
                           / (2.0 * spacing[d]);
                }
                now(i, j, k)[0] = old(i, j, k)[0] - dt * div + dt * planted;
            });
    }

    TEST(ChdrMonitors, ZeroNormalizationAndOrderingGuard) {
        auto fields = makeBoxFieldContainer({8, 8, 8}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2);
        Monitors monitors(fields.getJ());
        EXPECT_THROW(monitors.afterSolve(solver), std::logic_error);
        monitors.seedPreviousSource(fields.getJ());
        const auto report = monitors.advance(solver);
        EXPECT_EQ(report.divMagnetic_m.cell_m, 0.0);
        EXPECT_EQ(report.divRawMagnetic_m.global_m, 0.0);
        EXPECT_EQ(report.gauge_m.cell_m, 0.0);
        EXPECT_TRUE(report.continuityAvailable_m);
        EXPECT_EQ(report.continuity_m.cell_m, 0.0);
        for (const auto& value : report.faraday_m)
            EXPECT_EQ(value.cell_m, 0.0);
    }

    /// Exercise the retained-source wiring with independently constructed conserved charge.
    void checkLifecycleContinuity(EBReconstructionConvention ebReconstructionConvention) {
        auto fields = makeBoxFieldContainer({12, 10, 8}, {1, 1, 1});
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  ebReconstructionConvention);
        AJFourFieldsType previous;
        previous.initialize(fields.getMesh(), fields.getFL());
        fillNoise(previous, 211);
        Monitors monitors(fields.getJ());
        monitors.seedPreviousSource(previous);
        const double dt = solver.getDt(), weight = solver.levelContract().gaugeWeight();
        for (double planted : {0.0, 0.125, -0.125}) {
            fillNoise(fields.getJ(), 107);
            setConservingCharge(fields.getJ(), previous, dt, weight, planted);
            const auto report = monitors.advance(solver);
            EXPECT_TRUE(report.continuityAvailable_m);
            if (planted == 0.0) {
                EXPECT_LE(report.continuity_m.cell_m, 1e-12);
            }
            EXPECT_NEAR(report.continuity_m.maxAbs_m, std::abs(planted),
                        1e-12 * report.continuity_m.maxScale_m);
            EXPECT_NEAR(report.continuity_m.l2_m, std::abs(planted) * std::sqrt(12.0 * 10.0 * 8.0),
                        1e-12 * report.continuity_m.scaleL2_m);
            EXPECT_LE(chdr::solver::test::expectedContinuityMismatch(monitors.consumedSource(),
                                                                     previous, dt, weight, planted)
                          .cell_m,
                      1e-12);
            Kokkos::deep_copy(previous.getView(), monitors.consumedSource().getView());
        }
    }

    TEST(ChdrMonitors, LifecycleConservingAndPlantedContinuity) {
        for (auto ebReconstructionConvention :
             {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
              EBReconstructionConvention::R2Prime})
            checkLifecycleContinuity(ebReconstructionConvention);
    }

    TEST(ChdrMonitors, SparseAndReloadedSourcesAreUnavailable) {
        auto fields = makeBoxFieldContainer({8, 8, 8}, {1, 1, 1});
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2);
        Monitors monitors(fields.getJ());
        EXPECT_FALSE(monitors.advance(solver).continuityAvailable_m);
        EXPECT_TRUE(monitors.advance(solver).continuityAvailable_m);
        solver.solve();
        EXPECT_FALSE(monitors.advance(solver).continuityAvailable_m);
        EXPECT_TRUE(monitors.advance(solver).continuityAvailable_m);
        solver.initialize();
        EXPECT_THROW(monitors.evaluate(solver), std::logic_error);
        EXPECT_FALSE(monitors.advance(solver).continuityAvailable_m);
        solver.initialize();
        solver.solve();
        EXPECT_FALSE(monitors.advance(solver).continuityAvailable_m);
        solver.solve();
        monitors.seedPreviousSource(fields.getJ());
        EXPECT_TRUE(monitors.advance(solver).continuityAvailable_m);
    }

}  // namespace
