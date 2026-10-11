/** @file test_staggered.cpp
 * @brief Staggered difference composition, component positions and Standard staggered-R2
 * identities.
 * @details Cartesian identities specialize [clemens2001FiniteIntegration,
 * Eqs. (18)-(24)]. CompactGaugeWave defines the signed space/time differences
 * used to construct the discrete counterpart of a purely longitudinal potential
 * with zero E and B [chew2014generalizedgauge, section 7, Eqs. (76)-(80)] (deduction).
 * GaugeMaximum uses the centred scalar-gradient magnitude as the common reference
 * scale for collocated and staggered gauge-wave tests; runGaugeSteps accumulates
 * the sampled error and reference maxima across each specified run.
 *
 * TODO: REVIEW.
 */
#include <Kokkos_MathematicalConstants.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <limits>

#include "MonitorTestSupport.h"
#include "Monitors.h"
#include "SpatialDimension.h"
#include "StaggeredOperators.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_staggered.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;
    constexpr double TWO_PI = 2.0 * Kokkos::numbers::pi_v<double>;

    void fillStaggeredNoise(AJFourFieldsType& field, unsigned seed) {
        const auto view  = field.getView();
        const auto local = field.getLayout().getLocalNDIndex();
        const int firstX = local[0].first() - 1, firstY = local[1].first() - 1,
                  firstZ = local[2].first() - 1;
        Kokkos::parallel_for(
            "Compact random Nyquist fixture", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                const int x = i + firstX, y = j + firstY, z = k + firstZ;
                for (unsigned c = 0; c < 4; ++c) {
                    const unsigned hash = unsigned(x) * 1103515245u + unsigned(y) * 12345u
                                          + unsigned(z) * 2654435761u + c * 374761393u + seed;
                    view(i, j, k)[c] = double(hash % 1021) / 1021.0 + ((x + y + z) & 1 ? -0.2 : 0.2)
                                       + (x & 1 ? -0.3 : 0.3);
                }
            });
        field.fillHalo();
    }

    /// Compare two staggered difference passes with the separate Standard Laplacian on grid-scale
    /// data.
    template <class SolverType>
    void checkCompactComposition(SolverType& solver) {
        const auto view    = solver.A_n.getView();
        const auto spacing = solver.A_n.get_mesh().getMeshSpacing();
        monitorDetail::Workspace gradient(solver.sourceField());
        gradient.fill(KOKKOS_LAMBDA(int i, int j, int k, unsigned d) {
            return monitorDetail::derivative(monitorDetail::Component{view, 0u}, i, j, k, d,
                                             spacing[d], 1);
        });
        const auto x = gradient.component(0), y = gradient.component(1), z = gradient.component(2);
        const auto composition = reduceResidual(
            solver.A_n, KOKKOS_LAMBDA(int i, int j, int k) {
                return monitorDetail::derivative(x, i, j, k, 0, spacing[0], -1)
                       + monitorDetail::derivative(y, i, j, k, 1, spacing[1], -1)
                       + monitorDetail::derivative(z, i, j, k, 2, spacing[2], -1)
                       - monitorDetail::laplace<Stencil::Standard>(
                           monitorDetail::Component{view, 0u}, i, j, k, spacing);
            });
        EXPECT_LE(composition.cell_m, 1e-12);
    }

    void checkStaggeredReport(const MonitorReport& result) {
        EXPECT_LE(result.divMagnetic_m.cell_m, 1e-12);
        EXPECT_LE(result.divRawMagnetic_m.cell_m, 1e-12);
        EXPECT_LE(result.gauss_m.cell_m, 1e-12);
        for (const auto& component : result.faraday_m)
            EXPECT_LE(component.cell_m, 1e-12);
    }

    TEST(ChdrStaggered, RandomCompositionAndIdentities) {
        auto fields = makePeriodicFieldContainer({12, 10, 8}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2,
                                                  SpatialDiscretization::Staggered);
        fillStaggeredNoise(solver.A_n, 19);
        fillStaggeredNoise(solver.A_nm1, 71);
        checkCompactComposition(solver);
        Monitors monitors(fields.getJ());
        for (unsigned step = 0; step < 6; ++step) {
            fillStaggeredNoise(fields.getJ(), 117 + step);
            checkStaggeredReport(monitors.advance(solver));
        }
    }

    struct LinearSamples {
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double time,
                                                 unsigned c) const {
            return x + 2.0 * y + 3.0 * z + 4.0 * time + 10.0 * c;
        }
    };

    /// Independent position/time oracle for the fixed eight-cell sampling fixture.
    struct LinearSampleError {
        AJFourFieldsType::view_type potential_m, source_m;
        ippl::Vector<int, SpatialDim> first_m;
        double dt_m;
        KOKKOS_INLINE_FUNCTION monitorDetail::Expression operator()(int i, int j, int k) const {
            monitorDetail::Expression error;
            for (unsigned c = 0; c < 4; ++c) {
                const double x        = (i + first_m[0] + 0.5 + (c == 1 ? 0.5 : 0.0)) / 8.0;
                const double y        = (j + first_m[1] + 0.5 + (c == 2 ? 0.5 : 0.0)) / 8.0;
                const double z        = (k + first_m[2] + 0.5 + (c == 3 ? 0.5 : 0.0)) / 8.0;
                const double expected = LinearSamples{}(x, y, z, c == 0 ? 0.5 * dt_m : 0.0, c);
                error.value_m += Kokkos::abs(potential_m(i, j, k)[c] - expected)
                                 + Kokkos::abs(source_m(i, j, k)[c] - expected);
                error.scale_m += 2.0 * Kokkos::abs(expected);
            }
            return error;
        }
    };

    template <class SolverType>
    void checkSampledCoordinates(SolverType& solver) {
        const auto local = solver.A_n.getLayout().getLocalNDIndex();
        const ippl::Vector<int, SpatialDim> first(local[0].first() - 1, local[1].first() - 1,
                                                  local[2].first() - 1);
        const LinearSampleError error{solver.A_n.getView(), solver.sourceField().getView(), first,
                                      solver.getDt()};
        EXPECT_LE(reduceResidual(solver.A_n, error).cell_m, 1e-15);
    }

    void checkOutputOffsets(LevelContract contract) {
        for (unsigned d = 0; d < SpatialDim; ++d) {
            const auto electric = contract.positionOffset(Quantity::Electric, d);
            const auto magnetic = contract.positionOffset(Quantity::Magnetic, d);
            for (unsigned a = 0; a < SpatialDim; ++a) {
                EXPECT_EQ(electric[a], d == a ? 0.5 : 0.0);
                EXPECT_EQ(magnetic[a], d == a ? 0.0 : 0.5);
            }
        }
    }

    TEST(ChdrStaggered, SpatialOffsetsReachHistoryAndSource) {
        auto fields = makePeriodicFieldContainer({8, 8, 8}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2,
                                                  SpatialDiscretization::Staggered);
        loadHistory(solver, LinearSamples{});
        writeSource(solver, LinearSamples{}, 0);
        checkSampledCoordinates(solver);
        checkOutputOffsets(solver.levelContract());
        EXPECT_THROW((ChdrSolver<StandardSolverPeriodic>(
                         fields.getJ(), fields.getE(), fields.getB(),
                         EBReconstructionConvention::R1, SpatialDiscretization::Staggered)),
                     IpplException);
        EXPECT_THROW((ChdrSolver<StandardSolverPeriodic>(
                         fields.getJ(), fields.getE(), fields.getB(),
                         EBReconstructionConvention::R2Prime, SpatialDiscretization::Staggered)),
                     IpplException);
    }

    /// Explicit backward face differences provide the source oracle independently of the monitor.
    void setCompactConservingCharge(AJFourFieldsType& source, AJFourFieldsType& previous, double dt,
                                    double planted) {
        const auto now = source.getView(), old = previous.getView();
        const auto spacing = source.get_mesh().getMeshSpacing();
        Kokkos::parallel_for(
            "Compact conserving charge", source.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                double divergence = 0.0;
                for (unsigned d = 0; d < SpatialDim; ++d) {
                    const int di = d == 0, dj = d == 1, dk = d == 2;
                    divergence +=
                        (now(i, j, k)[d + 1] - now(i - di, j - dj, k - dk)[d + 1]) / spacing[d];
                }
                now(i, j, k)[0] = old(i, j, k)[0] - dt * divergence + dt * planted;
            });
    }

    TEST(ChdrStaggered, ConservingAndPlantedContinuity) {
        auto fields = makePeriodicFieldContainer({12, 10, 8}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        AJFourFieldsType previous;
        previous.initialize(fields.getMesh(), fields.getFL());
        constexpr double DT = 0.0125;
        for (double planted : {0.0, 0.125}) {
            fillStaggeredNoise(fields.getJ(), 107);
            fillStaggeredNoise(previous, 211);
            setCompactConservingCharge(fields.getJ(), previous, DT, planted);
            const auto result = Monitors::continuity(fields.getJ(), previous, DT, 1.0,
                                                     SpatialDiscretization::Staggered);
            EXPECT_LE(
                chdr::solver::test::expectedContinuityMismatch(
                    fields.getJ(), previous, DT, 1.0, planted, SpatialDiscretization::Staggered)
                    .cell_m,
                1e-12);
            if (planted == 0.0) {
                EXPECT_LE(result.cell_m, 1e-12);
            }
            EXPECT_NEAR(result.maxAbs_m, std::abs(planted), 1e-12 * result.maxScale_m);
        }
    }

    TEST(ChdrStaggered, NonStandardReadingRefused) {
        auto fields = makePeriodicFieldContainer({8, 8, 12}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        EXPECT_THROW((ChdrSolver<NonStandardSolverPeriodic>(
                         fields.getJ(), fields.getE(), fields.getB(),
                         EBReconstructionConvention::R2, SpatialDiscretization::Staggered)),
                     IpplException);
    }

}  // namespace
