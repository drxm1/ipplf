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
        auto fields = makeBoxFieldContainer({12, 10, 8}, {1.0, 1.0, 1.0}, 0);
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
        auto fields = makeBoxFieldContainer({8, 8, 8}, {1.0, 1.0, 1.0}, 0);
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
        auto fields = makeBoxFieldContainer({12, 10, 8}, {1.0, 1.0, 1.0}, 0);
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
        auto fields = makeBoxFieldContainer({8, 8, 12}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        EXPECT_THROW((ChdrSolver<NonStandardSolverPeriodic>(
                         fields.getJ(), fields.getE(), fields.getB(),
                         EBReconstructionConvention::R2, SpatialDiscretization::Staggered)),
                     IpplException);
    }

    /** @brief Gauge potentials generated from a scalar discrete wave. This is the version used in
     * the staggered discretization.
     * @details Let \f$\chi[k]\f$ be the scalar wave sampled at cell centres at
     * \f$t_k\f$. For EBReconstructionConvention::R2, the fixture uses
     * \f[ A[k]=D_+\chi[k],\qquad
     * \phi[k]=-(\chi[k+1]-\chi[k])/dt. \f]
     * A is sampled on its component faces. The scalar forward time difference
     * has the half-step time assigned to phi by LevelContract. EBReconstructionConvention
     * reads the previous scalar slot, whose gradient cancels the negative
     * vector-potential time difference in E; the staggered curl of A vanishes.
     * This is a discrete deduction from [fallahi2020mithra20fullwavesimulation,
     * Eqs. (3.8)-(3.9)] and [clemens2001FiniteIntegration, Eqs. (18)-(24)].
     */
    struct CompactGaugeWave {
        ippl::Vector<double, SpatialDim> waveNumber_m, compact_m;
        double omega_m, temporal_m;
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double time,
                                                 unsigned c) const {
            const double phase =
                waveNumber_m[0] * x + waveNumber_m[1] * y + waveNumber_m[2] * z - omega_m * time;
            return -(c == 0 ? temporal_m : compact_m[c - 1]) * Kokkos::sin(phase);
        }
    };

    template <class SolverType>
    CompactGaugeWave makeGaugeWave(SolverType& solver,
                                   ippl::Vector<double, SpatialDim> waveNumber) {
        const auto spacing = solver.A_n.get_mesh().getMeshSpacing();
        const double omega =
            discreteFrequency(waveNumber, spacing, solver.getDt(), Stencil::Standard);
        ippl::Vector<double, SpatialDim> compact;
        for (unsigned d = 0; d < SpatialDim; ++d)
            compact[d] = staggeredDerivativeSymbol(waveNumber[d], spacing[d]);
        return {waveNumber, compact, omega, timeSymbol(omega, solver.getDt())};
    }

    /// Reduce sampled E/B error and the centered-gradient reference scale together.
    struct GaugeMaximum {
        EMFieldsType::view_type electric_m, magnetic_m;
        CompactGaugeWave wave_m;
        ippl::Vector<double, SpatialDim> spacing_m, centered_m;
        ippl::Vector<int, SpatialDim> first_m;
        double time_m                          = 0.0;
        static constexpr double INFINITY_VALUE = std::numeric_limits<double>::infinity();
        KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& maximum,
                                               double& scale) const {
            const double phase = wave_m.waveNumber_m[0] * (i + first_m[0] + 0.5) * spacing_m[0]
                                 + wave_m.waveNumber_m[1] * (j + first_m[1] + 0.5) * spacing_m[1]
                                 + wave_m.waveNumber_m[2] * (k + first_m[2] + 0.5) * spacing_m[2]
                                 - wave_m.omega_m * time_m;
            for (unsigned d = 0; d < SpatialDim; ++d) {
                const double electric = electric_m(i, j, k)[d], magnetic = magnetic_m(i, j, k)[d];
                maximum = Kokkos::max(
                    maximum, Kokkos::isfinite(electric) ? Kokkos::abs(electric) : INFINITY_VALUE);
                maximum = Kokkos::max(
                    maximum, Kokkos::isfinite(magnetic) ? Kokkos::abs(magnetic) : INFINITY_VALUE);
                scale = Kokkos::max(
                    scale, Kokkos::abs(centered_m[d] * wave_m.temporal_m * Kokkos::cos(phase)));
            }
        }
    };

    template <class SolverType>
    GaugeMaximum gaugeMaximum(SolverType& solver, CompactGaugeWave wave) {
        const auto spacing = solver.A_n.get_mesh().getMeshSpacing();
        const auto local   = solver.A_n.getLayout().getLocalNDIndex();
        ippl::Vector<double, SpatialDim> centered;
        for (unsigned d = 0; d < SpatialDim; ++d)
            centered[d] = centeredSymbol(wave.waveNumber_m[d], spacing[d]);
        return {solver.electricField().getView(),
                solver.magneticField().getView(),
                wave,
                spacing,
                centered,
                ippl::Vector<int, SpatialDim>(local[0].first() - 1, local[1].first() - 1,
                                              local[2].first() - 1)};
    }

    template <class SolverType>
    double runGaugeSteps(SolverType& solver, CompactGaugeWave wave, unsigned steps) {
        auto sample  = gaugeMaximum(solver, wave);
        double error = 0.0, reference = 0.0;
        Monitors monitors(solver.sourceField());
        for (unsigned step = 0; step < steps; ++step) {
            chdr::solver::test::advanceWithFinalMonitor(solver, monitors, step + 1 == steps);
            sample.time_m    = solver.labelTime(Quantity::Electric, solver.completedSolves());
            double stepError = 0.0, stepReference = 0.0;
            Kokkos::parallel_reduce("R2s gauge sampled maximum", solver.A_n.getFieldRangePolicy(),
                                    sample, Kokkos::Max<double>(stepError),
                                    Kokkos::Max<double>(stepReference));
            error     = std::max(error, stepError);
            reference = std::max(reference, stepReference);
        }
        MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, solver.A_n.getLayout().comm);
        MPI_Allreduce(MPI_IN_PLACE, &reference, 1, MPI_DOUBLE, MPI_MAX,
                      solver.A_n.getLayout().comm);
        return error / reference;
    }

    void checkGaugeWave(std::array<int, SpatialDim> cells, std::array<double, SpatialDim> lengths,
                        ippl::Vector<double, SpatialDim> waveNumber, unsigned fixedSteps) {
        auto fields = makeBoxFieldContainer(cells, lengths, 0);
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2,
                                                  SpatialDiscretization::Staggered);
        const auto wave = makeGaugeWave(solver, waveNumber);
        loadHistory(solver, wave);
        const unsigned steps =
            fixedSteps ? fixedSteps
                       : unsigned(std::ceil(10.0 * TWO_PI / (wave.omega_m * solver.getDt())));
        const double error = runGaugeSteps(solver, wave, steps);
        EXPECT_LE(error, 1e-12);
        if (ippl::Comm->rank() == 0)
            std::cout << "R2s gauge max E/B = " << error << '\n';
    }

    TEST(ChdrStaggered, GaugeAxial) {
        checkGaugeWave({8, 8, 64}, {0.125, 0.125, 1.0}, {0.0, 0.0, TWO_PI}, 128);
    }
    TEST(ChdrStaggered, GaugeObliqueSmall) {
        checkGaugeWave({8, 8, 64}, {0.125, 0.125, 1.0}, {8 * TWO_PI, 8 * TWO_PI, TWO_PI}, 128);
    }
    TEST(ChdrStaggered, GaugeObliqueTenPeriods) {
        checkGaugeWave({32, 32, 32}, {1.0, 1.0, 1.0}, {TWO_PI, 2 * TWO_PI, 2 * TWO_PI}, 0);
    }
}  // namespace
