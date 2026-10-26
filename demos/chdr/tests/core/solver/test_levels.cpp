/** @file test_levels.cpp
 * @brief Compare gauge-wave reconstructions with zero fields and the exact R1 electric oracle.
 * @details GaugeWave supplies the exact discrete potentials and the R1 closed
 * form, derived from [christlieb2024GaugeConserving, section 3.1, proof of
 * Theorem 3.5, p. 8] and [fallahi2020mithra20fullwavesimulation,
 * Eqs. (3.17), (3.24)]. A6 labels these gauge-cancellation records.
 * No continuum approximation replaces the oracle here.
 */
#include <gtest/gtest.h>
#include <iostream>
#include <limits>

#include "CoreTestSupport.h"
#include "MonitorTestSupport.h"
#include "Waves.h"
#include "core/base/SpatialDimension.h"

static_assert(chdr::solver::SpatialDim == 3,
              "test_levels.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;
    using namespace chdr::solver::coreTest;

    /// @brief Value-only sampling metadata; it never owns or relocates a field container.
    struct GaugeGeometry {
        ippl::Vector<double, SpatialDim> spacing_m;
        ippl::Vector<int, SpatialDim> first_m;
        explicit GaugeGeometry(AJFourFieldsType& field)
            : spacing_m(field.get_mesh().getMeshSpacing()) {
            const auto local = field.getLayout().getLocalNDIndex();
            first_m          = {local[0].first(), local[1].first(), local[2].first()};
        }
    };

    /** @brief Maximum owned electric deviation from zero or the exact R1 closed form.
     * @details Coordinates and the completed integer solve count retain the original
     * sampling expression. Each component rejects nonfinite data before taking a
     * maximum; the MPI maximum is deferred until all requested steps are complete.
     */
    template <class SolverType>
    double electricError(SolverType& solver, GaugeWave wave, GaugeGeometry geometry) {
        const auto e                          = solver.electricField().getView();
        const auto h                          = geometry.spacing_m;
        const auto first                      = geometry.first_m;
        const double t                        = solver.completedSolves() * solver.getDt();
        const auto ebReconstructionConvention = solver.ebReconstructionConvention();
        const double infinity                 = std::numeric_limits<double>::infinity();
        return localMaximum(
            solver.A_n, KOKKOS_LAMBDA(int i, int j, int k) {
                const double x = (i + first[0] - 0.5) * h[0], y = (j + first[1] - 0.5) * h[1],
                             z = (k + first[2] - 0.5) * h[2];
                double maximum = 0;
                for (unsigned c = 0; c < SpatialDim; ++c) {
                    const double reference =
                        ebReconstructionConvention == EBReconstructionConvention::R1
                            ? wave.r1Electric(x, y, z, t, c)
                            : 0;
                    const double value = Kokkos::abs(e(i, j, k)[c] - reference);
                    maximum = Kokkos::max(maximum, Kokkos::isfinite(value) ? value : infinity);
                }
                return maximum;
            });
    }

    /// @brief Gauge-wave magnetic oracle is zero; inspect each component before reducing.
    template <class SolverType>
    double magneticError(SolverType& solver) {
        const auto b          = solver.magneticField().getView();
        const double infinity = std::numeric_limits<double>::infinity();
        return localMaximum(
            solver.A_n, KOKKOS_LAMBDA(int i, int j, int k) {
                double maximum = 0;
                for (unsigned c = 0; c < SpatialDim; ++c) {
                    const double value = Kokkos::abs(b(i, j, k)[c]);
                    maximum = Kokkos::max(maximum, Kokkos::isfinite(value) ? value : infinity);
                }
                return maximum;
            });
    }

    /** @brief Sample the absolute scalar-gradient reference at the phi slot used by E.
     * @details R2 uses previous phi; the other kernels use current phi. This is the
     * actual maximum of the sampled \f$|D_c\phi|\f$, rather than its continuous
     * amplitude envelope. The expression retains its spatial phase and scalar time label
     * explicitly. Source for labels: LevelContract; gradient symbol: Operators.h.
     */
    template <class SolverType>
    double gradientScale(SolverType& solver, GaugeWave wave, GaugeGeometry geometry,
                         double scalar) {
        const auto h     = geometry.spacing_m;
        const auto first = geometry.first_m;
        const double phiTime =
            solver.labelTime(Quantity::Phi, solver.completedSolves())
            - (solver.ebReconstructionConvention() == EBReconstructionConvention::R2
                   ? solver.getDt()
                   : 0);
        return localMaximum(
            solver.A_n, KOKKOS_LAMBDA(int i, int j, int k) {
                const double phase = wave.waveNumber_m[0] * (i + first[0] - 0.5) * h[0]
                                     + wave.waveNumber_m[1] * (j + first[1] - 0.5) * h[1]
                                     + wave.waveNumber_m[2] * (k + first[2] - 0.5) * h[2]
                                     - wave.omega_m * phiTime;
                double maximum = 0;
                for (unsigned c = 0; c < SpatialDim; ++c)
                    maximum = Kokkos::max(maximum,
                                          Kokkos::abs(scalar * wave.g_m[c] * Kokkos::cos(phase)));
                return maximum;
            });
    }

    /** @brief Record the axial endpoint's signed global Az Nyquist projection.
     * @details Only the explicitly declared \f$S=1\f$ record calls this helper. It adds
     * no amplitude threshold. The global z parity is independent of rank seams;
     * MPI_SUM operates on the original owned-cell projection every 64 steps.
     */
    template <class SolverType>
    void recordNyquist(SolverType& solver, GaugeGeometry geometry,
                       std::array<int, SpatialDim> cells, int step) {
        const auto potential = solver.A_n.getView();
        const auto first     = geometry.first_m;
        double nyquist       = 0;
        Kokkos::parallel_reduce(
            solver.A_n.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& sum) {
                sum += ((k + first[2] - 1) % 2 ? -1.0 : 1.0) * potential(i, j, k)[3];
            },
            nyquist);
        MPI_Allreduce(MPI_IN_PLACE, &nyquist, 1, MPI_DOUBLE, MPI_SUM, solver.A_n.getLayout().comm);
        if (ippl::Comm->rank() == 0)
            std::cout << "A6 axial S=1 record n=" << step
                      << " reconstruction=" << static_cast<int>(solver.ebReconstructionConvention())
                      << " Az_Nyquist=" << nyquist / (cells[0] * cells[1] * cells[2]) << '\n';
    }

    /// @brief Maxima retained over all solves, before the three final MPI reductions.
    struct GaugeErrors {
        double electric_m = 0, magnetic_m = 0, scale_m = 0;
        void observe(double electric, double magnetic, double scale) {
            electric_m = std::max(electric_m, electric);
            magnetic_m = std::max(magnetic_m, magnetic);
            scale_m    = std::max(scale_m, scale);
        }
    };

    /// @brief Advance the unchanged solver and retain each step's three independent observations.
    template <class SolverType>
    GaugeErrors evolveGauge(SolverType& solver, GaugeWave wave, std::array<int, SpatialDim> cells,
                            int steps, bool record) {
        const GaugeGeometry geometry(solver.A_n);
        const double scalar = solver.ebReconstructionConvention() == EBReconstructionConvention::R1
                                  ? wave.omega_m
                                  : timeSymbol(wave.omega_m, solver.getDt());
        GaugeErrors errors;
        Monitors monitors(solver.A_n);
        for (int step = 1; step <= steps; ++step) {
            chdr::solver::test::advanceWithFinalMonitor(solver, monitors, step == steps);
            const double electric = electricError(solver, wave, geometry);
            const double magnetic = magneticError(solver);
            const double scale    = gradientScale(solver, wave, geometry, scalar);
            errors.observe(electric, magnetic, scale);
            if (record && step % 64 == 0)
                recordNyquist(solver, geometry, cells, step);
        }
        return errors;
    }

    /// @brief Apply normalized accuracy limits, or report the designated boundary case only.
    template <class SolverType>
    void expectGauge(SolverType& solver, GaugeErrors errors, int steps, bool record) {
        errors.electric_m = globalMaximum(solver.A_n, errors.electric_m);
        errors.magnetic_m = globalMaximum(solver.A_n, errors.magnetic_m);
        errors.scale_m    = globalMaximum(solver.A_n, errors.scale_m);
        EXPECT_TRUE(std::isfinite(errors.electric_m));
        EXPECT_TRUE(std::isfinite(errors.magnetic_m));
        ASSERT_TRUE(std::isfinite(errors.scale_m));
        ASSERT_GT(errors.scale_m, 0.0);
        if (ippl::Comm->rank() == 0)
            std::cout << "A6 reconstruction="
                      << static_cast<int>(solver.ebReconstructionConvention()) << " steps=" << steps
                      << " E=" << errors.electric_m / errors.scale_m
                      << " B=" << errors.magnetic_m / errors.scale_m << " record=" << record
                      << '\n';
        if (!record) {
            EXPECT_LE(errors.electric_m / errors.scale_m, 1e-12);
            EXPECT_LE(errors.magnetic_m / errors.scale_m, 1e-12);
        }
    }

    /// @brief Keep the owning fields stationary while the solver borrows their mesh/layout.
    template <class PotentialSolverBase>
    void oneGaugeCase(std::array<int, SpatialDim> cells, std::array<double, SpatialDim> lengths,
                      ippl::Vector<double, SpatialDim> waveNumber, int steps, double safety,
                      bool record, EBReconstructionConvention ebReconstructionConvention) {
        auto fields = makeBoxFieldContainer(cells, lengths);
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               ebReconstructionConvention);
        if (safety > 0)
            solver.setSafetyFactor(safety);
        const GaugeWave wave(waveNumber, fields.getMesh().getMeshSpacing(), solver.getDt(),
                             solver.stencil(), ebReconstructionConvention);
        wave.requirePeriodic(lengths);
        loadHistory(solver, wave);
        if (steps == 0)
            steps = static_cast<int>(
                std::ceil(20 * Kokkos::numbers::pi_v<double> / (wave.omega_m * solver.getDt())));
        const auto errors = evolveGauge(solver, wave, cells, steps, record);
        expectGauge(solver, errors, steps, record);
    }

    /// @brief Reuse exactly the same frozen grid/mode/run length for each reconstruction.
    template <class PotentialSolverBase>
    void gaugeCase(std::array<int, SpatialDim> cells, std::array<double, SpatialDim> lengths,
                   ippl::Vector<double, SpatialDim> waveNumber, int steps, double safety = 0,
                   bool record = false) {
        for (auto ebReconstructionConvention :
             {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
              EBReconstructionConvention::R2Prime})
            oneGaugeCase<PotentialSolverBase>(cells, lengths, waveNumber, steps, safety, record,
                                              ebReconstructionConvention);
    }

    TEST(Levels, EarlyAxial) {
        gaugeCase<StandardSolverPeriodic>(
            {8, 8, 64}, {0.125, 0.125, 1},
            ippl::Vector<double, SpatialDim>(0, 0, 2 * Kokkos::numbers::pi_v<double>), 128);
    }
    TEST(Levels, EarlyOblique) {
        gaugeCase<StandardSolverPeriodic>(
            {8, 8, 64}, {0.125, 0.125, 1},
            2 * Kokkos::numbers::pi_v<double> * ippl::Vector<double, SpatialDim>(8, 8, 1), 128);
    }
    TEST(GaugeWave, ObliqueStandard) {
        gaugeCase<StandardSolverPeriodic>(
            {32, 32, 32}, {1, 1, 1},
            2 * Kokkos::numbers::pi_v<double> * ippl::Vector<double, SpatialDim>(1, 2, 2), 0);
    }
    TEST(GaugeWave, ObliqueNonStandard) {
        gaugeCase<NonStandardSolverPeriodic>(
            {32, 32, 48}, {1, 1, 1},
            2 * Kokkos::numbers::pi_v<double> * ippl::Vector<double, SpatialDim>(1, 2, 2), 0);
    }
    TEST(GaugeWave, AxialNonStandardStable) {
        gaugeCase<NonStandardSolverPeriodic>(
            {8, 8, 64}, {0.1875, 0.1875, 1},
            ippl::Vector<double, SpatialDim>(0, 0, 2 * Kokkos::numbers::pi_v<double>), 640, 0.98);
    }
    TEST(GaugeWave, AxialNonStandardBoundaryRecord) {
        gaugeCase<NonStandardSolverPeriodic>(
            {8, 8, 64}, {0.1875, 0.1875, 1},
            ippl::Vector<double, SpatialDim>(0, 0, 2 * Kokkos::numbers::pi_v<double>), 640, 1.0,
            true);
    }
}  // namespace
