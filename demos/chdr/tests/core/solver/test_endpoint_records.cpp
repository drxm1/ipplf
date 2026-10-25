/** @file test_endpoint_records.cpp
 * @brief Record the stability-bound endpoint on even and odd meshes without an amplitude pass rule.
 * @details The even Standard corner and NonStandard axial mode have a double
 * root at minus one. The exact sequence from current amplitude one and previous
 * zero is \f$a_n=(-1)^n(n+1)\f$. Deduction from
 * [fallahi2020mithra20fullwavesimulation, Eqs. (3.17), (3.24)].
 * The record also follows the host symbol using its actual rounded step.
 *
 * TODO: REVIEW.
 */
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>

#include "MonitorTestSupport.h"
#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "testing/Util.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_endpoint_records.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;
    using chdr::testing::util::preparePotentialHalo;

    /** @brief Fixed spatial mode metadata; all coordinates are global integer indices.
     * @details Integer-index sampling gives nonzero amplitude at the even-grid
     * checkerboard. Choosing phase by global integer index retains that amplitude; the
     * periodic wave number still has an integer number of cycles on each axis.
     */
    struct EndpointMode {
        ippl::Vector<double, SpatialDim> theta_m, waveNumber_m;
        ippl::Vector<int, SpatialDim> first_m;
        KOKKOS_INLINE_FUNCTION double operator()(int i, int j, int k) const {
            return Kokkos::cos(theta_m[0] * (i + first_m[0] - 1) + theta_m[1] * (j + first_m[1] - 1)
                               + theta_m[2] * (k + first_m[2] - 1));
        }
    };

    /// @brief Select the highest periodic corner/axial mode on the supplied even or odd grid.
    template <class SolverType>
    EndpointMode endpointMode(SolverType& solver, std::array<int, SpatialDim> cells) {
        EndpointMode mode;
        const auto h     = solver.A_n.get_mesh().getMeshSpacing();
        const auto local = solver.A_n.getLayout().getLocalNDIndex();
        mode.first_m     = {local[0].first(), local[1].first(), local[2].first()};
        for (unsigned d = 0; d < SpatialDim; ++d) {
            const int index =
                SolverType::stencil() == Stencil::NonStandard && d < 2 ? 0 : cells[d] / 2;
            mode.theta_m[d]      = 2 * Kokkos::numbers::pi_v<double> * index / cells[d];
            mode.waveNumber_m[d] = mode.theta_m[d] / h[d];
        }
        return mode;
    }

    /// @brief Set current scalar amplitude one and previous amplitude zero, then prepare halos.
    template <class SolverType>
    void initializeEndpoint(SolverType& solver, EndpointMode shape) {
        const auto current = solver.A_n.getView();
        Kokkos::parallel_for(
            solver.A_n.getFieldRangePolicy(), KOKKOS_LAMBDA(int i, int j, int k) {
                current(i, j, k)    = 0.0;
                current(i, j, k)[0] = shape(i, j, k);
            });
        solver.A_nm1 = 0.0;
        preparePotentialHalo(solver.A_n);
        preparePotentialHalo(solver.A_nm1);
    }

    /// @brief Global squared norm of the fixed mode, with the original owned-cell sum order.
    template <class SolverType>
    double endpointNorm(SolverType& solver, EndpointMode shape) {
        double norm = 0;
        Kokkos::parallel_reduce(
            solver.A_n.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& sum) {
                const double value = shape(i, j, k);
                sum += value * value;
            },
            norm);
        MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, solver.A_n.getLayout().comm);
        return norm;
    }

    /// @brief Project the evolved scalar field onto the unchanged mode and combine all ranks.
    template <class SolverType>
    double endpointProjection(SolverType& solver, EndpointMode shape) {
        const auto current = solver.A_n.getView();
        double projection  = 0;
        Kokkos::parallel_reduce(
            solver.A_n.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& sum) {
                sum += current(i, j, k)[0] * shape(i, j, k);
            },
            projection);
        MPI_Allreduce(MPI_IN_PLACE, &projection, 1, MPI_DOUBLE, MPI_SUM,
                      solver.A_n.getLayout().comm);
        return projection;
    }

    /// @brief Describe the exact mesh and effective factor beside the scalar evolution record.
    template <class SolverType>
    void describeEndpoint(SolverType& solver, std::array<int, SpatialDim> cells) {
        if (ippl::Comm->rank() != 0)
            return;
        std::cout << "endpoint stencil="
                  << (SolverType::stencil() == Stencil::Standard ? "Standard" : "NonStandard")
                  << " grid=" << cells[0] << 'x' << cells[1] << 'x' << cells[2]
                  << " safety=" << solver.safetyFactor() << "\nn,amplitude,symbol_prediction\n";
    }

    /// @brief Report a finite observation without imposing an unfrozen amplitude tolerance.
    void recordEndpointStep(int step, double projection, double norm, double predicted) {
        EXPECT_TRUE(std::isfinite(projection));
        if (ippl::Comm->rank() == 0)
            std::cout << std::setprecision(17) << step << ',' << projection / norm << ','
                      << predicted << '\n';
    }

    /** @brief Step the unchanged solver and the independent scalar symbol recurrence together.
     * @details Every observation precedes the next solve, including the initial
     * state. The recurrence's subtraction/multiplication order and 256 updates are
     * unchanged. Sources remain zero throughout this numerical stability record.
     */
    template <class SolverType>
    void evolveEndpoint(SolverType& solver, EndpointMode shape, double norm, double coefficient) {
        double previous = 0, predicted = 1;
        Monitors monitors(solver.A_n);
        for (int step = 0; step <= 256; ++step) {
            recordEndpointStep(step, endpointProjection(solver, shape), norm, predicted);
            if (step < 256) {
                chdr::solver::test::advanceWithFinalMonitor(solver, monitors, step == 255);
                const double next = coefficient * predicted - previous;
                previous          = predicted;
                predicted         = next;
            }
        }
    }

    /// @brief Keep the owning mesh/layout stationary while the endpoint experiment runs.
    template <class PotentialSolverBase>
    void endpoint(std::array<int, SpatialDim> cells) {
        auto fields = makeBoxFieldContainer(cells, {1, 1, 1});
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               EBReconstructionConvention::R1);
        solver.setTimeStepUnchecked(solver.maxStableTimeStep());
        const auto shape = endpointMode(solver, cells);
        initializeEndpoint(solver, shape);
        const double norm = endpointNorm(solver, shape);
        const double coefficient =
            2
            + solver.getDt() * solver.getDt()
                  * laplacianSymbol(shape.waveNumber_m, fields.getMesh().getMeshSpacing(),
                                    solver.stencil());
        describeEndpoint(solver, cells);
        evolveEndpoint(solver, shape, norm, coefficient);
    }

    TEST(TimeStepRecords, EvenStandard) {
        endpoint<StandardSolverPeriodic>({16, 16, 16});
    }
    TEST(TimeStepRecords, OddStandard) {
        endpoint<StandardSolverPeriodic>({15, 15, 15});
    }
    TEST(TimeStepRecords, EvenNonStandard) {
        endpoint<NonStandardSolverPeriodic>({16, 16, 24});
    }
    TEST(TimeStepRecords, OddNonStandard) {
        endpoint<NonStandardSolverPeriodic>({15, 15, 23});
    }
}  // namespace
