/** @file MonitorTestSupport.h
 * @brief Shared final-transition checks for evolved periodic test cases.
 */
#pragma once
#include <gtest/gtest.h>

#include "core/base/SpatialDimension.h"
#include "opt-in/diagnostics/Monitors.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "MonitorTestSupport.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    /** @brief Check divergence of the reconstructed B without advancing time.
     * @details Use the output values with the expanded curl's absolute scale.
     * The outer difference direction follows the centre/face or edge/node
     * locations; commuting periodic differences give zero divergence of curl
     * [clemens2001FiniteIntegration, Eqs. (18)-(24)]. This tests reconstruction-only
     * cases while retaining the same normalization as Monitors.
     * @tparam SolverType Solver exposing both four-potential histories, B and LevelContract.
     * @param solver State whose magnetic output has already been reconstructed.
     * @param work Destination for each B component and its uncancelled curl scale.
     */
    template <class SolverType>
    void fillReconstructedMagneticWorkspace(SolverType& solver, monitorDetail::Workspace& work) {
        const auto now = solver.A_n.getView(), old = solver.A_nm1.getView();
        const auto output = solver.magneticField().getView();
        const auto h      = solver.A_n.get_mesh().getMeshSpacing();
        const int direction =
            solver.spatialDiscretization() == SpatialDiscretization::Staggered ? 1 : 0;
        const double weight =
            solver.ebReconstructionConvention() == EBReconstructionConvention::R1 ? 1 : 0.5;
        work.fill(KOKKOS_LAMBDA(int i, int j, int k, unsigned c) {
            const auto bound =
                monitorDetail::curl(now, 1, c, i, j, k, h, direction) * weight
                + monitorDetail::curl(old, 1, c, i, j, k, h, direction) * (1 - weight);
            return monitorDetail::Expression{output(i, j, k)[c], bound.scale_m};
        });
    }

    /// Check the no-solve magnetic identity after the expanded-scale workspace is current.
    template <class SolverType>
    void requireReconstructedDivB(SolverType& solver) {
        monitorDetail::Workspace work(solver.A_n);
        fillReconstructedMagneticWorkspace(solver, work);
        const auto h = solver.A_n.get_mesh().getMeshSpacing();
        const int direction =
            solver.spatialDiscretization() == SpatialDiscretization::Staggered ? 1 : 0;
        const auto x = work.component(0), y = work.component(1), z = work.component(2);
        const auto result = reduceResidual(
            solver.A_n, KOKKOS_LAMBDA(int i, int j, int k) {
                return monitorDetail::derivative(x, i, j, k, 0, h[0], direction)
                       + monitorDetail::derivative(y, i, j, k, 1, h[1], direction)
                       + monitorDetail::derivative(z, i, j, k, 2, h[2], direction);
            });
        EXPECT_LE(result.cell_m, 1e-12);
    }

    /** @brief Check a prescribed continuity residual at every owned cell.
     * @details The source fixture constructs charge with independent neighbour
     * differences. This subtracts its prescribed value before normalization,
     * so a local downward error cannot hide behind the global maximum or norm.
     * The continuity expression and time weights are documented in Monitors::gauge.
     * @param currentSource Source consumed by the transition under test: charge at
     * component 0 and current density at components 1 through SpatialDim.
     * @param previousSource Source consumed one timestep earlier, with the same layout.
     * @param dt Positive timestep separating the two source histories.
     * @param weight Current-source interpolation weight from LevelContract.
     * @param expected Prescribed continuity residual to subtract at every owned cell.
     * @param spatialDiscretization Source-component positions and difference directions.
     * @pre Current and previous source halos are valid for the selected differences.
     */
    inline ResidualMetrics expectedContinuityMismatch(
        const AJFourFieldsType& currentSource, const AJFourFieldsType& previousSource, double dt,
        double weight, double expected,
        SpatialDiscretization spatialDiscretization = SpatialDiscretization::Collocated) {
        const auto now = currentSource.getView(), old = previousSource.getView();
        const auto h        = currentSource.get_mesh().getMeshSpacing();
        const int direction = spatialDiscretization == SpatialDiscretization::Staggered ? -1 : 0;
        return reduceResidual(
            currentSource, KOKKOS_LAMBDA(int i, int j, int k) {
                return (monitorDetail::Expression(now(i, j, k)[0])
                        - monitorDetail::Expression(old(i, j, k)[0]))
                           * (1.0 / dt)
                       + monitorDetail::divergence(now, 1, i, j, k, h, direction) * weight
                       + monitorDetail::divergence(old, 1, i, j, k, h, direction) * (1.0 - weight)
                       - monitorDetail::Expression(expected);
            });
    }

    /** @brief Require the magnetic, Faraday and Gauss formula identities.
     * @details Monitors documents the discrete formulas and their sources;
     * each cellwise ratio uses the uncancelled absolute evaluation as scale.
     * These checks concern algebraic identities, including for unstable runs
     * and manufactured sources that do not conserve charge.
     */
    inline void requireFinalIdentities(const MonitorReport& report) {
        EXPECT_LE(report.divMagnetic_m.cell_m, 1e-12);
        EXPECT_LE(report.divRawMagnetic_m.cell_m, 1e-12);
        for (const auto& residual : report.faraday_m)
            EXPECT_LE(residual.cell_m, 1e-12);
        EXPECT_LE(report.gauss_m.cell_m, 1e-12);
    }

    /// Bracket the final requested solve without changing any earlier transition.
    template <class SolverType>
    void advanceWithFinalMonitor(SolverType& solver, Monitors& monitors, bool final) {
        if (final)
            requireFinalIdentities(monitors.advance(solver));
        else
            solver.solve();
    }
}  // namespace chdr::solver::test
