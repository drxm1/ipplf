/** @file MonitorRecords.h
 * @brief Gauss-formula checks of evolved field data.
 * @details GaussIdentity compares the measured field defect
 * with predictions obtained by substituting the potential wave update into
 * the reconstructed fields. The formula is a fully discrete deduction from
 * [christlieb2024GaugeConserving, Theorem 3.5] and
 * [fallahi2020mithra20fullwavesimulation, Eqs. (3.6)-(3.9)]. LevelContract
 * documents the component sample times and positions. A physical defect
 * and a mismatch against its discrete formula are distinct reported quantities.
 * Formula expressions retain their absolute evaluations throughout both stencil
 * passes; MonitorExpressions.h explains the normalization and MPI reductions.
 * Kernels run on the field execution space with device views. Workspace and E
 * halos use IPPL communication; residual norms reduce scalars on the field
 * communicator without constructing host field mirrors. Formula validity assumes
 * the inherited vacuum wave update and the source retained before that update;
 * a derived material step must supply its own corresponding formulas.
 * Centered differences have nonconstant null modes: vanishing derivative-based
 * residuals alone do not exclude grid-scale vector potentials (deduction from
 * [trefethen1996FiniteDifferenceSpectral, section 5.3, p. 203]).
 */
#pragma once

#include "core/base/SpatialDimension.h"
#include "core/solver/ChdrSolver.h"
#include "opt-in/diagnostics/MonitorExpressions.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "MonitorRecords.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    /// Distinguish the physical Gauss defect from mismatch against its discrete formula.
    struct GaussReport {
        /** @brief Metrics of the measured Gauss-law defect.
         * @details monitorDetail::GaussDefect subtracts the retained consumed charge
         * density from the discrete divergence of the electric field selected by
         * evaluateGauss. Its documented time-level pairing applies to this report.
         * The defect can be nonzero even when the potential update satisfies its
         * discrete identity; see monitorDetail::GaussIdentity and its deduction from
         * [christlieb2024GaugeConserving, Theorem 3.5].
         */
        ResidualMetrics defect_m;

        /** @brief Metrics of measured Gauss defect minus its discrete prediction.
         * @details monitorDetail::GaussIdentity subtracts its reconstruction-specific
         * prediction from monitorDetail::GaussDefect before reduceResidual computes
         * these metrics. The prediction and its time/stencil choices are documented
         * there as a deduction from [christlieb2024GaugeConserving, Theorem 3.5].
         * A small identity mismatch indicates agreement with that formula; it does
         * not imply that the measured Gauss defect reported by defect_m is small.
         */
        ResidualMetrics identity_m;
    };
    namespace monitorDetail {
        /// Four-component potential or source view; the field's role determines the units.
        using FourComponentView = AJFourFieldsType::view_type;
        using ElectricFieldView = EMFieldsType::view_type;
        using CachedComponent   = StoredComponent<FourComponentView>;

        /** @brief Device-visible history and operator metadata for one completed solve.
         * @details now_m is the solver's A_n slot after the monitored solve; old_m
         * is A_nm1 after that solve. oldest_m is A_nm1 copied before the solve,
         * so these are three consecutive potential histories, newest first.
         * Each element contains scalar potential in component zero and Cartesian
         * vector potential in components 1 through SpatialDim. LevelContract supplies
         * the component-specific times and positions. The three scalar values belong
         * to actual evolved levels, including for delayed R2Prime Gauss checks.
         * Directions are zero for centered differences, or +1/-1 for the compact
         * forward/backward pair [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(20)].
         */
        struct RecordHistory {
            FourComponentView now_m, old_m, oldest_m;
            ippl::Vector<double, SpatialDim> spacing_m;
            double inverseDt_m;
            bool compact_m, prime_m;
            int forward_m, backward_m;
        };

        /** @brief Borrow the three potential histories needed by the Gauss formula.
         * @tparam PotentialSolverBase StandardSolverPeriodic or NonStandardSolverPeriodic,
         * as required by ChdrSolver, which supplies the history and level interface.
         * @param solver The ChdrSolver immediately after the monitored solve.
         * @param oldestA Independent A_nm1 snapshot captured before that solve.
         * @pre All three potential histories have valid halos on the same layout.
         */
        template <class PotentialSolverBase>
        RecordHistory recordHistory(const ChdrSolver<PotentialSolverBase>& solver,
                                    const AJFourFieldsType& oldestA) {
            const bool compact = solver.spatialDiscretization() == SpatialDiscretization::Staggered;
            const int forward  = compact ? 1 : 0;
            return {solver.A_n.getView(),
                    solver.A_nm1.getView(),
                    oldestA.getView(),
                    solver.A_n.get_mesh().getMeshSpacing(),
                    1.0 / solver.getDt(),
                    compact,
                    solver.ebReconstructionConvention() == EBReconstructionConvention::R2Prime,
                    forward,
                    -forward};
        }

        /** @brief Time derivative of the auxiliary gauge used by GaussIdentity.
         * @details The backward form pairs each scalar difference with the newer
         * vector slot. R2Prime pairs it with the preceding vector slot, so the
         * needed forward scalar difference is available one solve later.
         * R1 also uses the backward form for its Gauss prediction; its separately
         * reported Monitors::gauge uses the mean vector potential instead.
         * These pairings follow by substituting the field reconstruction into
         * the potential update [christlieb2024GaugeConserving, Theorem 3.5],
         * with the slot convention in RecordHistory (deduction).
         */
        KOKKOS_INLINE_FUNCTION Expression gaugeDifference(RecordHistory h, int i, int j, int k) {
            const auto now =
                (Expression(h.now_m(i, j, k)[0]) - Expression(h.old_m(i, j, k)[0])) * h.inverseDt_m
                + divergence(h.prime_m ? h.old_m : h.now_m, 1, i, j, k, h.spacing_m, h.backward_m);
            const auto old = (Expression(h.old_m(i, j, k)[0]) - Expression(h.oldest_m(i, j, k)[0]))
                                 * h.inverseDt_m
                             + divergence(h.prime_m ? h.oldest_m : h.old_m, 1, i, j, k, h.spacing_m,
                                          h.backward_m);
            return (now - old) * h.inverseDt_m;
        }

        /** @brief Evaluate electric divergence minus the retained consumed charge density.
         * @details electric_m stores E in components 0 through SpatialDim-1;
         * source_m stores charge density in component zero. The consumed source is
         * independent of subsequent caller writes; GaussIdentity documents its pairing.
         */
        struct GaussDefect {
            ElectricFieldView electric_m;
            FourComponentView source_m;
            ippl::Vector<double, SpatialDim> spacing_m;
            int backward_m;
            KOKKOS_INLINE_FUNCTION Expression operator()(int i, int j, int k) const {
                return divergence(electric_m, 0, i, j, k, spacing_m, backward_m)
                       - Expression(source_m(i, j, k)[0]);
            }
        };

        /** @brief Subtract the complete formula from the measured Gauss defect.
         * @details For collocated R2 at completed solve k, scalar-wave substitution gives
         * \f[ D_c\cdot E[k]-\rho_s[k-1]
         * =L_h\phi[k-1]-D_c\cdot D_c\phi[k-1]
         * -(\Gamma_0[k]-\Gamma_0[k-1])/dt. \f]
         * Here the source \f$\rho_s[j]\f$ advances potential slot \f$j\f$ to \f$j+1\f$,
         * with \f$\Gamma_0[k]=(\phi[k]-\phi[k-1])/dt+D_c\cdot A[k]\f$.
         * This is a deduction from
         * [christlieb2024GaugeConserving, Theorem 3.5] using the wave operator
         * documented in Operators.h. R1 reads its scalar gradient from the current
         * slot, so its two double-gradient terms combine to \f$D_c\cdot D_c\phi[k]\f$;
         * the wave-operator term still uses the previous scalar slot.
         * The compact reading omits the spatial mismatch because backward
         * divergence of forward gradient equals the Standard Laplacian
         * (Cartesian deduction from
         * [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(23)] by substitution).
         * @tparam Scheme Inherited solver's spatial wave stencil, fixed before kernel launch.
         */
        template <Stencil Scheme>
        struct GaussIdentity {
            RecordHistory history_m;
            GaussDefect defect_m;
            CachedComponent x_m, y_m, z_m;
            KOKKOS_INLINE_FUNCTION Expression operator()(int i, int j, int k) const {
                const auto h = history_m;
                const auto doubleGradient =
                    derivative(x_m, i, j, k, 0, h.spacing_m[0], h.backward_m)
                    + derivative(y_m, i, j, k, 1, h.spacing_m[1], h.backward_m)
                    + derivative(z_m, i, j, k, 2, h.spacing_m[2], h.backward_m);
                const auto mismatch =
                    h.compact_m ? Expression{}
                                : laplace<Scheme>(Component{h.old_m, 0u}, i, j, k, h.spacing_m)
                                      - doubleGradient;
                const auto prediction = mismatch - gaugeDifference(h, i, j, k);
                return defect_m(i, j, k) - prediction;
            }
        };

        /** @brief Cache a scalar gradient and its absolute evaluation before divergence.
         * @param scalarPotentialView Four-component potential view; only component zero
         * is differentiated. The reconstruction convention selects its history slot.
         * @param h Supplies mesh spacings and the selected difference direction.
         * @param workspace Receives the Cartesian gradient in its leading components.
         */
        inline void prepareGaussGradient(FourComponentView scalarPotentialView, RecordHistory h,
                                         Workspace& workspace) {
            workspace.fill(KOKKOS_LAMBDA(int i, int j, int k, unsigned d) {
                return derivative(Component{scalarPotentialView, 0u}, i, j, k, d, h.spacing_m[d],
                                  h.forward_m);
            });
        }

    }  // namespace monitorDetail

    /** @brief Compare the measured Gauss defect with its reconstruction-specific prediction.
     * @tparam SolverType Solver exposing ChdrSolver potential histories, fields and stencil.
     * @param solver The ChdrSolver immediately after the monitored solve; its electric halo may
     * be filled.
     * @param oldestA Four-component A_nm1 snapshot from before the monitored solve;
     * component zero is scalar potential and the remaining components are vector potential.
     * @param consumedSource Retained source used by that solve, with charge density in
     * component zero and Cartesian current density in components 1 through SpatialDim.
     * @param previousElectric Electric output before the solve, used by R2Prime;
     * its halo is refreshed when selected.
     * @param workspace Matching temporary storage for gradient values and absolute scales.
     * @details R2 and R1 use the new output E. R2Prime uses previousElectric:
     * its forward gauge requires the newly evolved scalar slot. consumedSource
     * always holds the source of this completed solve, unaffected by subsequent
     * caller writes. Staggered R2 uses \f$-\delta_t G\f$ alone; see GaussIdentity.
     * @pre The potential histories have valid halos; all layout ranks participate.
     * @details This function refreshes the selected E halo before differentiating it
     * and exchanges the gradient workspace between stencil passes.
     */
    template <class SolverType>
    GaussReport evaluateGauss(SolverType& solver, const AJFourFieldsType& oldestA,
                              const AJFourFieldsType& consumedSource,
                              EMFieldsType& previousElectric, monitorDetail::Workspace& workspace) {
        const auto h   = monitorDetail::recordHistory(solver, oldestA);
        auto& electric = h.prime_m ? previousElectric : solver.electricField();
        electric.fillHalo();
        const auto scalar = solver.ebReconstructionConvention() == EBReconstructionConvention::R1
                                ? h.now_m
                                : h.old_m;
        monitorDetail::prepareGaussGradient(scalar, h, workspace);
        const monitorDetail::GaussDefect defect{electric.getView(), consumedSource.getView(),
                                                h.spacing_m, h.backward_m};
        const monitorDetail::GaussIdentity<SolverType::stencil()> identity{
            h, defect, workspace.component(0), workspace.component(1), workspace.component(2)};
        return {reduceResidual(solver.A_n, defect), reduceResidual(solver.A_n, identity)};
    }

}  // namespace chdr::solver
