/** @file Monitors.h
 * @brief Optional residual diagnostics with explicit history and source ownership.
 * @details Each spatial derivative is a fixed linear combination of periodic
 * index shifts. Those shifts commute, so the gradient/curl/divergence
 * cancellations in the potential definitions carry over to these stencils.
 * Gauge and continuity diagnostics pair scalar time differences with the
 * corresponding vector/source samples. evaluateGauss also
 * accounts for the reconstruction's scalar slot and spatial wave operator
 * (deduction from [fallahi2020mithra20fullwavesimulation, Eqs. (3.6)-(3.9)] and
 * [christlieb2024GaugeConserving, Lemma 3.3 and Theorem 3.5]). Monitoring changes
 * neither the wave update nor its potential histories. All field kernels use
 * device views; field halos and reduction scalars cross MPI. No full-field
 * host mirror is needed by the residual diagnostics.
 */
#pragma once

#include <array>
#include <stdexcept>

#include "ChdrSolver.h"
#include "MonitorExpressions.h"
#include "MonitorRecords.h"
#include "SpatialDimension.h"

static_assert(chdr::solver::SpatialDim == 3,
              "Monitors.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    /** @brief Diagnostics of one completed solve and its retained input source.
     * @details gauss_m measures the formula mismatch; gaussDefect_m measures the
     * physical defect itself. R2Prime's Gauss entries refer to the preceding
     * output E because its forward scalar difference needs this completed solve.
     * The availability flag prevents an absent first continuity sample from
     * being interpreted as a measured zero.
     */
    struct MonitorReport {
        ResidualMetrics divMagnetic_m, divRawMagnetic_m, gauge_m, continuity_m;
        std::array<ResidualMetrics, SpatialDim> faraday_m;
        ResidualMetrics gauss_m, gaussDefect_m;
        bool continuityAvailable_m = false;
    };

    /** @brief Bracket a solve with snapshots needed by the residual identities.
     * @details beforeSolve retains the older potential, previous E for R2Prime, and
     * source about to be consumed. afterSolve obtains the two newer potential
     * levels from the solver. Consequently a later caller write to the live
     * source cannot change the completed-step Gauss result.
     * Callers may start monitoring at any solve count; reconstruct initial E
     * only at count zero. Every rank in the layout participates in the same order.
     * The object allocates 23 doubles per cell including halos: three four-component
     * histories, one three-component electric history, and two four-component work
     * arrays. A consecutive collocated R2 sample performs four full-view copies,
     * 14 halo fills and nine residual reductions (two MPI collectives each).
     * Monitoring is optional and may be sparse; continuity is unavailable after
     * a skipped solve or history reload unless its preceding source is explicitly seeded.
     * Formula checks assume the inherited vacuum step and sources written before
     * solve(); a material update needs its own corresponding identities.
     */
    class Monitors {
    public:
        /** @brief Allocate matching one-halo work arrays.
         * @param prototype Allocated field whose mesh/layout owners outlive this monitor.
         */
        explicit Monitors(const AJFourFieldsType& prototype)
            : workspace_m(checkedPeriodic(prototype)) {
            for (auto* field : {&oldest_m, &consumedSource_m, &previousSource_m}) {
                field->initialize(prototype.get_mesh(), prototype.getLayout());
                *field = 0.0;
            }
            previousElectric_m.initialize(prototype.get_mesh(), prototype.getLayout());
        }

        /** @brief Reject open layouts before allocating derivative workspaces.
         * Boundary residuals need separate operators and normalization; the
         * current reports describe the periodic interior identities only.
         */
        static const AJFourFieldsType& checkedPeriodic(const AJFourFieldsType& field) {
            if (!field.getLayout().isAllPeriodic_m)
                throw std::invalid_argument("Monitors require an all-periodic layout");
            return field;
        }

        /** @brief Declare the source preceding the next monitored step.
         * @param source Four-component source consumed by the immediately preceding solve:
         * charge density in component zero and Cartesian current density in components
         * 1 through SpatialDim. Its halo is refreshed before taking the snapshot.
         * @details This explicit one-use seed is the caller's assertion of adjacency.
         * Automatic reuse thereafter requires the same solver, history revision,
         * time step and consecutive solve count. Seed again after a sparse gap.
         */
        void seedPreviousSource(AJFourFieldsType& source) {
            source.fillHalo();
            Kokkos::deep_copy(previousSource_m.getView(), source.getView());
            havePreviousSource_m = true;
            explicitlySeeded_m   = true;
        }

        /** @brief Retain the state before the source is consumed by one solve.
         * @tparam SolverType Solver exposing ChdrSolver's history, field and level interface.
         * @param solver Shell whose next solve is to be monitored.
         * @pre The caller has written the current source and initialized potential halos.
         * @throws std::logic_error If another beforeSolve is still awaiting afterSolve,
         * or if the solver has non-unit material permittivity.
         * @details These gauge/Gauss identities assume vacuum. The homogeneous medium has
         * \f$\nabla\cdot\mathbf A+\varepsilon\partial_t\phi=0\f$
         * and \f$\nabla\cdot(\varepsilon\mathbf E)=\rho\f$
         * [chew2014generalizedgauge, section 2, Eqs. (4), (11)]. A material solver at
         * unit permittivity is accepted because these reduce to the vacuum forms.
         * Refusal happens before capturing inputs or marking this monitor prepared;
         * it rejects this diagnostic, not the solver's ability to advance a medium.
         * Continuity is evaluated here, using this and the preceding consumed
         * source at the matched time weight supplied by LevelContract::gaugeWeight.
         */
        template <class SolverType>
        void beforeSolve(SolverType& solver) {
            if (solver.hasMaterial() && solver.homogeneousPermittivity() != 1.0)
                throw std::logic_error("Monitors require vacuum or unit material permittivity");
            if (prepared_m)
                throw std::logic_error("Monitors::beforeSolve called twice without afterSolve");
            captureInputs(solver);
            preparedAt_m                   = solver.completedSolves();
            preparedSolver_m               = &solver;
            preparedRevision_m             = solver.historyRevision();
            preparedDt_m                   = solver.getDt();
            prepared_m                     = true;
            report_m.continuityAvailable_m = previousSourceMatches(solver);
            report_m.continuity_m          = {};
            if (report_m.continuityAvailable_m)
                report_m.continuity_m = continuity(
                    consumedSource_m, previousSource_m, solver.getDt(),
                    solver.levelContract().gaugeWeight(), solver.spatialDiscretization());
            explicitlySeeded_m = false;
        }

        /** @brief Evaluate the completed step and retain its source for continuity.
         * @tparam SolverType Same solver type passed to beforeSolve.
         * @param solver Same shell, after exactly the one bracketed solve.
         * @return Updated records whose availability flags describe retained history.
         */
        template <class SolverType>
        const MonitorReport& afterSolve(SolverType& solver) {
            if (!prepared_m || !retainedStepMatches(solver))
                throw std::logic_error("Monitor histories do not bracket one completed solve");
            evaluate(solver);
            Kokkos::deep_copy(previousSource_m.getView(), consumedSource_m.getView());
            havePreviousSource_m = true;
            previousAt_m         = preparedAt_m;
            previousSolver_m     = preparedSolver_m;
            previousRevision_m   = preparedRevision_m;
            previousDt_m         = preparedDt_m;
            prepared_m           = false;
            return report_m;
        }

        /// Convenience wrapper preserving the solver's original solve operation.
        template <class SolverType>
        const MonitorReport& advance(SolverType& solver) {
            beforeSolve(solver);
            solver.solve();
            return afterSolve(solver);
        }

        /** @brief Re-evaluate the retained step, including after a live-source overwrite.
         * @pre The solver has completed exactly the solve bracketed by beforeSolve.
         * @details Magnetic, Faraday, gauge and Gauss reductions retain
         * their fixed collective order. Temporary derivative values and their
         * absolute evaluations exchange halos together before the outer stencil.
         */
        template <class SolverType>
        const MonitorReport& evaluate(SolverType& solver) {
            if (!retainedStepMatches(solver))
                throw std::logic_error("No retained monitor history for this solve");
            magneticDiagnostics(solver);
            faradayDiagnostics(solver);
            gaugeAndFieldEquations(solver);
            return report_m;
        }

        /** @brief Gauge at the reconstruction's matched vector-potential time weight.
         * @param currentField Current four-component history: scalar potential in component
         * zero, Cartesian vector potential in components 1 through SpatialDim.
         * @param previousField The preceding history, one dt earlier for each component.
         * Both fields use the same geometry and component placement; their halos are filled.
         * @param dt Positive interval between corresponding components of the two histories.
         * @param weight Current-vector weight supplied by LevelContract::gaugeWeight.
         * @param spatialDiscretization Selects centered or backward staggered divergence.
         * @details \f$G_w=(\phi_n-\phi_{n-1})/dt+
         * D_c\cdot(wA_n+(1-w)A_{n-1})\f$, with weights 1, 1/2, and 0 for R2,
         * R1, and R2Prime. R2 pairs the scalar difference with current A, R1
         * with the mean of both A slots, and R2Prime with previous A. These
         * pairings align the vector sample with the scalar-difference time
         * (deduction from [christlieb2024GaugeConserving, Eq. (29)] and the
         * component times documented in LevelContract). For the current/previous
         * slots after completed solve \f$k\f$, the report times are \f$t_k\f$ for R2,
         * \f$t_{k-1/2}\f$ for R1 and \f$t_{k-1}\f$ for R2Prime, with \f$t_k=k\,dt\f$.
         * R1's mean gauge differs from the
         * auxiliary backward gauge in its Gauss formula. Staggered R2 uses
         * backward face-to-centre divergence [ryu2016potentialfdtd, Appendix A,
         * Eq. (23)]. Both fields' halos are synchronized before reduction.
         */
        static ResidualMetrics gauge(
            AJFourFieldsType& currentField, AJFourFieldsType& previousField, double dt,
            double weight,
            SpatialDiscretization spatialDiscretization = SpatialDiscretization::Collocated) {
            currentField.fillHalo();
            previousField.fillHalo();
            const auto current = currentField.getView(), previous = previousField.getView();
            const auto spacing = currentField.get_mesh().getMeshSpacing();
            const int backward = spatialDiscretization == SpatialDiscretization::Staggered ? -1 : 0;
            return reduceResidual(
                currentField, KOKKOS_LAMBDA(int i, int j, int k) {
                    return (monitorDetail::Expression(current(i, j, k)[0])
                            - monitorDetail::Expression(previous(i, j, k)[0]))
                               * (1.0 / dt)
                           + monitorDetail::divergence(current, 1, i, j, k, spacing, backward)
                                 * weight
                           + monitorDetail::divergence(previous, 1, i, j, k, spacing, backward)
                                 * (1.0 - weight);
                });
        }

        /** @brief Evaluate charge continuity using the same time pairing as the gauge.
         * @param source Source of the current step: charge density in component zero,
         * Cartesian current density in components 1 through SpatialDim.
         * @param previousSource Source of the immediately preceding step on the same layout.
         * Both source halos are refreshed by gauge.
         * @param dt Positive time spacing between corresponding source components.
         * @param weight Current-density time weight from LevelContract::gaugeWeight.
         * @param spatialDiscretization Selects centered or backward staggered divergence.
         * @details The residual is
         * \f$(\rho_n-\rho_{n-1})/\Delta t+
         * D_h\cdot(w\mathbf J_n+(1-w)\mathbf J_{n-1})\f$,
         * the source counterpart of gauge's matched expression (deduction from
         * [christlieb2024GaugeConserving, Lemma 3.3] with LevelContract's time pairing).
         */
        static ResidualMetrics continuity(
            AJFourFieldsType& source, AJFourFieldsType& previousSource, double dt, double weight,
            SpatialDiscretization spatialDiscretization = SpatialDiscretization::Collocated) {
            return gauge(source, previousSource, dt, weight, spatialDiscretization);
        }

        /// @return Most recently evaluated records and their availability flags.
        const MonitorReport& report() const { return report_m; }
        /// @return Retained input source of the bracketed solve; independent of live source writes.
        AJFourFieldsType& consumedSource() { return consumedSource_m; }

    private:
        template <class SolverType>
        bool previousSourceMatches(const SolverType& solver) const {
            return havePreviousSource_m
                   && (explicitlySeeded_m
                       || (previousSolver_m == &solver
                           && previousRevision_m == solver.historyRevision()
                           && previousAt_m == solver.completedSolves() - 1
                           && previousDt_m == solver.getDt()));
        }

        template <class SolverType>
        bool retainedStepMatches(const SolverType& solver) const {
            return preparedAt_m >= 0 && preparedSolver_m == &solver
                   && preparedRevision_m == solver.historyRevision()
                   && preparedDt_m == solver.getDt()
                   && solver.completedSolves() == preparedAt_m + 1;
        }

        /// Deep copies preserve histories across timeShift and later caller source writes.
        template <class SolverType>
        void captureInputs(SolverType& solver) {
            if (solver.completedSolves() == 0)
                solver.reconstructEB();
            solver.sourceField().fillHalo();
            Kokkos::deep_copy(oldest_m.getView(), solver.A_nm1.getView());
            Kokkos::deep_copy(previousElectric_m.getView(), solver.electricField().getView());
            Kokkos::deep_copy(consumedSource_m.getView(), solver.sourceField().getView());
        }

    public:
        // NVCC extended lambdas require their enclosing member function to be public.
        // These workspace helpers are implementation details; use advance/beforeSolve/afterSolve.
        /** @brief Cache output B and its uncancelled scale, or compute the raw curl.
         * @param solver Shell supplying valid potential halos and reconstructed B.
         * @param raw Select current raw curl rather than measured output B.
         * @details Internal workspace operation; advance controls its call order.
         */
        template <class SolverType>
        void fillMagneticWorkspace(SolverType& solver, bool raw) {
            const auto current = solver.A_n.getView(), previous = solver.A_nm1.getView();
            const auto output  = solver.magneticField().getView();
            const auto spacing = solver.A_n.get_mesh().getMeshSpacing();
            const int forward =
                solver.spatialDiscretization() == SpatialDiscretization::Staggered ? 1 : 0;
            const double weight =
                solver.ebReconstructionConvention() == EBReconstructionConvention::R1 ? 1.0 : 0.5;
            workspace_m.fill(KOKKOS_LAMBDA(int i, int j, int k, unsigned d) {
                const auto now = monitorDetail::curl(current, 1, d, i, j, k, spacing, forward);
                if (raw)
                    return now;
                const auto bound = now * weight
                                   + monitorDetail::curl(previous, 1, d, i, j, k, spacing, forward)
                                         * (1.0 - weight);
                return monitorDetail::Expression{output(i, j, k)[d], bound.scale_m};
            });
        }

    private:
        template <class SolverType>
        void magneticDiagnostics(SolverType& solver) {
            const int direction =
                solver.spatialDiscretization() == SpatialDiscretization::Staggered ? 1 : 0;
            fillMagneticWorkspace(solver, false);
            report_m.divMagnetic_m = divergenceOfWorkspace(solver.A_n, direction);
            fillMagneticWorkspace(solver, true);
            report_m.divRawMagnetic_m = divergenceOfWorkspace(solver.A_n, direction);
        }

    public:
        /** @brief Cache measured E while expanding its potential-based absolute scale.
         * @param solver Shell supplying potential histories and reconstructed E.
         * @details Internal workspace operation; advance controls its call order.
         */
        template <class SolverType>
        void fillElectricWorkspace(SolverType& solver) {
            const auto current = solver.A_n.getView(), previous = solver.A_nm1.getView();
            const auto scalar =
                solver.ebReconstructionConvention() == EBReconstructionConvention::R2 ? previous
                                                                                      : current;
            const auto output  = solver.electricField().getView();
            const auto spacing = solver.A_n.get_mesh().getMeshSpacing();
            const int forward =
                solver.spatialDiscretization() == SpatialDiscretization::Staggered ? 1 : 0;
            const double inverseDt = 1.0 / solver.getDt();
            workspace_m.fill(KOKKOS_LAMBDA(int i, int j, int k, unsigned d) {
                const auto bound = (monitorDetail::Expression(current(i, j, k)[d + 1])
                                    - monitorDetail::Expression(previous(i, j, k)[d + 1]))
                                       * inverseDt
                                   + monitorDetail::derivative(monitorDetail::Component{scalar, 0u},
                                                               i, j, k, d, spacing[d], forward);
                return monitorDetail::Expression{output(i, j, k)[d], bound.scale_m};
            });
        }

        /** @brief Evaluate Faraday using consecutive raw curls and the current output E.
         * @param solver Shell supplying valid current and previous potential halos.
         * @param d Cartesian residual component, from zero through two.
         * @return Reduction of the selected residual against its absolute scale.
         * @pre fillElectricWorkspace has prepared this monitor's workspace.
         * @details The residual is
         * \f$(B_{\mathrm{raw}}[k]-B_{\mathrm{raw}}[k-1])/dt+\operatorname{curl}_h E[k]\f$.
         * The vector-potential time difference cancels the raw-curl difference.
         * The remaining scalar gradient has zero discrete curl because the
         * coordinate differences commute (deduction from
         * [fallahi2020mithra20fullwavesimulation, Eqs. (3.8)-(3.9)]).
         */
        template <class SolverType>
        ResidualMetrics faradayComponent(SolverType& solver, unsigned d) {
            const auto current = solver.A_n.getView(), previous = solver.A_nm1.getView();
            const auto spacing = solver.A_n.get_mesh().getMeshSpacing();
            const int forward =
                solver.spatialDiscretization() == SpatialDiscretization::Staggered ? 1 : 0;
            const double inverseDt = 1.0 / solver.getDt();
            const unsigned a = (d + 1) % SpatialDim, b = (d + 2) % SpatialDim;
            const auto electricA = workspace_m.component(a), electricB = workspace_m.component(b);
            return reduceResidual(
                solver.A_n, KOKKOS_LAMBDA(int i, int j, int k) {
                    return (monitorDetail::curl(current, 1, d, i, j, k, spacing, forward)
                            - monitorDetail::curl(previous, 1, d, i, j, k, spacing, forward))
                               * inverseDt
                           + monitorDetail::derivative(electricB, i, j, k, a, spacing[a], forward)
                           - monitorDetail::derivative(electricA, i, j, k, b, spacing[b], forward);
                });
        }

    private:
        template <class SolverType>
        void faradayDiagnostics(SolverType& solver) {
            fillElectricWorkspace(solver);
            for (unsigned d = 0; d < SpatialDim; ++d)
                report_m.faraday_m[d] = faradayComponent(solver, d);
        }

        /// Keep formula checks separate from physical-defect magnitudes and record availability.
        template <class SolverType>
        void gaugeAndFieldEquations(SolverType& solver) {
            report_m.gauge_m =
                gauge(solver.A_n, solver.A_nm1, solver.getDt(),
                      solver.levelContract().gaugeWeight(), solver.spatialDiscretization());
            const auto gauss =
                evaluateGauss(solver, oldest_m, consumedSource_m, previousElectric_m, workspace_m);
            report_m.gauss_m       = gauss.identity_m;
            report_m.gaussDefect_m = gauss.defect_m;
        }

    public:
        /** @brief Reduce the divergence of the current prepared workspace.
         * @param prototype Field supplying owned-cell range and communicator.
         * @param direction Zero for centered, positive/negative for forward/backward differences.
         * @return Divergence residual and its propagated absolute scale.
         * @pre A preceding fill operation synchronized all workspace halos.
         */
        ResidualMetrics divergenceOfWorkspace(const AJFourFieldsType& prototype,
                                              int direction = 0) {
            const auto x = workspace_m.component(0), y = workspace_m.component(1),
                       z       = workspace_m.component(2);
            const auto spacing = prototype.get_mesh().getMeshSpacing();
            return reduceResidual(
                prototype, KOKKOS_LAMBDA(int i, int j, int k) {
                    return monitorDetail::derivative(x, i, j, k, 0, spacing[0], direction)
                           + monitorDetail::derivative(y, i, j, k, 1, spacing[1], direction)
                           + monitorDetail::derivative(z, i, j, k, 2, spacing[2], direction);
                });
        }

    private:
        monitorDetail::Workspace workspace_m;
        AJFourFieldsType oldest_m, consumedSource_m, previousSource_m;
        EMFieldsType previousElectric_m;
        MonitorReport report_m;
        long long preparedAt_m                = -1;
        long long previousAt_m                = -1;
        const void* preparedSolver_m          = nullptr;
        const void* previousSolver_m          = nullptr;
        unsigned long long preparedRevision_m = 0, previousRevision_m = 0;
        double preparedDt_m = 0.0, previousDt_m = 0.0;
        bool explicitlySeeded_m = false;
        bool prepared_m = false, havePreviousSource_m = false;
    };
}  // namespace chdr::solver
