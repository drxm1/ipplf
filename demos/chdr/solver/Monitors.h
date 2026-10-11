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
            : workspace_m(prototype) {
            for (auto* field : {&oldest_m, &consumedSource_m, &previousSource_m}) {
                field->initialize(prototype.get_mesh(), prototype.getLayout());
                *field = 0.0;
            }
            previousElectric_m.initialize(prototype.get_mesh(), prototype.getLayout());
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
