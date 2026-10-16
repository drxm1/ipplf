/** @file test_ranks.cpp
 * @brief Record gauge waves and monitor data crossing MPI decomposition seams.
 * @details Global-index snapshots and independent reductions are compared by
 * the CTest fixture. Monitor residuals are normalized before comparison; a
 * near-zero residual is never used as its own relative-error denominator.
 *
 * TODO: REVIEW.
 */
#include <gtest/gtest.h>

#include "BoxSetup.h"
#include "ChdrSolver.h"
#include "MonitorTestSupport.h"
#include "Monitors.h"
#include "RankRecords.h"
#include "SpatialDimension.h"
#include "Waves.h"

static_assert(chdr::solver::SpatialDim == 3,
              "test_ranks.cpp requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    /// Reproducible data varying in all directions, with substantial short-scale content.
    struct RoughData {
        bool nyquist_m = false;
        int nz_m       = 16;
        KOKKOS_FUNCTION double operator()(double x, double y, double z, double time,
                                          unsigned c) const {
            const double phase    = 127 * x + 193 * y + 263 * z + 0.71 * c;
            const double original = Kokkos::sin(phase) + 0.37 * Kokkos::cos(3 * phase - 2 * time);
            if (!nyquist_m)
                return original;
            const double pi = Kokkos::numbers::pi_v<double>;
            const double xy = Kokkos::sin(16 * pi * x) * Kokkos::sin(16 * pi * y);
            return original + 0.3 * xy * Kokkos::sin(nz_m * pi * z)
                   + 0.2 * xy * Kokkos::sin(2 * pi * z);
        }
    };

    /// Seed the previous source before exercising three complete monitor lifecycles.
    template <class SolverType>
    void evolveMonitorSeamCase(SolverType& solver, Monitors& monitors, RoughData data) {
        writeSource(solver, data, -1);
        monitors.seedPreviousSource(solver.sourceField());
        for (int step = 0; step < 3; ++step) {
            writeSource(solver, data, step);
            monitors.advance(solver);
        }
    }

    /// Retain absolute maxima and Euclidean norms independently of normalized residual ratios.
    inline std::vector<double> residualNorms(const MonitorReport& report) {
        std::vector<double> norms;
        for (const auto& metric :
             {report.divMagnetic_m, report.divRawMagnetic_m, report.gauge_m, report.continuity_m,
              report.gauss_m, report.faraday_m[0], report.faraday_m[1], report.faraday_m[2]}) {
            norms.push_back(metric.l2_m);
            norms.push_back(metric.scaleL2_m);
            norms.push_back(metric.maxAbs_m);
            norms.push_back(metric.maxScale_m);
        }
        return norms;
    }

    /// Record global residual ratios, with the absolute-evaluation scales from Monitors.
    inline void recordResiduals(const std::string& name, const MonitorReport& report) {
        std::vector<double> residuals = {report.divMagnetic_m.global_m,
                                         report.divRawMagnetic_m.global_m, report.gauge_m.global_m,
                                         report.continuity_m.global_m, report.gauss_m.global_m};
        for (int d = 0; d < static_cast<int>(SpatialDim); ++d)
            residuals.push_back(report.faraday_m[d].global_m);
        requireFinalIdentities(report);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            std::ofstream out(recordDirectory() / (name + "-residuals.json"));
            out << nlohmann::json{{"normalized_residuals", residuals},
                                  {"residual_norms", residualNorms(report)}}
                       .dump()
                << '\n';
        }
    }

    template <class PotentialSolverBase>
    void monitorRankRecord(SpatialDiscretization spatialDiscretization, bool nyquist = false) {
        constexpr bool NS = std::is_same_v<PotentialSolverBase, NonStandardSolverPeriodic>;
        auto fields       = makeBoxFieldContainer({16, 16, NS ? 24 : 16}, {1, 1, 1}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               EBReconstructionConvention::R2,
                                               spatialDiscretization);
        const RoughData data{nyquist, NS ? 24 : 16};
        loadHistory(solver, data);
        solver.evaluateInitialFields();
        Monitors monitors(solver.A_n);
        evolveMonitorSeamCase(solver, monitors, data);
        const auto name = std::string(nyquist ? "monitors-nyquist-" : "monitors-")
                          + (NS ? "NS" : "Std") + "-"
                          + std::to_string(static_cast<int>(spatialDiscretization));
        recordState(name, solver);
        recordResiduals(name, monitors.report());
    }

    /** @brief Evolve the discrete gauge wave over at least ten periods.
     * @details Histories use the exact discrete frequency; no forcing is applied.
     * The discrete wave and reconstruction are defined in Waves.h, from
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.17), (3.24)] and
     * [christlieb2024GaugeConserving, section 3.1, proof of Theorem 3.5].
     * T2 names these evolved discrete-gauge records.
     */
    template <class SolverType>
    void evolveGaugeRankCase(SolverType& solver, const GaugeWave& wave) {
        loadHistory(solver, wave);
        solver.evaluateInitialFields();
        const double pi = Kokkos::numbers::pi_v<double>;
        const int steps = std::ceil(20 * pi / wave.omega_m / solver.getDt());
        Monitors monitors(solver.A_n);
        for (int step = 0; step < steps; ++step)
            advanceWithFinalMonitor(solver, monitors, step == steps - 1);
    }

    template <class PotentialSolverBase>
    void gaugeRankRecord() {
        constexpr bool NS = std::is_same_v<PotentialSolverBase, NonStandardSolverPeriodic>;
        auto fields       = makeBoxFieldContainer({32, 32, NS ? 48 : 32}, {1, 1, 1}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               EBReconstructionConvention::R2);
        const double pi = Kokkos::numbers::pi_v<double>;
        const GaugeWave wave({2 * pi, 4 * pi, 4 * pi}, fields.getHr(), solver.getDt(),
                             ChdrSolver<PotentialSolverBase>::stencil(),
                             EBReconstructionConvention::R2);
        evolveGaugeRankCase(solver, wave);
        recordState(NS ? "T2-NS" : "T2-Std", solver);
    }

    TEST(RankCases, GaugeWaves) {
        gaugeRankRecord<StandardSolverPeriodic>();
        gaugeRankRecord<NonStandardSolverPeriodic>();
    }
    TEST(RankCases, MonitorSeams) {
        monitorRankRecord<StandardSolverPeriodic>(SpatialDiscretization::Collocated);
        monitorRankRecord<NonStandardSolverPeriodic>(SpatialDiscretization::Collocated);
        monitorRankRecord<StandardSolverPeriodic>(SpatialDiscretization::Staggered);
    }
    TEST(RankCases, NyquistMonitorSeams) {
        monitorRankRecord<StandardSolverPeriodic>(SpatialDiscretization::Collocated, true);
        monitorRankRecord<NonStandardSolverPeriodic>(SpatialDiscretization::Collocated, true);
        monitorRankRecord<StandardSolverPeriodic>(SpatialDiscretization::Staggered, true);
    }
}  // namespace chdr::solver::test
