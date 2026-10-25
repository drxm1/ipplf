/** @file test_gauge_growth.cpp
 * @brief Compare the first and last period of an evolved gauge residual.
 * @details The discrete gauge obeys a wave equation driven by the continuity
 * defect: [christlieb2024GaugeConserving, Lemma 3.3,
 * Eq. (38)]. This test measures the norm from actual IPPL arrays, not a modal
 * predictor, over at least ten periods of the discrete gauge and mixed Lorenz waves.
 * B5 labels the collocated first/last-period growth record; T2 is the discrete
 * gauge mode and T3 the continuum mixed wave. Staggered histories are records
 * without a ratio criterion because their gauge is at round-off.
 *
 * TODO: REVIEW.
 */
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <iostream>

#include "MonitorTestSupport.h"
#include "RankRecords.h"
#include "Waves.h"
#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "core/solver/ChdrSolver.h"
#include "opt-in/diagnostics/Monitors.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_gauge_growth.cpp requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    struct GaugeEnvelope {
        double first_m = 0, last_m = 0;
        double initial_m = 0;
        std::vector<double> history_m;
        int steps_m = 0;
        /// First and last windows are defined in physical time, not rounded period counts.
        void observe(double time, double period, double norm) {
            if (time <= period)
                first_m = std::max(first_m, norm);
            if (time >= 9 * period)
                last_m = std::max(last_m, norm);
        }
    };

    /// Observe the report wiring at both ends and the same standalone norm between them.
    template <class SolverType>
    ResidualMetrics gaugeObservation(SolverType& solver, Monitors& monitors, bool bracket) {
        if (bracket) {
            requireFinalIdentities(monitors.advance(solver));
            return monitors.report().gauge_m;
        }
        solver.solve();
        return Monitors::gauge(solver.A_n, solver.A_nm1, solver.getDt(), 1.0,
                               solver.spatialDiscretization());
    }

    template <class SolverType>
    GaugeEnvelope gaugeEnvelope(SolverType& solver, double period) {
        GaugeEnvelope result;
        result.steps_m = std::ceil(10 * period / solver.getDt());
        Monitors monitors(solver.A_n);
        for (int n = 1; n <= result.steps_m; ++n) {
            const auto gauge = gaugeObservation(solver, monitors, n == 1 || n == result.steps_m);
            EXPECT_TRUE(std::isfinite(gauge.l2_m));
            if (n == 1)
                result.initial_m = gauge.l2_m;
            result.history_m.push_back(gauge.l2_m);
            result.observe(n * solver.getDt(), period, gauge.l2_m);
        }
        return result;
    }

    /** @brief Closed norm of the collocated gauge mode, independent of monitor code.
     * @details Deduction: the gauge amplitude is \f$K_L^2-\sum_d g_d^2\f$;
     * GaugeWave documents its derivation from [christlieb2024GaugeConserving,
     * Eq. (29)] and [fallahi2020mithra20fullwavesimulation, Eqs. (3.17),(3.24)].
     * A non-Nyquist integer mode has cosine-square sum equal to half the cell
     * count by finite Fourier orthogonality [heinzel2002SpectrumSpectralDensity,
     * section 5, Eqs. (3)-(5)]. This fixture has neither a zero nor Nyquist mode.
     */
    template <class SolverType>
    double exactGaugeNorm(SolverType& solver, const GaugeWave& wave) {
        const auto h     = solver.A_n.get_mesh().getMeshSpacing();
        double amplitude = -laplacianSymbol(wave.waveNumber_m, h, SolverType::stencil());
        for (unsigned d = 0; d < SpatialDim; ++d)
            amplitude -= wave.g_m[d] * wave.g_m[d];
        const auto domain  = solver.A_n.getLayout().getDomain();
        const double cells = double(domain[0].length()) * domain[1].length() * domain[2].length();
        return std::abs(amplitude) * std::sqrt(cells / 2);
    }

    /// Archive every sampled norm; staggered histories deliberately have no ratio pass rule.
    template <class SolverType>
    void recordGaugeHistory(SolverType& solver, const GaugeEnvelope& result, bool gaugeWave) {
        if (ippl::Comm->rank() != 0)
            return;
        const std::string name = "gauge-history-" + std::to_string(int(SolverType::stencil())) + "-"
                                 + std::to_string(int(solver.spatialDiscretization()))
                                 + (gaugeWave ? "-discrete" : "-mixed");
        std::filesystem::create_directories(recordDirectory());
        std::ofstream out(recordDirectory() / (name + ".jsonl"));
        out << nlohmann::json{{"dt", solver.getDt()},
                              {"initial", result.initial_m},
                              {"first_period_max", result.first_m},
                              {"last_period_max", result.last_m},
                              {"norms", result.history_m},
                              {"staggered_record_only",
                               solver.spatialDiscretization() == SpatialDiscretization::Staggered}}
                   .dump()
            << '\n';
    }

    inline void requireGaugeEnvelope(const GaugeEnvelope& envelope, bool nonstandard,
                                     bool gaugeWave) {
        ASSERT_GT(envelope.first_m, 0);
        ASSERT_TRUE(std::isfinite(envelope.last_m));
        std::cout << "B5 stencil=" << (nonstandard ? "NS" : "Std")
                  << " family=" << (gaugeWave ? "T2" : "T3") << " first=" << envelope.first_m
                  << " last=" << envelope.last_m << " ratio=" << envelope.last_m / envelope.first_m
                  << " steps=" << envelope.steps_m << '\n';
        EXPECT_LE(envelope.last_m / envelope.first_m, 1.05);
    }

    /// Construct the discrete wave from the solver's actual step and component positions.
    template <class SolverType>
    GaugeWave growthWave(SolverType& solver) {
        const double pi = Kokkos::numbers::pi_v<double>;
        const auto h    = solver.A_n.get_mesh().getMeshSpacing();
        GaugeWave discrete({2 * pi, 4 * pi, 4 * pi}, h, solver.getDt(), SolverType::stencil(),
                           EBReconstructionConvention::R2);
        if (solver.spatialDiscretization() == SpatialDiscretization::Staggered)
            for (unsigned d = 0; d < SpatialDim; ++d)
                discrete.g_m[d] = staggeredDerivativeSymbol(discrete.waveNumber_m[d], h[d]);
        return discrete;
    }

    /// Frozen growth criterion plus the independent two-sided discrete-wave development check.
    template <class SolverType>
    void checkGrowth(SolverType& solver, const GaugeWave& discrete, const GaugeEnvelope& result,
                     bool gaugeWave) {
        if (solver.spatialDiscretization() == SpatialDiscretization::Collocated) {
            requireGaugeEnvelope(result, SolverType::stencil() == Stencil::NonStandard, gaugeWave);
            if (gaugeWave) {
                EXPECT_LE(std::abs(result.initial_m / exactGaugeNorm(solver, discrete) - 1), 1e-12);
                EXPECT_LE(std::abs(result.last_m / result.first_m - 1), 1e-12);
            }
        }
        recordGaugeHistory(solver, result, gaugeWave);
    }

    template <class PotentialSolverBase>
    void gaugeGrowth(bool gaugeWave, SpatialDiscretization spatialDiscretization =
                                         SpatialDiscretization::Collocated) {
        constexpr bool NS = std::is_same_v<PotentialSolverBase, NonStandardSolverPeriodic>;
        auto fields       = makeBoxFieldContainer({32, 32, NS ? 48 : 32}, {1, 1, 1}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               EBReconstructionConvention::R2,
                                               spatialDiscretization);
        const auto discrete = growthWave(solver);
        const auto mixed    = makeMixedWave();
        if (gaugeWave)
            loadHistory(solver, discrete);
        else
            loadHistory(solver, mixed);
        solver.evaluateInitialFields();
        const double period =
            2 * Kokkos::numbers::pi_v<double> / (gaugeWave ? discrete.omega_m : mixed.omega_m);
        checkGrowth(solver, discrete, gaugeEnvelope(solver, period), gaugeWave);
    }

    template <class PotentialSolverBase>
    void bothGaugeGrowthCases() {
        gaugeGrowth<PotentialSolverBase>(true);
        gaugeGrowth<PotentialSolverBase>(false);
    }
    TEST(LongGaugeGrowth, Standard) {
        bothGaugeGrowthCases<StandardSolverPeriodic>();
    }
    TEST(LongGaugeGrowth, NonStandard) {
        bothGaugeGrowthCases<NonStandardSolverPeriodic>();
    }
    TEST(LongGaugeGrowth, StaggeredRecord) {
        gaugeGrowth<StandardSolverPeriodic>(true, SpatialDiscretization::Staggered);
    }
}  // namespace chdr::solver::test
