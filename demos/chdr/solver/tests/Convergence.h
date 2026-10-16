/** @file Convergence.h
 * @brief Execute source-free and sourced Lorenz-wave refinement ladders.
 * @details WaveNorms.h defines continuum comparisons. For vacuum cases, Monitors
 * verifies the independent discrete identities on the last evolved step. Medium
 * cases use WaveNorms.h's material Gauss and gauge residuals instead. RankRecords.h
 * records owned fields for decomposition comparisons; approximation errors
 * remain separate JSONL records.
 * Helpers parameterized by PotentialSolverBase select the inherited Standard or
 * NonStandard potential update. SolverType instead denotes the complete ChdrSolver
 * interface, including its potential histories, source field and E/B outputs.
 *
 * TODO: REVIEW.
 */
#pragma once
#include <gtest/gtest.h>
#include <iostream>

#include "BoxSetup.h"
#include "ChdrSolver.h"
#include "CoreTestSupport.h"
#include "Monitors.h"
#include "RankRecords.h"
#include "SerialReference.h"
#include "SpatialDimension.h"
#include "WaveNorms.h"
#include "Waves.h"

static_assert(chdr::solver::SpatialDim == 3,
              "Convergence.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    /** @brief Select one analytic wave and the solver representation used to evolve it.
     * @details epsilon_m equal to zero selects the owner-free vacuum path.
     * For source-free medium waves, direction_m selects axial, face-diagonal or
     * body-diagonal propagation (0, 1 or 2). gamma_m selects a transverse wave
     * (0) or one with an additional longitudinal potential (1); makeMediumWave
     * supplies the corresponding scalar potential through the independent oracle.
     * Sourced and held cases have their own fixed wave vectors.
     */
    struct WaveOptions {
        bool sourced_m;
        EBReconstructionConvention ebReconstructionConvention_m;
        SpatialDiscretization spatialDiscretization_m;
        bool holdout_m, recordOnly_m;
        double epsilon_m = 0;
        int direction_m  = 0;
        double gamma_m   = 1;
    };

    /** @brief Return ordinary or anisotropic mesh counts for one refinement level.
     * @details All boxes have unit lengths. The NonStandard mesh has
     * \f$N_z=3n/2\f$ and \f$N_x=N_y=n\f$, giving transverse spacings
     * \f$1.5h_z\f$. The independent anisotropic family instead uses
     * \f$(N_x,N_y,N_z)=(5n/4,25n/16,n)\f$, so its spacing proportions are
     * \f$h_x:h_y:h_z=1:0.8:1.25\f$. These explicit geometric inputs stay fixed
     * throughout each refinement family; the coarsest anisotropic y count is odd.
     */
    template <class PotentialSolverBase>
    std::array<int, SpatialDim> waveCells(int n, bool holdout) {
        constexpr bool NS = std::is_same_v<PotentialSolverBase, NonStandardSolverPeriodic>;
        return holdout ? std::array<int, SpatialDim>{5 * n / 4, 25 * n / 16, n}
                       : std::array<int, SpatialDim>{n, n, NS ? 3 * n / 2 : n};
    }

    /** @brief Load continuum samples at the shell's labels and reconstruct initial fields.
     * @details Select the medium step before sampling either initial history slot.
     * makeMediumWave supplies continuum frequency and amplitudes independently
     * of the discrete update. Their evolution error is therefore retained in
     * this convergence check; the dispersion test instead initializes a discrete mode.
     */
    template <class SolverType>
    LorenzWave initializeWave(SolverType& solver, const WaveOptions& options) {
        if (options.holdout_m)
            solver.setSafetyFactor(0.9);
        auto wave       = options.sourced_m ? makeSourcedWave() : makeMixedWave();
        const double pi = Kokkos::numbers::pi_v<double>;
        if (options.holdout_m)
            wave = makeSourcedWave({4 * pi, -2 * pi, 6 * pi}, 4 * pi);
        if (options.epsilon_m > 0) {
            solver.setSafetyFactor(options.holdout_m ? 0.9 : std::sqrt(3.0) / 2);
            const ippl::Vector<double, SpatialDim> mode =
                options.holdout_m   ? ippl::Vector<double, SpatialDim>(-1, 2, 3)
                : options.sourced_m ? ippl::Vector<double, SpatialDim>(1, 2, 2)
                                    : ippl::Vector<double, SpatialDim>(1, options.direction_m > 0,
                                                                       options.direction_m > 1);
            wave = makeMediumWave(2 * pi * mode, options.epsilon_m, options.gamma_m,
                                  options.sourced_m);
        }
        wave.requirePeriodic({1, 1, 1});
        loadHistory(solver, wave);
        solver.evaluateInitialFields();
        return wave;
    }

    /** @brief Advance the actual IPPL solver to the nearest step of one period.
     * @details The comparison uses the time actually reached, not a nominal
     * period substituted afterward. Sources are written before each solve.
     * Monitoring the final transition retains its consumed source and previous
     * history without allocating diagnostics for every intermediate step.
     * A null monitors pointer selects the medium path, where the vacuum identity
     * formulas do not apply. waveErrors still measures material continuum residuals.
     */
    template <class SolverType>
    void advanceWave(SolverType& solver, const LorenzWave& wave, Monitors* monitors) {
        const double pi = Kokkos::numbers::pi_v<double>;
        const int steps = std::lround(2 * pi / wave.omega_m / solver.getDt());
        for (int k = 0; k < steps; ++k) {
            writeSource(solver, wave.sources(), k);
            if (k == steps - 1 && monitors)
                monitors->advance(solver);
            else
                solver.solve();
        }
    }

    /** @brief Check stored identities unless the caller requests a record-only case.
     * @details The staggered holdout uses recordOnly and skips these assertions.
     * Otherwise the normalized magnetic-divergence, Faraday and Gauss-formula
     * mismatches must each be at most \f$10^{-12}\f$.
     */
    inline void requireWaveIdentities(const MonitorReport& report, bool recordOnly) {
        if (recordOnly)
            return;
        EXPECT_LE(report.divMagnetic_m.cell_m, 1e-12);
        EXPECT_LE(report.divRawMagnetic_m.cell_m, 1e-12);
        for (int d = 0; d < static_cast<int>(SpatialDim); ++d)
            EXPECT_LE(report.faraday_m[d].cell_m, 1e-12);
        EXPECT_LE(report.gauss_m.cell_m, 1e-12);
    }

    /** @brief Stable identifiers for source-free, sourced and held wave records.
     * @details T3 means the mixed source-free wave, T4 the sourced Lorenz wave,
     * and H1 the independently held anisotropic case. The remaining fields are
     * stencil, reconstruction enum, reading enum and refinement count.
     */
    template <class SolverType>
    std::string waveCaseName(int n, const WaveOptions& options) {
        if (options.epsilon_m > 0)
            return std::string(options.holdout_m   ? "MediumHeld"
                               : options.sourced_m ? "MediumSource"
                                                   : "MediumWave")
                   + "-" + std::to_string(options.direction_m) + "-"
                   + std::to_string(int(options.gamma_m)) + "-"
                   + std::to_string(int(options.spatialDiscretization_m)) + "-" + std::to_string(n);
        return std::string(options.holdout_m   ? "H1"
                           : options.sourced_m ? "T4"
                                               : "T3")
               + "-" + (SolverType::stencil() == Stencil::NonStandard ? "NS" : "Std") + "-"
               + std::to_string(static_cast<int>(options.ebReconstructionConvention_m)) + "-"
               + std::to_string(static_cast<int>(options.spatialDiscretization_m)) + "-"
               + std::to_string(n);
    }

    /// Keep continuum errors at every mesh, separately from the field-comparison JSON inputs.
    template <class SolverType>
    void writeWaveErrors(const std::string& name, const WaveErrors& errors, SolverType& solver) {
        std::filesystem::create_directories(recordDirectory());
        std::ofstream out(recordDirectory() / (name + "-errors.jsonl"));
        out << nlohmann::json{{"electric", errors.electric_m},
                              {"magnetic", errors.magnetic_m},
                              {"gauss", errors.gauss_m},
                              {"gauge", errors.gauge_m},
                              {"epsilon", solver.homogeneousPermittivity()},
                              {"formulation", solver.formulationName()},
                              {"alpha", solver.formulationParameters().alpha_m}}
                   .dump()
            << '\n';
    }

    /// Keep final-step identity magnitudes beside, but separate from, continuum errors.
    inline void writeWaveIdentities(const std::string& name, const MonitorReport& report) {
        if (ippl::Comm->rank() != 0)
            return;
        std::ofstream out(recordDirectory() / (name + "-identities.jsonl"));
        std::vector<double> faraday, faradayGlobal;
        for (const auto& component : report.faraday_m) {
            faraday.push_back(component.cell_m);
            faradayGlobal.push_back(component.global_m);
        }
        out << nlohmann::json{{"div_B_cell", report.divMagnetic_m.cell_m},
                              {"div_raw_B_cell", report.divRawMagnetic_m.cell_m},
                              {"faraday_cell", faraday},
                              {"div_B_global", report.divMagnetic_m.global_m},
                              {"div_raw_B_global", report.divRawMagnetic_m.global_m},
                              {"faraday_global", faradayGlobal},
                              {"gauss_formula_global", report.gauss_m.global_m},
                              {"gauss_formula_cell", report.gauss_m.cell_m},
                              {"gauge_l2", report.gauge_m.l2_m}}
                   .dump()
            << '\n';
    }

    /** @brief Record the driven-mode amplitude prediction beside the near-field factor.
     * @details Substitute a Fourier particular solution into MITHRA's discrete
     * wave equation (deduction from [fallahi2020mithra20fullwavesimulation,
     * Eqs. (3.17), (3.24)]). The first number is the exact discrete particular
     * solution's relative amplitude error at this finite frequency. The second
     * is the near-light-speed longitudinal limit, not an expected equal number:
     * as the temporal and longitudinal symbols cancel, the transverse operator
     * is multiplied by Sz and its driven amplitude by 1/Sz (deduction).
     * Measured errors also include the homogeneous response to the initial slots.
     */
    template <class SolverType>
    void writeSourceAmplitude(const std::string& name, SolverType& solver, const LorenzWave& wave) {
        const auto h = solver.A_n.get_mesh().getMeshSpacing(), k = wave.waveNumber_m;
        const double lambda = -laplacianSymbol(k, h, Stencil::NonStandard);
        const double w      = timeSymbol(wave.omega_m, solver.getDt());
        const double alpha =
            0.25 * (1 + 0.02 / (std::pow(h[2] / h[0], 2) + std::pow(h[2] / h[1], 2)));
        const double sz = 1 - 4 * alpha * std::pow(std::sin(k[2] * h[2] / 2), 2);
        std::ofstream out(recordDirectory() / (name + "-source-amplitude.jsonl"));
        out << nlohmann::json{{"prediction_particular_amplitude_error",
                               wave.sourceFactor_m / (lambda - w * w) - 1},
                              {"comparison_inverse_Sz_minus_one", 1 / sz - 1}}
                   .dump()
            << '\n';
    }

    template <class SolverType>
    void recordWave(int n, const WaveOptions& options, SolverType& solver, const LorenzWave& wave,
                    const WaveErrors& errors) {
        const auto name = waveCaseName<SolverType>(n, options);
        std::cout << name << " E=" << errors.electric_m << " B=" << errors.magnetic_m
                  << " Gauss=" << errors.gauss_m << " G/divA=" << errors.gauge_m << '\n';
        if (n <= 32)
            recordState(name, solver);
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank != 0)
            return;
        writeWaveErrors(name, errors, solver);
        if (SolverType::stencil() == Stencil::NonStandard && options.sourced_m)
            writeSourceAmplitude(name, solver, wave);
    }

    /** @brief Run one mesh while keeping the field owner at a stable address.
     * @details Helpers receive references to this owner and solver; neither is
     * copied or moved after allocation. The source, mesh and layout therefore
     * outlive all kernels and monitor work fields. The independent anisotropic case is registered
     * separately as a disabled test and never entered by a development ladder.
     * For a medium holdout, n is a refinement label rather than an axis count:
     * n=16 gives 15 by 18 by 20 cells, doubled at each subsequent level.
     * Its y decomposition differs from the development cases' x decomposition.
     */
    template <class PotentialSolverBase>
    WaveErrors runWave(
        int n, bool sourced, EBReconstructionConvention ebReconstructionConvention,
        SpatialDiscretization spatialDiscretization = SpatialDiscretization::Collocated,
        bool holdout = false, bool recordOnly = false, double eps = 0, int direction = 0,
        double gamma = 1) {
        const WaveOptions options{sourced,
                                  ebReconstructionConvention,
                                  spatialDiscretization,
                                  holdout,
                                  recordOnly,
                                  eps,
                                  direction,
                                  gamma};
        const auto cells = eps > 0 && holdout
                               ? std::array<int, 3>{15 * n / 16, 18 * n / 16, 20 * n / 16}
                               : waveCells<PotentialSolverBase>(n, holdout);
        auto fields      = makeBoxFieldContainer(cells, {1, 1, 1}, eps > 0 && holdout ? 1 : 0);
        initializeFieldStorage(fields);
        auto handle = coreTest::makeTestSolver<PotentialSolverBase>(
            fields.getJ(), fields.getE(), fields.getB(), ebReconstructionConvention,
            spatialDiscretization, eps);
        auto& solver    = *handle;
        const auto wave = initializeWave(solver, options);
        std::unique_ptr<Monitors> monitors;
        if (eps == 0)
            monitors = std::make_unique<Monitors>(solver.A_n);
        advanceWave(solver, wave, monitors.get());
        // Repeat the 32-cell medium cases on MPI_COMM_SELF to isolate decomposition
        // effects from continuum error. Held cases request this comparison at every level.
        if (holdout || (eps > 0 && n == 32))
            verifySerialReference(solver, wave, int(solver.completedSolves()),
                                  spatialDiscretization, !recordOnly);
        if (monitors)
            requireWaveIdentities(monitors->report(), recordOnly);
        const auto errors = waveErrors(solver, wave);
        recordWave(n, options, solver, wave, errors);
        if (monitors)
            writeWaveIdentities(waveCaseName<ChdrSolver<PotentialSolverBase>>(n, options),
                                monitors->report());
        return errors;
    }

    /** @brief Estimate mesh-halving order with \c std::log2(coarse/fine).
     * @details For an error dominated by one power of mesh spacing, the base-two
     * logarithm of the coarse/fine ratio estimates that power
     * [salari2000CodeVerificationMms, section 4.2].
     * WaveNorms.h defines the errors; requireOrders selects the checked pairs.
     */
    inline double order(double coarse, double fine) {
        return std::log2(coarse / fine);
    }

    inline void recordOrders(const WaveErrors& a, const WaveErrors& b, std::size_t pair) {
        for (double value : {b.electric_m, b.magnetic_m, b.gauss_m, b.gauge_m})
            EXPECT_TRUE(std::isfinite(value));
        std::cout << "order pair=" << (16 << (pair - 1)) << ',' << (16 << pair)
                  << " E=" << order(a.electric_m, b.electric_m)
                  << " B=" << order(a.magnetic_m, b.magnetic_m)
                  << " Gauss=" << order(a.gauss_m, b.gauss_m)
                  << " G=" << order(a.gauge_m, b.gauge_m) << '\n';
    }

    /// R1 is a first-order electric/Gauss control on the finest pair only.
    inline void requireFirstOrderControl(const WaveErrors& a, const WaveErrors& b,
                                         std::size_t pair) {
        if (pair != 3)
            return;
        EXPECT_GE(order(a.electric_m, b.electric_m), 0.8);
        EXPECT_LE(order(a.electric_m, b.electric_m), 1.2);
        EXPECT_GE(order(a.gauss_m, b.gauss_m), 0.8);
        EXPECT_LE(order(a.gauss_m, b.gauss_m), 1.2);
    }

    /** @brief Report refinement orders and enforce the specified pair/control criteria.
     * @details Every call checks that the fine-mesh error values are finite and
     * reports the observed orders. recordOnly adds no order-bound assertions.
     * On waveLadder's 32-to-64 and 64-to-128 pairs, magnetic-field order must be
     * at least 1.8. R2 and R2Prime require the same electric/Gauss order; R1 instead
     * requires electric/Gauss order from 0.8 through 1.2 on 64-to-128 only.
     * In a medium, Gauss order is asserted only when the fine-mesh Gauss error
     * exceeds \f$10^{-12}\f$, using the same floor as the gauge-order check.
     * This floor skips only the residual's order assertion; electric and magnetic
     * field orders remain required on both checked pairs.
     * For the transverse cases with equal nonzero wave-number components on a
     * cubic grid, the discrete Gauss residual vanishes in exact arithmetic:
     * the divergence weights are proportional to the wave number and the
     * polarization is transverse (Deduction from the Gauss constraint
     * [chew2014generalizedgauge, Eq. (4)]). Its measured round-off residual has
     * no defined truncation-error order (Deduction from the order definition
     * [salari2000CodeVerificationMms, section 4.2]).
     * For sourced non-holdout cases, gauge order must also be at least 1.8 when
     * the fine-mesh relative gauge error exceeds \f$10^{-12}\f$.
     * These are end-of-period comparisons, not temporal-envelope bounds. A driven
     * mode plus homogeneous gauge modes can interfere at individual instants
     * (deduction from [christlieb2024GaugeConserving, Lemma 3.3, Eq. (38)]).
     */
    inline void requireOrders(const WaveErrors& a, const WaveErrors& b, std::size_t pair,
                              const WaveOptions& options) {
        recordOrders(a, b, pair);
        if (pair < 2 || options.recordOnly_m)
            return;
        if (options.ebReconstructionConvention_m == EBReconstructionConvention::R1)
            requireFirstOrderControl(a, b, pair);
        else {
            EXPECT_GE(order(a.electric_m, b.electric_m), 1.8);
            if (options.epsilon_m == 0 || b.gauss_m > 1e-12) {
                EXPECT_GE(order(a.gauss_m, b.gauss_m), 1.8);
            }
        }
        EXPECT_GE(order(a.magnetic_m, b.magnetic_m), 1.8);
        if (options.sourced_m && !options.holdout_m && b.gauge_m > 1e-12) {
            EXPECT_GE(order(a.gauge_m, b.gauge_m), 1.8);
        }
    }

    template <class PotentialSolverBase>
    void waveLadder(bool sourced, EBReconstructionConvention ebReconstructionConvention,
                    bool longRun,
                    SpatialDiscretization spatialDiscretization = SpatialDiscretization::Collocated,
                    bool holdout = false, bool recordOnly = false, double eps = 0,
                    int direction = 0, double gamma = 1) {
        const WaveOptions options{sourced,
                                  ebReconstructionConvention,
                                  spatialDiscretization,
                                  holdout,
                                  recordOnly,
                                  eps,
                                  direction,
                                  gamma};
        std::vector<WaveErrors> errors;
        for (int n : {16, 32, 64, 128}) {
            if (n == 128 && !longRun)
                break;
            errors.push_back(runWave<PotentialSolverBase>(n, sourced, ebReconstructionConvention,
                                                          spatialDiscretization, holdout,
                                                          recordOnly, eps, direction, gamma));
        }
        for (std::size_t pair = 1; pair < errors.size(); ++pair)
            requireOrders(errors[pair - 1], errors[pair], pair, options);
    }
}  // namespace chdr::solver::test
