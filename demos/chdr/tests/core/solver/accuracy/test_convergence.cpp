/** @file test_convergence.cpp
 * @brief Check refinement of the source-free and sourced Lorenz-wave families.
 * @details Each Convergence test runs one refinement ladder of Convergence.h for one
 * stencil, one wave family and one reconstruction: collocated, default step, one period,
 * split along x, Standard on \f$N^3\f$ or NonStandard on \f$N\times N\times 3N/2\f$ cells of
 * the unit box, \f$N=16,32,64\f$. The reference is the continuum wave of LorenzWave.h: the
 * source-free mixed wave or the sourced Lorenz-gauge wave, both with
 * \f$\mathbf k=2\pi(1,2,2)\f$. On the pair (32, 64), the observed order
 * \f$\log_2(e_{32}/e_{64})\f$ of the relative \f$L^2\f$ error of B must be at least 1.8. For
 * R2 and R2Prime, the E and Gauss orders must also be at least 1.8; for the sourced wave, the
 * gauge order too when its fine error exceeds \f$10^{-12}\f$. R1 has no E or Gauss bound
 * here; its first-order control runs on the pair (64, 128) in test_convergence_long.cpp. The
 * final identities must hold at \f$10^{-12}\f$ on every level, and the states at
 * \f$N=16\f$ and 32 are written as rank records for chdr.solver.ranks. CTest runs the tests as
 * chdr.solver.convergence.{Standard,NonStandard}.{Mixed,Sourced}.{R1,R2,R2Prime}.{serial,mpi2}.
 *
 * WaveNorms (chdr.solver.wave_norms.{serial,mpi2}) compares the gauge ratio of waveErrors
 * with an independent Fourier calculation, without a solve.
 *
 * The medium suites use R2, \f$\varepsilon=2.13\f$, safety factor \f$\sqrt3/2\f$, Standard
 * \f$N^3\f$ cells and both readings. MediumConvergence evolves six source-free plane waves at
 * speed \f$1/n\f$; MediumSource evolves one driven Lorenz-gauge wave. The E and B orders
 * must be at least 1.8 on (32, 64), and also on (64, 128) in the Long suites; the Gauss
 * order, and for MediumSource the gauge order, only when its fine error exceeds
 * \f$10^{-12}\f$. At \f$N=32\f$, each ladder repeats the run on MPI_COMM_SELF and compares
 * the fields: cells to \f$10^{-13}\f$, statistics to \f$10^{-12}\f$ (SerialReference.h).
 * CTest entries, each serial and mpi2: chdr.solver.medium.Convergence.{Collocated,Staggered},
 * chdr.solver.medium.Source, chdr.solver.medium.long.Convergence.{Collocated,Staggered} and
 * chdr.solver.medium.long.Source.{Collocated,Staggered}. MediumAnisotropic.Fields runs as
 * chdr.solver.medium.long.anisotropic.Fields.mpi2.
 */
#include <complex>

#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "tests/core/solver/accuracy/Convergence.h"
#include "tests/core/solver/accuracy/LorenzWave.h"

static_assert(
    chdr::SpatialDim == 3,
    "test_convergence.cpp requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    /** @brief Fourier oracle for the gauge norm of continuum-sampled histories.
     * @details Substitute the analytic mode into the matched time quotient and
     * centred divergence [christlieb2024GaugeConserving, Eq. (29)]. Orthogonality
     * cancels the common cosine norm, leaving a complex-amplitude ratio
     * [heinzel2002SpectrumSpectralDensity, section 5, Eqs. (3)-(5)]; this is a
     * deduction for these sampled histories, independent of WaveNorms reductions.
     */
    template <class SolverType>
    double loadedGaugeRatio(SolverType& solver, const LorenzWave& wave) {
        using Complex = std::complex<double>;
        const Complex imaginary(0, 1);
        Complex divergence = 0;
        const auto h       = solver.A_n.get_mesh().getMeshSpacing();
        for (unsigned d = 0; d < SpatialDim; ++d)
            divergence += imaginary * (std::sin(wave.waveNumber_m[d] * h[d]) / h[d])
                          * Complex(wave.realA_m[d], wave.imagA_m[d]);
        const auto current  = std::polar(1.0, -wave.omega_m * solver.labelTime(Quantity::Phi, 0));
        const auto previous = std::polar(1.0, -wave.omega_m * solver.labelTime(Quantity::Phi, -1));
        const double weight = solver.levelContract().gaugeWeight();
        const Complex gauge =
            Complex(wave.phiAmplitude(false), wave.phiAmplitude(true)) * (current - previous)
                / solver.getDt()
            + divergence * (weight + (1 - weight) * std::polar(1.0, wave.omega_m * solver.getDt()));
        return std::abs(gauge) / std::abs(divergence);
    }

    /// Check R1, R2, R2Prime against an independent Fourier calculation, on 16 by 16 by 16
    /// cells with sourced-wave histories and no solve:
    ///   - The predicted gauge ratio must exceed 1e-4, so that the comparison is not trivial.
    ///   - The relative difference of measured and predicted ratio must be at most 1e-12.
    TEST(WaveNorms, NonzeroGaugeMatchesIndependentMode) {
        auto fields = makeBoxFieldContainer({16, 16, 16}, {1, 1, 1}, 0);
        initializeFieldStorage(fields);
        for (auto ebReconstructionConvention :
             {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
              EBReconstructionConvention::R2Prime}) {
            ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                      ebReconstructionConvention);
            const auto wave = makeSourcedWave();
            loadHistory(solver, wave);
            solver.evaluateInitialFields();
            const double predicted = loadedGaugeRatio(solver, wave);
            ASSERT_GT(predicted, 1e-4);
            EXPECT_LE(std::abs(waveErrors(solver, wave).gauge_m / predicted - 1), 1e-12);
        }
    }
    /// Standard, mixed wave, R1: B order at least 1.8 on (32, 64).
    TEST(Convergence, StandardMixedR1) {
        waveLadder<StandardSolverPeriodic>(false, EBReconstructionConvention::R1, false);
    }
    /// Standard, mixed wave, R2: E, B and Gauss orders at least 1.8 on (32, 64).
    TEST(Convergence, StandardMixedR2) {
        waveLadder<StandardSolverPeriodic>(false, EBReconstructionConvention::R2, false);
    }
    /// Standard, mixed wave, R2Prime: E, B and Gauss orders at least 1.8 on (32, 64).
    TEST(Convergence, StandardMixedR2Prime) {
        waveLadder<StandardSolverPeriodic>(false, EBReconstructionConvention::R2Prime, false);
    }
    /// Standard, sourced wave, R1: B and gauge orders at least 1.8 on (32, 64).
    TEST(Convergence, StandardSourcedR1) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R1, false);
    }
    /// Standard, sourced wave, R2: E, B, Gauss and gauge orders at least 1.8 on (32, 64).
    TEST(Convergence, StandardSourcedR2) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R2, false);
    }
    /// Standard, sourced wave, R2Prime: E, B, Gauss and gauge orders at least 1.8 on (32, 64).
    TEST(Convergence, StandardSourcedR2Prime) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R2Prime, false);
    }
    /// NonStandard, mixed wave, R1: B order at least 1.8 on (32, 64).
    TEST(Convergence, NonStandardMixedR1) {
        waveLadder<NonStandardSolverPeriodic>(false, EBReconstructionConvention::R1, false);
    }
    /// NonStandard, mixed wave, R2: E, B and Gauss orders at least 1.8 on (32, 64).
    TEST(Convergence, NonStandardMixedR2) {
        waveLadder<NonStandardSolverPeriodic>(false, EBReconstructionConvention::R2, false);
    }
    /// NonStandard, mixed wave, R2Prime: E, B and Gauss orders at least 1.8 on (32, 64).
    TEST(Convergence, NonStandardMixedR2Prime) {
        waveLadder<NonStandardSolverPeriodic>(false, EBReconstructionConvention::R2Prime, false);
    }
    /// NonStandard, sourced wave, R1: B and gauge orders at least 1.8 on (32, 64).
    TEST(Convergence, NonStandardSourcedR1) {
        waveLadder<NonStandardSolverPeriodic>(true, EBReconstructionConvention::R1, false);
    }
    /// NonStandard, sourced wave, R2: E, B, Gauss and gauge orders at least 1.8 on (32, 64).
    TEST(Convergence, NonStandardSourcedR2) {
        waveLadder<NonStandardSolverPeriodic>(true, EBReconstructionConvention::R2, false);
    }
    /// NonStandard, sourced wave, R2Prime: E, B, Gauss and gauge orders at least 1.8 on
    /// (32, 64).
    TEST(Convergence, NonStandardSourcedR2Prime) {
        waveLadder<NonStandardSolverPeriodic>(true, EBReconstructionConvention::R2Prime, false);
    }

    /** @brief Refine each medium wave without changing its physical direction or polarization.
     * @details Source-free cases span axial, face-diagonal and body-diagonal
     * propagation, each with transverse-only and mixed potentials. Sourced cases
     * use one fixed driven wave from initializeWave. longRun includes the finest
     * 128-cell level, so both required refinement pairs are checked.
     * Collocated and staggered readings share the wave but compare their outputs
     * at the different component positions prescribed by LevelContract.
     */
    void mediumLadders(bool longRun, bool sourced, SpatialDiscretization reading) {
        for (int direction = 0; direction < (sourced ? 1 : 3); ++direction)
            for (int gamma = sourced ? 1 : 0; gamma <= 1; ++gamma)
                waveLadder<StandardSolverPeriodic>(sourced, EBReconstructionConvention::R2, longRun,
                                                   reading, false, false, 2.13, direction, gamma);
    }
    /// Six plane waves at \f$\varepsilon=2.13\f$, collocated: E and B orders at least 1.8 on
    /// (32, 64).
    TEST(MediumConvergence, Collocated) {
        mediumLadders(false, false, SpatialDiscretization::Collocated);
    }
    /// As MediumConvergence.Collocated in the staggered reading.
    TEST(MediumConvergence, Staggered) {
        mediumLadders(false, false, SpatialDiscretization::Staggered);
    }
    /// The collocated plane-wave ladders up to \f$N=128\f$: orders at least 1.8 on (32, 64)
    /// and (64, 128).
    TEST(MediumConvergenceLong, Collocated) {
        mediumLadders(true, false, SpatialDiscretization::Collocated);
    }
    /// As MediumConvergenceLong.Collocated in the staggered reading.
    TEST(MediumConvergenceLong, Staggered) {
        mediumLadders(true, false, SpatialDiscretization::Staggered);
    }
    // This ladder samples analytic physical charge/current at each step; it does
    // not manufacture a discrete forcing as test_manufactured.cpp does. The final
    // field errors therefore include source discretization and accumulated evolution.
    /// Driven wave at \f$\varepsilon=2.13\f$, collocated: E and B orders at least 1.8 on
    /// (32, 64), Gauss and gauge orders too when their fine errors exceed \f$10^{-12}\f$.
    TEST(MediumSource, Collocated) {
        mediumLadders(false, true, SpatialDiscretization::Collocated);
    }
    /// As MediumSource.Collocated in the staggered reading.
    TEST(MediumSource, Staggered) {
        mediumLadders(false, true, SpatialDiscretization::Staggered);
    }
    /// The collocated driven-wave ladder up to \f$N=128\f$, with the rules on both pairs.
    TEST(MediumSourceLong, Collocated) {
        mediumLadders(true, true, SpatialDiscretization::Collocated);
    }
    /// As MediumSourceLong.Collocated in the staggered reading.
    TEST(MediumSourceLong, Staggered) {
        mediumLadders(true, true, SpatialDiscretization::Staggered);
    }
    /** @brief Independent anisotropic refinement with both component placements.
     * @details The anisotropicCase flag in runWave selects an anisotropic mesh family
     * and y decomposition. Each level compares evolved E/B with continuum
     * samples at their labels and repeats the evolution on MPI_COMM_SELF.
     * The separate anisotropic phase test initializes a discrete single-root mode;
     * these field runs retain continuum initial histories.
     */
    TEST(MediumAnisotropic, Fields) {
        for (auto reading : {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered})
            waveLadder<StandardSolverPeriodic>(false, EBReconstructionConvention::R2, true, reading,
                                               true, false, 3.7);
    }
}  // namespace chdr::solver::test
