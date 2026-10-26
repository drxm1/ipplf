/** @file test_convergence.cpp
 * @brief Check refinement of the source-free and sourced Lorenz-wave families.
 */
#include <complex>

#include "Convergence.h"
#include "core/base/SpatialDimension.h"

static_assert(
    chdr::solver::SpatialDim == 3,
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

    /// Check R1, R2, R2Prime against an independent Fourier calculation:
    ///   - Predicted gauge ration must exceed 1e-4, avoiding a trivial zero comparison.
    ///   - Measured/reference relative discrepancy must be <= 1e-12.
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
    TEST(Convergence, StandardMixedR1) {
        waveLadder<StandardSolverPeriodic>(false, EBReconstructionConvention::R1, false);
    }
    TEST(Convergence, StandardMixedR2) {
        waveLadder<StandardSolverPeriodic>(false, EBReconstructionConvention::R2, false);
    }
    TEST(Convergence, StandardMixedR2Prime) {
        waveLadder<StandardSolverPeriodic>(false, EBReconstructionConvention::R2Prime, false);
    }
    TEST(Convergence, StandardSourcedR1) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R1, false);
    }
    TEST(Convergence, StandardSourcedR2) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R2, false);
    }
    TEST(Convergence, StandardSourcedR2Prime) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R2Prime, false);
    }
    TEST(Convergence, NonStandardMixedR1) {
        waveLadder<NonStandardSolverPeriodic>(false, EBReconstructionConvention::R1, false);
    }
    TEST(Convergence, NonStandardMixedR2) {
        waveLadder<NonStandardSolverPeriodic>(false, EBReconstructionConvention::R2, false);
    }
    TEST(Convergence, NonStandardMixedR2Prime) {
        waveLadder<NonStandardSolverPeriodic>(false, EBReconstructionConvention::R2Prime, false);
    }
    TEST(Convergence, NonStandardSourcedR1) {
        waveLadder<NonStandardSolverPeriodic>(true, EBReconstructionConvention::R1, false);
    }
    TEST(Convergence, NonStandardSourcedR2) {
        waveLadder<NonStandardSolverPeriodic>(true, EBReconstructionConvention::R2, false);
    }
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
    TEST(MediumConvergence, Collocated) {
        mediumLadders(false, false, SpatialDiscretization::Collocated);
    }
    TEST(MediumConvergence, Staggered) {
        mediumLadders(false, false, SpatialDiscretization::Staggered);
    }
    TEST(MediumConvergenceLong, Collocated) {
        mediumLadders(true, false, SpatialDiscretization::Collocated);
    }
    TEST(MediumConvergenceLong, Staggered) {
        mediumLadders(true, false, SpatialDiscretization::Staggered);
    }
    // This ladder samples analytic physical charge/current at each step; it does
    // not manufacture a discrete forcing as test_manufactured.cpp does. The final
    // field errors therefore include source discretization and accumulated evolution.
    TEST(MediumSource, Collocated) {
        mediumLadders(false, true, SpatialDiscretization::Collocated);
    }
    TEST(MediumSource, Staggered) {
        mediumLadders(false, true, SpatialDiscretization::Staggered);
    }
    TEST(MediumSourceLong, Collocated) {
        mediumLadders(true, true, SpatialDiscretization::Collocated);
    }
    TEST(MediumSourceLong, Staggered) {
        mediumLadders(true, true, SpatialDiscretization::Staggered);
    }
    /** @brief Independent anisotropic refinement with both component placements.
     * @details The holdout flag in runWave selects an anisotropic mesh family
     * and y decomposition. Each level compares evolved E/B with continuum
     * samples at their labels and repeats the evolution on MPI_COMM_SELF.
     * The separate held phase test initializes a discrete single-root mode;
     * these field runs retain continuum initial histories. Disabled by default
     * and absent from the development CTest filters.
     */
    TEST(MediumHoldout, DISABLED_Fields) {
        for (auto reading : {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered})
            waveLadder<StandardSolverPeriodic>(false, EBReconstructionConvention::R2, true, reading,
                                               true, false, 3.7);
    }
}  // namespace chdr::solver::test
