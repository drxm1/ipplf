/** @file test_convergence.cpp
 * @brief Check refinement of the source-free and sourced Lorenz-wave families.
 *
 * TODO: REVIEW.
 */
#include <complex>

#include "Convergence.h"
#include "SpatialDimension.h"

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
        auto fields = makePeriodicFieldContainer({16, 16, 16}, {1, 1, 1}, 0);
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
}  // namespace chdr::solver::test
