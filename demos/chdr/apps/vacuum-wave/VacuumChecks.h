/** @file VacuumChecks.h
 * @brief Fixed benchmark acceptance, separate from measurement and solver parameters.
 * @details The public README documents our selected thresholds. References start from
 * Fallahi, MITHRA 2.0, arXiv:2009.13645v1, Eqs. (3.6)-(3.20); see VacuumDiagnostics.h.
 */
#pragma once

#include <cmath>
#include <stdexcept>

#include "VacuumDiagnostics.h"

namespace chdr::vacuum {
    // Fixed benchmark gates documented in README.md; these do not tune the solver.
    constexpr double DiscreteTolerance     = 1e-10;
    constexpr double ZeroTolerance         = 1e-12;
    constexpr double MinimumPhaseAmplitude = 0.5;
    constexpr double PhaseSpeedTolerance   = 0.002;
    constexpr double CoarseFieldTolerance  = 0.01;
    constexpr double FineFieldTolerance    = 0.003;

    /// @brief Apply the benchmark's amplitude, direction and mean-speed policy.
    inline bool acceptedPhase(const chdr::vacuum::Diagnostics& d,
                              const chdr::vacuum::PhaseSample& phase) {
        if (!phase.valid_m || !(phase.amplitude_m > MinimumPhaseAmplitude)
            || phase.hasSpeed_m != (d.step_m > 0))
            return false;
        return !phase.hasSpeed_m
               || (phase.increment_m > -Kokkos::numbers::pi_v<double> && phase.increment_m < 0.0
                   && std::abs(phase.speed_m - 1.0) <= PhaseSpeedTolerance);
    }

    /// @brief Enforce the fixed acceptance contract documented in README.md.
    inline void validateStep(const chdr::vacuum::Diagnostics& d,
                             const chdr::vacuum::PhaseSample& phase, bool checkpoint, bool fine) {
        if (!d.finite_m || d.discreteError_m > DiscreteTolerance || d.zeroError_m > ZeroTolerance
            || !acceptedPhase(d, phase))
            throw std::runtime_error("Discrete, finite-value, zero-channel or phase check failed");
        const double tolerance = fine ? FineFieldTolerance : CoarseFieldTolerance;
        if (checkpoint && (d.relativeE_m > tolerance || d.relativeB_m > tolerance))
            throw std::runtime_error("Continuum field error exceeded its fixed tolerance");
    }

}  // namespace chdr::vacuum
