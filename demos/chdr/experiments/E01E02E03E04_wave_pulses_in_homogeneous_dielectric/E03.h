/** @file E03.h
 * @brief Entry point and physical parameters for the E03 spherical pulse.
 */
#pragma once

#include <array>

#include "core/input/RunInput.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03AnalyticContinuumReference.h"

namespace chdr::experiments::e03 {
    namespace reference =
        E01E02E03E04_wave_pulses_in_homogeneous_dielectric::E03_analytic_continuum_reference;
    using E03Reference = reference::E03AnalyticContinuumReferenceSampler;

    inline constexpr double E03Epsilon   = 2.13;
    inline constexpr double E03BoxCentre = 6.0;
    inline constexpr std::array E03Levels{96, 144, 216, 324};

    /// Load the continuum histories and measure their evolution through the E03 targets.
    void runE03(const input::ChdrYamlRunInput& request);
}  // namespace chdr::experiments::e03
