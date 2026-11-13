/** @file E02.h
 * @brief Short Mur-exit experiment; E02.typ defines residual norms and the electric probe.
 */
#pragma once

#include <array>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Protocol.h"

namespace chdr::experiments::e01e02 {
    inline constexpr std::array E02CellsPerPulseLevels{32, 64, 128, 256, 512, 1024};
    inline constexpr std::array E02Targets{
        PulseTarget{"travel_1w", 1.0, true, false}, PulseTarget{"travel_2w", 2.0, true, true},
        PulseTarget{"travel_4w", 4.0, false, true}, PulseTarget{"travel_6w", 6.0, false, true},
        PulseTarget{"final", 7.5, false, true}};

    /// E02 fixes eight transverse cells and saves relative errors only before boundary contact.
    inline constexpr PulseProtocol E02{"E02", 3, false, 8, E02CellsPerPulseLevels, E02Targets, 5};
}  // namespace chdr::experiments::e01e02
