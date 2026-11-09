/** @file E01.h
 * @brief Periodic long-propagation experiment; E01.typ defines the physical protocol.
 */
#pragma once

#include <array>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Protocol.h"

namespace chdr::experiments::e01e02 {
    inline constexpr std::array E01CellsPerPulseLevels{32, 64, 128, 256, 512, 1024, 2048, 4096};
    inline constexpr std::array E01Targets{
        PulseTarget{"travel_2w", 2.0, true, false}, PulseTarget{"quarter_1", 10.0, true, false},
        PulseTarget{"quarter_2", 20.0, true, false}, PulseTarget{"quarter_3", 30.0, true, false},
        PulseTarget{"final", 40.0, true, false}};

    /// E01 retains its configurable transverse count and original quarter-run frame cadence.
    inline constexpr PulseProtocol E01{"E01", 22, true, 0, E01CellsPerPulseLevels, E01Targets, 4};
}  // namespace chdr::experiments::e01e02
