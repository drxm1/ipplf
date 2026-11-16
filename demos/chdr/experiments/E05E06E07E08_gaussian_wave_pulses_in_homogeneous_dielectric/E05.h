/** @file E05.h
 * @brief E01 with the Gaussian profile; E05.typ defines the changes.
 */
#pragma once

#include <array>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01.h"

namespace chdr::experiments::e05e08 {
    /// Name of the executable that runs E05 to E08, written to every run record.
    inline constexpr const char* GaussianPulsesApplication =
        "experiments-gaussian-wave-pulses-in-homogeneous-dielectric";

    /// The six levels of E02. E01 needed 2048 and 4096 only for its compact profile.
    inline constexpr std::array E05CellsPerPulseLevels{32, 64, 128, 256, 512, 1024};

    /// E01's geometry, targets and selectable transverse mesh, with the shorter ladder.
    inline constexpr e01e02::PulseProtocol E05{"E05",
                                               22,
                                               true,
                                               0,
                                               E05CellsPerPulseLevels,
                                               e01e02::E01Targets,
                                               E05CellsPerPulseLevels.size(),
                                               false,
                                               GaussianPulsesApplication};
}  // namespace chdr::experiments::e05e08
