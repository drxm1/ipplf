/** @file E06.h
 * @brief E02 with the Gaussian profile; E06.typ defines the changes.
 */
#pragma once

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E02.h"
#include "experiments/E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/E05.h"

namespace chdr::experiments::e05e08 {
    /** @brief E02's box, Mur faces, targets and fixed transverse mesh on the E05 ladder.
     * @details E06.typ drops E02's comparison with the periodic run, so the study bundles
     * no E05 study.
     */
    inline constexpr e01e02::PulseProtocol E06{"E06",
                                               3,
                                               false,
                                               8,
                                               E05CellsPerPulseLevels,
                                               e01e02::E02Targets,
                                               E05CellsPerPulseLevels.size(),
                                               false,
                                               GaussianPulsesApplication};
}  // namespace chdr::experiments::e05e08
