/** @file E09.h
 * @brief E05 and E06 in vacuum with three solver configurations; E09.typ and E10.typ.
 */
#pragma once

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E02.h"
#include "experiments/E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/E05.h"

namespace chdr::experiments::e09e14 {
    /// Name of the executable that runs E09 to E14, written to every run record.
    inline constexpr const char* SolverComparisonApplication =
        "experiments-comparison-with-ippl-standard-solver";

    /// Four transverse cells per axis instead of eight: E05 showed transverse independence.
    inline constexpr int VacuumPlaneTransverseCells = 4;

    /// E05 in vacuum; the input selects the solver configuration.
    inline constexpr e01e02::PulseProtocol E09{"E09",
                                               22,
                                               true,
                                               VacuumPlaneTransverseCells,
                                               e05e08::E05CellsPerPulseLevels,
                                               e01e02::E01Targets,
                                               e05e08::E05CellsPerPulseLevels.size(),
                                               false,
                                               SolverComparisonApplication,
                                               1.0,
                                               true};

    /// E06 in vacuum: Mur faces with speed one; the input selects the solver configuration.
    inline constexpr e01e02::PulseProtocol E10{"E10",
                                               3,
                                               false,
                                               VacuumPlaneTransverseCells,
                                               e05e08::E05CellsPerPulseLevels,
                                               e01e02::E02Targets,
                                               e05e08::E05CellsPerPulseLevels.size(),
                                               false,
                                               SolverComparisonApplication,
                                               1.0,
                                               true};
}  // namespace chdr::experiments::e09e14
