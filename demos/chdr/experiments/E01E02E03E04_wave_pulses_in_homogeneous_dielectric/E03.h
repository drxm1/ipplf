/** @file E03.h
 * @brief Physical parameters of the E03 and E04 spherical pulse.
 */
#pragma once

#include <array>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03AnalyticContinuumReference.h"

namespace chdr::experiments::e03 {
    namespace reference =
        E01E02E03E04_wave_pulses_in_homogeneous_dielectric::E03_analytic_continuum_reference;
    /// The compact dipole pulse of E03 and E04. The run path accepts other samplers of this form.
    using E03Reference = reference::E03AnalyticContinuumReferenceSampler<>;

    inline constexpr double E03Epsilon   = 2.13;
    inline constexpr double E03BoxCentre = 6.0;
    inline constexpr std::array E03Levels{96, 144, 216, 324};

    /** @brief Names that a spherical pulse run writes to its run record.
     * @details The measurement identifiers in the CSV files stay those of E03.typ and E04.typ,
     * because the same definitions apply to every experiment that uses this run path.
     */
    struct SphericalExperiment {
        /// Experiment name, for example "E03".
        const char* name_m;
        /// Name of the executable that runs the experiment.
        const char* application_m;
    };
}  // namespace chdr::experiments::e03
