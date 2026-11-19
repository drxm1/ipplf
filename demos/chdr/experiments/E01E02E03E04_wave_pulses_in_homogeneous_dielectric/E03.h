/** @file E03.h
 * @brief Physical parameters of the E03 and E04 spherical pulse.
 */
#pragma once

#include <array>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03AnalyticContinuumReference.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/SolverConfiguration.h"

namespace chdr::experiments::e03 {
    namespace reference =
        E01E02E03E04_wave_pulses_in_homogeneous_dielectric::E03_analytic_continuum_reference;
    /// The compact dipole pulse of E03 and E04. The run path accepts other samplers of this form.
    using E03Reference = reference::E03AnalyticContinuumReferenceSampler<>;

    inline constexpr double E03Epsilon   = 2.13;
    inline constexpr double E03BoxCentre = 6.0;
    inline constexpr std::array E03Levels{96, 144, 216, 324};

    /** @brief The experiment-specific choices of a spherical pulse run.
     * @details The measurement identifiers in the CSV files stay those of E03.typ and E04.typ,
     * because the same definitions apply to every experiment that uses this run path.
     * The reference sampler must use the speed \f$1/\sqrt{\varepsilon_r}\f$ of this medium.
     */
    struct SphericalExperiment {
        /// Experiment name, for example "E03".
        const char* name_m;
        /// Name of the executable that runs the experiment.
        const char* application_m;
        /// Relative permittivity of the homogeneous medium: E03Epsilon, or 1 in vacuum.
        double relativePermittivity_m;
        /// O1VS in E03 and E04; other experiments may select it in the input.
        SolverConfiguration solver_m;
    };
}  // namespace chdr::experiments::e03
