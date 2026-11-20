/** @file E01E02Protocol.h
 * @brief Experiment-local choices consumed by the shared homogeneous pulse implementation.
 */
#pragma once

#include <span>

namespace chdr::experiments::e01e02 {
    /// Name of the executable that runs E01 to E04, written to every run record.
    inline constexpr const char* WavePulsesApplication =
        "experiments-wave-pulses-in-homogeneous-dielectric";

    /** @brief One physical observation target with its requested products.
     * @details Profiles are always saved. Relative continuum errors are meaningful only
     * while the reference norm is nonzero. Snapshot selection uses the field time label.
     */
    struct PulseTarget {
        const char* name_m;
        double travelWidths_m;
        bool continuumErrors_m;
        bool fieldSnapshot_m;
    };

    /** @brief Keep the experiments explicit without duplicating their numerical path.
     * @details The spans refer to the static arrays of the experiment headers, for example
     * E01.h or E02.h. The cell-count factor follows from the fixed domain length divided by
     * the full pulse width. A zero transverse count permits an explicitly selected transverse
     * mesh, as in E01. Only E02's study bundles the profiles of a completed E01 study for its
     * pre-contact comparison. With selectableSolver_m the input selects the solver
     * configuration (`parameters.solver_configuration`); otherwise the run uses O1VS.
     */
    struct PulseProtocol {
        const char* name_m;
        int longitudinalCellsPerPulse_m;
        bool periodicZ_m;
        int transverseCells_m;
        std::span<const int> resolutions_m;
        std::span<const PulseTarget> targets_m;
        unsigned minimumStudyLevels_m;
        bool bundlesE01Reference_m;
        const char* application_m;
        /// Relative permittivity of the homogeneous medium: 2.13 in E01 and E02, 1 in vacuum.
        double relativePermittivity_m;
        /// Whether the input selects the solver configuration; otherwise the run uses O1VS.
        bool selectableSolver_m;
    };
}  // namespace chdr::experiments::e01e02
