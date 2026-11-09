/** @file E01E02Protocol.h
 * @brief Experiment-local choices consumed by the shared homogeneous pulse implementation.
 */
#pragma once

#include <span>

namespace chdr::experiments::e01e02 {
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

    /** @brief Keep the two experiments explicit without duplicating their numerical path.
     * @details The spans refer to the static arrays in E01.h or E02.h. The cell-count
     * factor follows from the fixed domain length divided by the full pulse width.
     * A zero transverse count permits E01's explicitly selected transverse mesh.
     */
    struct PulseProtocol {
        const char* name_m;
        int longitudinalCellsPerPulse_m;
        bool periodicZ_m;
        int transverseCells_m;
        std::span<const int> resolutions_m;
        std::span<const PulseTarget> targets_m;
        unsigned minimumStudyLevels_m;
    };
}  // namespace chdr::experiments::e01e02
