/** @file PlanePresets.h
 * @brief Explicit regression layouts using the common plane-row representation.
 */
#pragma once
#include "apps/planar-interface/PlaneSetup.h"

namespace chdr::apps::planarInterface {
    /** @brief Complete the recorded presets without changing their layer or gate integers.
     * @details Their original pulse centre is probeCells + q/sigma beyond the face,
     * while TE probes are an additional half cell away. Keep those distinct positions
     * when reproducing the recorded coefficients. Generic requests use the conservative
     * layout calculation in resolvePlaneRow instead.
     */
    inline PlaneRow completePlanePreset(PlaneRow p) {
        completePlaneSpectrum(p);
        const double side     = p.incidentFromVacuum_m ? 1 : -1;
        const double face     = p.dielectricCells_m * p.h_m[0];
        const double offset   = p.polarization_m == Polarization::TE ? 0.5 : 1;
        p.envelopeHalfWidth_m = p.tailSigmas_m / p.sigmaK_m;
        p.pulseCentre_m       = face + side * (3 * p.h_m[0] + p.envelopeHalfWidth_m);
        p.incidentProbe_m     = face + side * (4 - offset) * p.h_m[0];
        p.otherProbe_m        = face - side * (4 - offset) * p.h_m[0];
        p.incidentProbeIndex_m =
            static_cast<int>(std::lround(p.incidentProbe_m / p.h_m[0] - offset));
        p.otherProbeIndex_m = static_cast<int>(std::lround(p.otherProbe_m / p.h_m[0] - offset));
        p.requestedAngleDegrees_m       = p.actualAngleDegrees_m;
        p.requestedCellsPerWavelength_m = p.cellsPerWavelength_m;
        return p;
    }

    /// Recorded normal TE pulse: 59 modes, cubic cells, 20 cells per dielectric wavelength.
    inline PlaneRow normalPlaneRow() {
        PlaneRow p{.cellsX_m          = 812,
                   .dielectricCells_m = 211,
                   .dt_m              = 0.9 / std::sqrt(3.0),
                   .omega0_m          = 2 * Kokkos::numbers::pi / (20 * std::sqrt(2.13)),
                   .gateIncident_m    = 1175,
                   .gateOther_m       = 1178,
                   .steps_m           = 1180};
        return completePlanePreset(p);
    }

    /// Recorded TM Brewster pulse: 88 modes; the tangent period is now along z.
    inline PlaneRow brewsterPlaneRow() {
        PlaneRow p{.polarization_m    = Polarization::TM,
                   .cellsX_m          = 3581,
                   .cellsZ_m          = 18,
                   .dielectricCells_m = 1484,
                   .h_m               = {1, 2, 2},
                   .kz_m              = 2 * Kokkos::numbers::pi / (18 * 2),
                   .dt_m              = 0.9 / std::sqrt(1.5),
                   .tailSigmas_m      = 7,
                   .cutoffFraction_m  = 0.5,
                   .gateIncident_m    = 6759,
                   .gateOther_m       = 6755,
                   .steps_m           = 6761};
        p.omega0_m = p.kz_m * std::sqrt(1 + 2.13) / std::sqrt(2.13);
        return completePlanePreset(p);
    }

    /// Recorded dielectric-side TM pulse: 78 modes, with both transverse axes permuted.
    inline PlaneRow dielectricPlaneRow() {
        PlaneRow p{.polarization_m       = Polarization::TM,
                   .incidentFromVacuum_m = false,
                   .cellsX_m             = 1828,
                   .cellsZ_m             = 52,
                   .dielectricCells_m    = 949,
                   .h_m                  = {1, 0.85, 1.25},
                   .epsilon_m            = 3.4,
                   .kz_m                 = 2 * Kokkos::numbers::pi / (52 * 1.25),
                   .dt_m             = 0.9 / std::sqrt(1 + 1 / (1.25 * 1.25) + 1 / (0.85 * 0.85)),
                   .tailSigmas_m     = 7,
                   .cutoffFraction_m = 0.3,
                   .gateIncident_m   = 3787,
                   .gateOther_m      = 3785,
                   .steps_m          = 3789};
        p.omega0_m = p.kz_m / (std::sqrt(3.4) * std::sin(Kokkos::numbers::pi / 10));
        return completePlanePreset(p);
    }
}  // namespace chdr::apps::planarInterface
