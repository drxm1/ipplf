/** @file PlaneRow.h
 * @brief Requested and resolved data for one aligned planar-interface measurement.
 */
#pragma once
#include <array>
#include <optional>
#include <vector>

namespace chdr::apps::planarInterface {
    /// Measured fixed Cartesian component: electric y for TE, magnetic y for TM.
    enum class Polarization {
        TE,
        TM
    };

    /** @brief Physical and numerical choices before choosing periodic cells and probe gates.
     * @details Spacings are in code length units. Angle is the continuum incidence angle
     * relative to the x-normal material interface. At fixed cell shape, integer transverse cells
     * adjust frequency slightly from the requested dielectric cells per wavelength. Alternatively,
     * a transverse period fixes that frequency and determines the angle. The two input
     * forms are exclusive [oskooi2013sources, section 4.5, Eq. (4.30)].
     * An optional physical interface position fixes the dielectric layer; the requested
     * coordinate is retained even when cell-centre sampling represents another boundary.
     */
    struct PlaneRequest {
        Polarization polarization_m = Polarization::TE;
        bool incidentFromVacuum_m   = true;
        double epsilon_m = 2.13, courant_m = 0.9, cellsPerWavelength_m = 20;
        std::array<double, 3> h_m{1, 1, 1};
        std::optional<double> angleDegrees_m, transversePeriod_m, cutoffFraction_m, facePosition_m;
        double tailSigmas_m = 8, separationFactor_m = 1;
        int probeCells_m = 3, timingSamples_m = 513, gatePaddingSteps_m = 0;
        int frequencyBins_m = 3;
        double bandSigmas_m = 3, incidentFloor_m = 1e-8;
    };

    /** @brief Complete input shared by initialization, observation and run provenance.
     * @details The dielectric occupies cells below dielectricCells_m, with a second
     * interface at the periodic x join. Propagation lies in x-z, with two periodic y
     * cells. Probe coordinates include the actual electric/magnetic component offset.
     * facePosition_m is the requested upper material interface; dielectricCells_m times
     * h_m[0] is the cell face representing the sampled material jump. The requested
     * interface is the coefficient phase origin. The identifier facePosition_m retains
     * the input/record name; it does not denote an arbitrary mesh cell face.
     * Gates include completed states zero through gate-1; their field label is
     * \f$(k-1/2)\Delta t\f$ [ryu2016potentialfdtd, Appendix A]. No derived member
     * initializers are used: setup fills all dependent quantities together. PlanePresets.h
     * documents the recorded configurations, their pulse assumptions and gate lengths.
     */
    struct PlaneRow {
        Polarization polarization_m = Polarization::TE;
        bool incidentFromVacuum_m   = true;
        int cellsX_m = 0, cellsY_m = 2, cellsZ_m = 2, dielectricCells_m = 0;
        std::array<double, 3> h_m{1, 1, 1};
        double epsilon_m = 2.13, courant_m = 0.9, kz_m = 0, dt_m = 0;
        double omega0_m = 0, kx0_m = 0, sigmaK_m = 0, cutoffK_m = 0;
        double tailSigmas_m = 8, cutoffFraction_m = 0, envelopeHalfWidth_m = 0;
        double pulseCentre_m = 0, incidentProbe_m = 0, otherProbe_m = 0;
        double facePosition_m    = 0;
        int incidentProbeIndex_m = 0, otherProbeIndex_m = 0;
        int gateIncident_m = 0, gateOther_m = 0, steps_m = 0;
        double requestedAngleDegrees_m = 0, actualAngleDegrees_m = 0;
        double requestedCellsPerWavelength_m = 20, cellsPerWavelength_m = 0;
        int timingSamples_m = 0, gatePaddingSteps_m = 0;
        double separationFactor_m = 1;
        double bandSigmas_m = 3, incidentFloor_m = 1e-8;
        std::vector<double> omega_m{};
    };
}  // namespace chdr::apps::planarInterface
