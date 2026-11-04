/** @file PlanePresets.h
 * @brief The three recorded plane-pulse layouts and their shared measurement conventions.
 * @details normalPlaneRow(), brewsterPlaneRow() and dielectricPlaneRow() return explicit
 * inputs to the common PlaneRow/PlanePulse path. With `case: plane`, parameters.preset
 * selects `normal`, `brewster` or `dielectric`. They are distinct
 * from the qualitative packet cases described in the application README.
 * All preset quantities use normalized solver units, with
 * \f$c=\varepsilon_0=\mu_0=1\f$ and relative permeability one.
 *
 * All three layouts use the x-z incidence plane, with y out of plane. The dielectric
 * fills cells with x indices from zero through dielectricCells_m - 1. Its measured
 * material interface is at \f$x_I=D h_x\f$, where \f$D\f$ is dielectricCells_m; another
 * interface lies at the periodic x join. Here the material interface coincides with a
 * cell face. A requested subcell interface in resolvePlaneRow need not do so.
 *
 * The initial potential is a Gaussian-weighted sum of periodic incident-medium modes.
 * Their positive normal magnitudes are \f$k_m=2\pi m/(N_xh_x)\f$, with weights
 * \f$\exp[-(k_m-k_{x0})^2/(2\sigma_k^2)]\f$ retained for
 * \f$|k_m-k_{x0}|<9\sigma_k\f$. Each mode uses the Standard lattice frequency
 * [fallahi2020mithra20fullwavesimulation, section 3.1.3, Eq. (3.17)]. The packet-width
 * choice is \f$\sigma_k=(1-f)(k_{x0}-k_{\mathrm{cut}})/q\f$, with \f$q\f$ = tailSigmas_m
 * and \f$f\f$ = cutoffFraction_m. The cutoff is zero for vacuum incidence; for dielectric
 * incidence it is the incident-medium normal root at the vacuum propagation cutoff.
 * These are the finite Gaussian packet conventions of completePlaneSpectrum; see
 * [trefethen1982GroupVelocity, section 2] for the group-velocity packet construction.
 *
 * The nominal envelope half-width is \f$w=q/\sigma_k\f$. The pulse centre is
 * \f$3h_x+w\f$ from the measured interface on the incident side. TE probes measure
 * \f$E_y\f$ at x-z cell-centre projections, at distance \f$3.5h_x\f$ on each side.
 * TM probes measure the two-level average of \f$H_y=B_y\f$ at x-z cell-corner
 * projections, at distance \f$3h_x\f$. These component locations follow the staggered
 * field construction [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(27)]. The Fourier
 * probes label both channels at \f$(k-1/2)\Delta t\f$; raw magnetic field files instead
 * use the integer potential time. A gate of length G includes completed states
 * k = 0 through G - 1, including the initialized sample. steps_m counts solves.
 *
 * The recorded layer lengths and gates are explicit integers, not recomputed arrival
 * estimates. Their original selection considered second-interface reflections,
 * periodic wraps and the incident packet's periodic image with lattice group speeds
 * [trefethen1982GroupVelocity, sections 1-2]. Such finite-band timing estimates do not
 * make the Gaussian compactly supported; changed configurations need control runs.
 * The README table collects the values. Full reference and extraction derivations are
 * in demos/chdr/experiments/interface-reference-cases/references.org.
 */
#pragma once
#include "apps/planar-interface/PlaneSetup.h"

namespace chdr::apps::planarInterface {
    /** @brief Complete the recorded presets without changing their layer or gate integers.
     * @details Preserve the pulse distance \f$3h_x+q/\sigma_k\f$ and the TE probe distance
     * \f$3.5h_x\f$ from the material interface as distinct quantities. The file-level
     * documentation defines the pulse, probes and sample times. Generic requests use
     * resolvePlaneRow to calculate their layout from sampled travel-time estimates.
     */
    inline PlaneRow completePlanePreset(PlaneRow p) {
        completePlaneSpectrum(p);
        const double side     = p.incidentFromVacuum_m ? 1 : -1;
        const double face     = p.dielectricCells_m * p.h_m[0];
        p.facePosition_m      = face;
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
        selectPlaneBins(p, 3, 3, 1e-8);
        return p;
    }

    /** @brief Normal-incidence TE preset from vacuum into relative permittivity 2.13.
     * @details Cubic unit cells give 20 cells per dielectric wavelength at
     * \f$\omega_0=2\pi/(20\sqrt{2.13})\f$. The Courant fraction is 0.9, q = 8 and f = 0.
     * The 812 by 2 by 2 box has D = 211 dielectric cells and 59 retained modes.
     * Incident/transmitted gates contain 1175/1178 samples; evolution has 1180 solves.
     * The frequency follows the nonmagnetic plane-wave relation
     * [ellingson2020ElectromagneticsVol2, sections 5.6-5.8].
     */
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

    /** @brief TM preset from vacuum at the central continuum Brewster angle.
     * @details Relative permittivity is 2.13; \f$\theta_B=\arctan\sqrt{2.13}\f$
     * [ellingson2020ElectromagneticsVol2, section 5.10, Eq. (5.253)]. The fixed tangent
     * wavenumber is \f$k_z=2\pi/(18\cdot2)\f$; only the central frequency is at Brewster
     * incidence [oskooi2013sources, section 4.5]. The 3581 by 2 by 18 box uses cell
     * spacings (1, 2, 2), D = 1484, Courant fraction 0.9, q = 7, f = 0.5 and 88 modes.
     * Incident/transmitted gates contain 6759/6755 samples; evolution has 6761 solves.
     */
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

    /** @brief TM preset from the dielectric at 18 degrees to the interface normal.
     * @details Relative permittivity is 3.4, \f$k_z=2\pi/(52\cdot1.25)\f$ and
     * \f$\omega_0=k_z/(\sqrt{3.4}\sin18^\circ)\f$, from tangential phase matching
     * [ellingson2020ElectromagneticsVol2, sections 5.7-5.8]. The 1828 by 2 by 52 box
     * uses spacings (1, 0.85, 1.25), D = 949, Courant fraction 0.9, q = 7, f = 0.3
     * and 78 modes. Incident/transmitted gates contain 3787/3785 samples; evolution
     * has 3789 solves. The y-z permutation preserves the original incidence-plane
     * spacing and coefficients while using the application's x-z convention.
     */
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
