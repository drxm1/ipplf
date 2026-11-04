/** @file PlaneSetup.h
 * @brief Resolve a periodic plane packet using numerical dispersion and travel times.
 */
#pragma once
#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "apps/planar-interface/PlaneRow.h"

namespace chdr::apps::planarInterface {
    /// Spatial difference symbol [fallahi2020mithra20fullwavesimulation, Eq. (3.17)].
    inline double planeSpatialSymbol(double k, double h) {
        return 2 * std::sin(k * h / 2) / h;
    }

    /// Incident and transmitted permittivities in the selected direction of travel.
    inline std::array<double, 2> planeMedia(const PlaneRow& p) {
        return p.incidentFromVacuum_m ? std::array{1.0, p.epsilon_m} : std::array{p.epsilon_m, 1.0};
    }

    /// Frequency of an incident lattice mode [fallahi2020mithra20fullwavesimulation, Eq. (3.17)].
    inline double planeFrequency(const PlaneRow& p, double normalK) {
        const double x    = planeSpatialSymbol(normalK, p.h_m[0]);
        const double z    = planeSpatialSymbol(p.kz_m, p.h_m[2]);
        const double sine = p.dt_m * std::sqrt((x * x + z * z) / planeMedia(p)[0]) / 2;
        if (!(sine >= 0 && sine < 1))
            throw std::invalid_argument("Plane frequency is outside the propagating time branch");
        return 2 / p.dt_m * std::asin(sine);
    }

    /** @brief Positive normal root, allowing the zero-frequency normal-incidence limit.
     * Deduction: solve the Standard dispersion for its x symbol
     * [fallahi2020mithra20fullwavesimulation, section 3.1.3, Eq. (3.17)].
     */
    inline double planeNormalK(const PlaneRow& p, double omega, double epsilon) {
        if (!std::isfinite(omega) || omega < 0 || omega * p.dt_m >= Kokkos::numbers::pi)
            throw std::invalid_argument(
                "Plane frequency must lie in the principal leapfrog branch");
        const double W  = 2 * std::sin(omega * p.dt_m / 2) / p.dt_m;
        const double Kz = planeSpatialSymbol(p.kz_m, p.h_m[2]);
        double sine2    = p.h_m[0] * p.h_m[0] * (epsilon * W * W - Kz * Kz) / 4;
        if (sine2 < 0 && sine2 > -32 * std::numeric_limits<double>::epsilon())
            sine2 = 0;
        if (!(sine2 >= 0 && sine2 < 1))
            throw std::invalid_argument(
                "Plane row requires propagating normal roots in both media");
        return 2 / p.h_m[0] * std::asin(std::sqrt(sine2));
    }

    /** @brief Normal group speed at fixed tangential wavenumber.
     * Deduction: differentiate \f$\varepsilon W^2=K_x^2+K_z^2\f$
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.17);
     * trefethen1982GroupVelocity, section 1, Eqs. (1.4)-(1.5)].
     */
    inline double planeGroupSpeed(const PlaneRow& p, double omega, double epsilon) {
        if (omega == 0 && p.kz_m == 0)
            return 1 / std::sqrt(epsilon);
        const double k = planeNormalK(p, omega, epsilon);
        const double W = 2 * std::sin(omega * p.dt_m / 2) / p.dt_m;
        const double speed =
            std::sin(k * p.h_m[0]) / (p.h_m[0] * epsilon * W * std::cos(omega * p.dt_m / 2));
        if (!(speed > 0) || !std::isfinite(speed))
            throw std::invalid_argument("Plane pulse band reaches a zero group-speed cutoff");
        return speed;
    }

    /// Reject non-finite or out-of-range requests before any distributed field allocation.
    inline void validatePlaneRequest(const PlaneRequest& r) {
        if (!std::isfinite(r.epsilon_m) || r.epsilon_m < 1 || !std::isfinite(r.courant_m)
            || r.courant_m <= 0 || r.courant_m >= 1 || !std::isfinite(r.cellsPerWavelength_m)
            || r.cellsPerWavelength_m <= 2 || !std::isfinite(r.tailSigmas_m) || r.tailSigmas_m <= 3
            || !std::isfinite(r.separationFactor_m) || r.separationFactor_m < 1)
            throw std::invalid_argument(
                "Invalid plane material, resolution, Courant or pulse settings");
        for (double h : r.h_m)
            if (!std::isfinite(h) || h <= 0)
                throw std::invalid_argument("Plane cell spacings must be finite and positive");
        if (r.angleDegrees_m && r.transversePeriod_m)
            throw std::invalid_argument(
                "angle_degrees and transverse_period are mutually exclusive");
        if (r.probeCells_m < 2 || r.timingSamples_m < 17 || r.gatePaddingSteps_m < 0)
            throw std::invalid_argument("Plane probes, timing samples or gate padding are invalid");
    }

    /// Check conversion before narrowing a computed cell or step count.
    inline int planeInteger(double value) {
        if (!std::isfinite(value) || value < 0 || value > std::numeric_limits<int>::max() / 2)
            throw std::invalid_argument(
                "Resolved plane cell or step count exceeds the supported range");
        return static_cast<int>(value);
    }

    /// Vacuum Standard CFL bound [fallahi2020mithra20fullwavesimulation, Eq. (3.20)].
    inline double planeTimeStep(const PlaneRequest& r) {
        double inverseSquareSum = 0;
        for (double h : r.h_m)
            inverseSquareSum += 1 / (h * h);
        return r.courant_m / std::sqrt(inverseSquareSum);
    }

    /// Copy independent choices and compute the central frequency from dielectric resolution.
    inline PlaneRow planeBaseRow(const PlaneRequest& r) {
        PlaneRow p;
        p.polarization_m                = r.polarization_m;
        p.incidentFromVacuum_m          = r.incidentFromVacuum_m;
        p.h_m                           = r.h_m;
        p.epsilon_m                     = r.epsilon_m;
        p.courant_m                     = r.courant_m;
        p.dt_m                          = planeTimeStep(r);
        p.tailSigmas_m                  = r.tailSigmas_m;
        p.timingSamples_m               = r.timingSamples_m;
        p.gatePaddingSteps_m            = r.gatePaddingSteps_m;
        p.separationFactor_m            = r.separationFactor_m;
        p.requestedCellsPerWavelength_m = r.cellsPerWavelength_m;
        p.omega0_m =
            2 * Kokkos::numbers::pi / (std::sqrt(p.epsilon_m) * r.cellsPerWavelength_m * p.h_m[0]);
        return p;
    }

    /** @brief Fix cell shape and choose one periodic tangent mode at the requested angle.
     * Deduction from the unit Bloch phase and continuum incidence relation
     * [oskooi2013sources, section 4.5, Eq. (4.30)]. Actual wavelength is recorded.
     */
    inline void resolvePlaneAngle(PlaneRow& p, double angleDegrees) {
        if (!std::isfinite(angleDegrees) || std::abs(angleDegrees) >= 90)
            throw std::invalid_argument(
                "Plane incidence angle must lie strictly between -90 and 90 degrees");
        p.requestedAngleDegrees_m = angleDegrees;
        if (angleDegrees == 0)
            return;
        const double sine    = std::sin(angleDegrees * Kokkos::numbers::pi / 180);
        const double targetK = std::sqrt(planeMedia(p)[0]) * p.omega0_m * sine;
        p.cellsZ_m           = std::max(
            2, planeInteger(std::round(2 * Kokkos::numbers::pi / (std::abs(targetK) * p.h_m[2]))));
        p.kz_m     = std::copysign(2 * Kokkos::numbers::pi / (p.cellsZ_m * p.h_m[2]), sine);
        p.omega0_m = p.kz_m / (std::sqrt(planeMedia(p)[0]) * sine);
    }

    /// A requested period must already fit an integer number of cells; frequency stays fixed.
    inline void resolvePlanePeriod(PlaneRow& p, double period) {
        const double count = period / p.h_m[2], rounded = std::round(count);
        if (!std::isfinite(period) || period <= 0 || rounded < 2
            || std::abs(count - rounded)
                   > 32 * std::numeric_limits<double>::epsilon() * std::max(1.0, count))
            throw std::invalid_argument(
                "transverse_period must be an integer number of at least two z cells");
        p.cellsZ_m = planeInteger(rounded);
        p.kz_m     = 2 * Kokkos::numbers::pi / period;
        p.requestedAngleDegrees_m =
            180 / Kokkos::numbers::pi
            * std::asin(p.kz_m / (std::sqrt(planeMedia(p)[0]) * p.omega0_m));
    }

    /** @brief Complete the propagating central mode and Gaussian normal spectrum.
     * The useful band stays above the common propagation cutoff; tailSigmas specifies
     * both its spectral extent and initial envelope extent. Dispersion and angle
     * relations: [fallahi2020mithra20fullwavesimulation, Eq. (3.17)] and
     * [oskooi2013sources, section 4.5]; Gaussian packet
     * construction: Trefethen section 2 [trefethen1982GroupVelocity].
     */
    inline void completePlaneSpectrum(PlaneRow& p) {
        const auto media  = planeMedia(p);
        const double sine = p.kz_m / (std::sqrt(media[0]) * p.omega0_m);
        if (!std::isfinite(sine) || std::abs(sine) >= 1
            || std::abs(p.kz_m) >= std::sqrt(media[1]) * p.omega0_m)
            throw std::invalid_argument("Plane centre must propagate in both continuum media");
        p.actualAngleDegrees_m = 180 / Kokkos::numbers::pi * std::asin(sine);
        p.cellsPerWavelength_m =
            2 * Kokkos::numbers::pi / (std::sqrt(p.epsilon_m) * p.omega0_m * p.h_m[0]);
        p.kx0_m = planeNormalK(p, p.omega0_m, media[0]);
        const double cutFrequency =
            2 / p.dt_m * std::asin(p.dt_m * std::abs(planeSpatialSymbol(p.kz_m, p.h_m[2])) / 2);
        p.cutoffK_m = p.incidentFromVacuum_m ? 0 : planeNormalK(p, cutFrequency, media[0]);
        p.sigmaK_m  = (1 - p.cutoffFraction_m) * (p.kx0_m - p.cutoffK_m) / p.tailSigmas_m;
        if (!(p.sigmaK_m > 0) || p.kx0_m + 9 * p.sigmaK_m >= Kokkos::numbers::pi / p.h_m[0])
            throw std::invalid_argument(
                "Plane spectrum is empty or reaches the normal Nyquist edge");
    }

    /** @brief Arrival estimates over a finite band using the same incident envelope duration.
     * @details Each path uses its incident- and transmitted-medium distance at each sampled
     * frequency. The initial half-duration is always width/v_inc, even after transmission.
     * This is a ray-envelope estimate, not an error bound; test denser sampling and larger
     * boxes before claiming isolation [trefethen1982GroupVelocity, sections 1-2].
     */
    struct PlaneTiming {
        std::vector<std::array<double, 2>> slowness_m;
        double width_m;

        explicit PlaneTiming(const PlaneRow& p)
            : width_m(p.envelopeHalfWidth_m) {
            const auto media = planeMedia(p);
            for (int i = 0; i < p.timingSamples_m; ++i) {
                const double fraction = 2.0 * i / (p.timingSamples_m - 1) - 1;
                const double k = std::max(0.0, p.kx0_m + fraction * p.tailSigmas_m * p.sigmaK_m);
                const double omega = planeFrequency(p, k);
                slowness_m.push_back({1 / planeGroupSpeed(p, omega, media[0]),
                                      1 / planeGroupSpeed(p, omega, media[1])});
            }
        }

        /// Earliest and latest envelope arrival for distances measured in code length units.
        std::array<double, 2> window(double incidentDistance, double otherDistance) const {
            std::array<double, 2> range{std::numeric_limits<double>::infinity(), 0};
            for (const auto& inverse : slowness_m) {
                const double centre = incidentDistance * inverse[0] + otherDistance * inverse[1];
                range[0]            = std::min(range[0], centre - width_m * inverse[0]);
                range[1]            = std::max(range[1], centre + width_m * inverse[0]);
            }
            return range;
        }
    };

    /// Find a safe integer layer length for a monotone arrival-time predicate.
    template <class Predicate>
    int planeLayerCells(int minimum, Predicate fits) {
        int low = minimum, high = minimum;
        while (!fits(high)) {
            if (high > std::numeric_limits<int>::max() / 4)
                throw std::invalid_argument("Plane layout cannot fit within supported cell counts");
            high *= 2;
        }
        while (high - low > 1) {
            const int middle            = low + (high - low) / 2;
            (fits(middle) ? high : low) = middle;
        }
        return fits(low) ? low : high;
    }

    /// Place physical component samples and the pulse from resolved layer lengths.
    inline void placePlanePacket(PlaneRow& p, int incidentCells, int otherCells, double distance) {
        p.cellsX_m          = incidentCells + otherCells;
        p.dielectricCells_m = p.incidentFromVacuum_m ? otherCells : incidentCells;
        const double side   = p.incidentFromVacuum_m ? 1 : -1;
        const double face   = p.dielectricCells_m * p.h_m[0];
        p.incidentProbe_m   = face + side * distance;
        p.otherProbe_m      = face - side * distance;
        p.pulseCentre_m     = p.incidentProbe_m + side * p.envelopeHalfWidth_m;
        const double offset = p.polarization_m == Polarization::TE ? 0.5 : 1;
        p.incidentProbeIndex_m =
            static_cast<int>(std::lround(p.incidentProbe_m / p.h_m[0] - offset));
        p.otherProbeIndex_m = static_cast<int>(std::lround(p.otherProbe_m / p.h_m[0] - offset));
    }

    /// Exclude a second-face return from both probes using numerical arrival estimates.
    inline int planeOtherLayer(const PlaneTiming& timing, const PlaneRow& p, double source,
                               double distance) {
        const double h    = p.h_m[0];
        const int minimum = planeInteger(std::ceil(distance / h)) + 2;
        const int cells   = planeLayerCells(minimum, [&](int n) {
            return timing.window(source, 2 * n * h - distance)[0] >= p.gateOther_m * p.dt_m
                   && timing.window(source + distance, 2 * n * h)[0] >= p.gateIncident_m * p.dt_m;
        });
        return planeInteger(std::ceil(cells * p.separationFactor_m));
    }

    /** @brief Exclude both complete periodic wraps and the analytic incident's periodic image.
     * Distances are added within each medium before taking the earliest arrival over the band
     * [trefethen1982GroupVelocity, sections 1-2; deduction applied to the listed ray paths].
     */
    inline int planeIncidentLayer(const PlaneTiming& timing, const PlaneRow& p, double source,
                                  double distance, int otherCells) {
        const double h = p.h_m[0], other = otherCells * h;
        const double endR = p.gateIncident_m * p.dt_m, endT = p.gateOther_m * p.dt_m;
        const int minimum = planeInteger(std::ceil((source + p.envelopeHalfWidth_m) / h)) + 4;
        const int cells   = planeLayerCells(minimum, [&](int n) {
            return timing.window(source + n * h - distance, other)[0] >= endR
                   && timing.window(source + n * h, other - distance)[0] >= endT
                   && timing.window(source - distance + n * h + other, 0)[0] >= endR
                   && timing.window(source + n * h, other + distance)[0] >= endT
                   && timing.window(source + n * h + distance, other)[0] >= endR;
        });
        return planeInteger(std::ceil(cells * p.separationFactor_m));
    }

    /** @brief Choose layers and gates using the completed spectrum and physical probe offsets.
     * The extra samples include the last half-step-labelled sample after the latest intended
     * arrival. Layer searches then exclude all listed return paths until after these gates.
     */
    inline void resolvePlaneLayout(PlaneRow& p, int probeCells) {
        const double h        = p.h_m[0];
        const double distance = (probeCells + (p.polarization_m == Polarization::TE ? 0.5 : 0)) * h;
        p.envelopeHalfWidth_m = std::ceil(p.tailSigmas_m / (p.sigmaK_m * h)) * h;
        const PlaneTiming timing(p);
        const double source = distance + p.envelopeHalfWidth_m;
        p.gateIncident_m = planeInteger(std::ceil(timing.window(source + distance, 0)[1] / p.dt_m)
                                        + 3.0 + p.gatePaddingSteps_m);
        p.gateOther_m    = planeInteger(std::ceil(timing.window(source, distance)[1] / p.dt_m) + 3.0
                                        + p.gatePaddingSteps_m);
        const int other  = planeOtherLayer(timing, p, source, distance);
        const int incident = planeIncidentLayer(timing, p, source, distance, other);
        placePlanePacket(p, incident, other, distance);
        p.steps_m = std::max(p.gateIncident_m, p.gateOther_m) + 1;
    }

    /// Resolve one complete x-normal row; changing a request always recomputes its dependent data.
    inline PlaneRow resolvePlaneRow(const PlaneRequest& r) {
        validatePlaneRequest(r);
        auto p = planeBaseRow(r);
        if (r.transversePeriod_m)
            resolvePlanePeriod(p, *r.transversePeriod_m);
        else
            resolvePlaneAngle(p, r.angleDegrees_m.value_or(0));
        p.cutoffFraction_m = r.cutoffFraction_m.value_or(p.kz_m == 0 ? 0 : 0.5);
        if (!std::isfinite(p.cutoffFraction_m) || p.cutoffFraction_m < 0 || p.cutoffFraction_m >= 1
            || (p.kz_m != 0 && p.cutoffFraction_m == 0))
            throw std::invalid_argument(
                "cutoff_fraction must be in [0,1), and positive for oblique incidence");
        completePlaneSpectrum(p);
        resolvePlaneLayout(p, r.probeCells_m);
        return p;
    }
}  // namespace chdr::apps::planarInterface
