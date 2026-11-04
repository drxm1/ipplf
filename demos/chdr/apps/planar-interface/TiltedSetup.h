/** @file TiltedSetup.h
 * @brief Exactly periodic tilted geometry, carrier selection and directional return estimates.
 */
#pragma once
#include "apps/planar-interface/TiltedDispersion.h"

namespace chdr::apps::planarInterface {
    /** @brief Four directional slownesses sampled over the useful normal spectrum.
     * @details The initial envelope duration is width divided by incoming group speed
     * for every path. This is a ray-envelope estimate, requiring an enlarged-box check
     * [trefethen1982GroupVelocity, sections 1-2], not a rigorous contamination bound.
     */
    struct TiltedTiming {
        std::vector<std::array<double, 4>> inverse_m;
        double width_m;
        explicit TiltedTiming(TiltedRow& r)
            : width_m(r.plane_m.envelopeHalfWidth_m) {
            r.floquetMargin_m = std::numeric_limits<double>::infinity();
            const auto& p     = r.plane_m;
            for (int i = 0; i < p.timingSamples_m; ++i) {
                const double k =
                    p.kx0_m + (2.0 * i / (p.timingSamples_m - 1) - 1) * p.tailSigmas_m * p.sigmaK_m;
                const double omega  = tiltedFrequency(r, k);
                const auto incoming = tiltedRoots(r, omega, 1),
                           other    = tiltedRoots(r, omega, p.epsilon_m);
                inverse_m.push_back({-1 / incoming[0].groupVelocity_m,
                                     1 / incoming[1].groupVelocity_m, -1 / other[0].groupVelocity_m,
                                     1 / other[1].groupVelocity_m});
                r.floquetMargin_m = std::min(r.floquetMargin_m, tiltedOrderMargin(r, omega));
            }
            if (!(r.floquetMargin_m > 0))
                throw std::invalid_argument("Tilted useful band permits another Floquet order");
        }

        /// Earliest/latest arrival, with distances ordered incident-, incident+, other-, other+.
        std::array<double, 2> window(std::array<double, 4> distances) const {
            std::array<double, 2> result{std::numeric_limits<double>::infinity(), 0};
            for (const auto& inverse : inverse_m) {
                const double centre =
                    std::inner_product(distances.begin(), distances.end(), inverse.begin(), 0.0);
                result[0] = std::min(result[0], centre - width_m * inverse[0]);
                result[1] = std::max(result[1], centre + width_m * inverse[0]);
            }
            return result;
        }
    };

    /// Bound actual component positions in both finite Cartesian probe boxes.
    inline std::array<double, 2> tiltedProbeDistances(const TiltedRow& r) {
        const double centre = r.probeCells_m * r.plane_m.h_m[0];
        const double extent =
            (r.probeHalfWidth_m + 1.0)
            * (r.cosine_m * r.plane_m.h_m[0] + std::abs(r.sine_m) * r.plane_m.h_m[2]);
        if (!(centre > extent + r.plane_m.h_m[0]))
            throw std::invalid_argument("Tilted probe boxes must lie wholly inside their medium");
        return {centre - extent, centre + extent};
    }

    /// Round a layer width up in primitive normal rows, keeping the same material period.
    inline double tiltedLayerWidth(const TiltedRow& r, double value) {
        return planeInteger(std::ceil(value / r.rowSpacing_m)) * r.rowSpacing_m;
    }

    /// Exclude the transmitted layer's second face from both observation gates.
    inline double tiltedOtherWidth(const TiltedRow& r, const TiltedTiming& timing,
                                   std::array<double, 2> d) {
        const auto& p     = r.plane_m;
        const int minimum = planeInteger(std::ceil(d[1] / r.rowSpacing_m)) + 2;
        const int rows    = planeLayerCells(minimum, [&](int n) {
            const double D = n * r.rowSpacing_m, S = r.sourceDistance_m;
            return timing.window({S, 0, D, D - d[1]})[0] >= p.gateOther_m * p.dt_m
                   && timing.window({S, d[0], D, D})[0] >= p.gateIncident_m * p.dt_m;
        });
        return tiltedLayerWidth(r, rows * r.rowSpacing_m * p.separationFactor_m);
    }

    /** @brief Exclude complete wraps in each direction and the analytic incident image.
     * @details Distances retain their directional medium speeds. These paths are the
     * normal-plane ray-envelope deduction from Trefethen sections 1-2; reflected roots
     * are not negated incident roots on the anisotropic rotated lattice.
     */
    inline double tiltedIncidentWidth(const TiltedRow& r, const TiltedTiming& timing,
                                      std::array<double, 2> d, double D) {
        const auto& p  = r.plane_m;
        const double S = r.sourceDistance_m, endR = p.gateIncident_m * p.dt_m,
                     endT = p.gateOther_m * p.dt_m;
        const int minimum =
            planeInteger(std::ceil((S + p.envelopeHalfWidth_m) / r.rowSpacing_m)) + 4;
        const int rows = planeLayerCells(minimum, [&](int n) {
            const double H = n * r.rowSpacing_m;
            return timing.window({S + H - d[1], 0, D, 0})[0] >= endR
                   && timing.window({S, H, 0, D - d[1]})[0] >= endT
                   && timing.window({S - d[1] + H + D, 0, 0, 0})[0] >= endR
                   && timing.window({S + H, 0, D + d[0], 0})[0] >= endT
                   && timing.window({S, H + d[0], 0, D})[0] >= endR;
        });
        return tiltedLayerWidth(r, rows * r.rowSpacing_m * p.separationFactor_m);
    }

    /// Resolve directional gates before choosing a compatible Cartesian period.
    inline double tiltedRequiredPeriod(TiltedRow& r) {
        auto& p               = r.plane_m;
        const auto distance   = tiltedProbeDistances(r);
        p.envelopeHalfWidth_m = tiltedLayerWidth(r, p.tailSigmas_m / p.sigmaK_m);
        r.sourceDistance_m    = distance[1] + p.envelopeHalfWidth_m;
        const TiltedTiming timing(r);
        const double S   = r.sourceDistance_m;
        p.gateIncident_m = planeInteger(std::ceil(timing.window({S, distance[1], 0, 0})[1] / p.dt_m)
                                        + 3 + p.gatePaddingSteps_m);
        p.gateOther_m    = planeInteger(std::ceil(timing.window({S, 0, distance[1], 0})[1] / p.dt_m)
                                        + 3 + p.gatePaddingSteps_m);
        r.dielectricWidth_m = tiltedOtherWidth(r, timing, distance);
        r.incidentWidth_m   = tiltedIncidentWidth(r, timing, distance, r.dielectricWidth_m);
        p.steps_m           = std::max(p.gateIncident_m, p.gateOther_m) + 1;
        return r.dielectricWidth_m + r.incidentWidth_m;
    }

    /// Derive the Gaussian width from the actual incoming root's distance to grazing.
    inline void completeTiltedSpectrum(TiltedRow& r) {
        auto& p              = r.plane_m;
        const double grazing = tiltedGrazingRoot(r);
        p.kx0_m              = tiltedRoots(r, p.omega0_m, 1)[0].normalK_m;
        p.sigmaK_m           = (1 - p.cutoffFraction_m) * (grazing - p.kx0_m) / p.tailSigmas_m;
        p.cutoffK_m          = grazing;
        if (!(p.sigmaK_m > 0))
            throw std::invalid_argument("Tilted spectrum has no incoming propagating band");
        p.actualAngleDegrees_m = 180 / Kokkos::numbers::pi * std::asin(r.kappa_m / p.omega0_m);
        p.cellsPerWavelength_m =
            2 * Kokkos::numbers::pi / (std::sqrt(p.epsilon_m) * p.omega0_m * p.h_m[0]);
    }

    /** @brief Exact Cartesian periodicity at fixed tangent mode.
     * @details For coprime \f$p,q\f$, use \f$N_x=Mp\f$, \f$N_z=Mq\f$ and
     * \f$h_z=cp h_x/(|s|q)\f$. Carriers have \f$b=j+\operatorname{sign}(s)a\f$,
     * hence fixed tangent wavenumber \f$2\pi j/L_t\f$.
     * This is a direct reciprocal-lattice deduction from Oskooi/Johnson Eq. (4.30).
     */
    inline void setTiltedPeriods(TiltedRow& r, int multiplier) {
        r.multiplier_m     = multiplier;
        r.plane_m.cellsX_m = planeInteger(double(multiplier) * r.primitiveX_m);
        r.plane_m.cellsZ_m = planeInteger(double(multiplier) * r.primitiveZ_m);
        r.periodRows_m     = planeInteger(double(multiplier) * r.primitiveX_m * r.primitiveZ_m);
        r.normalPeriod_m   = r.periodRows_m * r.rowSpacing_m;
        r.tangentPeriod_m  = multiplier * r.primitivePeriod_m;
    }

    /// Quantize the tangent carrier in the chosen period while retaining its physical angle.
    inline void selectTiltedTangentMode(TiltedRow& r, double requestedOmega, double sine) {
        r.tangentMode_m =
            sine == 0 ? 0
                      : static_cast<int>(std::round(requestedOmega * sine * r.tangentPeriod_m
                                                    / (2 * Kokkos::numbers::pi)));
        if (sine != 0 && r.tangentMode_m == 0)
            throw std::invalid_argument(
                "Tilted box has no tangent mode at the requested incidence");
        r.kappa_m          = 2 * Kokkos::numbers::pi * r.tangentMode_m / r.tangentPeriod_m;
        r.plane_m.omega0_m = sine == 0 ? requestedOmega : r.kappa_m / sine;
    }

    /// Quantize only the tangent mode; adjust actual omega to retain requested incidence.
    inline void resolveTiltedPeriod(TiltedRow& r, double requestedOmega) {
        const double sine = std::sin(r.plane_m.requestedAngleDegrees_m * Kokkos::numbers::pi / 180);
        r.kappa_m         = requestedOmega * sine;
        for (int iteration = 0; iteration < 24; ++iteration) {
            completeTiltedSpectrum(r);
            const double required = tiltedRequiredPeriod(r);
            const int M           = std::max(
                r.multiplier_m, planeInteger(std::ceil(
                                    required / (r.cosine_m * r.primitiveX_m * r.plane_m.h_m[0]))));
            if (M == r.multiplier_m && (sine == 0 || r.tangentMode_m != 0))
                return;
            setTiltedPeriods(r, M);
            selectTiltedTangentMode(r, requestedOmega, sine);
        }
        throw std::invalid_argument("Tilted periodic layout did not converge");
    }

    /// Reject channels without two continuum propagating modes before any field allocation.
    inline void qualifyTiltedFrequency(TiltedRow& r, double omega, double& margin) {
        if (!(omega > std::abs(r.kappa_m)))
            throw std::invalid_argument("Tilted carrier has no propagating continuum reference");
        tiltedRoots(r, omega, 1);
        tiltedRoots(r, omega, r.plane_m.epsilon_m);
        margin = std::min(margin, tiltedOrderMargin(r, omega));
        if (!(margin > 0))
            throw std::invalid_argument("Tilted retained spectrum permits another Floquet order");
    }

    /// Channel-weighted modal L1 tail accounting bounds the omitted homogeneous field amplitude.
    inline double tiltedModeFieldWeight(const TiltedRow& r, double omega, double weight) {
        const double dt = r.plane_m.dt_m, W = 2 * std::sin(omega * dt / 2) / dt;
        return weight
               * (r.plane_m.polarization_m == Polarization::TE ? W
                                                               : W * W * std::cos(omega * dt / 2));
    }

    /// Form one Cartesian-periodic candidate, excluding only nonprincipal z phases.
    inline std::optional<std::array<double, 5>> tiltedCarrier(const TiltedRow& r, int a) {
        const auto& p = r.plane_m;
        const int b   = r.tangentMode_m + (r.sine_m > 0 ? a : -a);
        if (std::abs(b) * 2 >= p.cellsZ_m)
            return {};
        const double x         = 2 * Kokkos::numbers::pi * a / (p.cellsX_m * p.h_m[0]);
        const double z         = 2 * Kokkos::numbers::pi * b / (p.cellsZ_m * p.h_m[2]);
        const double normal    = r.cosine_m * x + r.sine_m * z;
        const double deviation = (normal - p.kx0_m) / p.sigmaK_m;
        return std::array{x, z, tiltedFrequency(r, normal), std::exp(-0.5 * deviation * deviation),
                          normal};
    }

    /// Account for a principal candidate before retaining only its qualified incoming part.
    inline void accumulateTiltedMode(TiltedRow& r, const std::array<double, 5>& mode,
                                     double& retained, double& omitted) {
        const auto& p          = r.plane_m;
        const double deviation = (mode[4] - p.kx0_m) / p.sigmaK_m;
        const double field     = tiltedModeFieldWeight(r, mode[2], mode[3]);
        if (std::abs(deviation) >= p.tailSigmas_m || tiltedSpatialDerivative(r, mode[4]) >= 0) {
            omitted += field;
            return;
        }
        qualifyTiltedFrequency(r, mode[2], r.retainedMargin_m);
        r.modes_m.push_back(mode);
        retained += field;
    }

    /// Keep incoming exact Cartesian modes; account for all omitted principal-mode field weights.
    inline void selectTiltedModes(TiltedRow& r) {
        const auto& p   = r.plane_m;
        double retained = 0, omitted = 0;
        r.retainedMargin_m = std::numeric_limits<double>::infinity();
        for (int a = -(p.cellsX_m - 1) / 2; a <= (p.cellsX_m - 1) / 2; ++a) {
            const auto candidate = tiltedCarrier(r, a);
            if (!candidate)
                continue;
            accumulateTiltedMode(r, *candidate, retained, omitted);
        }
        if (r.modes_m.size() < 16 || !(retained > 0))
            throw std::invalid_argument("Tilted packet has insufficient retained periodic modes");
        r.discardedFieldWeightFraction_m = omitted / retained;
    }

    /// Use actual tilted dispersion for an arbitrary-density normal-wavenumber measurement band.
    inline void selectTiltedBins(TiltedRow& r, const PlaneRequest& input) {
        auto& p = r.plane_m;
        if (input.frequencyBins_m < 3 || input.frequencyBins_m % 2 == 0
            || !std::isfinite(input.bandSigmas_m) || !std::isfinite(input.incidentFloor_m)
            || input.bandSigmas_m <= 0 || input.bandSigmas_m > p.tailSigmas_m
            || input.incidentFloor_m <= 0 || input.incidentFloor_m >= 1)
            throw std::invalid_argument("Invalid tilted frequency band or incident floor");
        p.bandSigmas_m    = input.bandSigmas_m;
        p.incidentFloor_m = input.incidentFloor_m;
        p.omega_m         = {p.omega0_m};
        const int half    = (input.frequencyBins_m - 1) / 2;
        for (int i = -half; i <= half; ++i)
            if (i != 0)
                p.omega_m.push_back(
                    tiltedFrequency(r, p.kx0_m - i * p.bandSigmas_m * p.sigmaK_m / half));
        for (double omega : p.omega_m)
            qualifyTiltedFrequency(r, omega, r.floquetMargin_m);
    }
}  // namespace chdr::apps::planarInterface
