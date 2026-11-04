/** @file TiltedMeasurement.h
 * @brief Raw box spectra and known-wavevector projection at a tilted requested plane.
 */
#pragma once
#include "apps/planar-interface/TiltedPulse.h"

namespace chdr::apps::planarInterface {
    /// Three fitted face amplitudes and their independent spatial residuals for one frequency.
    struct TiltedBinResult {
        TiltedRoot incidentRoot_m{}, reflectedRoot_m{}, transmittedRoot_m{};
        std::complex<double> incident_m{}, reflected_m{}, transmitted_m{}, r_m{}, t_m{};
        double incidentResidual_m = 0, reflectedResidual_m = 0, transmittedResidual_m = 0;
        std::string reason_m = "not_evaluated";
        bool valid() const { return reason_m.empty(); }
    };

    /// Preserve point data alongside the derived fit; energy uses the existing invariant observer.
    struct TiltedResult {
        diagnostics::ProbeRecord incidentProbe_m, transmittedProbe_m;
        std::vector<std::complex<double>> analyticIncident_m;
        std::vector<TiltedBinResult> bins_m;
        double drift_m = 0, scalarEnergy_m = 0, absoluteRatio_m = 0, incidentPeak_m = 0;
    };

    /// Exact solver time step and periodic material use the same resolved row as the modes.
    inline InterfaceRun tiltedRun(const TiltedRow& r) {
        validateTiltedProbe(r, true);
        validateTiltedProbe(r, false);
        return InterfaceRun(tiltedRunBox(r), TiltedMaterialMap(r),
                            {.safetyFactor_m = r.plane_m.courant_m, .timeStep_m = r.plane_m.dt_m});
    }

    /// Box probe retains the actual y-component positions and half-step channel labels.
    inline PlaneFourierProbe tiltedProbe(InterfaceSolver& state, const TiltedRow& r,
                                         bool incident) {
        std::vector<double> frequencies;
        for (double omega : r.plane_m.omega_m)
            frequencies.push_back(omega * r.plane_m.dt_m / (2 * Kokkos::numbers::pi));
        const auto quantity =
            r.plane_m.polarization_m == Polarization::TM ? Quantity::Magnetic : Quantity::Electric;
        return PlaneFourierProbe(state, tiltedProbeBox(r, incident), {{quantity, 1}}, frequencies,
                                 0, incident ? r.plane_m.gateIncident_m : r.plane_m.gateOther_m);
    }

    /// Full-rate Fourier accumulation remains independent of sparse energy or field output.
    inline void advanceTiltedMeasurement(InterfaceRun& run, const TiltedRow& r,
                                         TiltedResult& result, PlaneFourierProbe& incident,
                                         PlaneFourierProbe& other, int cadence,
                                         const InterfaceRun::Observer& observer) {
        double baseline = 0;
        run.advance(r.plane_m.steps_m, [&](long long step, auto& state, const auto& owner) {
            incident.sample();
            other.sample();
            if (step % cadence == 0 || step == r.plane_m.steps_m)
                observePlaneEnergy(state, owner, step, result, baseline);
            if (observer)
                observer(step, state, owner);
        });
    }

    /// Evaluate the known normal/tangent plane-wave basis at a lifted component coordinate.
    inline std::complex<double> tiltedBasis(const TiltedRow& r, std::array<double, 3> point,
                                            double normalK) {
        const auto uv = tiltedCoordinates(r, point);
        return std::exp(std::complex<double>(0, normalK * uv[0] + r.kappa_m * uv[1]));
    }

    /** @brief Unweighted complex least-squares amplitude of a unit-modulus plane-wave basis.
     * Deduction: minimizing \f$\sum|f-Ab|^2\f$ gives \f$\operatorname{mean}(f\overline b)\f$
     * when \f$|b|=1\f$.
     * Phase conventions follow [schneider2010UnderstandingFdtd, section 5.8.2].
     * This fits one propagating order; evanescent contamination is checked by the residual
     * and a separate probe-distance change, not removed by the fit.
     */
    inline std::complex<double> fitTiltedAmplitude(const TiltedRow& r,
                                                   const diagnostics::ProbeRecord& probe,
                                                   const std::vector<std::complex<double>>& values,
                                                   double normalK) {
        std::complex<double> result = 0;
        for (size_t point = 0; point < values.size(); ++point)
            result +=
                values[point] * std::conj(tiltedBasis(r, probe.coordinates_m[point], normalK));
        return result / double(values.size());
    }

    /// RMS remaining field, normalized later by the fitted incoming amplitude.
    inline double tiltedFitError(const TiltedRow& r, const diagnostics::ProbeRecord& probe,
                                 const std::vector<std::complex<double>>& values, double normalK,
                                 std::complex<double> amplitude) {
        double square = 0;
        for (size_t point = 0; point < values.size(); ++point)
            square += std::norm(values[point]
                                - amplitude * tiltedBasis(r, probe.coordinates_m[point], normalK));
        return std::sqrt(square / values.size());
    }

    /// Gather one frequency without changing the probe's original point ordering.
    inline std::vector<std::complex<double>> tiltedProbeValues(
        const diagnostics::ProbeRecord& probe, size_t bin) {
        std::vector<std::complex<double>> result;
        for (size_t point = 0; point < probe.coordinates_m.size(); ++point)
            result.push_back(probe.value(point, 0, bin));
        return result;
    }

    /// Preserve analytic incident samples at every incident-side point for independent
    /// re-projection.
    inline std::vector<std::complex<double>> tiltedAnalyticValues(const TiltedRow& r,
                                                                  TiltedResult& result,
                                                                  size_t bin) {
        std::vector<std::complex<double>> values;
        for (const auto& coordinate : result.incidentProbe_m.coordinates_m)
            values.push_back(tiltedIncidentTransform(r, coordinate, r.plane_m.omega_m[bin],
                                                     r.plane_m.gateIncident_m));
        result.analyticIncident_m.insert(result.analyticIncident_m.end(), values.begin(),
                                         values.end());
        return values;
    }

    /// Solve each propagation direction independently before projecting a measured bin.
    inline TiltedBinResult tiltedBinRoots(const TiltedRow& r, size_t bin) {
        TiltedBinResult value;
        const auto incomingRoots = tiltedRoots(r, r.plane_m.omega_m[bin], 1);
        value.incidentRoot_m     = incomingRoots[0];
        value.reflectedRoot_m    = incomingRoots[1];
        value.transmittedRoot_m  = tiltedRoots(r, r.plane_m.omega_m[bin], r.plane_m.epsilon_m)[0];
        return value;
    }

    /// Subtract the analytic incident pointwise, before fitting a different outgoing basis.
    inline std::vector<std::complex<double>> tiltedReflectedValues(
        const TiltedResult& result, size_t bin, const std::vector<std::complex<double>>& incident) {
        auto reflected = tiltedProbeValues(result.incidentProbe_m, bin);
        for (size_t point = 0; point < incident.size(); ++point)
            reflected[point] -= incident[point];
        return reflected;
    }

    /// Return the fitted amplitude with its unnormalized spatial RMS error.
    inline std::pair<std::complex<double>, double> projectTiltedField(
        const TiltedRow& r, const diagnostics::ProbeRecord& probe,
        const std::vector<std::complex<double>>& values, double normalK) {
        const auto amplitude = fitTiltedAmplitude(r, probe, values, normalK);
        return {amplitude, tiltedFitError(r, probe, values, normalK, amplitude)};
    }

    /// Fit incoming, incident-subtracted reflected, and transmitted amplitudes independently.
    inline TiltedBinResult fitTiltedBin(const TiltedRow& r, TiltedResult& result, size_t bin) {
        auto value             = tiltedBinRoots(r, bin);
        const auto incident    = tiltedAnalyticValues(r, result, bin);
        const auto reflected   = tiltedReflectedValues(result, bin, incident);
        const auto transmitted = tiltedProbeValues(result.transmittedProbe_m, bin);
        std::tie(value.incident_m, value.incidentResidual_m) =
            projectTiltedField(r, result.incidentProbe_m, incident, value.incidentRoot_m.normalK_m);
        std::tie(value.reflected_m, value.reflectedResidual_m) = projectTiltedField(
            r, result.incidentProbe_m, reflected, value.reflectedRoot_m.normalK_m);
        std::tie(value.transmitted_m, value.transmittedResidual_m) = projectTiltedField(
            r, result.transmittedProbe_m, transmitted, value.transmittedRoot_m.normalK_m);
        return value;
    }

    /// Report the first failed finite-data, normalization or spatial-fit criterion.
    inline std::string tiltedFitReason(const TiltedRow& r, const TiltedBinResult& bin,
                                       double peak) {
        const double amplitude = std::abs(bin.incident_m);
        if (!finitePlaneSpectrum(bin.incident_m) || !finitePlaneSpectrum(bin.reflected_m)
            || !finitePlaneSpectrum(bin.transmitted_m) || !(amplitude > 0)
            || !std::isfinite(amplitude))
            return "nonfinite_or_zero_incident";
        if (amplitude < r.plane_m.incidentFloor_m * peak)
            return "weak_incident";
        for (double residual :
             {bin.incidentResidual_m, bin.reflectedResidual_m, bin.transmittedResidual_m}) {
            if (!std::isfinite(residual))
                return "nonfinite_fit_residual";
            if (residual > r.fitResidualLimit_m)
                return "spatial_fit_residual";
        }
        return {};
    }

    /// Invalid bins retain raw point data and fitted amplitudes but never plausible ratios.
    inline void qualifyTiltedFit(const TiltedRow& r, TiltedBinResult& bin, double peak) {
        const double amplitude = std::abs(bin.incident_m),
                     nan       = std::numeric_limits<double>::quiet_NaN();
        bin.r_m = bin.t_m = {nan, nan};
        if (amplitude > 0 && std::isfinite(amplitude)) {
            bin.incidentResidual_m /= amplitude;
            bin.reflectedResidual_m /= amplitude;
            bin.transmittedResidual_m /= amplitude;
        }
        bin.reason_m = tiltedFitReason(r, bin, peak);
        if (!bin.valid())
            return;
        bin.r_m = bin.reflected_m / bin.incident_m;
        bin.t_m = bin.transmitted_m / bin.incident_m;
        if (!finitePlaneSpectrum(bin.r_m) || !finitePlaneSpectrum(bin.t_m)) {
            bin.reason_m = "nonfinite_coefficients";
            bin.r_m = bin.t_m = {nan, nan};
        }
    }

    /// Gather each point spectrum once; all fits are root-only derived products of those raw data.
    inline void extractTiltedResult(const TiltedRow& r, PlaneFourierProbe& incident,
                                    PlaneFourierProbe& other, TiltedResult& result) {
        result.incidentProbe_m    = incident.gather();
        result.transmittedProbe_m = other.gather();
        if (ippl::Comm->rank() != 0)
            return;
        for (size_t bin = 0; bin < r.plane_m.omega_m.size(); ++bin) {
            result.bins_m.push_back(fitTiltedBin(r, result, bin));
            result.incidentPeak_m =
                std::max(result.incidentPeak_m, std::abs(result.bins_m.back().incident_m));
        }
        for (auto& bin : result.bins_m)
            qualifyTiltedFit(r, bin, result.incidentPeak_m);
    }

    /// One initialization and one evolution, with any field writer borrowing the caller-owned run.
    inline TiltedResult measureTiltedRun(
        InterfaceRun& run, const TiltedRow& r, int cadence = 10,
        InterfaceRun::Observer observer                                                    = {},
        std::function<void(InterfaceSolver&, const material::MaterialOwner&)> initialCheck = {}) {
        run.load(TiltedPulse(r));
        if (initialCheck)
            initialCheck(run.solver(), run.material());
        TiltedResult result;
        auto incident = tiltedProbe(run.solver(), r, true),
             other    = tiltedProbe(run.solver(), r, false);
        advanceTiltedMeasurement(run, r, result, incident, other, cadence, observer);
        extractTiltedResult(r, incident, other, result);
        return result;
    }

    /// Convenience for checks without external output borrowing the distributed field owners.
    inline TiltedResult measureTilted(const TiltedRow& r, int cadence = 10) {
        auto run = tiltedRun(r);
        return measureTiltedRun(run, r, cadence);
    }
}  // namespace chdr::apps::planarInterface
