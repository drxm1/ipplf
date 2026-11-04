/** @file TiltedInput.h
 * @brief Resolve a requested physical tilt without rounding that angle to a square diagonal.
 */
#pragma once
#include <tuple>

#include "apps/planar-interface/PlaneInput.h"
#include "apps/planar-interface/TiltedSetup.h"

namespace chdr::apps::planarInterface {
    /// Common physical checks allow a signed normal offset in the tilted geometry.
    inline void validateTiltedPhysicalRequest(const TiltedRequest& input) {
        auto base = input.plane_m;
        if (base.transversePeriod_m || !base.incidentFromVacuum_m)
            throw std::invalid_argument(
                "Tilted planes currently require vacuum incidence and angle_degrees input");
        if (base.facePosition_m && !std::isfinite(*base.facePosition_m))
            throw std::invalid_argument(
                "Tilted face_position must be a finite signed normal offset");
        base.facePosition_m.reset();
        validatePlaneRequest(base);
    }

    /// Primitive translations and pinned controls must be complete and mutually compatible.
    inline void validateTiltedPrimitiveRequest(const TiltedRequest& input) {
        if ((input.primitiveX_m == 0) != (input.primitiveZ_m == 0))
            throw std::invalid_argument("Supply both primitive tangent cell counts");
        if (input.primitiveX_m
            && (input.primitiveX_m < 1 || input.primitiveZ_m < 1 || input.primitiveX_m > 16
                || input.primitiveZ_m > 16
                || std::gcd(input.primitiveX_m, input.primitiveZ_m) != 1))
            throw std::invalid_argument(
                "Primitive tangent cell counts must be positive coprime integers <=16");
        if ((input.periodMultiplier_m > 0) != input.tangentMode_m.has_value()
            || (input.periodMultiplier_m && !input.primitiveX_m))
            throw std::invalid_argument(
                "Pinned tilted periods require period_cells_multiplier, tangent_mode and an "
                "explicit primitive pair");
    }

    /// Validate independent choices before the periodic compatibility search.
    inline void validateTiltedRequest(const TiltedRequest& input) {
        validateTiltedPhysicalRequest(input);
        validateTiltedPrimitiveRequest(input);
        const double angle = input.plane_m.angleDegrees_m.value_or(0);
        if (!std::isfinite(input.tiltDegrees_m) || std::abs(input.tiltDegrees_m) <= 0
            || std::abs(input.tiltDegrees_m) >= 90 || !std::isfinite(angle) || std::abs(angle) >= 90
            || input.probeHalfWidth_m < 0 || input.probeHalfWidth_m > 16
            || !(input.fitResidualLimit_m > 0 && input.fitResidualLimit_m < 1))
            throw std::invalid_argument("Invalid tilted orientation, probe or fit controls");
    }

    /// Complete controls that do not depend on the primitive metric adjustment.
    inline void completeTiltedControls(TiltedRow& r, const TiltedRequest& input) {
        r.plane_m.requestedAngleDegrees_m = input.plane_m.angleDegrees_m.value_or(0);
        r.plane_m.cutoffFraction_m        = input.plane_m.cutoffFraction_m.value_or(
            r.plane_m.requestedAngleDegrees_m == 0 ? 0 : 0.5);
        r.probeCells_m       = input.plane_m.probeCells_m;
        r.probeHalfWidth_m   = input.probeHalfWidth_m;
        r.fitResidualLimit_m = input.fitResidualLimit_m;
    }

    /// Use exactly the requested angle and choose hz to make the primitive translation periodic.
    inline TiltedRow tiltedGeometry(const TiltedRequest& input, int x, int z) {
        TiltedRow r;
        const auto& original = input.plane_m;
        auto base            = original;
        r.tiltDegrees_m      = input.tiltDegrees_m;
        r.facePosition_m     = original.facePosition_m.value_or(0);
        r.cosine_m           = std::cos(input.tiltDegrees_m * Kokkos::numbers::pi / 180);
        r.sine_m             = std::sin(input.tiltDegrees_m * Kokkos::numbers::pi / 180);
        r.primitiveX_m       = x;
        r.primitiveZ_m       = z;
        r.requestedShape_m   = {1, base.h_m[1] / base.h_m[0], base.h_m[2] / base.h_m[0]};
        base.h_m[2]          = r.cosine_m * x * base.h_m[0] / (std::abs(r.sine_m) * z);
        r.plane_m            = planeBaseRow(base);
        r.rowSpacing_m       = r.cosine_m * base.h_m[0] / z;
        r.primitivePeriod_m  = x * base.h_m[0] / std::abs(r.sine_m);
        completeTiltedControls(r, input);
        return r;
    }

    /// Pin compatible periods for controls that must retain actual frequency and kappa.
    inline void resolvePinnedTiltedPeriod(TiltedRow& r, const TiltedRequest& input) {
        setTiltedPeriods(r, input.periodMultiplier_m);
        r.tangentMode_m   = *input.tangentMode_m;
        const double sine = std::sin(r.plane_m.requestedAngleDegrees_m * Kokkos::numbers::pi / 180);
        if ((sine == 0 && r.tangentMode_m != 0) || (sine != 0 && r.tangentMode_m * sine <= 0))
            throw std::invalid_argument("Pinned tangent_mode is inconsistent with angle_degrees");
        r.kappa_m = 2 * Kokkos::numbers::pi * r.tangentMode_m / r.tangentPeriod_m;
        if (sine != 0)
            r.plane_m.omega0_m = r.kappa_m / sine;
        completeTiltedSpectrum(r);
        if (tiltedRequiredPeriod(r) > r.normalPeriod_m)
            throw std::invalid_argument(
                "Pinned tilted period cannot isolate the requested probes and gates");
    }

    /// Complete one primitive-pair candidate, including every retained carrier's admissibility.
    inline TiltedRow resolveTiltedCandidate(const TiltedRequest& input, int x, int z) {
        auto r = tiltedGeometry(input, x, z);
        if (!(r.plane_m.cutoffFraction_m >= 0 && r.plane_m.cutoffFraction_m < 1)
            || (r.plane_m.requestedAngleDegrees_m != 0 && r.plane_m.cutoffFraction_m == 0))
            throw std::invalid_argument(
                "Tilted oblique pulses require positive cutoff_fraction below one");
        if (input.periodMultiplier_m)
            resolvePinnedTiltedPeriod(r, input);
        else
            resolveTiltedPeriod(r, r.plane_m.omega0_m);
        const int dielectricRows = planeInteger(std::round(r.dielectricWidth_m / r.rowSpacing_m));
        r.incidentRows_m         = r.periodRows_m - dielectricRows;
        r.incidentWidth_m        = r.incidentRows_m * r.rowSpacing_m;
        selectTiltedModes(r);
        selectTiltedBins(r, input.plane_m);
        return r;
    }

    /// Search small primitive translations by the change from the requested cell shape.
    inline std::vector<std::tuple<double, int, int>> tiltedPrimitiveCandidates(
        const TiltedRequest& input) {
        std::vector<std::tuple<double, int, int>> result;
        const double alpha = std::abs(input.tiltDegrees_m) * Kokkos::numbers::pi / 180;
        for (int x = 1; x <= 6; ++x)
            for (int z = 1; z <= 6; ++z)
                if (std::gcd(x, z) == 1) {
                    const double shape = x / (z * std::tan(alpha));
                    result.emplace_back(
                        std::abs(std::log(shape * input.plane_m.h_m[0] / input.plane_m.h_m[2])), x,
                        z);
                }
        std::sort(result.begin(), result.end());
        return result;
    }

    /// Select the closest qualified shape; record the adjustment rather than changing the tilt.
    inline TiltedRow resolveTiltedRow(const TiltedRequest& input) {
        validateTiltedRequest(input);
        if (input.primitiveX_m)
            return resolveTiltedCandidate(input, input.primitiveX_m, input.primitiveZ_m);
        std::string reason;
        for (const auto& [change, x, z] : tiltedPrimitiveCandidates(input)) {
            try {
                return resolveTiltedCandidate(input, x, z);
            } catch (const std::invalid_argument& error) {
                reason = error.what();
            }
        }
        throw std::invalid_argument("No qualified small primitive tangent translation: " + reason);
    }

    /// Signed tangent indices support either incidence sign without narrowing silently.
    inline void readTiltedPeriod(const input::RunInput& request, TiltedRequest& result) {
        result.periodMultiplier_m = planeIntegerParameter(request, "period_cells_multiplier", 0);
        if (request.parameters_m.has_path("tangent_mode")) {
            const double value = planeParameter(request, "tangent_mode", 0);
            if (value != std::trunc(value))
                throw std::invalid_argument("tangent_mode requires a signed integer");
            const int magnitude  = planeInteger(std::abs(value));
            result.tangentMode_m = value < 0 ? -magnitude : magnitude;
        }
    }

    /// Read the tilted-only controls; the physical incidence angle remains angle_degrees.
    inline TiltedRequest readTiltedRequest(const input::RunInput& request) {
        TiltedRequest result;
        result.plane_m = readPlaneRequest(request);
        readTiltedPeriod(request, result);
        if (!request.parameters_m.has_path("probe_cells"))
            result.plane_m.probeCells_m = 12;
        if (!request.parameters_m.has_path("timing_samples"))
            result.plane_m.timingSamples_m = 129;
        result.tiltDegrees_m      = planeParameter(request, "interface_tilt_degrees", 0);
        result.primitiveX_m       = planeIntegerParameter(request, "primitive_cells_x", 0);
        result.primitiveZ_m       = planeIntegerParameter(request, "primitive_cells_z", 0);
        result.probeHalfWidth_m   = planeIntegerParameter(request, "probe_half_width_cells", 2);
        result.fitResidualLimit_m = planeParameter(request, "fit_residual_limit", 1e-3);
        return result;
    }

    /// Keep unsupported tilted physical requests explicit rather than silently ignoring them.
    inline TiltedRow tiltedRow(const input::RunInput& request) {
        input::validateParameters(request, {"polarization",
                                            "incidence",
                                            "epsilon",
                                            "courant",
                                            "cells_per_wavelength",
                                            "angle_degrees",
                                            "cell_size",
                                            "cell_shape_y",
                                            "cell_shape_z",
                                            "face_position",
                                            "tail_sigmas",
                                            "cutoff_fraction",
                                            "probe_cells",
                                            "timing_samples",
                                            "gate_padding_steps",
                                            "separation_factor",
                                            "frequency_bins",
                                            "band_sigmas",
                                            "incident_floor",
                                            "interface_tilt_degrees",
                                            "primitive_cells_x",
                                            "primitive_cells_z",
                                            "period_cells_multiplier",
                                            "tangent_mode",
                                            "probe_half_width_cells",
                                            "fit_residual_limit"});
        return resolveTiltedRow(readTiltedRequest(request));
    }
}  // namespace chdr::apps::planarInterface
