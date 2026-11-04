/** @file PlaneInput.h
 * @brief Strict application-local input for the common plane measurement path.
 */
#pragma once
#include <string>

#include "apps/planar-interface/PlanePresets.h"
#include "core/input/RunInput.h"

namespace chdr::apps::planarInterface {
    /// Read a finite scalar; strings and arrays are not silently converted into numbers.
    inline double planeParameter(const input::RunInput& input, const std::string& key,
                                 double fallback) {
        if (!input.parameters_m.has_path(key))
            return fallback;
        const auto node = input.parameters_m[key];
        if (!node.dtype().is_number() || node.number_of_elements() != 1)
            throw std::invalid_argument("Plane parameter '" + key + "' requires one number");
        const double value = node.to_float64();
        if (!std::isfinite(value))
            throw std::invalid_argument("Plane parameter '" + key + "' must be finite");
        return value;
    }

    /// Keep integer controls exact and reject narrowing overflow.
    inline int planeIntegerParameter(const input::RunInput& input, const std::string& key,
                                     int fallback) {
        const double value = planeParameter(input, key, fallback);
        if (std::trunc(value) != value)
            throw std::invalid_argument("Plane parameter '" + key + "' requires an integer");
        return planeInteger(value);
    }

    /// Keep enumerated values explicitly typed and report unsupported choices at the caller.
    inline std::string planeStringParameter(const input::RunInput& input, const std::string& key,
                                            const std::string& fallback) {
        if (!input.parameters_m.has_path(key))
            return fallback;
        if (!input.parameters_m[key].dtype().is_string())
            throw std::invalid_argument("Plane parameter '" + key + "' requires a string");
        return input.parameters_m[key].as_string();
    }

    /// Parse the two physical selectors without inferring polarization or incidence from a name.
    inline void readPlaneSelectors(const input::RunInput& input, PlaneRequest& r) {
        const auto polarization = planeStringParameter(input, "polarization", "TE");
        if (polarization != "TE" && polarization != "TM")
            throw std::invalid_argument("Plane polarization must be TE or TM");
        r.polarization_m     = polarization == "TE" ? Polarization::TE : Polarization::TM;
        const auto incidence = planeStringParameter(input, "incidence", "vacuum");
        if (incidence != "vacuum" && incidence != "dielectric")
            throw std::invalid_argument("Plane incidence must be vacuum or dielectric");
        r.incidentFromVacuum_m = incidence == "vacuum";
        if (input.parameters_m.has_path("angle_degrees"))
            r.angleDegrees_m = planeParameter(input, "angle_degrees", 0);
        if (input.parameters_m.has_path("transverse_period"))
            r.transversePeriod_m = planeParameter(input, "transverse_period", 0);
    }

    /// Read dimensionless shape, resolution and pulse controls; h_x is the code length unit.
    inline PlaneRequest readPlaneRequest(const input::RunInput& input) {
        PlaneRequest r;
        readPlaneSelectors(input, r);
        r.epsilon_m = planeParameter(input, "epsilon", r.epsilon_m);
        r.courant_m = planeParameter(input, "courant", r.courant_m);
        r.cellsPerWavelength_m =
            planeParameter(input, "cells_per_wavelength", r.cellsPerWavelength_m);
        r.h_m[1]       = planeParameter(input, "cell_shape_y", 1);
        r.h_m[2]       = planeParameter(input, "cell_shape_z", 1);
        r.tailSigmas_m = planeParameter(input, "tail_sigmas", r.tailSigmas_m);
        if (input.parameters_m.has_path("cutoff_fraction"))
            r.cutoffFraction_m = planeParameter(input, "cutoff_fraction", 0);
        r.probeCells_m       = planeIntegerParameter(input, "probe_cells", r.probeCells_m);
        r.timingSamples_m    = planeIntegerParameter(input, "timing_samples", r.timingSamples_m);
        r.gatePaddingSteps_m = planeIntegerParameter(input, "gate_padding_steps", 0);
        r.separationFactor_m = planeParameter(input, "separation_factor", 1);
        return r;
    }

    /// Explicit presets retain recorded regression inputs; their physics cannot be partly
    /// overridden.
    inline PlaneRow readPlanePreset(const input::RunInput& input) {
        input::validateParameters(input, {"preset"});
        const auto preset = planeStringParameter(input, "preset", "normal");
        if (preset == "normal")
            return normalPlaneRow();
        if (preset == "brewster")
            return brewsterPlaneRow();
        if (preset == "dielectric")
            return dielectricPlaneRow();
        throw std::invalid_argument("Plane preset must be normal, brewster or dielectric");
    }

    /// Empty input retains the recorded normal example; any physical request is resolved afresh.
    inline PlaneRow planeRow(const input::RunInput& input) {
        if (input.parameters_m.dtype().is_empty() || input.parameters_m.number_of_children() == 0
            || input.parameters_m.has_path("preset"))
            return readPlanePreset(input);
        input::validateParameters(
            input, {"polarization", "incidence", "epsilon", "courant", "cells_per_wavelength",
                    "angle_degrees", "transverse_period", "cell_shape_y", "cell_shape_z",
                    "tail_sigmas", "cutoff_fraction", "probe_cells", "timing_samples",
                    "gate_padding_steps", "separation_factor"});
        return resolvePlaneRow(readPlaneRequest(input));
    }
}  // namespace chdr::apps::planarInterface
