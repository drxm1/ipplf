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

    /// Read the numerical measurement band independently of the pulse's initialization spectrum.
    inline void readPlaneBand(const input::RunInput& input, PlaneRequest& r) {
        r.frequencyBins_m = planeIntegerParameter(input, "frequency_bins", r.frequencyBins_m);
        r.bandSigmas_m    = planeParameter(input, "band_sigmas", r.bandSigmas_m);
        r.incidentFloor_m = planeParameter(input, "incident_floor", r.incidentFloor_m);
    }

    /// Read cell_size and the material-interface coordinate face_position in code length units.
    inline void readPlaneGeometry(const input::RunInput& input, PlaneRequest& r) {
        r.h_m[0] = planeParameter(input, "cell_size", 1);
        r.h_m[1] = r.h_m[0] * planeParameter(input, "cell_shape_y", 1);
        r.h_m[2] = r.h_m[0] * planeParameter(input, "cell_shape_z", 1);
        if (input.parameters_m.has_path("face_position"))
            r.facePosition_m = planeParameter(input, "face_position", 0);
    }

    /// Read physical cell size and shape independently from pulse and resolution controls.
    inline PlaneRequest readPlaneRequest(const input::RunInput& input) {
        PlaneRequest r;
        readPlaneSelectors(input, r);
        readPlaneBand(input, r);
        readPlaneGeometry(input, r);
        r.epsilon_m = planeParameter(input, "epsilon", r.epsilon_m);
        r.courant_m = planeParameter(input, "courant", r.courant_m);
        r.cellsPerWavelength_m =
            planeParameter(input, "cells_per_wavelength", r.cellsPerWavelength_m);
        r.tailSigmas_m = planeParameter(input, "tail_sigmas", r.tailSigmas_m);
        if (input.parameters_m.has_path("cutoff_fraction"))
            r.cutoffFraction_m = planeParameter(input, "cutoff_fraction", 0);
        r.probeCells_m       = planeIntegerParameter(input, "probe_cells", r.probeCells_m);
        r.timingSamples_m    = planeIntegerParameter(input, "timing_samples", r.timingSamples_m);
        r.gatePaddingSteps_m = planeIntegerParameter(input, "gate_padding_steps", 0);
        r.separationFactor_m = planeParameter(input, "separation_factor", 1);
        return r;
    }

    /// The three named configurations are explicit inputs to the common row, not evolution paths.
    inline PlaneRow namedPlanePreset(const std::string& preset) {
        if (preset == "normal")
            return normalPlaneRow();
        if (preset == "brewster")
            return brewsterPlaneRow();
        if (preset == "dielectric")
            return dielectricPlaneRow();
        throw std::invalid_argument("Plane preset must be normal, brewster or dielectric");
    }

    /// Preserve preset physics while allowing a denser observation band without changing gates.
    inline PlaneRow readPlanePreset(const input::RunInput& input) {
        input::validateParameters(input,
                                  {"preset", "frequency_bins", "band_sigmas", "incident_floor"});
        auto p = namedPlanePreset(planeStringParameter(input, "preset", "normal"));
        PlaneRequest r;
        readPlaneBand(input, r);
        selectPlaneBins(p, r.frequencyBins_m, r.bandSigmas_m, r.incidentFloor_m);
        return p;
    }

    /// Observation-only requests retain the default physical preset instead of resolving a new box.
    inline bool onlyPlaneBandParameters(const input::RunInput& input) {
        for (conduit_index_t i = 0; i < input.parameters_m.number_of_children(); ++i) {
            const auto name = input.parameters_m.child(i).name();
            if (name != "frequency_bins" && name != "band_sigmas" && name != "incident_floor")
                return false;
        }
        return true;
    }

    /// Empty input retains the recorded normal example; any physical request is resolved afresh.
    inline PlaneRow planeRow(const input::RunInput& input) {
        if (input.parameters_m.has_path("preset") || onlyPlaneBandParameters(input))
            return readPlanePreset(input);
        input::validateParameters(
            input,
            {"polarization",          "incidence",      "epsilon",           "courant",
             "cells_per_wavelength",  "angle_degrees",  "transverse_period", "cell_size",
             "cell_shape_y",          "cell_shape_z",   "face_position",     "tail_sigmas",
             "cutoff_fraction",       "probe_cells",    "timing_samples",    "gate_padding_steps",
             "separation_factor",     "frequency_bins", "band_sigmas",       "incident_floor",
             "interface_tilt_degrees"});
        return resolvePlaneRow(readPlaneRequest(input));
    }
}  // namespace chdr::apps::planarInterface
