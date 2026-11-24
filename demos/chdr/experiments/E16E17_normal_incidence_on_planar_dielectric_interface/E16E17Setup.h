/** @file E16E17Setup.h
 * @brief Fixed material, mesh and observation schedule of E16.typ and E17.typ.
 */
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "core/input/RunInput.h"
#include "core/material/PlanarSlab.h"
#include "core/simulation/ChDRSimulationManager.h"
#include "experiments/E16E17_normal_incidence_on_planar_dielectric_interface/E16E17AnalyticContinuumReference.h"

namespace chdr::experiments::e16e17 {
    using Simulation = simulation::ChDRSimulationManager<>;
    using Solver     = Simulation::Solver;
    using Reference  = E16E17AnalyticContinuumReferenceSampler;

    /// Read one finite numeric scalar without accepting strings or arrays.
    inline double scalarParameter(const input::ChdrYamlRunInput& request, const char* key) {
        const auto node = request.parameters_m[key];
        if (!node.dtype().is_number() || node.number_of_elements() != 1)
            throw std::invalid_argument(std::string(key) + " must be a numeric scalar");
        const double value = node.to_float64();
        if (!std::isfinite(value))
            throw std::invalid_argument(std::string(key) + " must be finite");
        return value;
    }

    /// Select the five mesh levels specified for both directions of incidence.
    inline int cellsPerPulse(const input::ChdrYamlRunInput& request) {
        if (!request.parameters_m.has_path("cells_per_pulse"))
            throw std::invalid_argument("E16/E17 requires cells_per_pulse");
        const double value = scalarParameter(request, "cells_per_pulse");
        for (int allowed : {128, 256, 512, 1024, 2048})
            if (value == allowed)
                return allowed;
        throw std::invalid_argument("E16/E17 cells_per_pulse must be 128, 256, 512, 1024 or 2048");
    }

    /// Retain the prescribed three snapshots while allowing storage to be disabled for checks.
    inline void validateFieldOptions(const input::ChdrYamlRunInput& request) {
        if (!request.observations_m.contains("fields"))
            return;
        const auto& fields = request.observations_m.at("fields");
        if (fields.every_m || fields.first_m || fields.last_m || fields.includeFinal_m)
            throw std::invalid_argument("E16/E17 fields accepts only enabled, with fixed targets");
    }

    /// Movie profiles keep the initial and final states with an optional positive cadence.
    inline void validateMovieOptions(const input::ChdrYamlRunInput& request) {
        if (!request.observations_m.contains("movie_profiles"))
            return;
        const auto& movie = request.observations_m.at("movie_profiles");
        if (movie.first_m || movie.last_m || movie.includeFinal_m
            || (movie.every_m && *movie.every_m <= 0))
            throw std::invalid_argument(
                "E16/E17 movie_profiles accepts enabled and positive every");
    }

    /// Validate the named physical case before allocation or creation of output files.
    inline void validateInput(const input::ChdrYamlRunInput& request) {
        if ((request.caseName_m != "E16" && request.caseName_m != "E17")
            || request.units_m.system_m != "normalized" || request.outputDirectory_m.empty()
            || request.duration_m)
            throw std::invalid_argument(
                "E16/E17 requires its named case, normalized units and output");
        input::validateParameters(request, {"cells_per_pulse", "carrier_wavelength"});
        input::validateObservationNames(request, {"fields", "movie_profiles"});
        validateFieldOptions(request);
        validateMovieOptions(request);
    }

    /// Exchange the two media for E17, with the same incident profile and face position.
    inline Reference makeReference(const input::ChdrYamlRunInput& request) {
        const double wavelength = request.parameters_m.has_path("carrier_wavelength")
                                      ? scalarParameter(request, "carrier_wavelength")
                                      : 0.25;
        if (wavelength <= 0)
            throw std::invalid_argument("carrier_wavelength must be positive");
        const bool fromVacuum = request.caseName_m == "E16";
        return Reference(fromVacuum ? 1.0 : 2.13, fromVacuum ? 2.13 : 1.0, 2.5, 4.0, wavelength);
    }

    /// The transverse mesh is fixed and MPI partitions only the propagation direction.
    inline simulation::ChdrCartesianGridSettings makeGrid(int pulseCells) {
        return {{8, 8, 5 * pulseCells}, {0.125, 0.125, 2.0 / pulseCells}, 0, {true, true, true}, 2};
    }

    /// Sample the specified dielectric at cell centers, with its other face at the periodic wrap.
    inline auto materialMap(const Reference& reference) {
        return material::TwoMaterialMap<material::PlanarSlab>{
            {2, reference.zf_m, 10.0}, reference.epsilon2_m, reference.epsilon1_m};
    }

    /// One physical target selects separate potential and reconstructed-field states.
    struct Checkpoint {
        const char* name_m;
        double target_m;
        int potentialStep_m, fieldStep_m;
    };

    /// Invert LevelContract labels and choose the later state when two samples are equally near.
    inline int nearestStep(const Solver& solver, Quantity quantity, double target) {
        return static_cast<int>(
            std::floor((target - solver.labelTime(quantity, 0)) / solver.getDt() + 0.5));
    }

    /// Initial fields use the loaded state at -dt/2, as in the initial observation of E01.typ.
    inline std::array<Checkpoint, 3> checkpoints(const Solver& solver, const Reference& reference) {
        std::array<Checkpoint, 3> points{{{"initial", 0, 0, 0},
                                          {"midpoint", 1.5 / reference.v1_m, 0, 0},
                                          {"final", 3.0 / reference.v1_m, 0, 0}}};
        for (unsigned index = 1; index < points.size(); ++index) {
            auto& point           = points[index];
            point.potentialStep_m = nearestStep(solver, Quantity::VectorPotential, point.target_m);
            point.fieldStep_m     = nearestStep(solver, Quantity::Electric, point.target_m);
        }
        return points;
    }

    /// Require every run to include the full experiment's last required native sample.
    inline int selectedSteps(const input::ChdrYamlRunInput& request,
                             const std::array<Checkpoint, 3>& points) {
        const int full = std::max(points.back().potentialStep_m, points.back().fieldStep_m);
        if (request.steps_m && *request.steps_m != full)
            throw std::invalid_argument("E16/E17 run.steps must equal its prescribed final step");
        return full;
    }

    /// The full-domain arrays share the field observation indices.
    inline std::vector<int> snapshotSteps(const std::array<Checkpoint, 3>& points) {
        std::vector<int> steps;
        for (const auto& point : points)
            steps.push_back(point.fieldStep_m);
        return steps;
    }
}  // namespace chdr::experiments::e16e17
