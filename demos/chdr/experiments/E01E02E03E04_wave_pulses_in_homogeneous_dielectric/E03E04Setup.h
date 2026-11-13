/** @file E03E04Setup.h
 * @brief Shared physical configuration of the spherical pulse experiments.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "core/simulation/ChDRSimulationManager.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03.h"

namespace chdr::experiments::e03 {
    using SphericalSimulation = simulation::ChDRSimulationManager<>;

    /// Validate the cubic resolution before converting the YAML number to an index.
    inline int sphericalCells(const input::ChdrYamlRunInput& request) {
        if (request.parameters_m.has_path("cells_per_axis")) {
            const auto node = request.parameters_m["cells_per_axis"];
            if (node.dtype().is_number() && node.number_of_elements() == 1) {
                const double cells = node.to_float64();
                if (std::find(E03Levels.begin(), E03Levels.end(), cells) != E03Levels.end())
                    return static_cast<int>(cells);
            }
        }
        throw std::invalid_argument("cells_per_axis must be 96, 144, 216 or 324");
    }

    /// Represent the centered physical cube on the manager's translated mesh.
    inline simulation::ChdrCartesianGridSettings sphericalGrid(int n, bool periodic) {
        const double h = 2 * E03BoxCentre / n;
        return {{n, n, n}, {h, h, h}, 0, {periodic, periodic, periodic}, 2};
    }

    /** @brief Read the optional sharpness \f$\beta\f$ of the compact radial profile.
     * @details The default \f$\beta=1\f$ is the profile of E01. A larger value gives a pulse
     * with a narrower spectrum inside the same support; see the reference sampler.
     */
    inline double sphericalSharpness(const input::ChdrYamlRunInput& request) {
        if (!request.parameters_m.has_path("profile_sharpness"))
            return 1;
        const auto node = request.parameters_m["profile_sharpness"];
        if (node.dtype().is_number() && node.number_of_elements() == 1) {
            const double sharpness = node.to_float64();
            if (std::isfinite(sharpness) && sharpness > 0)
                return sharpness;
        }
        throw std::invalid_argument("profile_sharpness must be a positive number");
    }

    /// The dipole pulse of E03.typ, centred in the translated mesh coordinates.
    inline E03Reference sphericalReference(double sharpness = 1) {
        return {1 / std::sqrt(E03Epsilon), E03BoxCentre, sharpness};
    }

    /** @brief Select the same bulk update for periodic propagation and Mur exit.
     * @details O1VS advances in physical time. The Mur coefficient sees the dielectric
     * speed [chew2014generalizedgauge, Eq. (13), mur1981AbsorbingBoundaryFdtd, Eq. (16)].
     */
    inline SphericalSimulation makeSphericalSimulation(
        const simulation::ChdrCartesianGridSettings& grid, const E03Reference& ref) {
        const auto dielectric =
            KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
            marked = false;
            return E03Epsilon;
        };
        const auto speed =
            grid.periodic_m[0] ? solver::MurSpeed::matched() : solver::MurSpeed::physical(ref.v_m);
        return SphericalSimulation{
            grid, dielectric, {.murSpeed_m = speed, .timeStep_m = grid.spacing_m[0] / 2}};
    }

    /// Invert LevelContract's label map to select the nearest native observation.
    inline int nearestSphericalStep(SphericalSimulation::Solver& state, Quantity quantity,
                                    double target) {
        return static_cast<int>(
            std::floor((target - state.labelTime(quantity, 0)) / state.getDt() + 0.5));
    }

    /// The only optional observation is `fields`: full-volume frames for visualization.
    inline void validateSphericalRun(const input::ChdrYamlRunInput& request) {
        input::validateObservationNames(request, {"fields"});
        if (request.units_m.system_m != "normalized" || request.duration_m
            || request.outputDirectory_m.empty())
            throw std::invalid_argument(
                "Spherical pulse runs require normalized units, fixed duration and "
                "output.directory");
    }
}  // namespace chdr::experiments::e03
