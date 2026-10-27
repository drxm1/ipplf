/** @file MovieCases.h
 * @brief Reproducible open-box material layouts and initial data for field movies.
 */
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "AnnulusStart.h"
#include "TMAnnulusStart.h"
#include "core/material/RectangleXZ.h"
#include "core/simulation/SimulationConfig.h"

namespace chdr::apps::interfaceMovies {
    /// Cell-centred upper-left half-plane; ties on the diagonal belong to the dielectric.
    struct DiagonalXZ {
        KOKKOS_INLINE_FUNCTION bool contains(double x, double, double z) const { return x >= z; }
    };

    /// One movie input: z horizontal, x vertical, and two periodic y cells.
    struct MovieCase {
        std::string name_m;
        simulation::RunBox box_m;
        material::RectangleXZ region_m;
        AnnulusStart start_m;
        double epsilon_m = 2.13;
        int steps_m = 1424, frameEvery_m = 8;
        /// Select the half-plane and demo face solver; region_m is unused.
        bool diagonal_m = false;
        bool tm_m       = false;  ///< Initialize a discrete-curl TM annulus instead of TE A_y.
    };

    /// Use the actual box spacings when constructing either a driver or a test's TM sampler.
    inline TMAnnulusStart movieTMStart(const MovieCase& input) {
        return {input.start_m, input.box_m.spacing_m[0], input.box_m.spacing_m[2]};
    }

    /// Common mesh and pulse scales; the small configuration is a four-step smoke run only.
    inline MovieCase movieDefaults(const std::string& name, bool smoke) {
        const int n        = smoke ? 64 : 600;
        const double scale = n / 600.0;
        MovieCase result{name, {{n, 2, n}, {1, 1, 1}, 0, {false, true, false}, 0}, {}, {}};
        result.start_m.centreX_m    = 150 * scale;
        result.start_m.centreZ_m    = 300 * scale;
        result.start_m.radius_m     = 72 * scale;
        result.start_m.wavelength_m = 36 * scale;
        if (smoke)
            result.steps_m = 4;
        return result;
    }

    /// Finite approximation of a half-space, with a vacuum margin at every open face.
    inline MovieCase halfSpaceMovie(bool smoke) {
        auto result    = movieDefaults("half-space", smoke);
        const double n = result.box_m.cells_m[0], gap = result.start_m.wavelength_m;
        result.region_m = {0.5 * n, n - gap, gap, n - gap};
        return result;
    }

    /// Floating material brick with its right edge above the annulus centre.
    inline MovieCase brickMovie(bool smoke) {
        auto result    = movieDefaults("brick", smoke);
        const double n = result.box_m.cells_m[0], gap = result.start_m.wavelength_m;
        result.region_m = {0.5 * n, 0.8 * n, gap, 0.5 * n};
        result.steps_m  = smoke ? 4 : 1120;
        return result;
    }

    /// Finite-thickness slab, separated from all open faces by vacuum.
    inline MovieCase slabMovie(bool smoke) {
        auto result    = movieDefaults("slab", smoke);
        const double n = result.box_m.cells_m[0], gap = result.start_m.wavelength_m;
        result.region_m = {0.5 * n, 0.65 * n, gap, n - gap};
        result.steps_m  = smoke ? 4 : 1280;
        return result;
    }

    /// Half-plane reaching the upper x and lower z faces, with a vacuum pulse below the diagonal.
    inline MovieCase diagonalMovie(bool smoke) {
        auto result              = movieDefaults("diagonal", smoke);
        result.diagonal_m        = true;
        result.start_m.centreZ_m = 0.75 * result.box_m.cells_m[2];
        result.steps_m           = smoke ? 4 : 1672;
        return result;
    }

    /// Select geometry independently of polarization.
    inline MovieCase movieGeometry(const std::string& name, bool smoke) {
        if (name == "half-space")
            return halfSpaceMovie(smoke);
        if (name == "brick")
            return brickMovie(smoke);
        if (name == "slab")
            return slabMovie(smoke);
        if (name == "diagonal")
            return diagonalMovie(smoke);
        throw std::invalid_argument("Movie geometries: half-space, brick, slab, diagonal");
    }

    /// Add -tm to any geometry name; an unqualified name retains its TE start.
    inline MovieCase movieCase(const std::string& name, bool smoke = false) {
        const bool tm       = name.ends_with("-tm");
        const auto geometry = name.substr(0, name.size() - (tm ? 3 : 0));
        auto result         = movieGeometry(geometry, smoke);
        result.name_m       = name;
        result.tm_m         = tm;
        if (smoke && tm) {
            result.start_m.centreX_m += 1;  // Keep ghost-extension input layers zero.
            if (result.diagonal_m)
                result.start_m.centreZ_m -= 1;
        }
        return result;
    }

    /** @brief Check that the enclosing support disk stays clear of material and open faces.
     * @details The previous potential time is earlier, so its enclosing radius is smaller.
     * TE uses a one-cell face margin. TM leaves three layers zero so quadratic open-ghost
     * extrapolation cannot introduce an initial divergence (OpenFaceClosure). The distance
     * to the rectangle is computed from its nearest point. For diagonal TM, subtract
     * the cell half-width projected on the normal: the resulting half-plane contains
     * every sampled dielectric cell, giving a conservative staircase-clearance bound.
     */
    inline void checkMovieStart(const MovieCase& input, double currentTime) {
        const auto& p      = input.start_m;
        const auto& r      = input.region_m;
        const auto lengths = input.box_m.lengths();
        const double radius =
            p.outerRadius(currentTime)
            + (input.tm_m ? 0.5 * std::max(input.box_m.spacing_m[0], input.box_m.spacing_m[2]) : 0);
        const double dx = p.centreX_m - std::clamp(p.centreX_m, r.lowerX_m, r.upperX_m);
        const double dz = p.centreZ_m - std::clamp(p.centreZ_m, r.lowerZ_m, r.upperZ_m);
        const double cellReach =
            input.tm_m ? 0.5 * (input.box_m.spacing_m[0] + input.box_m.spacing_m[2]) : 0;
        const double distance = input.diagonal_m
                                    ? (p.centreZ_m - p.centreX_m - cellReach) / std::sqrt(2.0)
                                    : std::hypot(dx, dz);
        const double layers   = input.tm_m ? 3 : 1;
        if (p.centreX_m - radius <= layers * input.box_m.spacing_m[0]
            || p.centreX_m + radius >= lengths[0] - layers * input.box_m.spacing_m[0]
            || p.centreZ_m - radius <= layers * input.box_m.spacing_m[2]
            || p.centreZ_m + radius >= lengths[2] - layers * input.box_m.spacing_m[2]
            || distance <= radius)
            throw std::invalid_argument(
                "Both annulus histories must start in vacuum clear of faces");
    }
}  // namespace chdr::apps::interfaceMovies
