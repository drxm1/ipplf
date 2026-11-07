/** @file PacketCase.h
 * @brief Existing periodic packet cases expressed independently of input and output handling.
 */
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

#include "apps/planar-interface/PacketStart.h"
#include "core/material/PlanarSlab.h"
#include "core/simulation/SimulationConfig.h"

namespace chdr::apps::planarInterface {
    /// Numerical settings and initial data of one existing periodic packet case.
    struct PacketCase {
        std::string name_m;
        simulation::ChdrCartesianGridSettings box_m;
        simulation::ChdrSolverSettings solver_m;
        material::PlanarSlab slab_m;
        PacketStart start_m;
        int steps_m;
        int frameEvery_m = 10;
    };

    /// Retain the case's original pulse position, widths, direction and polarization.
    inline PacketStart packetStart(const std::string& name, int cells,
                                   const material::PlanarSlab& slab) {
        const bool brewster = name.starts_with("brewster");
        const double scale  = cells / 512.0;
        PacketStart start;
        start.slab_m            = slab;
        start.centreX_m         = 371 * scale;
        start.centreZ_m         = cells / (brewster ? 4.0 : 2.0);
        start.length_m          = cells;
        start.width_m           = 29 * scale;
        start.transverseWidth_m = 45 * scale;
        start.tm_m              = name == "brewster-tm";
        if (brewster) {
            start.directionX_m = -1 / std::sqrt(3.13);
            start.directionZ_m = std::sqrt(2.13 / 3.13);
        }
        return start;
    }

    /// Select an existing case; the small mesh retains the four-step smoke configuration.
    inline PacketCase packetCase(const std::string& name, bool smoke) {
        if (name != "half-space" && name != "slab" && name != "brewster-te"
            && name != "brewster-tm")
            throw std::invalid_argument("Packet cases: half-space, slab, brewster-te, brewster-tm");
        const int cells    = smoke ? 32 : 512;
        const double scale = cells / 512.0;
        const material::PlanarSlab slab{0, name == "slab" ? 156 * scale : 0, 230 * scale};
        return {name,
                {{cells, 2, cells}, {1, 1, 1}},
                {},
                slab,
                packetStart(name, cells, slab),
                smoke ? 4 : name.starts_with("brewster") ? 1000 : 560};
    }
}  // namespace chdr::apps::planarInterface
