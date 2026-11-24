/** @file E18E19Setup.h
 * @brief Mesh, material and observation steps for the oblique-incidence experiments.
 */
#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include "core/input/RunInput.h"
#include "core/material/PlanarSlab.h"
#include "core/simulation/ChDRSimulationManager.h"
#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19AnalyticContinuumReference.h"

namespace chdr::experiments::e18e19 {
    using Simulation = simulation::ChDRSimulationManager<>;
    using Solver     = Simulation::Solver;
    using Reference  = E18E19AnalyticContinuumReference;

    /// E18 fixes the duration to 28 tangential light-crossing times and dt to half a cell crossing.
    inline int finalStep(int cellsPerPeriod) {
        return 56 * cellsPerPeriod;
    }

    /// E18-S and E19-S select the solver time step through the safety factor.
    inline bool usesSafetyFactor(const input::ChdrYamlRunInput& request) {
        return request.caseName_m == "E18-S" || request.caseName_m == "E19-S";
    }

    /// Read the tangential resolution shared by the TE and TM refinement series.
    inline int cellsPerPeriod(const input::ChdrYamlRunInput& request) {
        if (request.parameters_m.has_path("cells_per_period")) {
            const auto cells = request.parameters_m["cells_per_period"];
            if (cells.dtype().is_number() && cells.number_of_elements() == 1)
                for (int allowed : {24, 48, 96, 192})
                    if (cells.to_float64() == allowed)
                        return allowed;
        }
        throw std::invalid_argument("E18/E19 cells_per_period must be 24, 48, 96 or 192");
    }

    /// Each input specifies one full study level, with the fixed duration and all observations.
    inline int validateInput(const input::ChdrYamlRunInput& request) {
        if ((request.caseName_m != "E18" && request.caseName_m != "E19"
             && !usesSafetyFactor(request))
            || request.units_m.system_m != "normalized" || request.outputDirectory_m.empty())
            throw std::invalid_argument(
                "E18/E19 requires its named case, normalized units and output");
        input::validateParameters(request, {"cells_per_period"});
        input::validateObservationNames(request, {});
        if (request.duration_m)
            throw std::invalid_argument("E18/E19 uses its prescribed final step, not run.duration");
        return cellsPerPeriod(request);
    }

    /// Cubic cells, two layers along y, and decomposition along the interface normal.
    inline simulation::ChdrCartesianGridSettings makeGrid(int cells) {
        const double spacing = 1.0 / cells;
        return {{cells, 2, static_cast<int>(Reference::BOX_LENGTH_Z) * cells},
                {spacing, spacing, spacing},
                0,
                {true, true, true},
                2};
    }

    /** @brief Select the time step prescribed by E18.typ, E19.typ, E18S.typ or E19S.typ.
     * @details The S cases delegate their stability bound to ChdrSolver through its checked
     * safety factor. The selected step is set before loading the potential histories.
     */
    inline simulation::ChdrSolverSettings solverSettings(
        const input::ChdrYamlRunInput& request, const simulation::ChdrCartesianGridSettings& grid) {
        if (usesSafetyFactor(request))
            return {.safetyFactor_m = 0.99};
        return {.timeStep_m = grid.spacing_m[0] / 2};
    }

    /** @brief Select the last reconstructed field sample prescribed by the experiment.
     * @details E18S.typ and E19S.typ use the E01 nearest-time rule at time 28.
     * Invert the constructed LevelContract through labelTime(), choosing the later sample
     * in a tie. E18 and E19 retain their original final step.
     */
    inline int selectedSteps(const input::ChdrYamlRunInput& request, const Solver& solver,
                             int cells) {
        const int steps =
            usesSafetyFactor(request)
                ? static_cast<int>(std::floor(
                      (28.0 - solver.labelTime(Quantity::Electric, 0)) / solver.getDt() + 0.5))
                : finalStep(cells);
        if (request.steps_m && *request.steps_m != steps)
            throw std::invalid_argument("E18/E19 run.steps must equal its prescribed final step");
        return steps;
    }

    /// Vacuum occupies the lower part of the periodic box, up to the specified cell face.
    inline auto materialMap() {
        return material::TwoMaterialMap<material::PlanarSlab>{
            {2, Reference::INTERFACE_Z, Reference::BOX_LENGTH_Z},
            Reference::DIELECTRIC_EPSILON,
            1.0};
    }

    /// Save the initial state and each successive quarter of the prescribed duration.
    inline std::vector<int> snapshotSteps(int steps) {
        return {0, steps / 4, steps / 2, 3 * steps / 4, steps};
    }
}  // namespace chdr::experiments::e18e19
