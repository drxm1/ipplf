#include "VacuumSetup.h"

#include <algorithm>
#include <stdexcept>

namespace chdr::apps::vacuumWave {
    /** @brief Describe the normalized periodic grid with MPI decomposition along z.
     * @param cellsZ Accepted z resolution: 64 (coarse) or 128 (fine).
     * @details Cells are (cellsZ/8, cellsZ/8, cellsZ) in the same normalized box.
     * The box is [0,1/8] x [0,1/8] x [0,1]; other resolutions are rejected.
     * The simulation manager allocates its fields with the common one-cell halos.
     */
    simulation::ChdrCartesianGridSettings vacuumBox(int cellsZ) {
        if (cellsZ != 64 && cellsZ != 128)
            throw std::invalid_argument("Supported z cell counts are 64 and 128");
        const double spacing = 1.0 / cellsZ;
        return {{cellsZ / 8, cellsZ / 8, cellsZ},
                {spacing, spacing, spacing},
                0,
                {true, true, true},
                2};
    }

    /** @brief Preserve the original Standard timestep and collocated field reconstruction.
     * @details StandardFDTDSolver::initialize selects half the smallest mesh spacing.
     * R1 calls its existing evaluate_EB reconstruction, and an empty formulation retains
     * its potential update. The manager validates the explicit timestep against the
     * Standard bound [fallahi2020mithra20fullwavesimulation, Eqs. (3.19)-(3.20)].
     * Keeping the exact step preserves the benchmark's quarter-period checkpoints.
     */
    simulation::ChdrSolverSettings vacuumSettings(
        const simulation::ChdrCartesianGridSettings& box) {
        simulation::ChdrSolverSettings settings;
        settings.reconstruction_m = EBReconstructionConvention::R1;
        settings.spatial_m        = SpatialDiscretization::Collocated;
        settings.formulation_m    = std::nullopt;
        settings.timeStep_m = *std::min_element(box.spacing_m.begin(), box.spacing_m.end()) / 2;
        return settings;
    }

}  // namespace chdr::apps::vacuumWave
