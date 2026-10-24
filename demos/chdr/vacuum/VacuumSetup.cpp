#include "VacuumSetup.h"

#include <array>
#include <stdexcept>

#include "VacuumWave.h"

namespace chdr::vacuum {
    /** @brief Construct the normalized periodic grid with MPI decomposition along z.
     * @return Owning mesh/layout container; field arrays are not yet allocated.
     * @pre IPPL is initialized; all ranks call this function.
     * @param cellsZ Accepted z resolution: 64 (coarse) or 128 (fine).
     * @details Cells are (cellsZ/8, cellsZ/8, cellsZ) in the same normalized box.
     * The box is [0,1/8] x [0,1/8] x [0,1]; other resolutions are rejected.
     * initializeFieldStorage later allocates arrays with default one-cell halos.
     */
    Fields makePeriodicFieldContainer(int cellsZ) {
        if (cellsZ != 64 && cellsZ != 128)
            throw std::invalid_argument("Supported z cell counts are 64 and 128");
        const std::array<int, Dim> cells{cellsZ / 8, cellsZ / 8, cellsZ};
        ippl::Vector<double, Dim> lower(0.0), upper(0.125, 0.125, 1.0), spacing;
        ippl::NDIndex<Dim> domain;
        for (unsigned d = 0; d < Dim; ++d) {
            domain[d]  = ippl::Index(cells[d]);
            spacing[d] = (upper[d] - lower[d]) / cells[d];
        }
        return Fields(spacing, lower, upper, {false, false, true}, domain, lower, true);
    }

    /** @brief Initialize both histories of the chosen wave in one owned-cell kernel.
     * @param solver Fresh periodic Standard solver on the agreed unit-length z box.
     * @pre All ranks participate; both histories share the same mesh and layout.
     * @details Fallahi, MITHRA 2.0, "FDTD for Wave Equation", Eqs. (3.14)-(3.15),
     * requires current and previous histories. We sample our exact wave at 0 and -dt;
     * referencePotential defines the normalized units.
     */
    void initializePotentials(VacuumSolver& solver) {
        const auto current    = solver.A_n.getView();
        const auto previous   = solver.A_nm1.getView();
        const auto& mesh      = solver.A_n.get_mesh();
        const double originZ  = mesh.getOrigin()[2];
        const double spacingZ = mesh.getMeshSpacing()[2];
        const double dt       = solver.getDt();
        const int firstZ      = solver.A_n.getLayout().getLocalNDIndex()[2].first();
        const int ghost       = solver.A_n.getNghost();
        Kokkos::parallel_for(
            "Vacuum-wave histories", solver.A_n.getFieldRangePolicy(),
            KOKKOS_LAMBDA(const int i, const int j, const int k) {
                const double z    = originZ + (k + firstZ - ghost + 0.5) * spacingZ;
                current(i, j, k)  = chdr::vacuum::referencePotential(z, 0.0);
                previous(i, j, k) = chdr::vacuum::referencePotential(z, -dt);
            });
        solver.prepareHistories();
    }

}  // namespace chdr::vacuum
