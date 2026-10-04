#include "VacuumSetup.h"

#include <algorithm>
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

    /** @brief Require at least two owned cells per axis on every rank.
     * @param layout Initialized, consistent 3D layout on an MPI intracommunicator.
     * @pre All ranks in layout.comm call this in matching collective order.
     * @throws std::runtime_error If any rank has fewer than two cells on an axis.
     * @details Checks extents for the default one-cell halo, not arbitrary stencils.
     * Any grid size/split axes; see HaloCells::exchangeBoundaries in src/Field/HaloCells.hpp.
     */
    void checkLocalDomain(const ippl::FieldLayout<Dim>& layout) {
        const auto& local    = layout.getLocalNDIndex();
        const int localValid = std::all_of(local.begin(), local.end(), [](const auto& axis) {
            return axis.length() >= 2;
        });
        int globalValid      = 0;
        MPI_Allreduce(&localValid, &globalValid, 1, MPI_INT, MPI_MIN, layout.comm);
        if (!globalValid)
            throw std::runtime_error(
                "Every rank must own at least two cells per local axis, required by IPPL halo "
                "exchange");
    }

    /** @brief Allocate E/B/J storage and zero the charge/current source.
     * @param fields Container owning the mesh, layout and resulting field storage.
     * @pre Fresh setup after checkLocalDomain; keep fields at this address afterward.
     * @details FELFieldContainer::initializeFields calls Field::initialize(mesh, layout)
     * without a halo argument, using its nghost=1 default (src/Field/Field.h).
     * Potential histories belong to the solver; evaluate_EB initializes E and B.
     */
    void initializeFieldStorage(Fields& fields) {
        fields.initializeFields();
        // Scalar assignment zeros every source component in owned and halo cells.
        fields.getJ() = 0.0;
    }

    /** @brief Finish pending writes and synchronize one potential's periodic halos.
     * @param potential Allocated field with one-cell halos and periodic layout/faces.
     * @pre All layout ranks call in matching order; IPPL/Kokkos are initialized.
     * @details Debug builds check local configuration; the fence completes kernel work.
     */
    void preparePotentialHalo(SourceField_t<double, Dim>& potential) {
        Kokkos::fence();
#ifndef NDEBUG
        if (!potential.getView().data())
            throw std::logic_error("Potential storage must be allocated before halo preparation");
        if (potential.getNghost() != 1 || !potential.getLayout().isAllPeriodic_m)
            throw std::logic_error("Expected one-cell halos on a periodic layout");
        for (unsigned face = 0; face < 2 * Dim; ++face) {
            const auto& bc = potential.getFieldBC()[face];
            if (!bc || bc->getBCType() != ippl::PERIODIC_FACE || bc->getFace() != face)
                throw std::logic_error("Expected a matching periodic condition on every face");
        }
#endif
        potential.fillHalo();
        potential.getFieldBC().apply(potential);
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
        preparePotentialHalo(solver.A_n);
        preparePotentialHalo(solver.A_nm1);
    }

}  // namespace chdr::vacuum
