/**
 * @file NonStandardFDTDSolver.hpp
 * @brief Impelmentation of the NonStandardFDTDSolver class functions.
 */

namespace ippl {

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    NonStandardFDTDSolver<EMField, SourceField, boundary_conditions>::NonStandardFDTDSolver(
        SourceField& source, EMField& E, EMField& B)
        requires(boundary_conditions == periodic || boundary_conditions == absorbing)
        : FDTDSolverBase<EMField, SourceField, boundary_conditions>(source, E, B) {
        validateMeshSpacing();
        initialize();
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    NonStandardFDTDSolver<EMField, SourceField, boundary_conditions>::NonStandardFDTDSolver(
        SourceField& source, EMField& E, EMField& B, const FDTDBoundary& boundary)
        requires(boundary_conditions == configured)
        : FDTDSolverBase<EMField, SourceField, boundary_conditions>(source, E, B, boundary) {
        validateMeshSpacing();
        initialize();
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void NonStandardFDTDSolver<EMField, SourceField, boundary_conditions>::validateMeshSpacing()
        const {
        auto hx = this->JN_mp->get_mesh().getMeshSpacing();
        if ((hx[2] / hx[0]) * (hx[2] / hx[0]) + (hx[2] / hx[1]) * (hx[2] / hx[1]) >= 1) {
            throw IpplException("NonStandardFDTDSolver Constructor",
                                "Dispersion-free CFL condition not satisfiable\n");
        }
    }

    /**
     * @brief Advances the simulation by one time step.
     *
     * This function updates the fields by one time step using non-standard FDTD method.
     * It calculates the new field values based on the current and previous field values,
     * as well as the source field. The boundary conditions are applied after the update.
     */
    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void NonStandardFDTDSolver<EMField, SourceField, boundary_conditions>::step() {
        const auto& ldom       = this->layout_mp->getLocalNDIndex();
        const int nghost       = this->A_n.getNghost();
        const auto aview       = this->A_n.getView();
        const auto anp1view    = this->A_np1.getView();
        const auto anm1view    = this->A_nm1.getView();
        const auto source_view = Maxwell<EMField, SourceField>::JN_mp->getView();

        const scalar calA =
            0.25
            * (1
               + 0.02
                     / (static_cast<scalar>(Kokkos::pow(this->hr_m[2] / this->hr_m[0], 2))
                        + static_cast<scalar>(Kokkos::pow(this->hr_m[2] / this->hr_m[1], 2))));
        nondispersive<scalar> ndisp{
            .a1 =
                2
                * (1
                   - (1 - 2 * calA) * static_cast<scalar>(Kokkos::pow(this->dt / this->hr_m[0], 2))
                   - (1 - 2 * calA) * static_cast<scalar>(Kokkos::pow(this->dt / this->hr_m[1], 2))
                   - static_cast<scalar>(Kokkos::pow(this->dt / this->hr_m[2], 2))),
            .a2 = static_cast<scalar>(Kokkos::pow(this->dt / this->hr_m[0], 2)),
            .a4 = static_cast<scalar>(Kokkos::pow(this->dt / this->hr_m[1], 2)),
            .a6 = static_cast<scalar>(Kokkos::pow(this->dt / this->hr_m[2], 2))
                  - 2 * calA * static_cast<scalar>(Kokkos::pow(this->dt / this->hr_m[0], 2))
                  - 2 * calA * static_cast<scalar>(Kokkos::pow(this->dt / this->hr_m[1], 2)),
            .a8 = static_cast<scalar>(Kokkos::pow(this->dt, 2))};
        const InteriorBounds interiorBounds(this->domain_m, this->boundary_m);
        Kokkos::parallel_for(
            ippl::getRangePolicy(aview, nghost), KOKKOS_LAMBDA(size_t i, size_t j, size_t k) {
                // global indices
                uint32_t ig = i + ldom.first()[0];
                uint32_t jg = j + ldom.first()[1];
                uint32_t kg = k + ldom.first()[2];
                // Mur owns only the selected outermost layers; periodic axes stay interior.
                if (interiorBounds.contains(ig, jg, kg)) {
                    anp1view(i, j, k) =
                        -anm1view(i, j, k) + ndisp.a1 * aview(i, j, k)
                        + ndisp.a2
                              * (calA * aview(i + 1, j, k - 1) + (1 - 2 * calA) * aview(i + 1, j, k)
                                 + calA * aview(i + 1, j, k + 1))
                        + ndisp.a2
                              * (calA * aview(i - 1, j, k - 1) + (1 - 2 * calA) * aview(i - 1, j, k)
                                 + calA * aview(i - 1, j, k + 1))
                        + ndisp.a4
                              * (calA * aview(i, j + 1, k - 1) + (1 - 2 * calA) * aview(i, j + 1, k)
                                 + calA * aview(i, j + 1, k + 1))
                        + ndisp.a4
                              * (calA * aview(i, j - 1, k - 1) + (1 - 2 * calA) * aview(i, j - 1, k)
                                 + calA * aview(i, j - 1, k + 1))
                        + ndisp.a6 * aview(i, j, k + 1) + ndisp.a6 * aview(i, j, k - 1)
                        + ndisp.a8 * source_view(i, j, k);
                }
            });
        Kokkos::fence();
        this->finishPotentialStep();
    }

    /**
     * @brief Initializes the solver.
     *
     * This function initializes the solver by setting up the layout, mesh, and fields.
     * It retrieves the mesh spacing, domain, and mesh size, and initializes the fields to zero.
     */
    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void NonStandardFDTDSolver<EMField, SourceField, boundary_conditions>::initialize() {
        if constexpr (boundary_conditions != periodic) {
            this->validateBoundaryFields(this->boundary_m);
        }
        // get layout and mesh
        this->layout_mp = &(this->JN_mp->getLayout());
        this->mesh_mp   = &(this->JN_mp->get_mesh());

        // get mesh spacing, domain, and mesh size
        this->hr_m     = this->mesh_mp->getMeshSpacing();
        this->dt       = this->hr_m[2];
        this->domain_m = this->layout_mp->getDomain();
        for (unsigned int i = 0; i < Dim; ++i)
            this->nr_m[i] = this->domain_m[i].length();

        // initialize fields
        this->A_nm1.initialize(*this->mesh_mp, *this->layout_mp);
        this->A_n.initialize(*this->mesh_mp, *this->layout_mp);
        this->A_np1.initialize(*this->mesh_mp, *this->layout_mp);

        this->A_nm1 = 0.0;
        this->A_n   = 0.0;
        this->A_np1 = 0.0;

        // set the mesh domain
        this->configurePotentialBoundaryConditions();
    };
}  // namespace ippl
