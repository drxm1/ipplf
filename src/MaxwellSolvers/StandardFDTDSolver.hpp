/**
 * @file StandardFDTDSolver.hpp
 * @brief Impelmentation of the StandardFDTDSolver class functions.
 */

namespace ippl {

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    StandardFDTDSolver<EMField, SourceField, boundary_conditions>::StandardFDTDSolver(
        SourceField& source, EMField& E, EMField& B)
        requires(boundary_conditions == periodic || boundary_conditions == absorbing)
        : FDTDSolverBase<EMField, SourceField, boundary_conditions>(source, E, B) {
        initialize();
    }

    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    StandardFDTDSolver<EMField, SourceField, boundary_conditions>::StandardFDTDSolver(
        SourceField& source, EMField& E, EMField& B, const FDTDBoundary& boundary)
        requires(boundary_conditions == configured)
        : FDTDSolverBase<EMField, SourceField, boundary_conditions>(source, E, B, boundary) {
        initialize();
    }

    /**
     * @brief Advances the simulation by one time step.
     *
     * This function updates the fields by one time step using the FDTD method.
     * It calculates the new field values based on the current and previous field values,
     * as well as the source field. The boundary conditions are applied after the update.
     */
    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    void StandardFDTDSolver<EMField, SourceField, boundary_conditions>::step() {
        const auto& ldom       = this->layout_mp->getLocalNDIndex();
        const int nghost       = this->A_n.getNghost();
        const auto aview       = this->A_n.getView();
        const auto anp1view    = this->A_np1.getView();
        const auto anm1view    = this->A_nm1.getView();
        const auto source_view = Maxwell<EMField, SourceField>::JN_mp->getView();

        // Calculate the coefficients for the nondispersive update
        const scalar a1 = scalar(2)
                          * (scalar(1) - Kokkos::pow(this->dt / this->hr_m[0], 2)
                             - Kokkos::pow(this->dt / this->hr_m[1], 2)
                             - Kokkos::pow(this->dt / this->hr_m[2], 2));
        const scalar a2 = Kokkos::pow(this->dt / this->hr_m[0], 2);
        const scalar a4 = Kokkos::pow(this->dt / this->hr_m[1], 2);
        const scalar a6 = Kokkos::pow(this->dt / this->hr_m[2], 2);
        const scalar a8 = Kokkos::pow(this->dt, 2);
        const InteriorBounds interiorBounds(this->domain_m, this->boundary_m);

        // Update the field values
        Kokkos::parallel_for(
            "Source field update", ippl::getRangePolicy(aview, nghost),
            KOKKOS_LAMBDA(const size_t i, const size_t j, const size_t k) {
                // global indices
                const uint32_t ig = i + ldom.first()[0];
                const uint32_t jg = j + ldom.first()[1];
                const uint32_t kg = k + ldom.first()[2];
                // Mur owns only the selected outermost layers; periodic axes stay interior.
                if (interiorBounds.contains(ig, jg, kg)) {
                    SourceVector_t interior = -anm1view(i, j, k) + a1 * aview(i, j, k)
                                              + a2 * (aview(i + 1, j, k) + aview(i - 1, j, k))
                                              + a4 * (aview(i, j + 1, k) + aview(i, j - 1, k))
                                              + a6 * (aview(i, j, k + 1) + aview(i, j, k - 1))
                                              + a8 * source_view(i, j, k);
                    anp1view(i, j, k) = interior;
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
    void StandardFDTDSolver<EMField, SourceField, boundary_conditions>::initialize() {
        if constexpr (boundary_conditions != periodic) {
            this->validateBoundaryFields(this->boundary_m);
        }
        // get layout and mesh
        this->layout_mp = &(this->JN_mp->getLayout());
        this->mesh_mp   = &(this->JN_mp->get_mesh());

        // get mesh spacing, timestep, domain, and mesh size
        this->hr_m = this->mesh_mp->getMeshSpacing();
        this->dt   = this->hr_m[0] / 2;
        for (unsigned int i = 0; i < Dim; ++i) {
            this->dt = std::min(this->dt, this->hr_m[i] / 2);
        }
        this->domain_m = this->layout_mp->getDomain();
        for (unsigned int i = 0; i < Dim; ++i)
            this->nr_m[i] = this->domain_m[i].length();

        // initialize fields
        this->A_nm1.initialize(*(this->mesh_mp), *(this->layout_mp));
        this->A_n.initialize(*(this->mesh_mp), *(this->layout_mp));
        this->A_np1.initialize(*(this->mesh_mp), *(this->layout_mp));

        // Initialize fields to zero
        this->A_nm1 = 0.0;
        this->A_n   = 0.0;
        this->A_np1 = 0.0;

        // set the mesh domain
        this->configurePotentialBoundaryConditions();
    }
}  // namespace ippl
