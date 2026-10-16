#ifndef IPPL_FDTD_H
#define IPPL_FDTD_H
#include <cstddef>
using std::size_t;
#include "Types/Vector.h"

#include "Utility/ViewUtils.h"

#include "FieldLayout/FieldLayout.h"
#include "MaxwellSolvers/AbsorbingBC.h"
#include "MaxwellSolvers/Maxwell.h"
#include "MaxwellSolvers/OpenFaceClosure.h"
#include "Meshes/UniformCartesian.h"
#include "Particle/ParticleBase.h"

namespace ippl {
    enum fdtd_bc {
        periodic,
        absorbing
    };

    /**
     * @class FDTDSolverBase
     * @brief Base class for FDTD solvers in the ippl library.
     *
     * @tparam EMField The type representing the electromagnetic field.
     * @tparam SourceField The type representing the source field.
     * @tparam boundary_conditions The boundary conditions to be applied (default is periodic).
     */
    template <typename EMField, typename SourceField, fdtd_bc boundary_conditions>
    class FDTDSolverBase : public Maxwell<EMField, SourceField> {
    public:
        constexpr static unsigned Dim = EMField::dim;
        using scalar                  = typename EMField::value_type::value_type;
        using Vector_t                = Vector<typename EMField::value_type::value_type, Dim>;
        using SourceVector_t          = typename SourceField::value_type;

        /**
         * @brief Constructor for the FDTDSolverBase class.
         *
         * @param source Reference to the source field.
         * @param E Reference to the electric field.
         * @param B Reference to the magnetic field.
         */
        FDTDSolverBase(SourceField& source, EMField& E, EMField& B);

        /**
         * @brief Solves the FDTD equations.
         */
        void solve() override;

        /**
         * @brief Sets periodic boundary conditions.
         */
        void setPeriodicBoundaryConditions();

        /** @brief Complete the ghosts after loading the two stored potential levels.
         * @details Exchange halos of A_n and A_nm1, then apply their periodic field conditions
         * or the absorbing-face extension. Owned samples and A_np1 remain unchanged; no Mur
         * update or time shift is performed. Call before advancing externally loaded histories.
         */
        void prepareHistories();

        /**
         * @brief Gets the time step size.
         *
         * @return The time step size.
         */
        scalar getDt() const { return dt; }

        SourceField A_n;    // Current field.
        SourceField A_np1;  // Field at the next time step.
        SourceField A_nm1;  // Field at the previous time step.

        /**
         * @brief Shifts the saved fields in time.
         */
        void timeShift();

        /**
         * @brief Steps the solver forward in time. This is a pure virtual function.
         */
        virtual void step() = 0;
        /**
         * @brief Evaluates the electric and magnetic fields.
         */
        void evaluate_EB();

        /**
         * @brief Initializes the solver. This is a pure virtual function.
         */
        virtual void initialize() = 0;

    protected:
        /**
         * @brief Complete the potential update and publish its final owned values.
         *
         * Absorbing faces, edges and corners must finish before halo exchange:
         * otherwise neighboring ranks retain the previous boundary values when
         * timeShift() copies the histories. Mur reads the old histories and
         * locally available new inward values; it needs no preceding exchange
         * when each local extent is at least two cells, as required by HaloCells.
         * This is the owner-to-ghost ordering described by Oskooi et al. (2010),
         * "Meep: A flexible free-software package for electromagnetic simulations
         * by the FDTD method", section 2.3.
         *
         * The periodic path retains its existing exchange and field-BC order.
         * Derived potential updates must call this once after their interior
         * kernel has completed and before shifting the histories.
         */
        virtual void finishPotentialStep();

        /**
         * @brief Applies the boundary conditions.
         */
        void applyBCs();

        typename SourceField::Mesh_t* mesh_mp;  // Pointer to the mesh.
        FieldLayout<Dim>* layout_mp;            // Pointer to the layout
        NDIndex<Dim> domain_m;                  // Domain of the mesh.
        Vector_t hr_m;                          // Mesh spacing.

        Vector<int, Dim> nr_m;  // Number of grid points in each direction.
        scalar dt;              // Time step size.
    };
}  // namespace ippl

#include "FDTDSolverBase.hpp"
#endif
