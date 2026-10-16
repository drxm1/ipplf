#ifndef IPPL_FDTD_H
#define IPPL_FDTD_H
#include <cstddef>
using std::size_t;
#include "Types/Vector.h"

#include "Utility/ViewUtils.h"

#include "FieldLayout/FieldLayout.h"
#include "MaxwellSolvers/AbsorbingBC.h"
#include "MaxwellSolvers/FDTDBoundary.h"
#include "MaxwellSolvers/Maxwell.h"
#include "MaxwellSolvers/OpenFaceClosure.h"
#include "Meshes/UniformCartesian.h"
#include "Particle/ParticleBase.h"

namespace ippl {
    /** @brief Construction-time boundary policy for the potential solvers.
     * periodic and absorbing select one kind on all faces. configured requires
     * an explicit FDTDBoundary constructor argument.
     */
    enum fdtd_bc {
        periodic,
        absorbing,
        configured
    };

    /**
     * @class FDTDSolverBase
     * @brief Base class for FDTD solvers in the ippl library.
     *
     * @tparam EMField The type representing the electromagnetic field.
     * @tparam SourceField The type representing the source field.
     * @tparam boundary_conditions Fixed all-face policy, or configured for an explicit description.
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
        FDTDSolverBase(SourceField& source, EMField& E, EMField& B)
            requires(boundary_conditions == periodic || boundary_conditions == absorbing);

        /** @brief Construct with an explicit, immutable face description.
         * @details Available only for configured. All ranks must supply the same
         * description. Validate it against the source layout and source/E/B ghost
         * widths before potential allocation or communication. The constructor
         * does not derive face kinds. The description must mark both faces of each
         * periodic layout axis periodic and both faces of each nonperiodic axis Mur;
         * otherwise construction throws. Only the Mur speeds are a free choice.
         * FDTDBoundary::fromLayout() builds such a description with one Mur speed.
         * @param source Source field supplying the layout and mesh.
         * @param E Electric field whose storage is validated.
         * @param B Magnetic field whose storage is validated.
         * @param boundary Face kinds and kernel-unit Mur speeds, copied into the solver.
         */
        FDTDSolverBase(SourceField& source, EMField& E, EMField& B, const FDTDBoundary& boundary)
            requires(boundary_conditions == configured);

        /**
         * @brief Solves the FDTD equations.
         */
        void solve() override;

        /**
         * @brief Sets periodic boundary conditions.
         */
        void setPeriodicBoundaryConditions();

        /// Return the construction-time description; initialize() preserves it.
        const FDTDBoundary& boundary() const { return boundary_m; }

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
         * The runtime open-face mask selects the path. With no open faces, retain
         * the existing exchange and field-BC order. Otherwise apply Mur, exchange
         * its final owned values, then extend the physical ghosts on open faces.
         * The last operation supplies the closure documented by fillOpenFaceGhosts;
         * periodic wraps must finish before that closure completes mixed corners.
         * Derived potential updates must call this once after their interior
         * kernel has completed and before shifting the histories.
         */
        virtual void finishPotentialStep();

        /**
         * @brief Applies the boundary conditions.
         */
        void applyBCs();

        /** @brief Validate the source layout and the source, E and B ghost-layer widths.
         * @details Pass the layout of JN_mp and the getNghost() values of JN_mp, En_mp
         * and Bn_mp to validateBoundary(), without MPI communication. These are the
         * fields supplied to the solver constructor; potential histories are checked separately.
         * @param boundary Face kinds and Mur speeds to check against the source layout.
         */
        void validateBoundaryFields(const FDTDBoundary& boundary) const;

        /// Restore the selected potential field conditions after initialization.
        void configurePotentialBoundaryConditions();

        const FDTDBoundary boundary_m;

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
