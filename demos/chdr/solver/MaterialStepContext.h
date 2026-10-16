/** @file MaterialStepContext.h
 * @brief Access to ChdrSolver fields and operations for one MaterialUpdate call.
 */
#pragma once

#include "Ippl.h"

#include "FELFieldContainer.hpp"
#include "SpatialDimension.h"

namespace chdr::solver {
    /// Four-component scalar/vector potential or charge/current storage type.
    /// Component zero stores \f$\phi\f$ or \f$\rho\f$; components 1 through SpatialDim
    /// store the Cartesian components of \f$\mathbf A\f$ or \f$\mathbf J\f$.
    using AJFourFieldsType = SourceField_t<double, SpatialDim>;

    /** @brief Give a MaterialUpdate access to ChdrSolver's fields and step operations.
     * @details ChdrSolver owns both a concrete implementation of this context
     * interface and a MaterialUpdate object. For one step, it passes its context
     * by reference to MaterialUpdate::advance(). The context forwards requests
     * to that ChdrSolver, so the material algorithm can access existing fields
     * and operations without depending directly on the concrete solver type.
     *
     * The MaterialUpdate borrows the context for that call: it must not delete
     * it or retain a reference for later use. The context itself remains owned
     * by ChdrSolver between calls. Its accessors lend references to existing
     * fields; they do not copy fields or transfer ownership. The update must not
     * replace their mesh/layout or retain these references for a later step.
     * Potential slots belong to ChdrSolver's inherited base. The physical source
     * field belongs to the caller that supplied it and must keep its original values.
     *
     * advance() writes next() and completes its halos and boundaries.
     * ChdrSolver::solve() then rotates histories, increments the solve count
     * and reconstructs E/B.
     * An update may use several device passes and halo exchanges. These virtual
     * calls run on the host; device kernels capture only views and scalars
     * [kokkos2026VirtualFunctions; kokkos2026ParallelDispatch].
     */
    class MaterialStepContext {
    public:
        /// Borrow A_nm1, the preceding history slot; it is an input to this update.
        virtual AJFourFieldsType& previous() = 0;
        /// Borrow A_n, the current history slot; it is an input to this update.
        virtual AJFourFieldsType& current() = 0;
        /// Borrow A_np1, where this update writes the next potential samples.
        virtual AJFourFieldsType& next() = 0;
        /// Borrow physical charge and current at ChdrSolver's declared source labels.
        virtual AJFourFieldsType& source() = 0;
        /// Return the physical time spacing used for histories and reconstruction.
        virtual double dt() const = 0;
        /** @brief Call the inherited step with temporary internal parameters.
         * @details ChdrSolver restores its time spacing and source binding on every
         * exit path. The inherited step itself finishes its halos and boundaries.
         * dtInternal and sourceInternal affect only this inherited call; dt() and
         * source() describe the caller's physical problem before and after it.
         * sourceInternal must remain alive until its launched device work finishes.
         */
        virtual void baseStep(double dtInternal, AJFourFieldsType& sourceInternal) = 0;
        /** @brief Complete next() after custom kernels at the solver's physical timestep.
         * @details All-periodic layouts exchange halos, then apply field conditions. With Mur
         * faces, use dt() and the stored boundary().speed values for Mur, exchange the completed
         * potentials, then fill exterior ghosts
         * [mur1981AbsorbingBoundaryFdtd, section IV, Eq. (16)]. A custom update with a different
         * internal timestep must account for that step in its boundary speeds; finishStep() does
         * not rescale them. baseStep() already completes its
         * boundaries while its temporary internal timestep is active; do not finish it again.
         */
        virtual void finishStep() = 0;

    protected:
        ~MaterialStepContext() = default;
    };

    /** @brief Compute one potential update for a material and provide its step bound.
     * @details ChdrSolver owns this algorithm object and passes its own
     * MaterialStepContext to advance(). The context is borrowed only for that
     * call, under the lifetime and field-access contract documented above.
     * A MaterialUpdate may own scratch fields; their borrowed mesh/layout
     * owners must outlive those fields. It must leave ChdrSolver's potential
     * history labels and solve count unchanged. ChdrSolver::solve() performs history
     * rotation and E/B reconstruction after the update returns.
     */
    class MaterialUpdate {
    public:
        virtual ~MaterialUpdate() = default;
        /// Produce the next potential slot, including valid halos and boundaries.
        virtual void advance(MaterialStepContext& context) = 0;
        /// Return the formulation's stability bound for the supplied mesh spacing.
        virtual double maxStableTimeStep(const ippl::Vector<double, SpatialDim>& spacing) const = 0;
    };
}  // namespace chdr::solver
