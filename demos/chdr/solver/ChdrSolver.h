/** @file ChdrSolver.h
 * @brief Periodic demo-level shell; wave updates remain in IPPL's base solvers.
 *
 * WHY ANOTHER SOLVER CLASS? WHY DO WE INHERIT FROM Standard / NonStandard FDTD SOLVERS?
 * This is done here because the regular FDTDSolverBase calls the following algorithm:
 * @code{.cpp}
 * void solve(...) {
 *     step();
 *     timeShift();
 *     evaluate_EB();
 * }
 * @endcode
 *
 * and evaluate_EB() is not virtual. We need custom behavior for the E/B field reconstruction,
 * basically for what evaluate_EB() does, since we also want to benchmark different such
 * reconstructions and therefore we override the virtual solve() function here to inject the
 * reconstruction in place of evaluate_EB().
 *
 * Furthermore we track the amount of steps that was already done for convenience.
 *
 * We also enable overriding the time step Delta_t to study the effect of different CFL S numbers.
 */
#pragma once

// Solver internals use the communication declarations provided by the IPPL umbrella.
#include "Ippl.h"

#include <type_traits>

#include "FELFieldContainer.hpp"
#include "LevelContract.h"
#include "MaxwellSolvers/NonStandardFDTDSolver.h"
#include "Operators.h"
#include "Reconstruction.h"
#include "SpatialDimension.h"
#include "StaggeredReconstruction.h"
#include "VacuumSetup.h"

static_assert(chdr::solver::SpatialDim == 3,
              "ChdrSolver.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {

    /// Caller-owned mesh, layout, source and output storage container.
    using ChDRFieldContainer = FELFieldContainer<double, SpatialDim>;

    /// Three-component electric or magnetic field storage type on an IPPL mesh.
    using EMFieldsType = VField_t<double, SpatialDim>;

    /// Four-component scalar/vector potential or charge/current storage type.
    /// Component zero stores \f$\phi\f$ or \f$\rho\f$; components 1 through SpatialDim
    /// store the Cartesian components of \f$\mathbf A\f$ or \f$\mathbf J\f$.
    using AJFourFieldsType = SourceField_t<double, SpatialDim>;

    /// Periodic Standard FDTD Solver, like in the vacuum example.
    /// Later on we will enable different boundary conditions,
    /// maybe Mur or PML_ADI.
    using StandardSolverPeriodic =
        ippl::StandardFDTDSolver<EMFieldsType, AJFourFieldsType, ippl::periodic>;

    /// Periodic longitudinally(===z)-smoothed NonStandard vacuum update.
    /// Later on we will enable different boundary conditions,
    /// maybe Mur or PML_ADI.
    using NonStandardSolverPeriodic =
        ippl::NonStandardFDTDSolver<EMFieldsType, AJFourFieldsType, ippl::periodic>;

    /** @brief Own the solve count and reconstruction choice, preserving base slot storage.
     * @tparam PotentialSolverBase StandardSolverPeriodic or NonStandardSolverPeriodic;
     * supplies the potential wave update, three history slots and caller-owned E/B bindings.
     * @pre Source, E and B are allocated on the same all-periodic, one-halo layout.
     * Every rank must own at least two cells per axis; construction checks this
     * collectively before the base solver exchanges potential halos.
     * @details In normalized units the component equations are
     * \f$(\partial_t^2-\nabla^2)\psi=s\f$, with \f$s=\rho\f$ for the scalar
     * potential and \f$s=J_d\f$ for each vector component. Thus
     * \f$\psi[k+1]=2\psi[k]-\psi[k-1]+\Delta t^2(L_h\psi[k]+s[k])\f$
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.6)-(3.7);
     * christlieb2024GaugeConserving, Eqs. (25)-(26)]. The printed MITHRA
     * Eqs. (3.15), (3.26) use a plus sign with its negatively defined source
     * \f$\zeta\f$; IPPL instead uses the positive source above. In Eq. (3.27),
     * the repeated x spacing in \f$\alpha'_6\f$ must include the y spacing;
     * IPPL's NonStandard update uses that corrected coefficient.
     * Halo communication remains in step(); reconstruction writes owned E/B cells.
     * Callers synchronize E/B halos before differentiating them.
     * The mesh, layout and source/output fields remain owned by the caller and must
     * keep their addresses for this solver's lifetime. Potential histories remain
     * the public three-slot storage supplied by the base. Selecting a reading changes
     * the samples and reconstruction, not the vacuum wave-update kernel.
     */
    template <class PotentialSolverBase>
    class ChdrSolver : public PotentialSolverBase {
        static_assert(std::is_same_v<PotentialSolverBase, StandardSolverPeriodic>
                      || std::is_same_v<PotentialSolverBase, NonStandardSolverPeriodic>);

    public:
        /** @brief Select the time reconstruction explicitly; compact reading requires Standard R2.
         * @param source Caller-owned charge/current field: component zero is charge density;
         * components 1 through SpatialDim are Cartesian current density components.
         * @param electric Caller-owned electric field (E) output, with Cartesian components 0
         * through SpatialDim-1 on the source mesh/layout.
         * @param magnetic Caller-owned magnetic field (B) output, with the same component ordering
         * and mesh/layout as electric.
         * @param ebReconstructionConvention Explicit temporal reading; use R2 for matched sourced
         * fields.
         * @param spatialDiscretization Spatial component placement, collocated unless
         * explicitly staggered.
         */
        ChdrSolver(
            AJFourFieldsType& source, EMFieldsType& electric, EMFieldsType& magnetic,
            const EBReconstructionConvention ebReconstructionConvention,
            const SpatialDiscretization spatialDiscretization = SpatialDiscretization::Collocated)
            : PotentialSolverBase(checkedSource(source, electric, magnetic), electric, magnetic)
            , contract_m{ebReconstructionConvention, spatialDiscretization} {
            // The staggered SpatialDiscretization currently works with R2 and the standard stencil
            // only.
            if (spatialDiscretization == SpatialDiscretization::Staggered
                && (ebReconstructionConvention != EBReconstructionConvention::R2
                    || stencil() != Stencil::Standard)) {
                throw IpplException("ChdrSolver",
                                    "Staggered reading requires R2 and the Standard stencil");
            }

            // Check that their layouts are periodic in every direction,
            // and that they have exactly one ghost-cell layer.
            requirePeriodic(source);
            requirePeriodic(electric);
            requirePeriodic(magnetic);
        }

        /** @brief Advance the potentials, relabel the slots, then reconstruct this solve's output.
         * @details The caller writes the source for the next step before this call.
         * Incrementing the count after timeShift keeps the slot and output labels in
         * LevelContract synchronized. No source history is implicitly resampled here.
         */
        void solve() override {
            this->step();
            this->timeShift();  // TODO: at a later point check if rotation wouldn't be more
                                // performant than the current copy impl
            ++completed_m;      // Keep track of the number of completed steps
            reconstructEB();    // generalization of evaluate_EB()
        }

        /** @brief Reset base storage and solve count while retaining an explicit step override.
         * @details Call before loading histories: the base initializer zeros its slots.
         * A previously set explicit timestep survives reinitialization so the loaded time labels
         * continue to use the step selected through this shell.
         * E/B are not reconstructed here; call evaluateInitialFields after loading
         * both potential histories before reading either output array.
         */
        void initialize() override {
            // Regular initialization as in the Base (Standard | NonStandard) solver
            PotentialSolverBase::initialize();

            // We allow different dt for studying their effect
            if (hasStepOverride_m) {
                this->dt = overriddenDt_m;
            }
            // We also keep track of the number of completed steps
            completed_m = 0;

            // Initialization zeroed the potential slots. Allow a new timestep to be selected
            // before loading histories, reconstructing E/B or advancing the solver.
            timeStepLocked_m = false;
            // Start a new history generation so monitors and probes can detect this reset.
            ++historyRevision_m;
        }

        /** @brief Reconstruct the loaded histories without stepping or changing their labels.
         * @pre Both potential slots have current periodic halos on their shared layout.
         * @details A_n is the current four-component potential slot; A_nm1 is its
         * predecessor, one time step earlier for each component. Component zero is the
         * scalar potential and components 1 through SpatialDim are the vector potential.
         * LevelContract specifies the scalar/vector time offsets for each convention.
         * The reconstructed fields follow
         * \f$\mathbf E=-\partial_t\mathbf A-\nabla\phi\f$ and
         * \f$\mathbf B=\nabla\times\mathbf A\f$
         * [fallahi2020mithra20fullwavesimulation, Eqs. (3.8)-(3.9)]; the selected
         * convention determines the discrete time differences and magnetic averaging.
         * @details loadHistory supplies those halos for analytic initialization; a normal
         * solve uses the base step's exchanges. R1 retains IPPL's shipped reconstruction;
         * R2 changes the scalar slot and averages B. R2Prime exercises the base electric
         * kernel with consistently shifted scalar samples (LevelContract sources).
         *
         * TODO: One might later check thoroughly whether it would make sense
         * to combine the currently sepparate Kokkos for loops for computing E in
         * reconstructElectric and computing B in reconstructMagnetic into one Kokkos loop. That
         * might maybe slightly improve performance.
         */
        void reconstructEB() {
            // Freeze the timestep even when reconstructing initial fields before the first solve.
            // TODO: Not ideal to have to call this at every step, maybe improve this at some point
            // in the future.
            timeStepLocked_m = true;
            if (spatialDiscretization() == SpatialDiscretization::Staggered) {
                reconstructStaggered(this->A_n, this->A_nm1, electricField(), magneticField(),
                                     this->getDt());
                return;
            }
            if (ebReconstructionConvention() == EBReconstructionConvention::R2)
                reconstructElectric(this->A_n, this->A_nm1, electricField(), this->getDt());
            else
                PotentialSolverBase::evaluate_EB();
            // R2Prime keeps the base E output but replaces its unaveraged B below.
            // Calling the combined base kernel therefore computes a B result that is discarded
            // (TODO: this is not optimal performance-wise, improve that later).
            // For R2 it is already better: first reconstructElectric is called to compute E then
            // reconstructMagnetic to compute B.
            if (ebReconstructionConvention() != EBReconstructionConvention::R1)
                reconstructMagnetic(this->A_n, this->A_nm1, magneticField(), true);
        }

        /** @brief Produce the unaveraged magnetic channel at the current potential time.
         * @param output Independent caller-owned destination for the Cartesian components
         * of \f$\mathbf B_{\mathrm{raw}}=\nabla_h\times\mathbf A_n\f$, using the
         * selected spatial curl [fallahi2020mithra20fullwavesimulation, Eq. (3.8)].
         * @pre Output is allocated on the same mesh/layout; current potential halos are valid.
         * @details This explicit destination keeps raw B separate from the output B used
         * by R2's common half-step observables [fallahi2020mithra20fullwavesimulation,
         * section 3.2.2, Eqs. (3.50)-(3.52)]. The current spatial reading chooses
         * its curl and component positions. Aliasing either bound output is refused.
         */
        void computeRawMagneticField(EMFieldsType& output) {
            if (output.getView().data() == electricField().getView().data()
                || output.getView().data() == magneticField().getView().data())
                throw IpplException("ChdrSolver::computeRawMagneticField",
                                    "Raw magnetic output must not alias E or B");
            requireMatching(output);
            if (spatialDiscretization() == SpatialDiscretization::Staggered)
                rawStaggeredMagnetic(this->A_n, output);
            else
                reconstructMagnetic(this->A_n, this->A_nm1, output, false);
        }
        /// @brief Reconstruct the loaded histories without advancing the solve count.
        void evaluateInitialFields() { reconstructEB(); }
        /// Preserve the shell's selected reconstruction when using the familiar base entry point.
        void evaluate_EB() { reconstructEB(); }
        /// @return Number of completed solve calls since initialization.
        long long completedSolves() const { return completed_m; }
        /// @return History generation, independent of the resettable completed-solve count.
        unsigned long long historyRevision() const { return historyRevision_m; }
        /// @return Selected spatial component placement.
        SpatialDiscretization spatialDiscretization() const {
            return contract_m.spatialDiscretization_m;
        }
        /// @return Selected temporal reconstruction.
        EBReconstructionConvention ebReconstructionConvention() const {
            return contract_m.ebReconstructionConvention_m;
        }
        /// @return Immutable shared position/time table for the selected reading.
        const LevelContract& levelContract() const { return contract_m; }
        /// Borrow the caller-owned output/source arrays; access does not exchange halos.
        EMFieldsType& electricField() { return *this->En_mp; }
        /// @return Borrowed caller-owned output B; its halos are not exchanged here.
        EMFieldsType& magneticField() { return *this->Bn_mp; }
        /// @return Borrowed charge/current array written before the next solve.
        AJFourFieldsType& sourceField() { return *this->JN_mp; }

        /// Bound field ownership is fixed at construction; rebuild the shell to change it.
        void setSources(AJFourFieldsType&) override {
            throw IpplException("ChdrSolver::setSources", "Construct a new shell to rebind fields");
        }
        /// Hide the unchecked base rebinding entry point for normal shell use.
        void setEMFields(EMFieldsType&, EMFieldsType&) {
            throw IpplException("ChdrSolver::setEMFields",
                                "Construct a new shell to rebind fields");
        }

        /// The class static_assert deliberately rejects unknown or absorbing base families.
        static constexpr Stencil stencil() {
            return std::is_same_v<PotentialSolverBase, StandardSolverPeriodic>
                       ? Stencil::Standard
                       : Stencil::NonStandard;
        }

    private:
        /// Validate collectively before the base constructor performs halo communication.
        static AJFourFieldsType& checkedSource(AJFourFieldsType& source,
                                               const EMFieldsType& electric,
                                               const EMFieldsType& magnetic) {
            requirePeriodic(source);
            requirePeriodic(electric);
            requirePeriodic(magnetic);
            if (&source.getLayout() != &electric.getLayout()
                || &source.getLayout() != &magnetic.getLayout()
                || &source.get_mesh() != &electric.get_mesh()
                || &source.get_mesh() != &magnetic.get_mesh())
                throw IpplException(
                    "ChdrSolver", "ChDRFieldContainer must share the same mesh and layout owners");
            chdr::vacuum::checkLocalDomain(source.getLayout());
            return source;
        }
        void requireMatching(const EMFieldsType& field) {
            requirePeriodic(field);
            if (&field.getLayout() != &sourceField().getLayout()
                || &field.get_mesh() != &sourceField().get_mesh())
                throw IpplException("ChdrSolver", "Output must share the source mesh and layout");
        }
        /// Inspect field metadata only; potential, source and E/B storage are all accepted.
        template <class FieldType>
        static void requirePeriodic(const FieldType& field) {
            if (!field.getLayout().isAllPeriodic_m || field.getNghost() != 1)
                throw IpplException("ChdrSolver",
                                    "Requires an all-periodic layout and one-cell halos");
        }
        LevelContract contract_m;
        long long completed_m                = 0;
        double overriddenDt_m                = 0;
        bool hasStepOverride_m               = false;
        bool timeStepLocked_m                = false;
        unsigned long long historyRevision_m = 0;
    };
}  // namespace chdr::solver
