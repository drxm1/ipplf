/** @file ChdrSolver.h
 * @brief Potential solver on top of IPPL's Standard and NonStandard FDTD solvers, with periodic
 * or Mur faces.
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

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include "MaxwellSolvers/NonStandardFDTDSolver.h"
#include "MaxwellSolvers/StandardFDTDSolver.h"
#include "core/base/FieldContainer.h"
#include "core/base/LevelContract.h"
#include "core/base/SpatialDimension.h"
#include "core/solver/formulation/MaterialFormulation.h"
#include "core/solver/formulation/MaterialStepContext.h"
#include "core/solver/operators/Operators.h"
#include "core/solver/operators/Reconstruction.h"
#include "core/solver/operators/StaggeredReconstruction.h"

static_assert(chdr::solver::SpatialDim == 3,
              "ChdrSolver.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {

    /** @brief Select the propagation speed assumed by the Mur boundary formulas.
     * @details physical() specifies a speed using the time coordinate of
     * ChdrSolver::getDt(), with vacuum propagation speed normalized to one.
     * HomogeneousUpdate calls the potential solver with the shorter timestep
     * \f$\Delta t_{\mathrm{update}}=\Delta t/n\f$, where
     * \f$n=\sqrt{\varepsilon_r}\f$; see HomogeneousUpdate and
     * [chew2014generalizedgauge, section 2, Eqs. (11)-(13)].
     * The Mur coefficients depend on speed times timestep
     * [mur1981AbsorbingBoundaryFdtd, section IV, Eq. (16)]. Therefore kernel()
     * returns \f$n v\f$ for an explicitly supplied speed \f$v\f$, so that
     * \f$(n v)\Delta t_{\mathrm{update}}=v\Delta t\f$ (deduction).
     * Without a homogeneous material update, the conversion uses \f$n=1\f$.
     *
     * matched() instead selects exactly one for the Mur formulas, equal to the
     * propagation speed in the reused potential-update equations. In the time
     * coordinate of getDt(), this corresponds to speed \f$1/n\f$. The name
     * denotes this choice in this API; it does not assert perfect absorption.
     */
    class MurSpeed {
    public:
        /// Select exactly one as the speed passed to the Mur coefficient constructors.
        static MurSpeed matched() { return MurSpeed(); }

        /** @brief Specify the boundary propagation speed before timestep rescaling.
         * @param speed Finite, positive speed in the time coordinate of ChdrSolver::getDt();
         * vacuum propagation speed is one. kernel() converts it for the potential update.
         */
        static MurSpeed physical(double speed) {
            if (!std::isfinite(speed) || speed <= 0)
                throw std::invalid_argument("Mur speed must be finite and positive");
            return MurSpeed(speed);
        }
        /** @brief Return the speed to pass to the Mur coefficient constructors.
         * @param refractiveIndex Ratio of ChdrSolver::getDt() to the timestep used
         * by the potential update: the homogeneous medium's refractive index, or one
         * without that rescaling. The caller supplies the validated material value.
         * @return Exactly one for matched(); otherwise the supplied speed multiplied
         * by refractiveIndex. The private zero sentinel selects the first case;
         * it is never passed to Mur as a zero propagation speed.
         */
        double kernel(double refractiveIndex) const {
            return physical_m == 0 ? 1.0 : physical_m * refractiveIndex;
        }

    private:
        explicit MurSpeed(double physical = 0)
            : physical_m(physical) {}
        double physical_m;
    };

    /// Configured Standard base; retain the alias name for existing ChdrSolver callers.
    using StandardSolverPeriodic =
        ippl::StandardFDTDSolver<EMFieldsType, AJFourFieldsType, ippl::configured>;

    /// Longitudinally (z) smoothed NonStandard update; ChdrSolver permits periodic faces only.
    using NonStandardSolverPeriodic =
        ippl::NonStandardFDTDSolver<EMFieldsType, AJFourFieldsType, ippl::configured>;

    /** @brief Own the solve count and reconstruction choice, preserving base slot storage.
     * @tparam PotentialSolverBase StandardSolverPeriodic or NonStandardSolverPeriodic;
     * supplies the potential wave update, three history slots and caller-owned E/B bindings.
     * @pre Source, E and B are allocated on the same one-halo layout.
     * Open faces use Mur conditions; NonStandard supports periodic layouts only.
     * Every rank must own at least two cells per axis. Each rank checks the stored
     * partition metadata before the base solver exchanges potential halos.
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
         * @param murSpeed Speed assumed by the Mur formulas. The default passes exactly
         * one, equal to the vacuum potential solver's propagation speed. An explicit
         * physical() speed is unchanged here because no material timestep rescaling occurs.
         */
        ChdrSolver(
            AJFourFieldsType& source, EMFieldsType& electric, EMFieldsType& magnetic,
            const EBReconstructionConvention ebReconstructionConvention,
            const SpatialDiscretization spatialDiscretization = SpatialDiscretization::Collocated,
            MurSpeed murSpeed                                 = MurSpeed::matched())
            : ChdrSolver(source, electric, magnetic, ebReconstructionConvention,
                         spatialDiscretization,
                         checkedBoundary(source, electric, magnetic, murSpeed.kernel(1.0))) {}

        /** @brief Bind a material formulation on the source mesh and layout.
         * @details Mesh and layout owners must outlive this ChdrSolver; the material
         * owner itself is borrowed only during construction. Validity uses its cached
         * statistics before the delegated constructor can enter a collective.
         * The source remains physical charge/current; the formulation applies its
         * factors internally. Reconstruction always uses the physical timestep.
         * The inherited default step is retained: on cubic cells it corresponds to
         * \f$S=(\sqrt3/2)/\sqrt\varepsilon\f$ for the homogeneous update
         * (deduction from [fallahi2020mithra20fullwavesimulation, Eq. (3.20)]).
         * Select setSafetyFactor explicitly to use a prescribed medium step.
         * The owner must use the exact source mesh/layout objects, not merely an
         * equivalent grid. The delegating argument validates the owner before reading
         * its permittivity, converts the speed and validates the final face description.
         * Base construction receives that description; workspace allocation follows.
         * HomogeneousUpdate copies the uniform permittivity, so it does not retain
         * the owner or a reference to its field. Source and E/B storage remain borrowed.
         * O1VSUpdate retains reference-counted coefficient views and uses the vacuum
         * timestep bound documented in O1VSUpdate::maxStableTimeStep.
         * @pre O1-VS with Mur faces requires alpha = 1 and assumes vacuum permittivity
         * in the outer layers for the standard boundary path. Face permittivity is not
         * checked; an experimental boundary override must supply its own documented
         * assumptions and speed treatment.
         * @param murSpeed Speed assumed by the Mur formulas. The default passes exactly
         * one to agree with the reused potential-update equations. For Homogeneous,
         * physical(v) passes v times the refractive index to compensate for the shorter
         * update timestep; MurSpeed documents the conversion and its sources.
         */
        ChdrSolver(AJFourFieldsType& source, EMFieldsType& electric, EMFieldsType& magnetic,
                   EBReconstructionConvention reconstruction, SpatialDiscretization spatial,
                   const material::MaterialOwner& owner, FormulationParameters parameters = {},
                   MurSpeed murSpeed = MurSpeed::matched())
            : ChdrSolver(source, electric, magnetic, reconstruction, spatial,
                         checkedMaterialBoundary(source, electric, magnetic, owner, parameters,
                                                 {stencil(), spatial, reconstruction}, murSpeed)) {
            update_m =
                makeMaterialUpdate(parameters, owner, {stencil(), spatial, reconstruction}, source);
            homogeneousPermittivity_m = owner.statistics().homogeneous()
                                            ? owner.homogeneousValue()
                                            : std::numeric_limits<double>::quiet_NaN();
            formulationParameters_m   = parameters;
        }

        /// The service adapter borrows this solver's address; bindings cannot be copied or moved.
        ChdrSolver(const ChdrSolver&)            = delete;
        ChdrSolver& operator=(const ChdrSolver&) = delete;
        ChdrSolver(ChdrSolver&&)                 = delete;
        ChdrSolver& operator=(ChdrSolver&&)      = delete;

        /** @brief Produce the next potential slot without relabelling or counting a solve.
         * @details Each selected update finishes its own halos and boundaries.
         * The owner-free path retains the inherited step unchanged.
         * Both paths write A_np1 from A_n and A_nm1. A direct step() call leaves
         * those slot names, completedSolves() and E/B outputs unchanged; solve()
         * subsequently rotates histories, increments the count and reconstructs E/B.
         */
        void step() override {
            if (update_m)
                update_m->advance(context_m);
            else
                PotentialSolverBase::step();
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
         * continue to use the step selected through this ChdrSolver.
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
         * @pre Both potential slots have current halos and closed physical ghosts.
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
        /// Preserve the reconstruction selected in this ChdrSolver when using the base entry point.
        void evaluate_EB() { reconstructEB(); }
        /// @return Number of completed solve calls since initialization.
        long long completedSolves() const { return completed_m; }
        /// @return History generation, independent of the resettable completed-solve count.
        unsigned long long historyRevision() const { return historyRevision_m; }
        /** @brief Declare a custom or analytic history load and freeze its time spacing.
         * @pre Call before writing both slots, after selecting the time step.
         * @details Reloading at count zero invalidates monitor/probe histories through
         * historyRevision. An evolved solver must be initialized before reloading.
         */
        void markHistoryLoaded() {
            if (completed_m != 0)
                throw IpplException("ChdrSolver::markHistoryLoaded",
                                    "Initialize before reloading histories");
            timeStepLocked_m = true;
            ++historyRevision_m;
        }
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
        /// Whether a material formulation is bound rather than the inherited vacuum step.
        bool hasMaterial() const { return bool(update_m); }
        /// Uniform relative permittivity, cached at construction; the owner-free
        /// vacuum route returns one. Reserved nonuniform metadata is NaN.
        double homogeneousPermittivity() const { return homogeneousPermittivity_m; }
        /// Parameters to record beside stencil, spatial placement and reconstruction.
        /// Without a material update these are defaults, not a selected formulation.
        const FormulationParameters& formulationParameters() const {
            return formulationParameters_m;
        }
        /// Stable formulation label for run records; distinguish the owner-free route.
        const char* formulationName() const {
            return hasMaterial()
                       ? chdr::solver::formulationName(formulationParameters_m.formulation_m)
                       : "Vacuum";
        }

        /** @brief Interpret an integer solve/source index through the shared label table.
         * @param quantity Output, potential or source channel whose label is requested.
         * @param completed Completed-solve count for outputs; source index for inputs.
         * @param phiIsZero Explicit declaration required for the unique R1 electric label.
         * @return Absolute channel time for the current fixed step.
         */
        double labelTime(Quantity quantity, long long completed, bool phiIsZero = false) const {
            return contract_m.labelTime(quantity, completed, this->getDt(), phiIsZero);
        }

        /// Borrow the caller-owned output/source arrays; access does not exchange halos.
        EMFieldsType& electricField() { return *this->En_mp; }
        /// @return Borrowed caller-owned output B; its halos are not exchanged here.
        EMFieldsType& magneticField() { return *this->Bn_mp; }
        /// @return Borrowed charge/current array written before the next solve.
        AJFourFieldsType& sourceField() { return *this->JN_mp; }

        /// Bound field ownership is fixed at construction; construct a new ChdrSolver to change it.
        void setSources(AJFourFieldsType&) override {
            throw IpplException("ChdrSolver::setSources",
                                "Construct a new ChdrSolver to rebind fields");
        }
        /// Hide the unchecked base rebinding entry point for normal ChdrSolver use.
        void setEMFields(EMFieldsType&, EMFieldsType&) {
            throw IpplException("ChdrSolver::setEMFields",
                                "Construct a new ChdrSolver to rebind fields");
        }

        /// The class static_assert restricts PotentialSolverBase to StandardSolverPeriodic
        /// or NonStandardSolverPeriodic; this function identifies which stencil it uses.
        static constexpr Stencil stencil() {
            return std::is_same_v<PotentialSolverBase, StandardSolverPeriodic>
                       ? Stencil::Standard
                       : Stencil::NonStandard;
        }
        /** @brief Return the selected update's stability bound for both step-selection methods.
         * @details The bound applies to the physical step returned by getDt(), not
         * HomogeneousUpdate's reduced internal step. Delegation makes safetyFactor()
         * and setSafetyFactor() use the same material-dependent normalization.
         * The owner-free route retains the stencil's vacuum bound.
         * @return Positive maximum stable time spacing for this update and mesh.
         */
        double maxStableTimeStep() const {
            if (update_m)
                return update_m->maxStableTimeStep(this->hr_m);
            return chdr::solver::maxStableTimeStep(this->hr_m, stencil());
        }
        /// @return Current physical time step divided by the selected update's stability bound.
        double safetyFactor() const { return this->getDt() / maxStableTimeStep(); }

        /** @brief Select an effective CFL factor before the first solve.
         * @param safety Positive fraction of maxStableTimeStep within the stated range.
         * @details Standard accepts \f$0<S<1\f$; NonStandard accepts
         * \f$0<S\leq1\f$. The equality exception is audited as a double-root record.
         * Sources: Operators.h; gedney2011IntroductionFdtd, section 3.5, Eq. (3.63).
         * The factor multiplies the effective bound, not the base constructor's default
         * step. Loading or reconstructing histories freezes the step until initialize().
         */
        void setSafetyFactor(double safety) {
            if (!std::isfinite(safety) || safety <= 0 || safety > 1
                || (stencil() == Stencil::Standard && safety == 1))
                throw IpplException("ChdrSolver::setSafetyFactor",
                                    "Safety factor outside allowed range");
            setTimeStepUnchecked(safety * maxStableTimeStep());
        }
        /** @brief Set a deliberate test/record step without imposing the stability interval.
         * @param step Finite positive time spacing, before loading, reconstruction or solving.
         * @details TimeStep.GrowsAboveBound uses an unstable step to check that the growth
         * assertion is sensitive. Finite positive steps are required, and the timestep lock still
         * applies. This also supports exact-period records; ordinary runs use setSafetyFactor.
         *
         * THIS IS MOSTLY USED FOR TESTS, AND SHOULD USUALLY NOT BE USED IN PRODUCTION OTHER THAN IN
         * "setSafetyFactor()".
         */
        void setTimeStepUnchecked(double step) {
            if (!std::isfinite(step) || step <= 0.0)
                throw IpplException("ChdrSolver::setTimeStepUnchecked",
                                    "Step must be finite and positive");
            if (completed_m != 0 || timeStepLocked_m)
                throw IpplException("ChdrSolver::setTimeStepUnchecked",
                                    "Initialize before changing loaded histories' time step");
            this->dt = overriddenDt_m = step;
            hasStepOverride_m         = true;
        }

    private:
        /** @brief Construct the base with validated face kinds and converted Mur speeds.
         * description stores the speeds used with the potential-update timestep,
         * not the unconverted speeds supplied through MurSpeed::physical().
         * Only the public constructors call this overload. Their argument helpers finish
         * storage, material and boundary validation before base initialization can
         * communicate. No derived member is read to select the boundary.
         */
        ChdrSolver(AJFourFieldsType& source, EMFieldsType& electric, EMFieldsType& magnetic,
                   EBReconstructionConvention reconstruction, SpatialDiscretization spatial,
                   const ippl::FDTDBoundary& description)
            : PotentialSolverBase(source, electric, magnetic, description)
            , contract_m{reconstruction, spatial} {
            if (spatial == SpatialDiscretization::Staggered
                && (reconstruction != EBReconstructionConvention::R2
                    || stencil() != Stencil::Standard))
                throw IpplException("ChdrSolver",
                                    "Staggered reading requires R2 and the Standard stencil");
        }

        /** @brief Validate material metadata, then convert the requested Mur speed.
         * Read cached metadata only; the owner need not outlive construction.
         * For Homogeneous, multiply an explicit speed by the refractive index before
         * passing it to the base solver; otherwise use a conversion factor of one.
         * MurSpeed::kernel explains how this compensates for timestep rescaling.
         * No not-yet-initialized ChdrSolver member supplies the conversion factor.
         */
        static ippl::FDTDBoundary checkedMaterialBoundary(
            AJFourFieldsType& source, const EMFieldsType& electric, const EMFieldsType& magnetic,
            const material::MaterialOwner& owner, const FormulationParameters& parameters,
            const SolverDiscretizationSelection& selection, MurSpeed speed) {
            checkedMaterialSource(source, owner, parameters, selection);
            const double index = parameters.formulation_m == MaterialFormulation::Homogeneous
                                     ? std::sqrt(owner.homogeneousValue())
                                     : 1.0;
            return checkedBoundary(source, electric, magnetic, speed.kernel(index));
        }

        /** @brief Check immutable material metadata before delegating to the base constructor.
         * @details Formulation refusal uses rank-identical cached statistics and
         * performs no collective or allocation. Field identity is a caller invariant.
         * The owner already reduced statistics during construction; this function
         * does not synchronize ranks. Parameters and selection must agree across
         * ranks, and each rank must bind its own matching mesh/layout objects.
         * checkedMaterialBoundary calls this before reading homogeneousValue() or
         * constructing the base, so rejected metadata never selects a boundary speed.
         */
        static AJFourFieldsType& checkedMaterialSource(
            AJFourFieldsType& source, const material::MaterialOwner& owner,
            const FormulationParameters& parameters,
            const SolverDiscretizationSelection& selection) {
            validateMaterialUpdate(parameters, owner.statistics(), selection, owner.layout());
            if (&source.get_mesh() != &owner.mesh() || &source.getLayout() != &owner.layout())
                throw IpplException("ChdrSolver", "Material must share the source mesh and layout");
            return source;
        }

        /** @brief Temporarily substitute the inherited kernel's time spacing and source.
         * @details The kernel captures coefficients and source views by value.
         * Restoring the host bindings therefore preserves physical reconstruction
         * labels even if the inherited step throws. The caller's source is never changed.
         * The local guard saves references to the base members and copies their old
         * values before either assignment. Its destructor restores both on normal
         * return and exception unwinding. This restores parameters, not partially
         * written field data; it is not a rollback of a failed numerical step.
         * Direct assignment to JN_mp is confined here because the public source
         * binding is fixed for the solver's lifetime. internalSource must outlive
         * all device work launched by StandardFDTDSolver::step; the material update
         * owns that scratch field. No history rotation or reconstruction happens here.
         */
        void baseStep(double internalDt, AJFourFieldsType& internalSource) {
            struct Restore {
                double& dt_m;
                AJFourFieldsType*& source_m;
                double savedDt_m;
                AJFourFieldsType* savedSource_m;
                ~Restore() {
                    dt_m     = savedDt_m;
                    source_m = savedSource_m;
                }
            } restore{this->dt, this->JN_mp, this->dt, this->JN_mp};
            // beware: here the timestep is rescaled GLOBALLY across the whole domain
            this->dt    = internalDt;
            this->JN_mp = &internalSource;
            PotentialSolverBase::step();
        }

        /** @brief Lend the existing potential slots and ChdrSolver operations to an update.
         * @details This member stores a reference to its enclosing ChdrSolver, not copies
         * of the fields. The ChdrSolver is noncopyable and nonmovable so that the reference
         * cannot silently point at an old instance. The update borrows this context
         * only during advance(); all dispatch through this adapter occurs on the host.
         */
        class StepContext final : public MaterialStepContext {
        public:
            explicit StepContext(ChdrSolver& solver)
                : solver_m(solver) {}
            AJFourFieldsType& previous() override { return solver_m.A_nm1; }
            AJFourFieldsType& current() override { return solver_m.A_n; }
            AJFourFieldsType& next() override { return solver_m.A_np1; }
            AJFourFieldsType& source() override { return solver_m.sourceField(); }
            double dt() const override { return solver_m.getDt(); }
            void baseStep(double step, AJFourFieldsType& source) override {
                solver_m.baseStep(step, source);
            }
            void finishStep() override { solver_m.finishPotentialStep(); }

        private:
            ChdrSolver& solver_m;
        };

        /// Validate rank-identical metadata before the base can perform halo communication.
        static AJFourFieldsType& checkedSource(AJFourFieldsType& source,
                                               const EMFieldsType& electric,
                                               const EMFieldsType& magnetic) {
            requireHalo(source);
            requireHalo(electric);
            requireHalo(magnetic);
            if (&source.getLayout() != &electric.getLayout()
                || &source.getLayout() != &magnetic.getLayout()
                || &source.get_mesh() != &electric.get_mesh()
                || &source.get_mesh() != &magnetic.get_mesh())
                throw IpplException(
                    "ChdrSolver", "ChDRFieldContainer must share the same mesh and layout owners");
            requireLocalOwnership(source.getLayout());
            const auto description = ippl::FDTDBoundary::fromLayout(source.getLayout());
            if (description.murMask() && stencil() == Stencil::NonStandard)
                throw IpplException("ChdrSolver", "NonStandard requires periodic faces");
            return source;
        }

        /** @brief Preserve ChdrSolver's partition-refusal exception without MPI communication.
         * Every rank has the complete partition table, so all ranks reject the
         * same undersized domain before any solver can exchange a halo.
         */
        static void requireLocalOwnership(const ippl::FieldLayout<SpatialDim>& layout) {
            const auto domains = layout.getHostLocalDomains();
            for (unsigned rank = 0; rank < domains.extent(0); ++rank)
                for (unsigned d = 0; d < SpatialDim; ++d)
                    if (domains(rank)[d].length() < 2)
                        throw std::runtime_error("ChdrSolver requires two owned cells per axis");
        }
        void requireMatching(const EMFieldsType& field) {
            requireHalo(field);
            if (&field.getLayout() != &sourceField().getLayout()
                || &field.get_mesh() != &sourceField().get_mesh())
                throw IpplException("ChdrSolver", "Output must share the source mesh and layout");
        }
        /// Inspect field metadata only; potential, source and E/B storage are all accepted.
        template <class FieldType>
        static void requireHalo(const FieldType& field) {
            if (field.getNghost() != 1)
                throw IpplException("ChdrSolver", "Requires one-cell halos");
        }
        /** @brief Validate storage and the converted Mur speed before base construction.
         * checkedSource preserves ChdrSolver's storage checks, including the partition
         * exception type. validateBoundary then checks the speed that Mur will receive,
         * global domain and face intersections without allocating or communicating.
         */
        static ippl::FDTDBoundary checkedBoundary(AJFourFieldsType& source,
                                                  const EMFieldsType& electric,
                                                  const EMFieldsType& magnetic, double speed) {
            checkedSource(source, electric, magnetic);
            const auto description = ippl::FDTDBoundary::fromLayout(source.getLayout(), speed);
            ippl::validateBoundary(description, source.getLayout(), source.getNghost(),
                                   electric.getNghost(), magnetic.getNghost());
            return description;
        }
        LevelContract contract_m;
        /// Own formulation workspace; a null pointer selects the inherited vacuum step.
        std::unique_ptr<MaterialUpdate> update_m;
        /// Borrow the stable address of this ChdrSolver without creating a new adapter per step.
        StepContext context_m{*this};
        /// Cached material metadata for queries; the formulation owns its numerical data.
        double homogeneousPermittivity_m = 1.0;
        FormulationParameters formulationParameters_m;
        long long completed_m                = 0;
        double overriddenDt_m                = 0;
        bool hasStepOverride_m               = false;
        bool timeStepLocked_m                = false;
        unsigned long long historyRevision_m = 0;
    };
}  // namespace chdr::solver
