/** @file LevelContract.h
 * @brief Define the space and time sampling contract for solver quantities.
 * @details LevelContract gives fields, potentials and sources a shared description
 * of their sample times and component positions for each reconstruction and spatial
 * discretization. Source/history loaders, field monitors and Fourier probes use this
 * description to interpret stored values consistently. It also supplies the
 * time-centering weight used by the gauge monitor.
 * Time labels are integer half-step offsets relative to the completed solve count;
 * spatial offsets are measured from cell centres in mesh-spacing units.
 */
#pragma once

#include <Kokkos_Core.hpp>

#include <array>
#include <stdexcept>

#include "HalfStepOffset.h"
#include "SpatialDimension.h"

namespace chdr::solver {
    /// Select spatial component placement and its associated discrete operators.
    /// Both choices use IPPL's regular field storage.
    enum class SpatialDiscretization {
        /// Field components all live in the cell-centres in space.
        Collocated,
        /// Field components can live displaced in space by half steps of the cell sizes.
        Staggered
    };

    /**
     * @brief Select the time convention and reconstruction of electric E and
     * magnetic B from the scalar potential phi and vector potential A.
     *
     * @details Fields are reconstructed from the potentials through
     * \f[
     * \begin{aligned}
     * \vec{E} &= -\partial_t\vec{A} - \nabla\phi, \\
     * \vec{B} &= \nabla\times\vec{A}.
     * \end{aligned}
     * \f]
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.8)-(3.9)].
     * The electric field has two contributions.
     * The difference between current and previous vector potentials approximates the first
     * contribution halfway between their sample times. The scalar-potential gradient must represent
     * that same time for the two contributions to form a consistently timed electric-field sample.
     *
     * Below, current A represents \f$t_k=k\Delta t\f$ and previous A represents \f$t_{k-1}\f$.
     * Their difference is centred at \f$t_{k-1/2}\f$. "Current" and "previous" identify storage
     * slots. The physical times assigned to the scalar-potential slots depend on this convention.
     *
     * R2 and R2Prime also average the magnetic fields obtained from current and previous A, giving
     * magnetic output at the same half-step as E.
     *
     * For the collocated R2 reading, the reconstruction is (DERIVATION)
     * \f[
     * \begin{aligned}
     * \vec{E}^{\,k-1/2}
     *   &= -\frac{\vec{A}[k]-\vec{A}[k-1]}{\Delta t}
     *      -\nabla_h\phi[k-1], \\
     * \vec{B}^{\,k-1/2}
     *   &= \frac{1}{2}\left(
     *        \nabla_h\times\vec{A}[k]
     *        +\nabla_h\times\vec{A}[k-1]\right).
     * \end{aligned}
     * \f]
     * Square brackets identify storage indices; superscripts identify
     * physical sample times. In R2, \f$\phi[k-1]\f$ represents \f$t_{k-1/2}\f$.
     * Here, the spatial derivatives use centred differences.
     * This is a vector restatement of MITHRA's component formulas,
     * with its step index n replaced by k-1:
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.50)-(3.55)].
     * SpatialDiscretization::Staggered uses different spatial operators and component
     * positions with the same temporal convention.
     * The compatible potential/source equations from
     * [christlieb2024GaugeConserving, Eqs. (25)-(29)] are shown in LevelContract.
     * The storage-slot conventions below are deductions for this solver.
     */
    enum class EBReconstructionConvention {
        /// The same behavior as IPPL's existing `evaluate_EB()` reconstruction.
        /// Current phi represents \f$t_k\f$ and supplies the scalar gradient,
        /// while the A difference \f$t_{k-1/2}\f$.
        /// Thus, E generally combines contributions from different times.
        /// Magnetic output uses the current-time A and represents \f$t_k\f$.
        R1,

        /// A half-step convention for scalar potential and charge.
        /// The current phi is at \f$t_{k+1/2}\f$, so reconstruction reads the previous phi at
        /// \f$t_{k-1/2}\f$ to match the A difference. Magnetic output is averaged
        /// as described above.
        R2,

        /// Use an alternative convention for equivalence tests.
        /// The current phi represents \f$t_{k-1/2}\f$, so reconstruction reads the current phi to
        /// match the A difference. Scalar(phi)-history initialization and
        /// charge-source sampling must shift consistently relative to R2. Changing the field
        /// formula alone is not sufficient. Magnetic output is averaged as in R2.
        R2Prime
    };

    /// Identify the physical quantity whose time label or component position is requested.
    /// This enum selects metadata; it does not attach a physical type to field storage.
    enum class Quantity {
        Electric,
        Magnetic,
        RawMagnetic,
        Phi,
        VectorPotential,
        Charge,
        Current
    };

    /** @brief Component time labels, in half steps after completed solve k.
     *
     * @details In order \f$\vec{E}, \vec{B}, \vec{B}_\text{raw}, \phi, \vec{A}, \rho, \vec{J}\f$,
     * the half-step labels are
     * \f$(\{-1,0\},0,0,0,0,0,0)\f$ for `EBReconstructionConvention::R1`,
     * \f$(-1,-1,0,+1,0,+1,0)\f$ for `EBReconstructionConvention::R2`,
     * and \f$(-1,-1,0,-1,0,-1,0)\f$ for `EBReconstructionConvention::R2Prime`.
     *
     * @details The `phiIsZero`-flag selects the single-label `R1` electric case.
     * `R2` reads previous \f$phi\f$.
     * R2Prime reads current phi, with the corresponding history/source shift.
     * For general smooth data, define \f$t_m=(k-1/2)\Delta t\f$. The backward
     * A difference is centred at \f$t_m\f$, whereas R1's scalar gradient uses
     * \f$\phi[k]\f$ at \f$t_m+\Delta t/2\f$. Taylor expansion therefore gives
     * \f$E_{R1}-E(t_m)=-(\Delta t/2)\nabla\partial_t\phi(t_m)
     * +O(\Delta t^2+h^2)\f$, with \f$h=\max_d h_d\f$. This explains its first-order
     * control role for a time-dependent scalar gradient; setting phi to zero removes that term.
     *
     * This is a deduction from the field definition in
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.9)] and the centered stencil
     * [trefethen1996FiniteDifferenceSpectral, Eq. (4.1.7), p. 151]. R1's raw B retains its
     * whole-step label. Sources: fallahi2020mithra20fullwavesimulation, section 3.1.4;
     * christlieb2024GaugeConserving, Eqs. (25)-(29). Slot table is a deduction.
     *
     * The compatible potential/source time levels appear in the explicit
     * time-discrete system below.
     *
     * [christlieb2024GaugeConserving, Eqs. (25)-(27)]:
     * \f[
     * \begin{aligned}
     * \frac{\phi^{n+3/2}-2\phi^{n+1/2}+\phi^{n-1/2}}
     *      {c^2\Delta t^2}
     *   -\nabla^2\phi^{n+1/2}
     *   &= \frac{\rho^{n+1/2}}{\epsilon_0}, \\
     * \frac{\vec{A}^{n+1}-2\vec{A}^{n}+\vec{A}^{n-1}}
     *      {c^2\Delta t^2}
     *   -\nabla^2\vec{A}^{n}
     *   &= \mu_0\vec{J}^{n}, \\
     * \frac{\phi^{n+3/2}-\phi^{n+1/2}}{c^2\Delta t}
     *   +\nabla\cdot\vec{A}^{n+1}
     *   &= 0.
     * \end{aligned}
     * \f]
     *
     * The matching continuity equation [ibid., Eq. (28)] is:
     * \f[
     * \frac{\rho^{n+3/2}-\rho^{n+1/2}}{\Delta t}
     *   +\nabla\cdot\vec{J}^{n+1}=0.
     * \f]
     *
     * The gauge and continuity residuals [ibid., Eq. (29)] are:
     * \f[
     * \begin{aligned}
     * \epsilon_1^k
     *   &= \frac{\phi^{k+1/2}-\phi^{k-1/2}}{c^2\Delta t}
     *      +\nabla\cdot\vec{A}^{k}, \\
     * \epsilon_2^k
     *   &= \frac{\rho^{k+1/2}-\rho^{k-1/2}}{\Delta t}
     *      +\nabla\cdot\vec{J}^{k}.
     * \end{aligned}
     * \f]
     */
    struct LevelContract {
        EBReconstructionConvention ebReconstructionConvention_m = EBReconstructionConvention::R2;
        SpatialDiscretization spatialDiscretization_m           = SpatialDiscretization::Collocated;

        /** @brief Component displacement from its cell centre, in mesh-spacing units.
         * @param quantity Channel whose spatial placement is requested.
         * @param component Cartesian component index, from zero through two.
         * @return Device-callable array of dimensionless x/y/z offsets.
         * @details Staggered scalar/charge remain at centres; vector potential,
         * current and electric components are on their positive faces; magnetic
         * components are on the positive transverse edges.
         * Cartesian placement follows [clemens2001FiniteIntegration, section 1,
         * Figs. 1-3; ryu2016potentialfdtd, Appendix A, Eqs. (19)-(24)]; the
         * positive offsets are this storage convention's specialization (deduction).
         * @pre Vector component is 0, 1 or 2. Offsets are dimensionless fractions of
         * mesh spacing; callers add them to their cell-centre coordinate before sampling.
         * This small value-only function is safe to capture in a device sampling kernel.
         */
        KOKKOS_INLINE_FUNCTION Kokkos::Array<double, SpatialDim> positionOffset(
            Quantity quantity, unsigned component) const {
            Kokkos::Array<double, SpatialDim> offset{0.0, 0.0, 0.0};

            // In the collocated spatial discretization all Quantity::* live at cell centers.
            if (spatialDiscretization_m == SpatialDiscretization::Collocated
                // Even in the staggered spatial discretization, phi and rho still live at cell
                // center.
                || quantity == Quantity::Phi || quantity == Quantity::Charge)
                return offset;
            // Magnetic quantities are either marked as Quantity::Magnetic or Quantity::RawMagnetic.
            // Raw magnetic is not averaged.
            const bool magnetic =
                quantity == Quantity::Magnetic || quantity == Quantity::RawMagnetic;
            for (unsigned d = 0; d < SpatialDim; ++d) {
                offset[d] = (magnetic ? d != component : d == component) ? 0.5 : 0.0;
            }
            return offset;
        }

        /** @brief Return the sample-time offset in integer half timesteps relative to the step
         * index.
         * @return Offset m: -1 means half a timestep earlier, 0 the indexed time,
         * and +1 half a timestep later.
         * @details For index k, labelTime() converts this offset to
         * \f$t_{\mathrm{sample}}=(k+m/2)\Delta t\f$.
         * Outputs use the completed solve count as k. Sources use the source index
         * written before the next solve. This integer encoding represents the time
         * assignments derived in the class documentation from
         * [fallahi2020mithra20fullwavesimulation, section 3.1.4, Eqs. (3.50)-(3.55)].
         * R1 can reconstruct E with nonzero phi, but its two contributions generally
         * represent different times. This API therefore refuses a single E label
         * unless the caller explicitly declares zero scalar potential.
         * Identity monitors cannot validate externally supplied source timestamps:
         * a common charge/current shift can preserve their algebraic identities.
         * Analytic source comparisons must independently pin these labels.
         */
        HalfStepOffset halfStepLabel(Quantity quantity, bool phiIsZero = false) const {
            if (quantity == Quantity::Electric) {
                // Reject a single time label for mixed-time R1 electric output;
                // this does not prohibit reconstructing that output (see class documentation).
                if (ebReconstructionConvention_m == EBReconstructionConvention::R1 && !phiIsZero) {
                    throw std::invalid_argument(
                        "R1 electric field has no unique label with phi != 0");
                }

                // Minus one half-step labels E at \f$t_{k-1/2}\f$ for either
                // supported spatial discretization (see the class time-level table).
                return HalfStepOffset{-1};
            }
            if (ebReconstructionConvention_m != EBReconstructionConvention::R1) {
                if (quantity == Quantity::Magnetic)
                    return HalfStepOffset{-1};
                if (quantity == Quantity::Phi || quantity == Quantity::Charge)
                    return HalfStepOffset{
                        ebReconstructionConvention_m == EBReconstructionConvention::R2 ? 1 : -1};
            }
            return HalfStepOffset{0};
        }

        /** Convert an integer index directly to time; never advance a floating-point clock by
         * addition.
         *
         * @param completed k
         * */
        double labelTime(Quantity quantity, long long completed, double dt,
                         bool phiIsZero = false) const {
            return (static_cast<double>(completed)
                    + 0.5 * halfStepLabel(quantity, phiIsZero).count())
                   * dt;
        }

        /** @brief Weight current A to match the time of the scalar difference in the gauge monitor.
         * @return Weight w of current A; previous A has weight 1-w.
         * @details This supplies a temporal rule for the gauge monitor; it does not
         * evaluate a residual or modify the fields. In the solver's normalized units,
         * the monitor forms the discrete Lorenz-gauge residual
         * \f[
         * G_w = \frac{\phi[k]-\phi[k-1]}{\Delta t}
         *       + D_h\cdot\left(w\vec{A}[k]+(1-w)\vec{A}[k-1]\right).
         * \f]
         * Here k is the completed solve count. Brackets identify current [k]
         * and previous [k-1] storage slots, as elsewhere in this contract.
         * The vector slots represent \f$t_k\f$ and \f$t_{k-1}\f$ in all conventions;
         * the scalar slots have convention-dependent physical times.
         * Both terms in the residual must represent the same time:
         * - R1: the scalar difference is centred at \f$t_{k-1/2}\f$;
         *   use the average of current and previous A, so \f$w=1/2\f$.
         * - R2: the scalar difference is centred at \f$t_k\f$;
         *   use current A, so \f$w=1\f$.
         * - R2Prime: the scalar difference is centred at \f$t_{k-1}\f$;
         *   use previous A, so \f$w=0\f$.
         * These choices are a deduction from the class time-level table and
         * [christlieb2024GaugeConserving, Eqs. (25)-(29), p. 5], with the
         * gauge residual defined in Eq. (29). SpatialDiscretization independently
         * selects the spatial divergence. Matching the gauge terms in time does
         * not give R1's mixed-time electric output a single time label.
         */
        KOKKOS_INLINE_FUNCTION constexpr double gaugeWeight() const {
            return ebReconstructionConvention_m == EBReconstructionConvention::R1   ? 0.5
                   : ebReconstructionConvention_m == EBReconstructionConvention::R2 ? 1.0
                                                                                    : 0.0;
        }
    };
}  // namespace chdr::solver
