/** @file Reconstruction.h
 * @brief Reconstruct cell-centred fields from scalar and vector potential histories.
 * @details These kernels use SpatialDiscretization::Collocated and support
 * histories produced by either the Standard or NonStandard potential solver.
 * They read the supplied histories and write field outputs without advancing
 * the potentials. The caller supplies synchronized potential halos and histories
 * sampled at the times specified by LevelContract for the selected convention.
 */
#pragma once

#include "Operators.h"
#include "SpatialDimension.h"

namespace chdr::solver {
    /** @brief Reconstruct the collocated R2 electric field from two potential histories.
     * @tparam FourPotentialA Field with scalar-plus-vector potential values at each cell.
     * @tparam VectorFieldE Field with SpatialDim electric components at each cell.
     * @param curA Current potential slot; component 0 stores the scalar potential,
     * and components 1 through SpatialDim store the vector potential.
     * @param prevA Previous potential slot, with the same component layout.
     * @param E Electric-field output on the same mesh and local index layout.
     * @param dt Positive timestep separating the vector-potential histories.
     * @details The R2 reconstruction is
     * \f$\mathbf{E}^{k-1/2}=-(\mathbf{A}[k]-\mathbf{A}[k-1])/\Delta t-D_c\phi[k-1]\f$.
     * Square brackets identify storage slots; the superscript identifies sample time.
     * The previous scalar slot represents the same half-step time as the vector
     * difference; see LevelContract. This is the slot convention applied to
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.53)-(3.55)] (deduction).
     * Preserve the subtraction and multiplication order of FDTDSolverBase::evaluate_EB
     * for the R2/R2Prime bitwise comparison. Input views have const elements;
     * the kernel captures only views and scalars, never an owning solver.
     * @pre Both histories use collocated samples on the output layout; the previous
     * scalar potential has synchronized one-cell halos. Mesh spacings are positive.
     * @details Clear the full output view, then reconstruct owned cells. Output halos
     * remain zero until the caller refreshes them. The final fence completes the
     * field writes before return; this function performs no MPI exchange or host copy.
     */
    template <class FourPotentialA, class VectorFieldE>
    void reconstructElectric(const FourPotentialA& curA, const FourPotentialA& prevA,
                             VectorFieldE& E, double dt) {
        static_assert(SpatialDim == 3,
                      "reconstructElectric uses three-index access and explicit x/y/z components");
        using ConstPotentialView   = typename FourPotentialA::view_type::const_type;
        E                          = 0.0;  // TODO: Is this write really needed?
        const ConstPotentialView a = curA.getView(), aPrev = prevA.getView();
        const auto e = E.getView();
        const auto inverseSpacing =
            ippl::Vector<double, SpatialDim>(0.5) / curA.get_mesh().getMeshSpacing();
        const double inverseDt = 1.0 / dt;
        // Write owned cells only; neighbouring potential samples may lie in halos.
        Kokkos::parallel_for(
            "R2 electric reconstruction", curA.getFieldRangePolicy(),
            KOKKOS_LAMBDA(size_t i, size_t j, size_t k) {
                // Components 1–3 store the vector potential.
                // Difference the current and previous slots using the supplied timestep.
                const ippl::Vector<double, SpatialDim> dAdt{
                    (a(i, j, k)[1] - aPrev(i, j, k)[1]) * inverseDt,
                    (a(i, j, k)[2] - aPrev(i, j, k)[2]) * inverseDt,
                    (a(i, j, k)[3] - aPrev(i, j, k)[3]) * inverseDt,
                };
                // In R2, the previous scalar slot has the same half-step time
                // as the vector-potential time difference; see LevelContract.
                // Alternative implementation: difference all four components, then select phi (more
                // unnecessary operations).
                //
                // const ippl::Vector<double, 4> dx =
                //     (aPrev(i + 1, j, k) - aPrev(i - 1, j, k)) * inverseSpacing[0];
                // const ippl::Vector<double, 4> dy =
                //     (aPrev(i, j + 1, k) - aPrev(i, j - 1, k)) * inverseSpacing[1];
                // const ippl::Vector<double, 4> dz =
                //     (aPrev(i, j, k + 1) - aPrev(i, j, k - 1)) * inverseSpacing[2];
                // const ippl::Vector<double, SpatialDim> gradientPhi{dx[0], dy[0], dz[0]};

                // Only scalar component 0 contributes to the electric-field gradient.
                const double dx =
                    (aPrev(i + 1, j, k)[0] - aPrev(i - 1, j, k)[0]) * inverseSpacing[0];
                const double dy =
                    (aPrev(i, j + 1, k)[0] - aPrev(i, j - 1, k)[0]) * inverseSpacing[1];
                const double dz =
                    (aPrev(i, j, k + 1)[0] - aPrev(i, j, k - 1)[0]) * inverseSpacing[2];
                const ippl::Vector<double, SpatialDim> gradientPhi{dx, dy, dz};
                e(i, j, k) = -dAdt - gradientPhi;
            });
        // Complete the field writes before returning to the caller.
        // Completion semantics: Kokkos Core API, fence.
        // https://kokkos.org/kokkos-core-wiki/API/core/parallel-dispatch/fence.html
        Kokkos::fence();
    }

    /** @brief Reconstruct collocated magnetic output, optionally averaged over two histories.
     * @tparam FourPotentialA Field with scalar potential in component 0 and vector
     * potential in components 1 through SpatialDim, all sampled at cell centres.
     * @tparam VectorFieldB Field with SpatialDim magnetic components at each cell centre.
     * @param curA Current potential slot; vector components represent \f$t_n\f$.
     * @param prevA Previous potential slot; vector components represent \f$t_{n-1}\f$.
     * It is read only when average is true. Neither scalar component is read.
     * @param B Allocated magnetic output on the same mesh and local index layout.
     * @param average Select the mean of the two curls when true; otherwise use only curA.
     * @pre Current vector-potential halos are valid; averaging also needs previous halos.
     * @details The continuous relation is \f$\mathbf{B}=\nabla\times\mathbf{A}\f$
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.8)]. With \f$C_c\f$ denoting
     * the centred curl implemented by curl, the two outputs are
     * \f[\mathbf{B}_{\mathrm{raw}}^n=C_c\mathbf{A}[n],\qquad
     * \mathbf{B}^{n-1/2}=\frac{1}{2}
     *   \left(C_c\mathbf{A}[n]+C_c\mathbf{A}[n-1]\right).\f]
     * The averaging follows [ibid., Eqs. (3.50)-(3.52)]; curl documents the spatial
     * stencil and its source. Evaluate the two curls separately before addition to retain
     * the arithmetic order. Clear the full output view, write owned cells and fence
     * before return. Output halos remain zero; the caller performs any subsequent exchange.
     * The kernel captures read-only potential views and scalars, not an owning solver.
     */
    template <class FourPotentialA, class VectorFieldB>
    void reconstructMagnetic(const FourPotentialA& curA, const FourPotentialA& prevA,
                             VectorFieldB& B, bool average) {
        static_assert(SpatialDim == 3, "reconstructMagnetic uses three-index field access");
        using ConstPotentialView   = typename FourPotentialA::view_type::const_type;
        B                          = 0.0;
        const ConstPotentialView a = curA.getView(), aPrev = prevA.getView();
        const auto b = B.getView();
        const auto h = curA.get_mesh().getMeshSpacing();
        Kokkos::parallel_for(
            "Magnetic reconstruction", curA.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                const auto now = curl(a, i, j, k, h, 1);
                if (average)
                    b(i, j, k) = 0.5 * (now + curl(aPrev, i, j, k, h, 1));
                else
                    b(i, j, k) = now;
            });
        Kokkos::fence();
    }
}  // namespace chdr::solver
