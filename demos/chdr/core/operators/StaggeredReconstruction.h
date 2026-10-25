/** @file StaggeredReconstruction.h
 * @brief R2 reconstruction at the natural face/edge locations of the Standard reading.
 * @details Arrays are interpreted at their component offsets; no interpolation
 * or extra wave update occurs. E uses the previous scalar slot and the one-step
 * vector-potential difference [christlieb2024GaugeConserving, section 3.1,
 * proof of Theorem 3.5, p. 8].
 * The staggered curl follows [ryu2016potentialfdtd, Appendix A, Eq. (24)].
 */
#pragma once

#include "SpatialDimension.h"
#include "StaggeredOperators.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "StaggeredReconstruction.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    namespace staggeredDetail {
        /** @brief Evaluate staggered electric components on a cell's positive faces.
         * @tparam FourPotentialView Device-callable three-dimensional view whose elements store the
         * scalar potential in component 0 and vector-potential components in 1 through SpatialDim.
         * E.g. `view(i, j, k)[1] == A_x`.
         * @param curAView Current R2 potential slot: vector potential at \f$t_n\f$ and
         * scalar potential at \f$t_{n+1/2}\f$. The current scalar component is not read.
         * @param prevAView Previous R2 potential slot: vector potential at \f$t_{n-1}\f$
         * and scalar potential at \f$t_{n-1/2}\f$.
         * @param spacing Positive Cartesian mesh spacings, ordered x, y, z.
         * @param inverseDt Reciprocal \f$1/\Delta t\f$ of the positive timestep separating
         * the vector-potential histories.
         * @param i Local x array index, including the halo offset.
         * @param j Local y array index, including the halo offset.
         * @param k Local z array index, including the halo offset; not a temporal index.
         * @pre Both views use the same index layout and staggered component positions.
         * Vector components at the selected index and previous scalar samples at that
         * index and its positive-axis neighbours are valid; required halos are synchronized.
         * @return Electric components at \f$t_{n-1/2}\f$ on the selected cell's positive
         * x, y and z faces. The returned components occupy different spatial positions.
         * @details Scalar samples lie at cell centres. Vector component d+1 lies on the
         * positive d face of the same logical cell. The forward scalar difference uses
         * the two centres adjoining that face, matching the vector-potential position.
         * With n denoting the temporal index, the reconstruction is
         * \f$E_d^{n-1/2}=-(A_d[n]-A_d[n-1])/\Delta t-D_{+,d}\phi[n-1]\f$.
         * Square brackets identify storage slots; superscripts identify sample times.
         * The R2 slot mapping follows LevelContract and is a deduction from
         * [christlieb2024GaugeConserving, section 3.1, Eqs. (25)-(27), and the proof of
         * Theorem 3.5, p. 8]. The spatial placement specializes
         * [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(24)] to the positive-face
         * indexing convention. TODO: VERIFY.
         */
        template <class FourPotentialView>
        KOKKOS_INLINE_FUNCTION ippl::Vector<double, SpatialDim> electricAt(
            FourPotentialView curAView, FourPotentialView prevAView,
            ippl::Vector<double, SpatialDim> spacing, double inverseDt, int i, int j, int k) {
            ippl::Vector<double, SpatialDim> dAdt, gradPhi;
            for (unsigned d = 0; d < SpatialDim; ++d) {
                dAdt[d] = (curAView(i, j, k)[d + 1] - prevAView(i, j, k)[d + 1]) * inverseDt;
                gradPhi[d] =
                    staggeredDifference(Component{prevAView, 0u}, i, j, k, d, spacing[d], true);
            }
            return -dAdt - gradPhi;
        }
    }  // namespace staggeredDetail

    /** @brief Reconstruct the unaveraged magnetic field from the current vector potential.
     * @tparam FourPotentialA Field storing scalar potential in component 0 and vector
     * potential in components 1 through SpatialDim on the staggered mesh.
     * @tparam VectorFieldB Field storing the SpatialDim magnetic components on edges.
     * @param currentA Current potential history; its vector components represent \f$t_n\f$
     * and lie on positive cell faces. The scalar component is not read.
     * @param B Allocated magnetic output with the same mesh and local index layout.
     * @pre Vector-potential samples required by the forward curl have valid one-cell halos.
     * @details The magnetic-potential relation and its spatial discretization are
     * \f[\mathbf{B}=\nabla\times\mathbf{A},\qquad
     *     \mathbf{B}_{\mathrm{raw}}^n=C_+\mathbf{A}[n].\f]
     * The continuum relation has a positive sign
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.8)]. Here \f$C_+\f$ is the
     * forward face-to-edge curl implemented by staggeredCurl
     * [ryu2016potentialfdtd, Appendix A, Eq. (24)]. Component d lies on the edge
     * displaced half a cell along the other two axes from its logical cell centre.
     * Clear the complete output view, write owned cells, then fence before returning.
     * Output halos remain zero until the caller synchronizes them.
     */
    template <class FourPotentialA, class VectorFieldB>
    void rawStaggeredMagnetic(const FourPotentialA& currentA, VectorFieldB& B) {
        using ConstPotentialView              = typename FourPotentialA::view_type::const_type;
        const ConstPotentialView currentAView = currentA.getView();
        const auto magneticView               = B.getView();
        const auto spacing                    = currentA.get_mesh().getMeshSpacing();
        B                                     = 0.0;
        Kokkos::parallel_for(
            "ChDR compact raw magnetic field", currentA.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                magneticView(i, j, k) = staggeredCurl(currentAView, i, j, k, spacing, 1);
            });
        Kokkos::fence();
    }

    /** @brief Reconstruct staggered electric and magnetic output at the R2 half step.
     * @tparam FourPotentialA Field storing scalar potential in component 0 and vector
     * potential in components 1 through SpatialDim; positions follow LevelContract.
     * @tparam VectorFieldEB Three-component field storage shared by the E and B outputs.
     * @param currentA Current R2 slot: vector potential at \f$t_n\f$, scalar at
     * \f$t_{n+1/2}\f$. The current scalar component is not read.
     * @param previousA Previous R2 slot: vector potential at \f$t_{n-1}\f$, scalar at
     * \f$t_{n-1/2}\f$.
     * @param E Allocated electric output; component d lies on the positive d face.
     * @param B Allocated magnetic output; component d lies on an edge parallel to axis d,
     * displaced half a cell along the other two axes from the logical cell centre.
     * @param dt Positive timestep separating the vector-potential histories.
     * @pre All fields share a mesh and index layout; required potential halos are synchronized.
     * @details staggeredDetail::electricAt defines the electric reconstruction. Magnetic
     * reconstruction applies \f$\mathbf{B}=\nabla\times\mathbf{A}\f$
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.8)] to both vector histories:
     * \f[\mathbf{B}^{n-1/2}=\frac{1}{2}
     *       \left(C_+\mathbf{A}[n]+C_+\mathbf{A}[n-1]\right).\f]
     * The averaging follows [ibid., Eqs. (3.50)-(3.52)]; \f$C_+\f$ is the
     * staggeredCurl operator documented in rawStaggeredMagnetic. The two curls are
     * evaluated separately before addition, preserving their arithmetic order.
     * Both outputs have the half-step time label. Clear full output views, write owned
     * cells and fence before return; output halos remain zero. Captured read-only potential
     * views and scalars avoid capturing the owning solver. No potential update occurs here.
     */
    template <class FourPotentialA, class VectorFieldEB>
    void reconstructStaggered(const FourPotentialA& currentA, const FourPotentialA& previousA,
                              VectorFieldEB& E, VectorFieldEB& B, double dt) {
        using ConstPotentialView               = typename FourPotentialA::view_type::const_type;
        const ConstPotentialView currentAView  = currentA.getView(),
                                 previousAView = previousA.getView();
        const auto electricView = E.getView(), magneticView = B.getView();
        const auto spacing     = currentA.get_mesh().getMeshSpacing();
        const double inverseDt = 1.0 / dt;
        E                      = 0.0;
        B                      = 0.0;
        Kokkos::parallel_for(
            "ChDR compact R2 fields", currentA.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                electricView(i, j, k) = staggeredDetail::electricAt(currentAView, previousAView,
                                                                    spacing, inverseDt, i, j, k);
                magneticView(i, j, k) = (staggeredCurl(currentAView, i, j, k, spacing, 1)
                                         + staggeredCurl(previousAView, i, j, k, spacing, 1))
                                        * 0.5;
            });
        Kokkos::fence();
    }
}  // namespace chdr::solver
