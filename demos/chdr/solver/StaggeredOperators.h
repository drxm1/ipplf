/** @file StaggeredOperators.h
 * @brief Adjacent-sample differences and curl for the staggered spatial discretization.
 * @details Callers supply samples at the component positions documented by LevelContract;
 * these helpers evaluate differences without relocating data or checking its placement.
 * A scalar difference represents the derivative halfway between its input samples
 * [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(22)].
 * The incidence identities follow [clemens2001FiniteIntegration,
 * section 2, Eqs. (18)-(24)] and [teixeira2014LatticeMaxwell, section 4,
 * Eqs. (13)-(15)]. Here each Cartesian operator is a difference between two
 * adjacent samples. The coordinate shifts commute, giving the corresponding
 * discrete cancellations on a periodic mesh (deduction).
 */
#pragma once

#include "Operators.h"
#include "SpatialDimension.h"

namespace chdr::solver {
    /** @brief Difference to the adjacent half-cell location, or its absolute evaluation.
     * @tparam Absolute Select the sum of input magnitudes used for diagnostic scaling.
     * @param sample Scalar accessor callable as sample(i, j, k).
     * @param i Local x index, including any halo offset.
     * @param j Local y index, including any halo offset.
     * @param k Local z index, including any halo offset.
     * @param axis Differentiation direction: zero for x, one for y, two for z.
     * @param h Positive distance between adjacent input samples along the selected axis.
     * @param forward Use the positive neighbour when true, the negative neighbour otherwise.
     * @pre axis < SpatialDim and h > 0.
     * @pre The accessor is valid at (i, j, k) and the selected adjacent index.
     * For field-backed accessors, any required ghost-cell value is synchronized before this call.
     * @details Write \f$f_n=\mathtt{sample}(i,j,k)\f$, where n indexes only the selected
     * axis; the other two indices stay fixed. Its physical coordinate \f$x_n\f$ includes
     * any existing component offset. For Absolute=false, the input pairs and output
     * positions are [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(20)]:
     * \f[
     * \begin{aligned}
     * \text{forward:}\quad (D_+f)_{n+1/2}
     *   &= \frac{f_{n+1}-f_n}{h}, & x_{\mathrm{out}} &= x_n+\frac{h}{2},\\
     * \text{backward:}\quad (D_-f)_{n-1/2}
     *   &= \frac{f_n-f_{n-1}}{h}, & x_{\mathrm{out}} &= x_n-\frac{h}{2}.
     * \end{aligned}
     * \f]
     * The forward input pair is at \f$x_n,x_n+h\f$; the backward pair is at
     * \f$x_n-h,x_n\f$. The returned value has the stated midpoint location;
     * the helper does not return coordinates or interpolate a stored field.
     * Consequently \f$\sum_dD_{-,d}D_{+,d}=L_{Std}\f$ by direct substitution
     * (deduction from the difference/divergence definitions in the same source,
     * Appendix A, Eqs. (19)-(23)).
     */
    template <bool Absolute = false, class Sample>
    KOKKOS_INLINE_FUNCTION double staggeredDifference(Sample sample, int i, int j, int k,
                                                      unsigned axis, double h, bool forward) {
        static_assert(SpatialDim == 3, "staggeredDifference uses a three-index spatial accessor");
        const int di = axis == 0, dj = axis == 1, dk = axis == 2;
        const double plus  = forward ? sample(i + di, j + dj, k + dk) : sample(i, j, k);
        const double minus = forward ? sample(i, j, k) : sample(i - di, j - dj, k - dk);
        if constexpr (Absolute)
            return (Kokkos::abs(plus) + Kokkos::abs(minus)) / h;
        else
            return (plus - minus) / h;
    }

    /** @brief Curl from adjacent component samples at their natural staggered positions.
     * @param view Device view containing the selected three vector components.
     * @param i Local x index, including the halo offset.
     * @param j Local y index, including the halo offset.
     * @param k Local z index, including the halo offset.
     * @param spacing Positive Cartesian mesh spacings.
     * @param firstComponent First vector component: one for a four-potential, zero for E/B.
     * @param forward Forward face-to-edge differences when true; backward edge-to-face
     * differences otherwise. Staggered E/B reconstruction uses only the forward form.
     * @pre Every neighbour read is owned or has a synchronized one-cell halo.
     * @return Three cyclic curl components at the selected natural positions.
     * @details Component d of a forward result lies half a cell along each
     * transverse axis from the scalar centre; this is the Cartesian curl
     * specialization of [ryu2016potentialfdtd, Appendix A, Eq. (24)].
     * This read-only device helper performs no communication or data movement.
     */
    template <class View>
    KOKKOS_INLINE_FUNCTION ippl::Vector<double, SpatialDim> staggeredCurl(
        View view, int i, int j, int k, ippl::Vector<double, SpatialDim> spacing,
        unsigned firstComponent = 0, bool forward = true) {
        static_assert(SpatialDim == 3, "staggeredCurl requires three spatial dimensions");
        ippl::Vector<double, SpatialDim> result;
        for (unsigned d = 0; d < SpatialDim; ++d) {
            const unsigned a = (d + 1) % SpatialDim, b = (d + 2) % SpatialDim;
            result[d] = staggeredDifference(Component<View>{view, firstComponent + b}, i, j, k, a,
                                            spacing[a], forward)
                        - staggeredDifference(Component<View>{view, firstComponent + a}, i, j, k, b,
                                              spacing[b], forward);
        }
        return result;
    }

    /** @brief Signed Fourier factor at the natural half-cell location; phase comes from the offset.
     * @param waveNumber Signed physical wave number along the differentiated axis.
     * @param spacing Positive cell spacing along that axis.
     * @return Signed factor, retaining the wave number's orientation.
     * @details \f$K=2\sin(kh/2)/h\f$ follows by Fourier substitution into
     * [ryu2016potentialfdtd, Appendix A, Eqs. (19)-(20)] (deduction).
     */
    inline double staggeredDerivativeSymbol(double waveNumber, double spacing) {
        return waveSymbol(waveNumber, spacing);
    }
}  // namespace chdr::solver
