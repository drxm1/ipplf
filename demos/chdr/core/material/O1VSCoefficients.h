/** @file O1VSCoefficients.h
 * @brief Cell and positive-face coefficient rules for the variable-scalar gauge.
 */
#pragma once
#include "core/material/MaterialOwner.h"

namespace chdr::material {
    /** @brief Store \f$\varepsilon^{-2}\f$ at cell centres, without the gauge parameter.
     * @details Chew's \f$\chi=\alpha\varepsilon^2\mu\f$ with \f$\mu=1\f$ gives this factor;
     * the update supplies \f$1/\alpha\f$ [chew2014generalizedgauge, Eq. (18), p. 3].
     * Here epsilon is the stored cell-centred relative permittivity. Keeping alpha
     * out of the cached field allows the material coefficient to be reused for
     * different gauge parameters without changing the material map or cache key.
     */
    struct InverseSquare {
        static constexpr const char* Name = "inverse-square";
        static constexpr int ForwardAxis  = -1;
        KOKKOS_INLINE_FUNCTION double operator()(ScalarFieldConstView eps, int i, int j,
                                                 int k) const {
            return 1 / (eps(i, j, k) * eps(i, j, k));
        }
    };

    /** @brief Harmonic permittivity at the positive face normal to axis d.
     * @details The interface lies midway between centre samples, giving
     * \f$H(a,b)=2ab/(a+b)\f$ [eymard2000FiniteVolume, section 7.1, Example 7.1,
     * p. 22]. The output index denotes the face between that cell and its next
     * neighbour along d, not the negative face. ForwardAxis tells MaterialOwner
     * which lower physical ghost needs its boundary pair recomputed.
     *
     * Equal inputs return the original value directly, preserving homogeneous
     * bits without forming products and a quotient. Input halos must be valid
     * before evaluation; MaterialOwner supplies them. The caller is responsible
     * for selecting material values on which this coefficient is well defined.
     */
    template <unsigned d>
    struct HarmonicFace {
        static_assert(d < 3, "A face axis must be x, y or z");
        static constexpr const char* Name = d == 0   ? "harmonic-x"
                                            : d == 1 ? "harmonic-y"
                                                     : "harmonic-z";
        static constexpr int ForwardAxis  = d;
        KOKKOS_INLINE_FUNCTION double operator()(ScalarFieldConstView eps, int i, int j,
                                                 int k) const {
            CellIndex upper(i, j, k);
            ++upper[d];
            const double a = eps(i, j, k), b = eps(upper[0], upper[1], upper[2]);
            if (a == b)
                return a;
            const double product = 2 * a * b;
            // Reciprocal form avoids the overflowing product for admitted large inputs.
            return Kokkos::isfinite(product) ? product / (a + b) : 2 / (1 / a + 1 / b);
        }
    };
}  // namespace chdr::material
