/** @file O1VSOperator.h
 * @brief Device expressions for the variable-scalar generalized gauge.
 */
#pragma once
#include "core/formulation/MaterialStepContext.h"
#include "core/material/O1VSCoefficients.h"

namespace chdr::solver {
    /** @brief Read one potential level and immutable face/node coefficients.
     * @details Components 1--3 lie on positive faces; component zero lies at nodes.
     * With backward divergence G, this evaluates \f$g=G(\varepsilon A)\f$,
     * \f$\gamma=G(A)\f$ and \f$\chi=\alpha\varepsilon^2\f$
     * [chew2014generalizedgauge, Eqs. (15)--(18)].
     * Views own allocation references; mesh/layout and synchronized one-cell halos
     * must remain valid. Mixed-sign edge reads need no extra scratch or exchange.
     */
    struct O1VSOperator {
        AJFourFieldsType::view_type::const_type potential_m;
        Kokkos::Array<material::ScalarFieldConstView, 3> face_m;
        material::ScalarFieldConstView inverseSquare_m;
        ippl::Vector<double, 3> spacing_m;
        double inverseAlpha_m;

        O1VSOperator(const AJFourFieldsType& potential, const material::MaterialOwner& owner,
                     double alpha)
            : potential_m(potential.getView())
            , face_m{owner.derive(material::HarmonicFace<0>{}),
                     owner.derive(material::HarmonicFace<1>{}),
                     owner.derive(material::HarmonicFace<2>{})}
            , inverseSquare_m(owner.derive(material::InverseSquare{}))
            , spacing_m(owner.mesh().getMeshSpacing())
            , inverseAlpha_m(1 / alpha) {}

        /// Backward divergence, weighted by face permittivity when requested.
        KOKKOS_INLINE_FUNCTION double divergence(int i, int j, int k, bool weighted) const {
            double result = 0;
            for (unsigned d = 0; d < 3; ++d) {
                material::CellIndex p(i, j, k);
                --p[d];
                const double here = weighted ? face_m[d](i, j, k) : 1;
                const double back = weighted ? face_m[d](p[0], p[1], p[2]) : 1;
                result += (here * potential_m(i, j, k)[d + 1]
                           - back * potential_m(p[0], p[1], p[2])[d + 1])
                          / spacing_m[d];
            }
            return result;
        }
        KOKKOS_INLINE_FUNCTION double g(int i, int j, int k) const {
            return divergence(i, j, k, true);
        }
        KOKKOS_INLINE_FUNCTION double gamma(int i, int j, int k) const {
            return divergence(i, j, k, false);
        }
        /// Component zero uses \f$\chi^{-1}\f$; vector component \f$d+1\f$ uses
        /// \f$\varepsilon_d^{-1}\f$.
        KOKKOS_INLINE_FUNCTION double inverseMass(unsigned c, int i, int j, int k) const {
            return c == 0 ? inverseAlpha_m * inverseSquare_m(i, j, k) : 1 / face_m[c - 1](i, j, k);
        }
        /** @brief Evaluate \f$(C_h A)_d=D_d^+\gamma-\varepsilon_dD_d^+(g/\chi)\f$.
         * Deduction from Chew's generalized-gauge vector equation (17), with
         * forward gradients and backward divergence as in [ryu2016potentialfdtd, App. A].
         */
        KOKKOS_INLINE_FUNCTION double curlCorrection(unsigned d, int i, int j, int k) const {
            material::CellIndex p(i, j, k);
            ++p[d];
            return (gamma(p[0], p[1], p[2]) - gamma(i, j, k)
                    - face_m[d](i, j, k)
                          * (g(p[0], p[1], p[2]) * inverseMass(0, p[0], p[1], p[2])
                             - g(i, j, k) * inverseMass(0, i, j, k)))
                   / spacing_m[d];
        }
        /** @brief Flux correction \f$D_\Phi\phi=\sum_dD_d^-[(\varepsilon_d-1)D_d^+\phi]\f$.
         * This factors the variable-coefficient scalar equation (16) of Chew
         * into the Standard Laplacian and a correction that vanishes exactly in vacuum.
         */
        KOKKOS_INLINE_FUNCTION double scalarCorrection(int i, int j, int k) const {
            double result = 0;
            for (unsigned d = 0; d < 3; ++d) {
                material::CellIndex p(i, j, k), m(i, j, k);
                ++p[d];
                --m[d];
                result += ((face_m[d](i, j, k) - 1)
                               * (potential_m(p[0], p[1], p[2])[0] - potential_m(i, j, k)[0])
                           - (face_m[d](m[0], m[1], m[2]) - 1)
                                 * (potential_m(i, j, k)[0] - potential_m(m[0], m[1], m[2])[0]))
                          / (spacing_m[d] * spacing_m[d]);
            }
            return result;
        }
    };
}  // namespace chdr::solver
