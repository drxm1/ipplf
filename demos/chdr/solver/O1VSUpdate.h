/** @file O1VSUpdate.h
 * @brief Effective-source O1-VS update through the unchanged Standard kernel.
 */
#pragma once
#include "O1VSOperator.h"
#include "Operators.h"

namespace chdr::solver {
    /** @brief Add variable-material terms before the physical-time Standard step.
     * @details For each component with inverse mass m, the effective source is
     * \f$s=mJ+(m-1)L_c A-m C_hA\f$; for the scalar component the last term is
     * \f$+mD_\Phi\Phi\f$. Substitution into the Standard recurrence gives Chew's
     * equations (16)--(18) [chew2014generalizedgauge], discretized by forward/backward
     * differences [ryu2016potentialfdtd, App. A]. The caller retains physical sources.
     * One four-component workspace and one device pass precede baseStep's exchange.
     */
    class O1VSUpdate final : public MaterialUpdate {
    public:
        /// Derive immutable coefficients collectively; prototype supplies allocation geometry.
        O1VSUpdate(const material::MaterialOwner& owner, double alpha, AJFourFieldsType& prototype)
            : operator_m(prototype, owner, alpha) {
            effective_m.initialize(prototype.get_mesh(), prototype.getLayout());
        }
        /** @brief Write effective sources on the inherited kernel's interior mask.
         * @details All axes are periodic in the admitted formulation. Capturing views
         * by value keeps coefficients alive on the field execution space. Copying J
         * directly when \f$m=1\f$ and its correction is zero preserves signed-zero bits;
         * no arithmetic is added to the inherited vacuum expression.
         */
        void advance(MaterialStepContext& context) override {
            auto op           = operator_m;
            op.potential_m    = context.current().getView();
            const auto source = context.source().getView(), output = effective_m.getView();
            const auto& layout = context.current().getLayout();
            const auto first   = layout.getLocalNDIndex().first();
            const ippl::InteriorBounds bounds(layout.getDomain(),
                                              ippl::FDTDBoundary::fromLayout(layout));
            Kokkos::parallel_for(
                "O1-VS effective source", context.current().getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k) {
                    if (!bounds.contains(i + first[0], j + first[1], k + first[2]))
                        return;
                    for (unsigned c = 0; c < 4; ++c) {
                        const double mass = op.inverseMass(c, i, j, k);
                        const double correction =
                            (mass - 1) * laplacian(op, c, i, j, k)
                            + mass
                                  * (c == 0 ? op.scalarCorrection(i, j, k)
                                            : -op.curlCorrection(c - 1, i, j, k));
                        output(i, j, k)[c] = mass == 1 && correction == 0
                                                 ? source(i, j, k)[c]
                                                 : mass * source(i, j, k)[c] + correction;
                    }
                });
            context.baseStep(context.dt(), effective_m);
        }
        /** @brief Return the vacuum Standard bound for alpha and epsilon at least one.
         * The generalized-gauge stiffness is bounded by the vacuum symbol under
         * these restrictions (deduction from [kreiss2002SecondOrderWave, section 2,
         * text after Eq. (2.9)]). ChdrSolver excludes the stability endpoint.
         */
        double maxStableTimeStep(const ippl::Vector<double, 3>& spacing) const override {
            return chdr::solver::maxStableTimeStep(spacing, Stencil::Standard);
        }
        /// Centred component Laplacian [fallahi2020mithra20fullwavesimulation, section 3.1.2].
        KOKKOS_INLINE_FUNCTION static double laplacian(const O1VSOperator& op, unsigned c, int i,
                                                       int j, int k) {
            double result = 0;
            for (unsigned d = 0; d < 3; ++d) {
                material::CellIndex p(i, j, k), m(i, j, k);
                ++p[d];
                --m[d];
                result += (op.potential_m(p[0], p[1], p[2])[c] - 2 * op.potential_m(i, j, k)[c]
                           + op.potential_m(m[0], m[1], m[2])[c])
                          / (op.spacing_m[d] * op.spacing_m[d]);
            }
            return result;
        }

    private:
        O1VSOperator operator_m;
        AJFourFieldsType effective_m;
    };
}  // namespace chdr::solver
