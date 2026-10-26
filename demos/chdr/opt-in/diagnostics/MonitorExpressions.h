/** @file MonitorExpressions.h
 * @brief Device expression values with their absolute evaluations.
 * @details For the linear difference expressions used here, the absolute
 * evaluation replaces every subtraction by addition of input magnitudes and
 * every coefficient by its magnitude. Applying the absolute outer stencil to
 * the inner absolute evaluation keeps cancellations out of the scale while retaining
 * the original residual value. This is a deduction from the absolute error
 * scale in [jeannerod2013ImprovedErrorBounds, Eq. (1.6)]. Theorem 4.2 concerns
 * the paper's specified dot-product evaluation; it is not an error certificate
 * for these nested expressions when the compiler contracts operations to FMA.
 * Here the absolute evaluation defines a cancellation-free diagnostic scale.
 * Expression carries both values; Workspace exchanges both intermediate halos;
 * ResidualMetrics defines the resulting cellwise and global ratios.
 *
 * Through each calculation, two numbers are carried:
 * 1. value_m: the calculated result
 * 2. scale_m: a magnitude scale for interpreting the residual
 * Subtraction subtracts values but adds scales.
 */
#pragma once

#include <cmath>
#include <limits>

#include "core/base/FieldContainer.h"
#include "core/base/SpatialDimension.h"
#include "core/operators/Operators.h"
#include "core/operators/StaggeredOperators.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "MonitorExpressions.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    namespace monitorDetail {
        /// Carry the computed expression and its cancellation-free absolute evaluation.
        struct Expression {
            double value_m = 0.0, scale_m = 0.0;
            KOKKOS_INLINE_FUNCTION Expression() = default;
            KOKKOS_INLINE_FUNCTION Expression(double value, double scale)
                : value_m(value)
                , scale_m(scale) {}
            KOKKOS_INLINE_FUNCTION explicit Expression(double value)
                : value_m(value)
                , scale_m(Kokkos::abs(value)) {}
            KOKKOS_INLINE_FUNCTION Expression operator+(Expression other) const {
                return {value_m + other.value_m, scale_m + other.scale_m};
            }
            KOKKOS_INLINE_FUNCTION Expression operator-(Expression other) const {
                return {value_m - other.value_m, scale_m + other.scale_m};
            }
            KOKKOS_INLINE_FUNCTION Expression operator*(double factor) const {
                return {value_m * factor, scale_m * Kokkos::abs(factor)};
            }
        };

        /** @brief Read one field component together with its input magnitude.
         * @details The lookup matches chdr::solver::Component in Operators.h.
         * That accessor returns double; this one returns Expression, initializing
         * both the value and its absolute-magnitude scale for monitor calculations.
         * @tparam FieldView Three-dimensional view of component-indexed values, such as
         * potential, source or E/B storage. This accessor does not change view elements.
         */
        template <class FieldView>
        struct Component {
            FieldView view_m;
            unsigned component_m;
            KOKKOS_INLINE_FUNCTION Expression operator()(int i, int j, int k) const {
                return Expression(view_m(i, j, k)[component_m]);
            }
        };

        /** @brief Bind a component without assigning a physical meaning to its storage.
         * @tparam FieldView Three-dimensional view of component-indexed field values.
         * @param view Field samples, indexed by (i,j,k) and then by component.
         * @param index Component index within each stored value.
         */
        template <class FieldView>
        Component<FieldView> component(FieldView view, unsigned index) {
            return {view, index};
        }

        /// Select the value or absolute evaluation without a per-sample runtime choice.
        template <class ExpressionAccessor, bool Absolute>
        struct Part {
            ExpressionAccessor sample_m;
            KOKKOS_INLINE_FUNCTION double operator()(int i, int j, int k) const {
                const auto result = sample_m(i, j, k);
                if constexpr (Absolute)
                    return result.scale_m;
                else
                    return result.value_m;
            }
        };

        /** @brief Evaluate one spatial derivative and its absolute evaluation.
         * @tparam ExpressionAccessor Device-callable sample(i,j,k) returning Expression.
         * @param sample Scalar expression at local array indices, including halo offsets.
         * @param i First local spatial index.
         * @param j Second local spatial index.
         * @param k Third local spatial index.
         * @param axis Cartesian derivative direction, in [0, SpatialDim).
         * @param spacing Positive mesh spacing along axis.
         * @param direction Zero selects centered differences; positive selects forward
         * and negative selects backward differences. See centeredDifference and
         * staggeredDifference for the stencils and their literature sources.
         * @pre The selected stencil's input samples, including halos, are valid.
         */
        template <class ExpressionAccessor>
        KOKKOS_INLINE_FUNCTION Expression derivative(ExpressionAccessor sample, int i, int j, int k,
                                                     unsigned axis, double spacing,
                                                     int direction = 0) {
            if (direction != 0)
                return {staggeredDifference(Part<ExpressionAccessor, false>{sample}, i, j, k, axis,
                                            spacing, direction > 0),
                        staggeredDifference<true>(Part<ExpressionAccessor, true>{sample}, i, j, k,
                                                  axis, spacing, direction > 0)};
            return {
                centeredDifference(Part<ExpressionAccessor, false>{sample}, i, j, k, axis, spacing),
                centeredDifference<true>(Part<ExpressionAccessor, true>{sample}, i, j, k, axis,
                                         spacing)};
        }

        /// Retain the solver's stencil choice in both device-side expression evaluations.
        template <Stencil Scheme, class ExpressionAccessor>
        KOKKOS_INLINE_FUNCTION Expression laplace(ExpressionAccessor sample, int i, int j, int k,
                                                  ippl::Vector<double, SpatialDim> spacing) {
            return {
                laplacian<Scheme>(Part<ExpressionAccessor, false>{sample}, i, j, k, spacing),
                laplacian<Scheme, true>(Part<ExpressionAccessor, true>{sample}, i, j, k, spacing)};
        }

        /** @brief Sum the Cartesian derivatives of a stored vector's components.
         * @tparam FieldView Three-dimensional view of component-indexed field values.
         * @param view Three-dimensional view of component-indexed field values.
         * @param firstVectorComponent Index of the x component: one for four-component
         * potential/source storage, zero for E/B storage.
         * @param i First local spatial index, including the halo offset.
         * @param j Second local spatial index, including the halo offset.
         * @param k Third local spatial index, including the halo offset.
         * @param spacing Positive Cartesian mesh spacings.
         * @param direction Zero selects centered differences; positive selects forward
         * and negative selects backward differences.
         * @details The absolute evaluation propagates through the sum as described
         * by Expression; the selected spatial stencil is documented in Operators.h
         * and StaggeredOperators.h.
         */
        template <class FieldView>
        KOKKOS_INLINE_FUNCTION Expression divergence(FieldView view, unsigned firstVectorComponent,
                                                     int i, int j, int k,
                                                     ippl::Vector<double, SpatialDim> spacing,
                                                     int direction = 0) {
            Expression result;
            for (unsigned d = 0; d < SpatialDim; ++d)
                result = result
                         + derivative(Component<FieldView>{view, firstVectorComponent + d}, i, j, k,
                                      d, spacing[d], direction);
            return result;
        }

        /** @brief Evaluate one Cartesian curl component with its absolute evaluation.
         * @tparam FieldView Three-dimensional view of component-indexed field values.
         * @param view Three-dimensional view of component-indexed field values.
         * @param firstVectorComponent Index of the x component, as in divergence.
         * @param d Output curl component, in [0, SpatialDim).
         * @param i First local spatial index, including the halo offset.
         * @param j Second local spatial index, including the halo offset.
         * @param k Third local spatial index, including the halo offset.
         * @param spacing Positive Cartesian mesh spacings.
         * @param direction Zero selects centered differences; positive selects forward
         * and negative selects backward differences.
         * @details The cyclic component ordering follows the Cartesian curl
         * [schneider2010UnderstandingFdtd, section 2.5, Eq. (2.24), p. 16].
         * Centered or staggered differences are selected without changing that ordering;
         * their stencils and sources are documented in the operator helpers.
         */
        template <class FieldView>
        KOKKOS_INLINE_FUNCTION Expression curl(FieldView view, unsigned firstVectorComponent,
                                               unsigned d, int i, int j, int k,
                                               ippl::Vector<double, SpatialDim> spacing,
                                               int direction = 0) {
            const unsigned a = (d + 1) % SpatialDim, b = (d + 2) % SpatialDim;
            return derivative(Component<FieldView>{view, firstVectorComponent + b}, i, j, k, a,
                              spacing[a], direction)
                   - derivative(Component<FieldView>{view, firstVectorComponent + a}, i, j, k, b,
                                spacing[b], direction);
        }

        /** @brief Read a cached expression and its separately stored absolute evaluation.
         * @tparam FieldView Three-dimensional view of matching component-indexed arrays.
         * @details value_m and scale_m have identical layouts; component_m selects
         * the same cached residual component from both arrays.
         */
        template <class FieldView>
        struct StoredComponent {
            FieldView value_m, scale_m;
            unsigned component_m;
            KOKKOS_INLINE_FUNCTION Expression operator()(int i, int j, int k) const {
                return {value_m(i, j, k)[component_m], scale_m(i, j, k)[component_m]};
            }
        };

        /** @brief Reusable mesh storage for intermediate monitor values and their scales.
         * @details Two one-halo fields store the value and scale of each Expression
         * separately. Construction allocates these fields using the prototype's geometry;
         * it does not copy the prototype's values. The monitor retains this workspace
         * and reuses its allocations for successive intermediate calculations.
         *
         * fill() evaluates the supplied formula on owned cells, stores both parts,
         * waits for the writes to finish, then exchanges both fields' halos.
         * component() exposes a stored pair as an Expression accessor. A subsequent
         * derivative can therefore read neighbouring intermediate values and scales,
         * including those owned by another rank.
         *
         * For example, the magnetic monitor stores a raw curl and its scale here,
         * then evaluates a divergence residual and its propagated scale after the halo exchange.
         */
        struct Workspace {
            AJFourFieldsType value_m, scale_m;
            /** @brief Allocate residual storage using an existing field's geometry.
             * @param prototype Field whose mesh/layout owners outlive the workspace.
             * Its values are not read; the new fields retain pointers to its mesh and layout.
             */
            explicit Workspace(const AJFourFieldsType& prototype) {
                value_m.initialize(prototype.get_mesh(), prototype.getLayout());
                scale_m.initialize(prototype.get_mesh(), prototype.getLayout());
                value_m = 0.0;
                scale_m = 0.0;
            }
            auto component(unsigned index) const {
                return StoredComponent{value_m.getView(), scale_m.getView(), index};
            }
            /** @brief Evaluate owned components, then exchange value and scale halos.
             * @tparam ExpressionFormula Device-copyable formula(i,j,k,component) returning
             * Expression for each selected component.
             * @param formula Expression evaluated at local indices including halo offsets.
             * @param components Number of leading components to write, at most four.
             */
            template <class ExpressionFormula>
            void fill(ExpressionFormula formula, unsigned components = SpatialDim) {
                auto value = value_m.getView(), scale = scale_m.getView();
                Kokkos::parallel_for(
                    "ChDR monitor intermediate", value_m.getFieldRangePolicy(),
                    KOKKOS_LAMBDA(int i, int j, int k) {
                        for (unsigned d = 0; d < components; ++d) {
                            const auto result = formula(i, j, k, d);
                            value(i, j, k)[d] = result.value_m;
                            scale(i, j, k)[d] = result.scale_m;
                        }
                    });
                Kokkos::fence();
                value_m.fillHalo();
                scale_m.fillHalo();
            }
        };
    }  // namespace monitorDetail

    /** @brief Owned-cell residual maxima and unnormalized global Euclidean norms.
     * @details \f$r_{\mathrm{cell}}=\max |R|/S\f$ and \f$r_{\mathrm{global}}=\max|R|/\max S\f$;
     * \f$0/0=0\f$. Scale is the absolute evaluation, including inner stencils.
     * Reference scale: [jeannerod2013ImprovedErrorBounds, Eq. (1.6)]; the file
     * documentation distinguishes this diagnostic from a formal error bound.
     */
    struct ResidualMetrics {
        double cell_m = 0.0, global_m = 0.0, maxAbs_m = 0.0, maxScale_m = 0.0;
        double l2_m = 0.0, scaleL2_m = 0.0;
    };

    namespace monitorDetail {
        /// Local accumulators mirror the five Kokkos reductions, before any MPI call.
        struct LocalResiduals {
            double cell_m = 0.0, absolute_m = 0.0, scale_m = 0.0;
            double square_m = 0.0, scaleSquare_m = 0.0;
        };

        /** @brief Accumulate one owned cell without discarding a nonfinite residual.
         * @details Each maximum and sum uses the same expression evaluation.
         * Ghost cells never enter the range policy. The value and absolute scale
         * travel together, so a cancellation inside a derivative does not make
         * the denominator spuriously small; see Expression and ResidualMetrics.
         */
        template <class ExpressionFormula>
        struct ResidualReduction {
            ExpressionFormula formula_m;
            static constexpr double INFINITY_VALUE = std::numeric_limits<double>::infinity();
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& maxCell,
                                                   double& maxAbs, double& maxScale,
                                                   double& sumSquare,
                                                   double& sumScaleSquare) const {
                const auto result = formula_m(i, j, k);
                const bool finite =
                    Kokkos::isfinite(result.value_m) && Kokkos::isfinite(result.scale_m);
                const double magnitude = finite ? Kokkos::abs(result.value_m) : INFINITY_VALUE;
                const double bound     = finite ? result.scale_m : INFINITY_VALUE;
                const double ratio =
                    !finite ? INFINITY_VALUE
                            : (bound > 0.0 ? magnitude / bound
                                           : (magnitude == 0.0 ? 0.0 : INFINITY_VALUE));
                maxCell  = Kokkos::max(maxCell, ratio);
                maxAbs   = Kokkos::max(maxAbs, magnitude);
                maxScale = Kokkos::max(maxScale, bound);
                sumSquare += magnitude * magnitude;
                sumScaleSquare += bound * bound;
            }
        };

        /// Keep the device reduction in the field's execution space; no host mirror is made.
        template <class FieldType, class ExpressionFormula>
        LocalResiduals reduceOwned(const FieldType& field, ExpressionFormula formula) {
            LocalResiduals local;
            Kokkos::parallel_reduce(
                "ChDR residual reduction", field.getFieldRangePolicy(),
                ResidualReduction<ExpressionFormula>{formula}, Kokkos::Max<double>(local.cell_m),
                Kokkos::Max<double>(local.absolute_m), Kokkos::Max<double>(local.scale_m),
                Kokkos::Sum<double>(local.square_m), Kokkos::Sum<double>(local.scaleSquare_m));
            return local;
        }

        /** @brief Form global diagnostics on the field layout's communicator.
         * @details Maxima use MPI_MAX; squared norms use MPI_SUM. Communication
         * contains only plain doubles and follows the same collective order on
         * every participating rank, including an MPI_COMM_SELF witness.
         * Square roots are applied only after the global sums.
         */
        template <class FieldType>
        ResidualMetrics globalizeResidual(const FieldType& field, LocalResiduals local) {
            double localMax[3]{local.cell_m, local.absolute_m, local.scale_m}, globalMax[3];
            double localSum[2]{local.square_m, local.scaleSquare_m}, globalSum[2];
            MPI_Allreduce(localMax, globalMax, 3, MPI_DOUBLE, MPI_MAX, field.getLayout().comm);
            MPI_Allreduce(localSum, globalSum, 2, MPI_DOUBLE, MPI_SUM, field.getLayout().comm);
            const double ratio =
                std::isfinite(globalMax[1]) && std::isfinite(globalMax[2])
                    ? (globalMax[2] > 0.0 ? globalMax[1] / globalMax[2] : globalMax[0])
                    : std::numeric_limits<double>::infinity();
            return {
                globalMax[0],           ratio, globalMax[1], globalMax[2], std::sqrt(globalSum[0]),
                std::sqrt(globalSum[1])};
        }
    }  // namespace monitorDetail

    /** @brief Evaluate a residual on owned cells, then reduce on its layout.
     * @tparam FieldType Field supplying geometry and communication metadata; its stored
     * values need not be the values used by the residual.
     * @tparam ExpressionFormula Device-copyable formula(i,j,k) returning Expression.
     * @param field Supplies the owned-cell iteration range and layout communicator.
     * The formula determines which field data are read; this field is not modified.
     * @param formula Device-callable expression returning both value and absolute scale.
     * @pre All ranks belonging to field.getLayout().comm call in matching order.
     * @return Cellwise and global roundoff ratios and unnormalized Euclidean norms.
     * The normalization is the absolute evaluation, not the norm of the
     * already-cancelled top-level terms [jeannerod2013ImprovedErrorBounds,
     * Eq. (1.6)]. Propagating this scale through nested derivatives
     * is the deduction described in the file documentation.
     */
    template <class FieldType, class ExpressionFormula>
    ResidualMetrics reduceResidual(const FieldType& field, ExpressionFormula formula) {
        return monitorDetail::globalizeResidual(field, monitorDetail::reduceOwned(field, formula));
    }
}  // namespace chdr::solver
