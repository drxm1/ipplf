/** @file OwnedStatistics.h
 * @brief Stable owned-cell maxima and nonnegative sums for field diagnostics.
 *
 * Collects three rank-local summaries of a field:
 * 1. largest absolute value,
 * 2. sum of absolute values,
 * 3. sum of squared values
 * So this file provides statistics that are computed from rank-owned field cells.
 */
#pragma once

#include "Ippl.h"

#include <array>
#include <limits>
#include <utility>

#include "core/base/SpatialDimension.h"

static_assert(
    chdr::SpatialDim == 3,
    "OwnedStatistics.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    struct OwnedStatistics {
        double maxAbs_m, sumAbs_m, sumSquares_m;
        std::array<double, 3> values() const { return {maxAbs_m, sumAbs_m, sumSquares_m}; }
    };

    namespace statisticsDetail {
        /// Two statistics per owned cell, stored in the input field's device space.
        template <class FieldType>
        using PairView = Kokkos::View<double* [2], typename FieldType::view_type::device_type>;

        /** @brief Store absolute-value and square sums for each owned cell.
         * @tparam FieldView Three-dimensional view of component-indexable cell values.
         * @tparam StatisticPairs Two-column view: column 0 is the absolute sum;
         * column 1 is the sum of squares. Rows follow owned x/y/z index order.
         * @tparam COMPONENTS Number of stored components to read from each cell.
         */
        template <class FieldView, class StatisticPairs, unsigned COMPONENTS>
        struct CellStatistics {
            FieldView view_m;
            StatisticPairs pairs_m;
            std::size_t ny_m, nz_m;
            int ghost_m;
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k) const {
                const auto cell =
                    (std::size_t(i - ghost_m) * ny_m + j - ghost_m) * nz_m + k - ghost_m;
                double absolute = 0.0, squares = 0.0;
                for (unsigned c = 0; c < COMPONENTS; ++c) {
                    const double value = view_m(i, j, k)[c];
                    absolute += Kokkos::abs(value);
                    squares += value * value;
                }
                pairs_m(cell, 0) = absolute;
                pairs_m(cell, 1) = squares;
            }
        };

        /** @brief One leaf contains the two positive statistics of one owned cell.
         * @details Components retain their original order. Keeping each cell's at
         * most four terms together adds only a bounded local summation depth; the
         * subsequent tree handles the potentially millions of cells. Ghosts are
         * excluded by the field's owned range policy, not by their stored values.
         */
        template <class FieldType>
        PairView<FieldType> ownedPairs(const FieldType& field) {
            const auto local     = field.getLayout().getLocalNDIndex();
            const std::size_t ny = local[1].length(), nz = local[2].length();
            PairView<FieldType> pairs("owned statistic pairs", local[0].length() * ny * nz);
            const auto view               = field.getView();
            constexpr unsigned COMPONENTS = FieldType::value_type::dim;
            using Kernel = CellStatistics<decltype(view), PairView<FieldType>, COMPONENTS>;
            Kokkos::parallel_for("initialize owned statistic pairs", field.getFieldRangePolicy(),
                                 Kernel{view, pairs, ny, nz, field.getNghost()});
            return pairs;
        }

        /** @brief Sum adjacent rows; copy an unpaired final row unchanged.
         * @tparam StatisticPairs Two-column device view described by CellStatistics.
         * @param input Previous tree level; read the first count rows.
         * @param output Distinct destination with space for the next tree level.
         * @param count Number of valid input rows.
         * @param execution Execution instance ordering this level after earlier work.
         */
        template <class StatisticPairs>
        void pairwiseLevel(StatisticPairs input, StatisticPairs output, std::size_t count,
                           const typename StatisticPairs::execution_space& execution) {
            using Policy = Kokkos::RangePolicy<typename StatisticPairs::execution_space,
                                               Kokkos::IndexType<std::size_t>>;
            Kokkos::parallel_for(
                "pairwise statistic level", Policy(execution, 0, (count + 1) / 2),
                KOKKOS_LAMBDA(std::size_t pair) {
                    const auto left = 2 * pair, right = left + 1;
                    for (unsigned statistic = 0; statistic < 2; ++statistic)
                        output(pair, statistic) =
                            right < count ? input(left, statistic) + input(right, statistic)
                                          : input(left, statistic);
                });
        }

        /** @brief Sum leaves in an explicit balanced tree using two device buffers.
         * @details Jeannerod and Rump [jeannerod2013ImprovedErrorBounds, p. 7]
         * give the pairwise inner-product factor \f$\gamma_\ell\f$ with
         * \f$\ell=\lceil\log_2 n\rceil+1\f$; their model is Sec. 1, Eqs. (1.1),(1.3).
         * Here \f$n\f$ is the nonempty inner-product length, \f$u\f$ is unit roundoff, and
         * \f$\gamma_\ell=\ell u/(1-\ell u)\f$ requires \f$\ell u<1\f$.
         * Applying that depth argument to our fixed-size cell leaves adds a constant
         * local component depth (deduction). The path length therefore grows
         * logarithmically in cell count, rather than with a sequential field traversal.
         * The error-bound argument assumes round-to-nearest arithmetic without
         * overflow or underflow, as in the cited model. It concerns the local tree;
         * the final communicator reduction is performed separately by ownedStatistics().
         * Every level uses the same execution instance. Only the final pair is copied
         * to the host, after its kernels have completed; the fields are never moved.
         * The input is a disposable leaf buffer and may be overwritten by later levels.
         */
        template <class StatisticPairs>
        std::array<double, 2> pairwiseSums(
            StatisticPairs input, const typename StatisticPairs::execution_space& execution = {}) {
            std::size_t count = input.extent(0);
            if (count == 0)
                return {0.0, 0.0};
            StatisticPairs output("pairwise statistic scratch", (count + 1) / 2);
            while (count > 1) {
                pairwiseLevel(input, output, count, execution);
                std::swap(input, output);
                count = (count + 1) / 2;
            }
            execution.fence("statistic pair ready for host");
            std::array<double, 2> result{};
            Kokkos::deep_copy(result[0], Kokkos::subview(input, 0, 0));
            Kokkos::deep_copy(result[1], Kokkos::subview(input, 0, 1));
            return result;
        }

        /// Preserve the maximum reduction independently of the new summation tree.
        template <class FieldType>
        double localMaximum(const FieldType& field) {
            const auto view               = field.getView();
            constexpr unsigned COMPONENTS = FieldType::value_type::dim;
            double maximum                = 0.0;
            const double infinity         = std::numeric_limits<double>::infinity();
            Kokkos::parallel_reduce(
                "owned statistic maximum", field.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k, double& result) {
                    for (unsigned c = 0; c < COMPONENTS; ++c) {
                        const double value = view(i, j, k)[c];
                        result             = Kokkos::fmax(
                            result, Kokkos::isfinite(value) ? Kokkos::abs(value) : infinity);
                    }
                },
                Kokkos::Max<double>(maximum));
            return maximum;
        }
    }  // namespace statisticsDetail

    /** @brief Compute the global absolute maximum, absolute sum and sum of squares.
     * @pre Collective on the field communicator; a three-dimensional vector field
     * with at most four components and valid owned values.
     * @details Device leaves and binary levels use the field's execution/memory
     * space. The final MPI operations remain SUM for the two positive sums and
     * MAX for the absolute maximum, on plain doubles. These are raw array
     * statistics: there is no mesh-volume multiplier or sample-count normalization.
     * Gathering individual field values is a separate operation, so comparisons
     * can detect errors in the actual distributed reductions independently.
     * @tparam FieldType IPPL field with component-indexable values; this helper
     * accepts potential, source and electromagnetic fields without assigning time labels.
     * @param field Read-only field whose owned values are measured on every rank.
     */
    template <class FieldType>
    OwnedStatistics ownedStatistics(const FieldType& field) {
        static_assert(FieldType::dim == SpatialDim && FieldType::value_type::dim <= 4);
        const auto execution = field.getFieldRangePolicy().space();
        const auto sums =
            statisticsDetail::pairwiseSums(statisticsDetail::ownedPairs(field), execution);
        std::array<double, 3> values{statisticsDetail::localMaximum(field), sums[0], sums[1]};
        MPI_Allreduce(MPI_IN_PLACE, values.data() + 1, 2, MPI_DOUBLE, MPI_SUM,
                      field.getLayout().comm);
        MPI_Allreduce(MPI_IN_PLACE, values.data(), 1, MPI_DOUBLE, MPI_MAX, field.getLayout().comm);
        return {values[0], values[1], values[2]};
    }
}  // namespace chdr::solver::test
