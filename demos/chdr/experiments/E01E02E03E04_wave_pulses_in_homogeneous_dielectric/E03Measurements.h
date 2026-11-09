/** @file E03Measurements.h
 * @brief Native continuum errors and material Gauss/gauge residuals for E03.
 */
#pragma once

#include <Kokkos_Core.hpp>

#include <Kokkos_MathematicalFunctions.hpp>
#include <cmath>
#include <cstdint>
#include <mpi.h>
#include <stdexcept>

#include "Types/Vector.h"

#include "core/base/LevelContract.h"
#include "core/solver/operators/StaggeredOperators.h"

namespace chdr::experiments::e03 {
    /** @brief E03.meas01--05 absolute norms at one quantity's native time.
     * @details Reference normalization is formed after the MPI reduction. Scalar
     * zero-loading checks use errorLinf_m directly, without dividing by zero.
     */
    struct E03FieldNorms {
        double errorL2_m, errorLinf_m, referenceL2_m, referenceLinf_m;

        double relativeL2() const {
            if (referenceL2_m <= 0)
                throw std::runtime_error("E03 relative L2 error has zero reference norm");
            return errorL2_m / referenceL2_m;
        }
        double relativeLinf() const {
            if (referenceLinf_m <= 0)
                throw std::runtime_error("E03 relative maximum error has zero reference norm");
            return errorLinf_m / referenceLinf_m;
        }
    };

    /// Residual centres in the whole domain, a fixed central cube, or its complement.
    enum class ConstraintRegion {
        All,
        Interior,
        Boundary
    };

    /// Raw E03.meas06 norms, measured centres, and omitted centres with unavailable inputs.
    struct E03ResidualNorms {
        double l2_m, linf_m;
        std::uint64_t count_m, unavailable_m;
    };

    namespace e03MeasurementDetail {
        /// Local sums and maxima, before multiplying squared sums by the cell volume.
        struct FieldReduction {
            double errorSquares_m = 0, referenceSquares_m = 0;
            double errorMaximum_m = 0, referenceMaximum_m = 0;
            std::uint64_t nonfinite_m = 0;
        };

        /// Scalar reference at a native position. Vector components use Cartesian indices.
        template <class Reference>
        KOKKOS_INLINE_FUNCTION double referenceValue(const Reference& reference, Quantity quantity,
                                                     double x, double y, double z, double time,
                                                     unsigned component) {
            if (quantity == Quantity::Phi)
                return reference.referenceFourPotentialA(x, y, z, time, 0);
            if (quantity == Quantity::VectorPotential)
                return reference.referenceFourPotentialA(x, y, z, time, component + 1);
            if (quantity == Quantity::Electric)
                return reference.referenceE(x, y, z, time, component);
            return reference.referenceB(x, y, z, time, component);
        }

        /// Reject configurations whose physical placement differs from this experiment.
        inline void requireE03Contract(LevelContract contract) {
            if (contract.spatialDiscretization_m != SpatialDiscretization::Staggered
                || contract.ebReconstructionConvention_m != EBReconstructionConvention::R2)
                throw std::invalid_argument("E03 measurements require staggered R2 fields");
        }

        /** @brief Reduce field errors over owned samples at the contract's positions.
         * @details No field is interpolated or copied to host. The Kokkos reduction
         * executes in the field's execution space. Native offsets follow
         * [ryu2016potentialfdtd, Appendix A, Eqs. (19)--(24)] and LevelContract.
         */
        template <class Field, class Reference>
        FieldReduction localFieldNorms(Field& field, Reference reference, LevelContract contract,
                                       Quantity quantity, double time) {
            const auto view = field.getView();
            const auto h = field.get_mesh().getMeshSpacing(), origin = field.get_mesh().getOrigin();
            const auto first              = field.getLayout().getLocalNDIndex().first();
            const int ghost               = field.getNghost();
            const unsigned count          = quantity == Quantity::Phi ? 1 : 3;
            const unsigned firstComponent = quantity == Quantity::VectorPotential ? 1 : 0;
            FieldReduction local;
            Kokkos::parallel_reduce(
                "E03 native continuum errors", field.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k, double& errorSquares, double& referenceSquares,
                              double& errorMaximum, double& referenceMaximum,
                              std::uint64_t& nonfinite) {
                    const double x = origin[0] + (i + first[0] - ghost + 0.5) * h[0];
                    const double y = origin[1] + (j + first[1] - ghost + 0.5) * h[1];
                    const double z = origin[2] + (k + first[2] - ghost + 0.5) * h[2];
                    for (unsigned c = 0; c < count; ++c) {
                        const auto offset = contract.positionOffset(quantity, c);
                        const double exact =
                            referenceValue(reference, quantity, x + offset[0] * h[0],
                                           y + offset[1] * h[1], z + offset[2] * h[2], time, c);
                        const double error = view(i, j, k)[c + firstComponent] - exact;
                        if (!Kokkos::isfinite(exact) || !Kokkos::isfinite(error)) {
                            ++nonfinite;
                            continue;
                        }
                        errorSquares += error * error;
                        referenceSquares += exact * exact;
                        errorMaximum     = Kokkos::max(errorMaximum, Kokkos::abs(error));
                        referenceMaximum = Kokkos::max(referenceMaximum, Kokkos::abs(exact));
                    }
                },
                local.errorSquares_m, local.referenceSquares_m,
                Kokkos::Max<double>(local.errorMaximum_m),
                Kokkos::Max<double>(local.referenceMaximum_m), local.nonfinite_m);
            return local;
        }

        /** @brief Combine absolute norms before forming reference-normalized errors.
         * @details The volume-weighted norm and component maximum follow
         * [salari2000CodeVerificationMms, section 4.1, Eqs. (15), (18)].
         * E03 defines the additional division by its sampled reference norm.
         */
        inline E03FieldNorms globalFieldNorms(FieldReduction local, double cellVolume,
                                              MPI_Comm communicator) {
            double squares[]{local.errorSquares_m, local.referenceSquares_m};
            double peaks[]{local.errorMaximum_m, local.referenceMaximum_m};
            MPI_Allreduce(MPI_IN_PLACE, squares, 2, MPI_DOUBLE, MPI_SUM, communicator);
            MPI_Allreduce(MPI_IN_PLACE, peaks, 2, MPI_DOUBLE, MPI_MAX, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &local.nonfinite_m, 1, MPI_UINT64_T, MPI_SUM, communicator);
            if (local.nonfinite_m || !std::isfinite(squares[0]) || !std::isfinite(squares[1])
                || !std::isfinite(peaks[0]) || !std::isfinite(peaks[1]))
                throw std::runtime_error("E03 field measurement contains nonfinite data");
            return {std::sqrt(cellVolume * squares[0]), peaks[0],
                    std::sqrt(cellVolume * squares[1]), peaks[1]};
        }

        /// Select the fixed physical region using residual cell-centre coordinates.
        KOKKOS_INLINE_FUNCTION bool inRegion(ConstraintRegion region, double x, double y, double z,
                                             double interiorHalfWidth) {
            const bool interior = Kokkos::abs(x) <= interiorHalfWidth
                                  && Kokkos::abs(y) <= interiorHalfWidth
                                  && Kokkos::abs(z) <= interiorHalfWidth;
            return region == ConstraintRegion::All
                   || (region == ConstraintRegion::Interior ? interior : !interior);
        }

        /** @brief Backward face-to-centre divergence of an E or vector-potential view.
         * @details Uses the solver's staggeredDifference and component accessors
         * [ryu2016potentialfdtd, Appendix A, Eq. (23)]. Inputs must have valid halos.
         */
        template <class View>
        KOKKOS_INLINE_FUNCTION double divergence(View view, int i, int j, int k,
                                                 ippl::Vector<double, 3> h,
                                                 unsigned firstComponent) {
            double result = 0;
            for (unsigned d = 0; d < 3; ++d)
                result += solver::staggeredDifference(
                    solver::Component<View>{view, firstComponent + d}, i, j, k, d, h[d], false);
            return result;
        }

        /// Raw local residual statistics. Counts distinguish absent samples from zero values.
        struct ResidualReduction {
            double squares_m = 0, maximum_m = 0;
            std::uint64_t count_m = 0, unavailable_m = 0, nonfinite_m = 0;
        };

        /** @brief Evaluate one material constraint from numerical histories.
         * @details Gauss is epsilon times the backward divergence of R2 E. Gauge is
         * the backward divergence of current A plus epsilon times the scalar history
         * difference [chew2014generalizedgauge, Eqs. (4), (11), and
         * christlieb2024GaugeConserving, section 3.1, Eqs. (25)--(27)].
         * Gauss omits lower open-face centres requiring exterior E. Gauge retains
         * centres supplied by the solver's recorded potential-ghost closure.
         */
        template <bool Gauss, class Solver>
        ResidualReduction localConstraint(Solver& state, double epsilon, ConstraintRegion region,
                                          double interiorHalfWidth) {
            const auto current = state.A_n.getView(), previous = state.A_nm1.getView();
            const auto electric = state.electricField().getView();
            const auto h        = state.A_n.get_mesh().getMeshSpacing();
            const auto& layout  = state.A_n.getLayout();
            const auto first = layout.getLocalNDIndex().first(), lower = layout.getDomain().first();
            const auto& domain = layout.getDomain();
            const ippl::Vector<double, 3> middle{lower[0] + 0.5 * domain[0].length(),
                                                 lower[1] + 0.5 * domain[1].length(),
                                                 lower[2] + 0.5 * domain[2].length()};
            const ippl::Vector<int, 3> periodic{layout.isPeriodic(0), layout.isPeriodic(1),
                                                layout.isPeriodic(2)};
            const int ghost        = state.A_n.getNghost();
            const double inverseDt = 1.0 / state.getDt();
            ResidualReduction local;
            Kokkos::parallel_reduce(
                "E03 material constraint", state.A_n.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k, double& squares, double& maximum,
                              std::uint64_t& count, std::uint64_t& unavailable,
                              std::uint64_t& nonfinite) {
                    const ippl::Vector<int, 3> global{i + first[0] - ghost, j + first[1] - ghost,
                                                      k + first[2] - ghost};
                    if (!inRegion(region, (global[0] + 0.5 - middle[0]) * h[0],
                                  (global[1] + 0.5 - middle[1]) * h[1],
                                  (global[2] + 0.5 - middle[2]) * h[2], interiorHalfWidth))
                        return;
                    if constexpr (Gauss) {
                        for (unsigned d = 0; d < 3; ++d)
                            if (!periodic[d] && global[d] == lower[d]) {
                                ++unavailable;
                                return;
                            }
                    }
                    const double residual =
                        Gauss ? epsilon * divergence(electric, i, j, k, h, 0)
                              : divergence(current, i, j, k, h, 1)
                                    + epsilon * (current(i, j, k)[0] - previous(i, j, k)[0])
                                          * inverseDt;
                    ++count;
                    if (!Kokkos::isfinite(residual)) {
                        ++nonfinite;
                        return;
                    }
                    squares += residual * residual;
                    maximum = Kokkos::max(maximum, Kokkos::abs(residual));
                },
                local.squares_m, Kokkos::Max<double>(local.maximum_m), local.count_m,
                local.unavailable_m, local.nonfinite_m);
            return local;
        }

        /// Complete one residual's collective norm and retain omitted-centre counts.
        inline E03ResidualNorms globalConstraint(ResidualReduction local, double volume,
                                                 MPI_Comm communicator) {
            std::uint64_t counts[]{local.count_m, local.unavailable_m, local.nonfinite_m};
            MPI_Allreduce(MPI_IN_PLACE, &local.squares_m, 1, MPI_DOUBLE, MPI_SUM, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &local.maximum_m, 1, MPI_DOUBLE, MPI_MAX, communicator);
            MPI_Allreduce(MPI_IN_PLACE, counts, 3, MPI_UINT64_T, MPI_SUM, communicator);
            if (counts[2] || !std::isfinite(local.squares_m) || !std::isfinite(local.maximum_m))
                throw std::runtime_error("E03 constraint measurement contains nonfinite data");
            return {std::sqrt(volume * local.squares_m), counts[0] ? local.maximum_m : 0, counts[0],
                    counts[1]};
        }
    }  // namespace e03MeasurementDetail

    /** @brief E03.meas01--05 native vector/scalar error and reference norms.
     * @param field E, R2 B or four-potential storage appropriate to quantity.
     * @param time Actual solver labelTime for this quantity and observed state.
     * @details Phi reads slot zero. A reads slots one through three. The continuum
     * evaluator receives LevelContract's component positions. All ranks participate.
     */
    template <class Field, class Reference>
    E03FieldNorms measureE03Field(Field& field, const Reference& reference, LevelContract contract,
                                  Quantity quantity, double time) {
        e03MeasurementDetail::requireE03Contract(contract);
        if (quantity != Quantity::Phi && quantity != Quantity::VectorPotential
            && quantity != Quantity::Electric && quantity != Quantity::Magnetic)
            throw std::invalid_argument("Unsupported E03 measured quantity");
        const auto h = field.get_mesh().getMeshSpacing();
        const auto local =
            e03MeasurementDetail::localFieldNorms(field, reference, contract, quantity, time);
        return e03MeasurementDetail::globalFieldNorms(local, h[0] * h[1] * h[2],
                                                      field.getLayout().comm);
    }

    /** @brief E03.meas06 electric Gauss-law residual at the R2 half-time.
     * @details Exchange reconstructed E halos, then reduce unique owned residual centres.
     * Exterior E samples are unavailable and their dependent centres are counted separately.
     * epsilon is the homogeneous relative permittivity, with zero free charge.
     */
    template <class Solver>
    E03ResidualNorms measureE03Gauss(Solver& state, double epsilon,
                                     ConstraintRegion region  = ConstraintRegion::All,
                                     double interiorHalfWidth = 5.0) {
        e03MeasurementDetail::requireE03Contract(state.levelContract());
        state.electricField().fillHalo();
        const auto h = state.A_n.get_mesh().getMeshSpacing();
        const auto local =
            e03MeasurementDetail::localConstraint<true>(state, epsilon, region, interiorHalfWidth);
        return e03MeasurementDetail::globalConstraint(local, h[0] * h[1] * h[2],
                                                      state.A_n.getLayout().comm);
    }

    /** @brief E03.meas06 Lorenz-gauge residual at the current vector-potential time.
     * @details prepareHistories exchanges both potential slots and applies the solver's
     * physical ghost closure. Its scalar difference is centred at the current A time
     * by the staggered R2 contract [christlieb2024GaugeConserving, Eqs. (25)--(27)].
     */
    template <class Solver>
    E03ResidualNorms measureE03Gauge(Solver& state, double epsilon,
                                     ConstraintRegion region  = ConstraintRegion::All,
                                     double interiorHalfWidth = 5.0) {
        e03MeasurementDetail::requireE03Contract(state.levelContract());
        state.prepareHistories();
        const auto h = state.A_n.get_mesh().getMeshSpacing();
        const auto local =
            e03MeasurementDetail::localConstraint<false>(state, epsilon, region, interiorHalfWidth);
        return e03MeasurementDetail::globalConstraint(local, h[0] * h[1] * h[2],
                                                      state.A_n.getLayout().comm);
    }
}  // namespace chdr::experiments::e03
