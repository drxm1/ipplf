/** @file E01E02Measurements.h
 * @brief Global numerical norms and continuum errors for the shared transverse pulse.
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
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01AnalyticContinuumReference.h"

namespace chdr::experiments::e01e02 {
    namespace e01Reference = chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::
        E01_analytic_continuum_reference;
    using E01Reference = e01Reference::E01AnalyticContinuumReferenceSampler;

    /// Reference-normalized errors of all three Cartesian components over the whole domain.
    struct PulseErrorNorms {
        double relativeL2_m;
        double relativeLinf_m;
    };

    /** @brief Unweighted squared norm of a numerical vector field over unique owned samples.
     * @details All three components contribute, without exterior ghosts. The uniform cell
     * volume cancels in E02's final/initial norm ratio [salari2000CodeVerificationMms,
     * section 4.1, Eq. (15)]. This reduction does not sample an analytic reference.
     */
    template <class Field>
    double pulseSquaredNorm(Field& field) {
        const auto view         = field.getView();
        double squared          = 0;
        std::uint64_t nonfinite = 0;
        Kokkos::parallel_reduce(
            "Pulse numerical vector norm", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& sum, std::uint64_t& invalid) {
                for (unsigned c = 0; c < 3; ++c) {
                    const double value = view(i, j, k)[c];
                    if (!Kokkos::isfinite(value))
                        ++invalid;
                    else
                        sum += value * value;
                }
            },
            squared, nonfinite);
        MPI_Allreduce(MPI_IN_PLACE, &squared, 1, MPI_DOUBLE, MPI_SUM, field.getLayout().comm);
        MPI_Allreduce(MPI_IN_PLACE, &nonfinite, 1, MPI_UINT64_T, MPI_SUM, field.getLayout().comm);
        if (nonfinite || !std::isfinite(squared))
            throw std::runtime_error("Pulse norm contains nonfinite numerical data");
        return squared;
    }

    namespace pulseMeasurementDetail {
        /// Local reductions, combined across ranks before forming either error ratio.
        struct ErrorReduction {
            double squaredError_m     = 0;
            double squaredReference_m = 0;
            double maximumError_m     = 0;
            double maximumReference_m = 0;
            std::uint64_t nonfinite_m = 0;
        };

        /// Evaluate the requested Cartesian component; four-potential component zero is excluded.
        KOKKOS_INLINE_FUNCTION double referenceValue(const E01Reference& reference,
                                                     Quantity quantity, double x, double y,
                                                     double z, double time, unsigned component) {
            if (quantity == Quantity::VectorPotential)
                return reference.referenceFourPotentialA(x, y, z, time, component + 1);
            if (quantity == Quantity::Electric)
                return reference.referenceE(x, y, z, time, component);
            return reference.referenceB(x, y, z, time, component);
        }

        /** @brief Form the relative norms after global sums and separate global maxima.
         * @details The discrete norms follow [salari2000CodeVerificationMms, section 4.1,
         * Eqs. (15), (18)]. E01.meas01--04 in E01.typ specify normalization by the sampled
         * reference norm. Uniform cell volume cancels between numerator and denominator.
         * Nonfinite data or a zero reference norm invalidate the measurement on every rank.
         */
        inline PulseErrorNorms globalNorms(ErrorReduction local, MPI_Comm communicator) {
            double squares[]{local.squaredError_m, local.squaredReference_m};
            double peaks[]{local.maximumError_m, local.maximumReference_m};
            MPI_Allreduce(MPI_IN_PLACE, squares, 2, MPI_DOUBLE, MPI_SUM, communicator);
            MPI_Allreduce(MPI_IN_PLACE, peaks, 2, MPI_DOUBLE, MPI_MAX, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &local.nonfinite_m, 1, MPI_UINT64_T, MPI_SUM, communicator);
            if (local.nonfinite_m || !std::isfinite(squares[0]) || !std::isfinite(squares[1])
                || !std::isfinite(peaks[0]) || !std::isfinite(peaks[1]) || squares[1] <= 0
                || peaks[1] <= 0)
                throw std::runtime_error(
                    "Pulse error measurement has nonfinite data or zero reference norm");
            return {std::sqrt(squares[0] / squares[1]), peaks[0] / peaks[1]};
        }
    }  // namespace pulseMeasurementDetail

    /** @brief Measure one numerical vector field at its native positions and time.
     * @param field Current A_n, electricField() or magneticField(), matching quantity.
     * @param reference Continuum reference with the experiment's fixed physical parameters.
     * @param contract The solver's staggered R2 level contract.
     * @param quantity VectorPotential, Electric or Magnetic. Potential samples use slots 1--3;
     * field samples use slots 0--2. All three components contribute, including reference zeros.
     * @param time Actual labelTime(quantity, completedSolves()), not the requested target time.
     * @details Only owned cells contribute, once globally. The reference uses the same
     * component positions as the numerical field [ryu2016potentialfdtd, Appendix A,
     * Eqs. (19)--(24)]. In particular, \f$B_y\f$ has a longitudinal half-cell offset.
     * The caller supplies R2-averaged magneticField(), whose time matches E
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.50)--(3.55)]. No interpolation or
     * analytic replacement of reconstructed fields is performed. Views are used only
     * during this measurement; no numerical field or halo is modified.
     */
    template <class Field>
    PulseErrorNorms measurePulseFieldError(Field& field, const E01Reference& reference,
                                           LevelContract contract, Quantity quantity, double time) {
        const auto view = field.getView();
        const auto h = field.get_mesh().getMeshSpacing(), origin = field.get_mesh().getOrigin();
        const auto domain = field.getLayout().getLocalNDIndex();
        const ippl::Vector<int, 3> first{domain[0].first(), domain[1].first(), domain[2].first()};
        const int ghost                = field.getNghost();
        const unsigned componentOffset = quantity == Quantity::VectorPotential ? 1 : 0;
        pulseMeasurementDetail::ErrorReduction local;
        Kokkos::parallel_reduce(
            "Pulse continuum field error", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& errorSquares, double& referenceSquares,
                          double& maximumError, double& maximumReference,
                          std::uint64_t& nonfinite) {
                const double x = origin[0] + (i + first[0] - ghost + 0.5) * h[0];
                const double y = origin[1] + (j + first[1] - ghost + 0.5) * h[1];
                const double z = origin[2] + (k + first[2] - ghost + 0.5) * h[2];
                for (unsigned c = 0; c < 3; ++c) {
                    const auto offset  = contract.positionOffset(quantity, c);
                    const double exact = pulseMeasurementDetail::referenceValue(
                        reference, quantity, x + offset[0] * h[0], y + offset[1] * h[1],
                        z + offset[2] * h[2], time, c);
                    const double numerical = view(i, j, k)[c + componentOffset];
                    const double error     = numerical - exact;
                    if (!Kokkos::isfinite(numerical) || !Kokkos::isfinite(exact)
                        || !Kokkos::isfinite(error)) {
                        ++nonfinite;
                        continue;
                    }
                    errorSquares += error * error;
                    referenceSquares += exact * exact;
                    maximumError     = Kokkos::max(maximumError, Kokkos::abs(error));
                    maximumReference = Kokkos::max(maximumReference, Kokkos::abs(exact));
                }
            },
            local.squaredError_m, local.squaredReference_m,
            Kokkos::Max<double>(local.maximumError_m),
            Kokkos::Max<double>(local.maximumReference_m), local.nonfinite_m);
        return pulseMeasurementDetail::globalNorms(local, field.getLayout().comm);
    }
}  // namespace chdr::experiments::e01e02
