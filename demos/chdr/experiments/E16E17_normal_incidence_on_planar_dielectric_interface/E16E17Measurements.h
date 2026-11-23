/** @file E16E17Measurements.h
 * @brief Owned-cell field errors and energies on the two sides of the dielectric face.
 */
#pragma once

#include <Kokkos_Core.hpp>

#include <Kokkos_MathematicalFunctions.hpp>
#include <array>
#include <cmath>
#include <cstdint>
#include <mpi.h>
#include <stdexcept>

#include "Types/Vector.h"

#include "core/base/LevelContract.h"
#include "core/material/MaterialOwner.h"
#include "experiments/E16E17_normal_incidence_on_planar_dielectric_interface/E16E17AnalyticContinuumReference.h"

namespace chdr::experiments::e16e17 {
    namespace measurementDetail {
        /// Select a continuum component in the same numbering as the solver's vector fields.
        KOKKOS_INLINE_FUNCTION double referenceValue(
            const E16E17AnalyticContinuumReferenceSampler& reference, Quantity quantity, double x,
            double y, double z, double time, unsigned component) {
            if (quantity == Quantity::VectorPotential)
                return reference.referenceFourPotentialA(x, y, z, time, component + 1);
            if (quantity == Quantity::Electric)
                return reference.referenceE(x, y, z, time, component);
            return reference.referenceB(x, y, z, time, component);
        }

        /// Combine owned-cell sums and reject invalid values collectively, including overflow.
        template <std::size_t Count>
        std::array<double, Count> globalSums(std::array<double, Count> sums,
                                             std::uint64_t nonfinite, MPI_Comm communicator) {
            MPI_Allreduce(MPI_IN_PLACE, sums.data(), Count, MPI_DOUBLE, MPI_SUM, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &nonfinite, 1, MPI_UINT64_T, MPI_SUM, communicator);
            for (double value : sums)
                if (!std::isfinite(value))
                    ++nonfinite;
            if (nonfinite)
                throw std::runtime_error("E16/E17 measurement contains nonfinite field data");
            return sums;
        }
    }  // namespace measurementDetail

    /** @brief Measure E16.meas01 or E16.meas02, including their E17 counterparts.
     * @return Squared sums in the order error side 1, reference side 1, error side 2,
     * reference side 2. A zero reference norm is retained for the initial empty side.
     * @details The reference is evaluated at the component positions from LevelContract
     * and the supplied actual field time. This follows the sampled-error construction of
     * [salari2000CodeVerificationMms, section 4.1], restricted to each side in E16.typ.
     * Side membership uses the storing cell's center. In particular, a magnetic sample
     * on the face belongs to side 1, as specified in E16.typ. Every owned component
     * contributes once. The reduction uses the field execution space and copies only
     * scalar sums to the host for MPI. It leaves the field and its halos unchanged.
     */
    template <class Field>
    std::array<double, 4> measureSideFieldSums(
        Field& field, const E16E17AnalyticContinuumReferenceSampler& reference,
        LevelContract contract, Quantity quantity, double time, double interfaceZ) {
        const auto view = field.getView();
        const auto h = field.get_mesh().getMeshSpacing(), origin = field.get_mesh().getOrigin();
        const auto first = field.getLayout().getLocalNDIndex().first();
        const int ghost  = field.getNghost();
        std::array<double, 4> sums{};
        std::uint64_t nonfinite = 0;
        Kokkos::parallel_reduce(
            "E16/E17 side field errors", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& error1, double& reference1, double& error2,
                          double& reference2, std::uint64_t& invalid) {
                const double x = origin[0] + (i + first[0] - ghost + 0.5) * h[0];
                const double y = origin[1] + (j + first[1] - ghost + 0.5) * h[1];
                const double z = origin[2] + (k + first[2] - ghost + 0.5) * h[2];
                for (unsigned component = 0; component < 3; ++component) {
                    const auto offset  = contract.positionOffset(quantity, component);
                    const double exact = measurementDetail::referenceValue(
                        reference, quantity, x + offset[0] * h[0], y + offset[1] * h[1],
                        z + offset[2] * h[2], time, component);
                    const double value      = view(i, j, k)[component];
                    const double difference = value - exact;
                    if (!Kokkos::isfinite(value) || !Kokkos::isfinite(exact)
                        || !Kokkos::isfinite(difference)) {
                        ++invalid;
                        continue;
                    }
                    if (z < interfaceZ) {
                        error1 += difference * difference;
                        reference1 += exact * exact;
                    } else {
                        error2 += difference * difference;
                        reference2 += exact * exact;
                    }
                }
            },
            sums[0], sums[1], sums[2], sums[3], nonfinite);
        return measurementDetail::globalSums(sums, nonfinite, field.getLayout().comm);
    }

    /** @brief Measure the two sampled field energies of E16.meas03 and E17.
     * @details E16.typ defines the cell quadrature of
     * \f$(\epsilon_r |E|^2+|B|^2)/2\f$ with unit permeability, deduced from
     * [ellingson2020ElectromagneticsVol2, Eqs. (3.20), (3.22)]. The cell permittivity
     * weights all electric components, including components whose reference vanishes.
     * The caller supplies E and R2-averaged B at the same half-step and on the same
     * layout. Samples on the interface belong to their storing cell. Uniform cell
     * volume is applied after summation. Only owned cells enter the device reduction.
     */
    template <class Field>
    std::array<double, 2> measureSideEnergy(Field& electric, Field& magnetic,
                                            material::ScalarFieldConstView epsilon,
                                            double interfaceZ) {
        const auto electricView = electric.getView(), magneticView = magnetic.getView();
        const auto h         = electric.get_mesh().getMeshSpacing();
        const double originZ = electric.get_mesh().getOrigin()[2];
        const int firstZ     = electric.getLayout().getLocalNDIndex()[2].first();
        const int ghost      = electric.getNghost();
        std::array<double, 2> sums{};
        std::uint64_t nonfinite = 0;
        Kokkos::parallel_reduce(
            "E16/E17 side field energies", electric.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& energy1, double& energy2,
                          std::uint64_t& invalid) {
                const double z = originZ + (k + firstZ - ghost + 0.5) * h[2];
                double energy  = 0;
                for (unsigned component = 0; component < 3; ++component) {
                    const double e = electricView(i, j, k)[component];
                    const double b = magneticView(i, j, k)[component];
                    if (!Kokkos::isfinite(e) || !Kokkos::isfinite(b)
                        || !Kokkos::isfinite(epsilon(i, j, k)))
                        ++invalid;
                    else
                        energy += epsilon(i, j, k) * e * e + b * b;
                }
                if (z < interfaceZ)
                    energy1 += energy;
                else
                    energy2 += energy;
            },
            sums[0], sums[1], nonfinite);
        sums = measurementDetail::globalSums(sums, nonfinite, electric.getLayout().comm);
        for (double& energy : sums)
            energy *= 0.5 * h[0] * h[1] * h[2];
        return sums;
    }
}  // namespace chdr::experiments::e16e17
