/** @file E18E19Measurements.h
 * @brief Spatial Fourier sums, side energies and second-face field samples of E18 and E19.
 */
#pragma once

#include <Kokkos_Core.hpp>

#include <Kokkos_MathematicalFunctions.hpp>
#include <array>
#include <cmath>
#include <cstdint>
#include <mpi.h>
#include <stdexcept>

#include "core/base/LevelContract.h"
#include "core/material/MaterialOwner.h"
#include "core/material/O1VSCoefficients.h"

namespace chdr::experiments::e18e19 {
    namespace measurementDetail {
        /// Combine owned-cell sums before taking magnitudes, and reject invalid data on all ranks.
        template <std::size_t Count>
        std::array<double, Count> globalSums(std::array<double, Count> sums,
                                             std::uint64_t nonfinite, MPI_Comm communicator) {
            MPI_Allreduce(MPI_IN_PLACE, sums.data(), Count, MPI_DOUBLE, MPI_SUM, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &nonfinite, 1, MPI_UINT64_T, MPI_SUM, communicator);
            for (double value : sums)
                if (!std::isfinite(value))
                    ++nonfinite;
            if (nonfinite)
                throw std::runtime_error("E18/E19 measurement contains nonfinite field data");
            return sums;
        }

        /// Combine maxima and their validity counts without summing the maxima themselves.
        inline double globalMaximum(double peak, std::uint64_t nonfinite, MPI_Comm communicator) {
            MPI_Allreduce(MPI_IN_PLACE, &peak, 1, MPI_DOUBLE, MPI_MAX, communicator);
            MPI_Allreduce(MPI_IN_PLACE, &nonfinite, 1, MPI_UINT64_T, MPI_SUM, communicator);
            if (nonfinite || !std::isfinite(peak))
                throw std::runtime_error("E18/E19 peak contains nonfinite field data");
            return peak;
        }

        /** @brief Accumulate one complex Fourier sum from owned vacuum cells.
         * @details E18.meas01 and E19.meas01 define the negative complex exponential.
         * The selected field is component y of E or B_R2. Side membership uses the
         * storing cell center, while the phase uses LevelContract's component position.
         * Consequently, a magnetic sample on the interface is included in side 1.
         * The sign of normalWavenumber selects incident or reflected propagation.
         */
        template <class Field>
        std::array<double, 2> localFourierSum(Field& field, LevelContract contract,
                                              Quantity quantity, double kx, double normalWavenumber,
                                              double interfaceZ, std::uint64_t& nonfinite) {
            const auto view = field.getView();
            const auto h = field.get_mesh().getMeshSpacing(), origin = field.get_mesh().getOrigin();
            const auto first  = field.getLayout().getLocalNDIndex().first();
            const auto offset = contract.positionOffset(quantity, 1);
            const int ghost   = field.getNghost();
            std::array<double, 2> sums{};
            Kokkos::parallel_reduce(
                "E18/E19 local Fourier sum", field.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k, double& real, double& imaginary,
                              std::uint64_t& invalid) {
                    const double centerZ = origin[2] + (k + first[2] - ghost + 0.5) * h[2];
                    if (centerZ >= interfaceZ)
                        return;
                    const double x = origin[0] + (i + first[0] - ghost + 0.5 + offset[0]) * h[0];
                    const double z = centerZ + offset[2] * h[2];
                    const double value = view(i, j, k)[1];
                    const double phase = kx * x + normalWavenumber * z;
                    if (!Kokkos::isfinite(value) || !Kokkos::isfinite(phase)) {
                        ++invalid;
                        return;
                    }
                    real += value * Kokkos::cos(phase);
                    imaginary -= value * Kokkos::sin(phase);
                },
                sums[0], sums[1], nonfinite);
            return sums;
        }

        /// Use the same cached positive-face coefficients as O1VS for the TM energy.
        inline Kokkos::Array<material::ScalarFieldConstView, 3> energyPermittivity(
            const material::MaterialOwner& owner, bool tm) {
            if (tm)
                return {owner.derive(material::HarmonicFace<0>{}),
                        owner.derive(material::HarmonicFace<1>{}),
                        owner.derive(material::HarmonicFace<2>{})};
            return {owner.epsilon(), owner.epsilon(), owner.epsilon()};
        }
    }  // namespace measurementDetail

    /** @brief Return the twelve Fourier magnitudes defined by E18.meas01 or E19.meas01.
     * @details The caller supplies Ey with Quantity::Electric for TE and B_R2,y with
     * Quantity::Magnetic for TM. reflected=false forms I and reflected=true forms P.
     * Real and imaginary parts are summed over MPI ranks before taking the absolute
     * value. The field execution space performs every local reduction, and only
     * 24 scalar sums and the invalid-value count move to the host for MPI.
     * No field or halo is modified. Floating-point reduction order can vary with
     * the execution space and the MPI decomposition.
     */
    template <class Field>
    std::array<double, 12> measureFourierMagnitude(Field& field, LevelContract contract,
                                                   Quantity quantity, double kx,
                                                   const std::array<double, 12>& normalWavenumbers,
                                                   double interfaceZ, bool reflected) {
        std::array<double, 24> sums{};
        std::uint64_t nonfinite = 0;
        for (unsigned index = 0; index < normalWavenumbers.size(); ++index) {
            std::uint64_t invalid = 0;
            const double q        = (reflected ? -1 : 1) * normalWavenumbers[index];
            const auto local = measurementDetail::localFourierSum(field, contract, quantity, kx, q,
                                                                  interfaceZ, invalid);
            sums[2 * index]  = local[0];
            sums[2 * index + 1] = local[1];
            nonfinite += invalid;
        }
        sums = measurementDetail::globalSums(sums, nonfinite, field.getLayout().comm);
        std::array<double, 12> magnitudes{};
        for (unsigned index = 0; index < magnitudes.size(); ++index)
            magnitudes[index] = std::hypot(sums[2 * index], sums[2 * index + 1]);
        return magnitudes;
    }

    /** @brief Return the sampled energies on side 1 and side 2 for E18.meas02 or E19.meas02.
     * @details The density is \f$(\epsilon_r |E|^2+|B|^2)/2\f$ in the specified
     * units, from [ellingson2020ElectromagneticsVol2, Eqs. (3.20), (3.22)].
     * E18 uses cell permittivity for every component. E19 uses the harmonic
     * positive-face coefficient cached by the update, as prescribed in E19.meas02.
     * Each component belongs to its storing cell, including samples on the interface.
     * The caller supplies E and B_R2 at their common LevelContract time and layout.
     * Only owned cells enter the execution-space reduction. The host receives two sums.
     */
    template <class Field>
    std::array<double, 2> measureSideEnergy(Field& electric, Field& magnetic,
                                            const material::MaterialOwner& owner, bool tm,
                                            double interfaceZ) {
        const auto electricView = electric.getView(), magneticView = magnetic.getView();
        const auto epsilon   = measurementDetail::energyPermittivity(owner, tm);
        const auto h         = electric.get_mesh().getMeshSpacing();
        const double originZ = electric.get_mesh().getOrigin()[2];
        const int firstZ     = electric.getLayout().getLocalNDIndex()[2].first();
        const int ghost      = electric.getNghost();
        std::array<double, 2> sums{};
        std::uint64_t nonfinite = 0;
        Kokkos::parallel_reduce(
            "E18/E19 side field energies", electric.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& energy1, double& energy2,
                          std::uint64_t& invalid) {
                const double z = originZ + (k + firstZ - ghost + 0.5) * h[2];
                double energy  = 0;
                for (unsigned component = 0; component < 3; ++component) {
                    const double e            = electricView(i, j, k)[component];
                    const double b            = magneticView(i, j, k)[component];
                    const double permittivity = epsilon[component](i, j, k);
                    if (!Kokkos::isfinite(e) || !Kokkos::isfinite(b)
                        || !Kokkos::isfinite(permittivity))
                        ++invalid;
                    else
                        energy += permittivity * e * e + b * b;
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

    /** @brief Return the largest electric component over all cells or the two second-face layers.
     * @details E18.meas03, also used by E19, selects the storing cells k=0 and k=Nz-1.
     * The caller retains the initial whole-domain maximum and the maximum of these
     * layer samples over the energy observation times. The ratio is the specified b.
     * Halos are excluded, and nonfinite owned samples are rejected collectively.
     */
    template <class Field>
    double measureElectricPeak(Field& electric, bool secondFaceOnly) {
        const auto view    = electric.getView();
        const auto domainZ = electric.getLayout().getDomain()[2];
        const int firstZ   = electric.getLayout().getLocalNDIndex()[2].first();
        const int lowerZ = domainZ.first(), upperZ = domainZ.last(), ghost = electric.getNghost();
        double peak             = 0;
        std::uint64_t nonfinite = 0;
        Kokkos::parallel_reduce(
            "E18/E19 electric component peak", electric.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& maximum, std::uint64_t& invalid) {
                const int globalZ = k + firstZ - ghost;
                double candidate  = 0;
                if (!secondFaceOnly || globalZ == lowerZ || globalZ == upperZ)
                    for (unsigned component = 0; component < 3; ++component) {
                        const double value = view(i, j, k)[component];
                        if (!Kokkos::isfinite(value))
                            ++invalid;
                        else
                            candidate = Kokkos::fmax(candidate, Kokkos::abs(value));
                    }
                maximum = Kokkos::fmax(maximum, candidate);
            },
            Kokkos::Max<double>(peak), nonfinite);
        return measurementDetail::globalMaximum(peak, nonfinite, electric.getLayout().comm);
    }
}  // namespace chdr::experiments::e18e19
