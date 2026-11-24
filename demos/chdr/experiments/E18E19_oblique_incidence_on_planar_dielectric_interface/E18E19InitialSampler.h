/** @file E18E19InitialSampler.h
 * @brief Tabulate the oblique packet once per local z position and potential slot.
 */
#pragma once

#include <Kokkos_Core.hpp>

#include <Kokkos_MathematicalFunctions.hpp>

#include "core/base/LevelContract.h"
#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19AnalyticContinuumReference.h"

namespace chdr::experiments::e18e19 {
    /** @brief Device sampler for the two initial potential histories of E18 and E19.
     * @details The complex z amplitudes are evaluated on the host and copied once
     * to the default execution space. The sampling kernel then applies the common
     * tangential phase. This factorization is exact for the finite sum in E18-phases.
     * Each rank stores only its owned z positions. Spatial offsets and times come
     * from LevelContract, which also supplies them to solver::writeSamples.
     * The table is used only during history loading, before any solver advance.
     */
    struct E18E19InitialSampler {
        using Table = Kokkos::View<double****>;

        Table amplitudes_m;
        Kokkos::Array<double, 3> zOffsets_m;
        int firstGlobalZ_m;
        double originZ_m, spacingZ_m, currentTime_m, timeStep_m;

        /** @brief Evaluate one stored potential sample from its complex z amplitude.
         * @pre The position is an owned sample at the component's contract offset,
         * and time is either the current or the previous vector-potential label.
         * Scalar potential is exactly zero at either of its own contract times.
         */
        KOKKOS_INLINE_FUNCTION double operator()(double x, double, double z, double time,
                                                 unsigned component) const {
            if (component == 0)
                return 0;
            const unsigned cartesian = component - 1;
            const int index          = static_cast<int>(Kokkos::floor(
                (z - originZ_m) / spacingZ_m - firstGlobalZ_m - zOffsets_m[cartesian]));
            const int slot           = time < currentTime_m - 0.5 * timeStep_m ? 1 : 0;
            const double phase       = E18E19AnalyticContinuumReference::KX * x;
            return amplitudes_m(slot, index, cartesian, 0) * Kokkos::cos(phase)
                   - amplitudes_m(slot, index, cartesian, 1) * Kokkos::sin(phase);
        }
    };

    /// Fill the host mirror, then transfer its real and imaginary amplitudes to the device.
    inline void fillInitialAmplitudes(E18E19InitialSampler& result,
                                      const E18E19AnalyticContinuumReference& reference) {
        auto host = Kokkos::create_mirror_view(result.amplitudes_m);
        for (unsigned component = 0; component < 3; ++component) {
            for (int slot = 0; slot < 2; ++slot)
                for (std::size_t cell = 0; cell < result.amplitudes_m.extent(1); ++cell) {
                    const double z =
                        result.originZ_m
                        + (result.firstGlobalZ_m + cell + 0.5 + result.zOffsets_m[component])
                              * result.spacingZ_m;
                    const auto value = reference.potentialAmplitude(
                        z, result.currentTime_m - slot * result.timeStep_m, component);
                    host(slot, cell, component, 0) = value.real();
                    host(slot, cell, component, 1) = value.imag();
                }
        }
        Kokkos::deep_copy(result.amplitudes_m, host);
    }

    /** @brief Build both z tables at the contract's vector positions and time labels.
     * @details The host mirror holds real and imaginary parts separately, so no
     * host complex object is captured by a device kernel. The deep copy completes
     * before the returned sampler can be passed to the shared history loader.
     * The metadata overload also permits a direct table check without evolving fields.
     */
    inline E18E19InitialSampler makeInitialSampler(
        const E18E19AnalyticContinuumReference& reference, int localCellsZ, int firstGlobalZ,
        double originZ, double spacingZ, LevelContract contract, double timeStep) {
        E18E19InitialSampler result{
            E18E19InitialSampler::Table("E18/E19 initial z amplitudes", 2, localCellsZ, 3, 2),
            {},
            firstGlobalZ,
            originZ,
            spacingZ,
            contract.labelTime(Quantity::VectorPotential, 0, timeStep),
            timeStep};
        for (unsigned component = 0; component < 3; ++component)
            result.zOffsets_m[component] =
                contract.positionOffset(Quantity::VectorPotential, component)[2];
        fillInitialAmplitudes(result, reference);
        return result;
    }

    /// Extract the actual local mesh and selected solver contract before history loading.
    template <class Solver>
    E18E19InitialSampler makeInitialSampler(const E18E19AnalyticContinuumReference& reference,
                                            Solver& solver) {
        const auto local   = solver.A_n.getLayout().getLocalNDIndex();
        const auto spacing = solver.A_n.get_mesh().getMeshSpacing();
        const auto origin  = solver.A_n.get_mesh().getOrigin();
        return makeInitialSampler(reference, local[2].length(), local[2].first(), origin[2],
                                  spacing[2], solver.levelContract(), solver.getDt());
    }
}  // namespace chdr::experiments::e18e19
