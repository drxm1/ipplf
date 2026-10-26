/** @file Waves.h
 * @brief Analytic generators; samples use global positions and explicit labels.
 */
#pragma once

#include <array>
#include <stdexcept>

#include "PeriodicMode.h"
#include "core/base/LevelContract.h"
#include "core/base/SpatialDimension.h"
#include "core/operators/Operators.h"

static_assert(chdr::solver::SpatialDim == 3,
              "Waves.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    /** @brief Exact discrete gauge wave, with the scalar potential matched to its time difference.
     * @details Deduction: insert \f$A=D_c\chi\f$, \f$\phi=-\delta_t\chi\f$ in
     * [christlieb2024GaugeConserving, section 3.1, proof of Theorem 3.5], using
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.17), (3.24)]. R1 instead
     * samples \f$-\partial_t\chi\f$ at its whole-step label.
     * The centered spatial differences commute with the time difference, so the
     * matched R2/R2Prime electric terms cancel and the curl of A vanishes. This
     * does not impose a zero Lorenz residual: for R2 it is
     * \f$(D_c\cdot D_c-L_h)\chi\f$, since the wave recurrence uses \f$L_h\f$.
     * That residual follows by substitution into the gauge definition in
     * [christlieb2024GaugeConserving, Eq. (29)]; the centered-gradient composition
     * differs from the solver Laplacian (see Operators).
     */
    struct GaugeWave {
        ippl::Vector<double, SpatialDim> waveNumber_m, g_m;
        double omega_m, dt_m;
        EBReconstructionConvention ebReconstructionConvention_m;
        GaugeWave(ippl::Vector<double, SpatialDim> waveNumber, ippl::Vector<double, SpatialDim> h,
                  double dt, Stencil stencil, EBReconstructionConvention ebReconstructionConvention)
            : waveNumber_m(waveNumber)
            , omega_m(discreteFrequency(waveNumber, h, dt, stencil))
            , dt_m(dt)
            , ebReconstructionConvention_m(ebReconstructionConvention) {
            for (unsigned d = 0; d < SpatialDim; ++d)
                g_m[d] = centeredSymbol(waveNumber[d], h[d]);
        }
        /** @brief Sample a scalar or vector potential at a time supplied by LevelContract.
         * @details Components one through three use the centered spatial symbol.
         * Component zero uses the matching one-step time symbol for R2/R2Prime,
         * while R1 uses the continuum time derivative at its whole-step label.
         * The caller supplies the label; this functor never adds a second time shift.
         */
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double t,
                                                 unsigned c) const {
            const double phase =
                waveNumber_m[0] * x + waveNumber_m[1] * y + waveNumber_m[2] * z - omega_m * t;
            const double scalar = ebReconstructionConvention_m == EBReconstructionConvention::R1
                                      ? omega_m
                                      : 2 * Kokkos::sin(omega_m * dt_m / 2) / dt_m;
            return -(c == 0 ? scalar : g_m[c - 1]) * Kokkos::sin(phase);
        }
        /** @brief Exact discrete R1 electric oracle at its current potential time.
         * @details For \f$\theta=k\cdot x-\omega_h t\f$, the component is
         * \f$g_c[\omega_h\cos\theta-W\cos(\theta+\omega_h\Delta t/2)]\f$.
         * This follows by substituting the class potentials into
         * [fallahi2020mithra20fullwavesimulation, Eq. (3.9)], with the shipped
         * backward time difference. The gauge-wave test compares
         * this full expression, never only its leading first-order expansion.
         */
        KOKKOS_INLINE_FUNCTION double r1Electric(double x, double y, double z, double t,
                                                 unsigned c) const {
            const double phase =
                waveNumber_m[0] * x + waveNumber_m[1] * y + waveNumber_m[2] * z - omega_m * t;
            const double w = 2 * Kokkos::sin(omega_m * dt_m / 2) / dt_m;
            return g_m[c]
                   * (omega_m * Kokkos::cos(phase) - w * Kokkos::cos(phase + omega_m * dt_m / 2));
        }
        void requirePeriodic(std::array<double, SpatialDim> lengths) const {
            requirePeriodicMode(waveNumber_m, lengths);
        }
    };
}  // namespace chdr::solver

#include "LorenzWave.h"
