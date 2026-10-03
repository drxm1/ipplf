/** @file VacuumWave.h
 * @brief Normalized analytic four-potential used to initialize the periodic vacuum demo.
 */
#pragma once
#include <Kokkos_Core.hpp>

#include <Kokkos_MathematicalConstants.hpp>

#include "Types/Vector.h"

namespace chdr::vacuum {
    /** @brief Four-potential of the chosen +z vacuum wave (c=1, wavelength=1).
     * @param z Cell-center z coordinate in normalized length units.
     * @param time Physical sample time in normalized units.
     * @return (phi, Ax, Ay, Az); no allocation or field modification.
     * @details Ax = sin(2*pi*(z-time))/(2*pi); other components are zero.
     * This reference follows Fallahi, MITHRA 2.0 (arXiv:2009.13645v1),
     * "Wave Equation", Eqs. (3.6)-(3.9).
     */
    KOKKOS_INLINE_FUNCTION ippl::Vector<double, 4> referencePotential(double z, double time) {
        constexpr double WaveNumber = 2.0 * Kokkos::numbers::pi_v<double>;
        return ippl::Vector<double, 4>(0.0, Kokkos::sin(WaveNumber * (z - time)) / WaveNumber, 0.0,
                                       0.0);
    }
}  // namespace chdr::vacuum
