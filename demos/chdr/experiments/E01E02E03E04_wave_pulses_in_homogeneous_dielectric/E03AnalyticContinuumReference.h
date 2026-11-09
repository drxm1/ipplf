/** @file E03AnalyticContinuumReference.h
 * @brief Native continuum samples of the spherical TE and TM pulses in E03.typ.
 */
#pragma once

#include <Kokkos_Array.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01AnalyticContinuumReference.h"

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::
    E03_analytic_continuum_reference {
    enum class Polarization {
        TEr,
        TMr
    };
    inline constexpr double w  = E01_analytic_continuum_reference::w;
    inline constexpr double R0 = 2 * w;

    /** @brief Evaluate the central outgoing pulse in the mesh coordinate system.
     * @details E03.typ constructs the two potential pairs from a compact radial wave.
     * Chew's homogeneous equations and reconstruction give both fields
     * [chew2014generalizedgauge, Eqs. (5)--(6), (11)--(13)]. Only the central copy
     * contributes during E03. The same unwrapped pulse supplies E04's free-space reference.
     */
    struct E03AnalyticContinuumReferenceSampler {
        double v_m;
        Polarization polarization_m;
        /// Translation of the spherical centre along each Cartesian axis.
        double centre_m = 0;

        /// Use four-potential storage order: phi, Ax, Ay, Az.
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double time,
                                                 unsigned component) const {
            const auto radial = radialSample(x, y, z, time);
            if (component == 0)
                return polarization_m == Polarization::TEr ? 0
                                                           : -v_m * radial.a_m * projection(radial);
            const unsigned c = component - 1;
            return polarization_m == Polarization::TEr ? radial.a_m * normalCrossDipole(radial, c)
                                                       : -dipole(c) * radial.fPrimeOverR_m;
        }

        KOKKOS_INLINE_FUNCTION double referenceFourPotentialA(double x, double y, double z,
                                                              double time,
                                                              unsigned component) const {
            return (*this)(x, y, z, time, component);
        }

        /// The TM electric field is the TE magnetic field multiplied by the wave speed.
        KOKKOS_INLINE_FUNCTION double referenceE(double x, double y, double z, double time,
                                                 unsigned component) const {
            const auto radial = radialSample(x, y, z, time);
            return v_m
                   * (polarization_m == Polarization::TEr
                          ? radial.c_m * normalCrossDipole(radial, component)
                          : teMagnetic(radial, component));
        }

        /// Both magnetic fields use the positive sign of curl A.
        KOKKOS_INLINE_FUNCTION double referenceB(double x, double y, double z, double time,
                                                 unsigned component) const {
            const auto radial = radialSample(x, y, z, time);
            return polarization_m == Polarization::TEr
                       ? teMagnetic(radial, component)
                       : -radial.c_m * normalCrossDipole(radial, component);
        }

    private:
        /// E03's radial derivatives a, b, c and the terms needed by the two potentials.
        struct RadialSample {
            Kokkos::Array<double, 3> normal_m{};
            double a_m = 0, b_m = 0, c_m = 0, aOverR_m = 0, fPrimeOverR_m = 0;
        };

        /** @brief Differentiate the compact generator as specified in E03.typ.
         * @details With \f$\psi=f(r-R_0-vt)/r\f$, \f$a=\psi_r\f$ and
         * \f$b=\psi_{rr}\f$. Check support before dividing by radius. The allowed
         * initial histories and evolved times put the origin outside the pulse.
         */
        KOKKOS_INLINE_FUNCTION RadialSample radialSample(double x, double y, double z,
                                                         double time) const {
            x -= centre_m;
            y -= centre_m;
            z -= centre_m;
            const double r = Kokkos::sqrt(x * x + y * y + z * z);
            const double s = r - R0 - v_m * time;
            if (r == 0 || Kokkos::abs(s) >= w)
                return {};
            const double u = s / w, q = 1 - u * u, q2 = q * q;
            const double f      = w * w * E01_analytic_continuum_reference::p(s);
            const double fPrime = -2 * u * f / (w * q2);
            const double fSecond =
                f / (w * w) * (4 * u * u / (q2 * q2) - 2 / q2 - 8 * u * u / (q2 * q));
            const double a = (fPrime - f / r) / r, c = (fSecond - fPrime / r) / r;
            const double b = (fSecond - 2 * fPrime / r + 2 * f / (r * r)) / r;
            return {{x / r, y / r, z / r}, a, b, c, a / r, fPrime / r};
        }

        KOKKOS_INLINE_FUNCTION static double dipole(unsigned c) {
            return (c + 1) / Kokkos::sqrt(14.0);
        }

        KOKKOS_INLINE_FUNCTION static double projection(const RadialSample& radial) {
            return radial.normal_m[0] * dipole(0) + radial.normal_m[1] * dipole(1)
                   + radial.normal_m[2] * dipole(2);
        }

        KOKKOS_INLINE_FUNCTION static double normalCrossDipole(const RadialSample& radial,
                                                               unsigned c) {
            const unsigned j = (c + 1) % 3, k = (c + 2) % 3;
            return radial.normal_m[j] * dipole(k) - radial.normal_m[k] * dipole(j);
        }

        /// Apply curl(grad psi cross d) as Hess(psi) d minus d times the Laplacian.
        KOKKOS_INLINE_FUNCTION static double teMagnetic(const RadialSample& radial, unsigned c) {
            return (radial.b_m - radial.aOverR_m) * radial.normal_m[c] * projection(radial)
                   - (radial.b_m + radial.aOverR_m) * dipole(c);
        }
    };
}  // namespace
   // chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::E03_analytic_continuum_reference
