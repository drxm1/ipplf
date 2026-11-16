/** @file E03AnalyticContinuumReference.h
 * @brief Native continuum samples of the spherical dipole pulse in E03.typ.
 */
#pragma once

#include <Kokkos_Array.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01AnalyticContinuumReference.h"

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::
    E03_analytic_continuum_reference {
    inline constexpr double w  = E01_analytic_continuum_reference::w;
    inline constexpr double R0 = 2 * w;

    /// The radial generator \f$f=w^2p\f$ and its first two derivatives at one value of \f$s\f$.
    struct RadialGenerator {
        double f_m = 0, fPrime_m = 0, fSecond_m = 0;
    };

    /** @brief The compact radial profile of E03.typ with sharpness \f$\beta>0\f$.
     * @details The profile is \f$p_\beta(s)=A_0\exp\bigl(\beta\,(1-1/(1-(s/w)^2))\bigr)\f$
     * for \f$|s|<w\f$ and zero otherwise. \f$\beta=1\f$ is the profile of E01.
     * A larger \f$\beta\f$ keeps the support and the smoothness, and approaches the
     * Gaussian \f$A_0\exp(-\beta s^2/w^2)\f$ near the center. Its flanks are less steep
     * relative to its width, so its spectrum is narrower and a given mesh resolves it
     * better (the factor in the exponent equals \f$-\beta u^2/(1-u^2)\f$, \f$u=s/w\f$).
     * A radial profile supplies `contains(s)`, `generator(s)` and `record(node)`.
     */
    struct E03SharpenedCompactProfile {
        double sharpness_m = 1;

        /// The generator is zero outside the support; the caller then skips all divisions.
        KOKKOS_INLINE_FUNCTION bool contains(double s) const { return Kokkos::abs(s) < w; }

        /** @brief Differentiate the compact generator as specified in E03.typ.
         * @details With \f$u=s/w\f$, \f$q=1-u^2\f$ and \f$f=w^2p_\beta\f$, the chain rule gives
         * \f$f'=-2\beta u f/(wq^2)\f$ and
         * \f$f''=f\,(4\beta^2u^2/q^4-2\beta/q^2-8\beta u^2/q^3)/w^2\f$.
         * @pre contains(s).
         */
        KOKKOS_INLINE_FUNCTION RadialGenerator generator(double s) const {
            const double u = s / w, q = 1 - u * u, q2 = q * q, beta = sharpness_m;
            const double f       = w * w * sharpenedProfile(s);
            const double fPrime  = beta * (-2 * u * f / (w * q2));
            const double fSecond = f / (w * w)
                                   * (beta * beta * 4 * u * u / (q2 * q2) - beta * 2 / q2
                                      - beta * 8 * u * u / (q2 * q));
            return {f, fPrime, fSecond};
        }

        /// Record the sharpness \f$\beta\f$ of the run (by reference; a copied node is detached).
        template <class Node>
        void record(Node& node) const {
            node["profile_sharpness"] = sharpness_m;
        }

    private:
        /** @brief Evaluate \f$p_\beta(s)\f$ inside its support.
         * @details For \f$\beta=1\f$ the expression and its rounding equal those of the
         * E01 profile, so that results with the default sharpness do not change.
         * @pre \f$|s|<w\f$.
         */
        KOKKOS_INLINE_FUNCTION double sharpenedProfile(double s) const {
            return E01_analytic_continuum_reference::A0
                   * Kokkos::exp(sharpness_m * (1 - 1 / (1 - Kokkos::pow(s / w, 2))));
        }
    };

    /** @brief Evaluate the central outgoing pulse in the mesh coordinate system.
     * @details E03.typ constructs the potentials from a radial wave
     * \f$\psi=f(r-R_0-vt)/r\f$ and a fixed unit vector \f$\mathbf d\f$:
     * \f[
     *   \mathbf A=\mathbf d\,\partial_t\psi/v,\qquad \phi=-v\,\mathbf d\cdot\nabla\psi .
     * \f]
     * This is the pulse of an electric dipole along \f$\mathbf d\f$ (E03.typ calls it the
     * TM\f$_r\f$ case: the magnetic field has no radial component). Both potentials are
     * nonzero, so the pulse exercises the scalar potential, its gradient in the electric
     * field and the Lorenz condition. Its dual with \f$\phi=0\f$ (TE\f$_r\f$) has the same
     * fields with \f$\mathbf E\f$ and \f$\mathbf B\f$ exchanged and adds no code path; it is
     * not implemented. Chew's homogeneous equations and reconstruction give the fields
     * [chew2014generalizedgauge, Eqs. (5)--(6), (11)--(13)]. Only the central copy
     * contributes during E03. The same unwrapped pulse supplies E04's free-space reference.
     * @tparam Profile Radial profile, as E03SharpenedCompactProfile.
     */
    template <class Profile = E03SharpenedCompactProfile>
    struct E03AnalyticContinuumReferenceSampler {
        double v_m;
        /// Translation of the spherical centre along each Cartesian axis.
        double centre_m = 0;
        /// Radial profile \f$p\f$ of the generator \f$f=w^2p\f$.
        Profile profile_m{};

        /// Use four-potential storage order: phi, Ax, Ay, Az.
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double z, double time,
                                                 unsigned component) const {
            const auto radial = radialSample(x, y, z, time);
            if (component == 0)
                return -v_m * radial.a_m * projection(radial);
            return -dipole(component - 1) * radial.fPrimeOverR_m;
        }

        KOKKOS_INLINE_FUNCTION double referenceFourPotentialA(double x, double y, double z,
                                                              double time,
                                                              unsigned component) const {
            return (*this)(x, y, z, time, component);
        }

        /// \f$\mathbf E=v\,[(b-a/r)\,\mathbf n(\mathbf n\cdot\mathbf d)-(b+a/r)\,\mathbf d]\f$.
        KOKKOS_INLINE_FUNCTION double referenceE(double x, double y, double z, double time,
                                                 unsigned component) const {
            return v_m * hessianTerm(radialSample(x, y, z, time), component);
        }

        /// \f$\mathbf B=-c\,(\mathbf n\times\mathbf d)\f$ with \f$c=f''/r-f'/r^2\f$.
        KOKKOS_INLINE_FUNCTION double referenceB(double x, double y, double z, double time,
                                                 unsigned component) const {
            const auto radial = radialSample(x, y, z, time);
            return -radial.c_m * normalCrossDipole(radial, component);
        }

        /// Write the center, the speed and the profile parameters to the run record.
        template <class Node>
        void record(Node node) const {
            node["centre_in_mesh_coordinates"] = centre_m;
            node["speed"]                      = v_m;
            profile_m.record(node);
        }

    private:
        /// E03's radial derivatives a, b, c and the terms needed by the potentials.
        struct RadialSample {
            Kokkos::Array<double, 3> normal_m{};
            double a_m = 0, b_m = 0, c_m = 0, aOverR_m = 0, fPrimeOverR_m = 0;
        };

        /** @brief Differentiate the radial wave as specified in E03.typ.
         * @details With \f$\psi=f(r-R_0-vt)/r\f$, \f$a=\psi_r\f$ and
         * \f$b=\psi_{rr}\f$. Check the profile's support before dividing by radius. The
         * allowed initial histories and evolved times put the origin outside a compact pulse.
         */
        KOKKOS_INLINE_FUNCTION RadialSample radialSample(double x, double y, double z,
                                                         double time) const {
            x -= centre_m;
            y -= centre_m;
            z -= centre_m;
            const double r = Kokkos::sqrt(x * x + y * y + z * z);
            const double s = r - R0 - v_m * time;
            if (r == 0 || !profile_m.contains(s))
                return {};
            const auto g   = profile_m.generator(s);
            const double a = (g.fPrime_m - g.f_m / r) / r, c = (g.fSecond_m - g.fPrime_m / r) / r;
            const double b = (g.fSecond_m - 2 * g.fPrime_m / r + 2 * g.f_m / (r * r)) / r;
            return {{x / r, y / r, z / r}, a, b, c, a / r, g.fPrime_m / r};
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

        /// Component c of Hess(psi) d minus d times the Laplacian of psi.
        KOKKOS_INLINE_FUNCTION static double hessianTerm(const RadialSample& radial, unsigned c) {
            return (radial.b_m - radial.aOverR_m) * radial.normal_m[c] * projection(radial)
                   - (radial.b_m + radial.aOverR_m) * dipole(c);
        }
    };
}  // namespace
   // chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::E03_analytic_continuum_reference
