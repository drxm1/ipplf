/** @file MediumOracles.h
 * @brief Expected values of the homogeneous-dielectric tests: step bound, discrete symbol,
 * source scaling, plane and driven waves, constraint factors and the dispersion record.
 * @details Normalized units \f$c = \varepsilon_0 = \mu_0 = 1\f$, \f$\mu_r = 1\f$, relative
 * permittivity \f$\varepsilon \ge 1\f$, index \f$n = \sqrt\varepsilon\f$. With the Lorenz gauge
 * \f$\nabla\cdot\mathbf A = -\varepsilon\,\partial_t\phi\f$ the potentials in a homogeneous medium
 * obey \f$\nabla^2\phi - \varepsilon\,\partial_t^2\phi = -\rho/\varepsilon\f$ and
 * \f$\nabla^2\mathbf A - \varepsilon\,\partial_t^2\mathbf A = -\mathbf J\f$
 * [chew2014generalizedgauge, section 2, Eqs. (11)-(13), p. 3]. Divided by \f$-\varepsilon\f$
 * they are four uncoupled wave equations with speed \f$1/n\f$,
 * \f$\partial_t^2\psi - \varepsilon^{-1}\nabla^2\psi = s_\varepsilon\f$ (Deduction), and the
 * scheme under test is the vacuum update [fallahi2020mithra20fullwavesimulation,
 * section 3.1.1, Eqs. (3.6)-(3.7), p. 12] applied to them:
 * \f[ \psi^{k+1} = 2\psi^k - \psi^{k-1} + \Delta t^2\bigl(\varepsilon^{-1}L_h\psi^k +
 * s_\varepsilon^k\bigr). \f] These host functions are written independently of the solver and its
 * helpers; the tests take their expected values from here. Vectors are Cartesian \f$(x, y, z)\f$.
 */
#pragma once
#include "Ippl.h"

#include <cmath>

#include "core/base/SpatialDimension.h"

static_assert(chdr::SpatialDim == 3,
              "MediumOracles.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    inline double dotProduct(ippl::Vector<double, SpatialDim> a,
                             ippl::Vector<double, SpatialDim> b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
    inline ippl::Vector<double, SpatialDim> crossProduct(ippl::Vector<double, SpatialDim> a,
                                                         ippl::Vector<double, SpatialDim> b) {
        return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    }

    /** @brief Largest stable step \f$\Delta t_{\max} = \sqrt{\varepsilon}\,(\sum_d
     * h_d^{-2})^{-1/2}\f$.
     * @details Deduction by the mode substitution of mediumPhasePerStep: every mode has a real
     * phase when \f$(\Delta t^2/4\varepsilon)K^2 \le 1\f$, and the largest \f$K^2\f$ is
     * \f$\sum_d 4/h_d^2\f$. This is the vacuum bound
     * [fallahi2020mithra20fullwavesimulation, section 3.1.3, Eqs. (3.19)-(3.20), p. 14] with
     * \f$c \to 1/n\f$, the rule for a medium [gedney2011IntroductionFdtd, section 3.5,
     * Eqs. (3.60), (3.62), p. 51]. The bound is strict: at equality the two roots of the
     * leapfrog recurrence coincide [kreiss2002SecondOrderWave, section 2, p. 1944]. Evaluated
     * as one square root of the quotient, deliberately not as a product of two roots. Serves
     * the step-bound check: the step of ChdrSolver equals the safety factor times this value.
     */
    inline double mediumMaxStableTimeStep(double eps, ippl::Vector<double, SpatialDim> h) {
        return std::sqrt(eps / (1 / (h[0] * h[0]) + 1 / (h[1] * h[1]) + 1 / (h[2] * h[2])));
    }

    /** @brief Spatial symbol \f$K_d = (2/h_d)\sin(k_dh_d/2)\f$ of one second difference.
     * @details Inserting \f$e^{\mathrm i k_d x}\f$ in \f$(f_{j+1} - 2f_j + f_{j-1})/h_d^2\f$ gives
     * \f$-(4/h_d^2)\sin^2(k_dh_d/2)\f$ (Deduction), the terms of
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.17), p. 14].
     */
    inline double discreteWaveNumber(double waveNumber, double h) {
        return 2 * std::sin(waveNumber * h / 2) / h;
    }
    /// @brief \f$K^2 = \sum_d K_d^2\f$, the Laplacian symbol with its sign removed.
    inline double discreteWaveNumberSquared(ippl::Vector<double, SpatialDim> waveNumber,
                                            ippl::Vector<double, SpatialDim> h) {
        double sum = 0;
        for (unsigned d = 0; d < SpatialDim; ++d) {
            const double symbol = discreteWaveNumber(waveNumber[d], h[d]);
            sum += symbol * symbol;
        }
        return sum;
    }

    /** @brief Phase per step \f$\vartheta\f$ of a single periodic mode of the medium update.
     * @details Inserting \f$e^{\mathrm i(\mathbf k\cdot\mathbf x - \vartheta k)}\f$ in the
     * source-free update gives \f$4\sin^2(\vartheta/2) = (\Delta t^2/\varepsilon)K^2\f$, so
     * \f[ \vartheta = 2\arcsin\Bigl(\frac{\Delta t}{2}\sqrt{K^2/\varepsilon}\Bigr), \f]
     * the dispersion relation [fallahi2020mithra20fullwavesimulation, section 3.1.3,
     * Eq. (3.17), p. 14] with \f$c \to 1/n\f$ (Deduction); the frequency of the mode is
     * \f$\vartheta/\Delta t\f$. Serves the eigenmode check: the measured cumulative phase per step
     * against this value, with the step the solver actually uses.
     */
    inline double mediumPhasePerStep(ippl::Vector<double, SpatialDim> waveNumber,
                                     ippl::Vector<double, SpatialDim> h, double dt, double eps) {
        return 2 * std::asin(dt / 2 * std::sqrt(discreteWaveNumberSquared(waveNumber, h) / eps));
    }

    /** @brief Factor from the scaled forcing of the update to the physical source ChdrSolver
     * accepts: \f$\rho = \varepsilon^2 s_{\varepsilon,0}\f$, \f$J_d =
     * \varepsilon\,s_{\varepsilon,d}\f$.
     * @details Dividing the medium equations by \f$-\varepsilon\f$ gives
     * \f$s_\varepsilon = (\rho/\varepsilon^2, \mathbf J/\varepsilon)\f$ in the positive source sign
     * of the update, \f$\partial_t^2\psi - \nabla^2\psi = (\rho, \mathbf J)\f$ at \f$\varepsilon =
     * 1\f$ (Deduction from [chew2014generalizedgauge, Eqs. (12)-(13)] with \f$\mu = 1\f$): the
     * charge is divided once by Gauss's law \f$\nabla\cdot(\varepsilon\mathbf E) = \rho\f$ and once
     * by the normalization of \f$\varepsilon\,\partial_t^2\phi\f$; the current enters Ampere's law
     * without a factor. Serves the manufactured-solution check, which computes \f$s_\varepsilon\f$
     * discretely and writes the physical sources.
     * @param component 0 for the charge, 1 to 3 for the current components.
     */
    inline double physicalSourceFactor(unsigned component, double eps) {
        return component == 0 ? eps * eps : eps;
    }

    /** @brief Frequency \f$\omega = \lvert\mathbf k\rvert/n\f$ of a source-free plane wave.
     * @details The transverse wave "has the velocity of light, or \f$1/\sqrt{\mu\varepsilon}\f$"
     * [chew2014generalizedgauge, section 2, text after Eq. (18), p. 3]; the source-free medium
     * equation gives \f$\omega^2 = \lvert\mathbf k\rvert^2/\varepsilon\f$ for every component.
     */
    inline double mediumFrequency(ippl::Vector<double, SpatialDim> waveNumber, double eps) {
        return std::sqrt(dotProduct(waveNumber, waveNumber) / eps);
    }

    /** @brief Scalar amplitude \f$\hat\phi = \mathbf k\cdot\hat{\mathbf a}/(\varepsilon\,\omega)\f$
     * of the wave \f$\mathbf A = \mathrm{Re}[\hat{\mathbf a}\,e^{\mathrm i(\mathbf k\cdot\mathbf x
     * - \omega t)}]\f$.
     * @details From the medium gauge \f$\nabla\cdot\mathbf A + \varepsilon\,\partial_t\phi = 0\f$
     * [chew2014generalizedgauge, Eq. (11), p. 3, \f$\mu = 1\f$]:
     * \f$\mathrm i\mathbf k\cdot\hat{\mathbf a} - \mathrm i\varepsilon\omega\hat\phi = 0\f$
     * (Deduction). Applied to the real and the imaginary part of \f$\hat{\mathbf a}\f$ separately;
     * for a driven wave pass its frequency \f$\Omega\f$. Serves the plane-wave and driven-wave
     * convergence checks.
     */
    inline double mediumScalarAmplitude(ippl::Vector<double, SpatialDim> waveNumber,
                                        ippl::Vector<double, SpatialDim> amplitudePart, double eps,
                                        double omega) {
        return dotProduct(waveNumber, amplitudePart) / (eps * omega);
    }

    /// @brief Real and imaginary parts of a complex vector amplitude.
    struct ComplexAmplitude {
        ippl::Vector<double, SpatialDim> real_m, imag_m;
    };

    /** @brief Amplitude \f$\hat{\mathbf a} = \mathbf e_1 + 0.5\mathrm i\,\mathbf e_2 + \mathrm
     * i\gamma\hat{\mathbf k}\f$ with a right-handed basis \f$(\mathbf e_1, \mathbf e_2,
     * \hat{\mathbf k})\f$ for any \f$\mathbf k\f$.
     * @details \f$\mathbf e_1 = \mathbf k\times\mathbf e_m/\lvert\mathbf k\times\mathbf
     * e_m\rvert\f$ with
     * \f$\mathbf e_m\f$ the coordinate axis least aligned with \f$\mathbf k\f$ (smallest
     * \f$\lvert k_m\rvert\f$, the first on ties), \f$\mathbf e_2 = \mathbf k\times\mathbf
     * e_1/\lvert\mathbf k\rvert\f$. Deduction: \f$\lvert\mathbf k\times\mathbf e_m\rvert^2 =
     * \lvert\mathbf k\rvert^2 - k_m^2 \ge
     * \tfrac23\lvert\mathbf k\rvert^2 > 0\f$, so the basis never degenerates, also for \f$\mathbf
     * k\f$ along an axis. The transverse part is a manufactured elliptic polarization; \f$\gamma =
     * 0\f$ gives \f$\phi = 0\f$, \f$\gamma = 1\f$ adds a longitudinal part, which with its
     * \f$\hat\phi\f$ is pure gauge and carries no field when \f$\varepsilon\omega^2 = \lvert\mathbf
     * k\rvert^2\f$ [chew2014generalizedgauge, section 7, Eqs. (76)-(80) and the text after Eq.
     * (80), p. 11].
     */
    inline ComplexAmplitude mediumWaveAmplitude(ippl::Vector<double, SpatialDim> waveNumber,
                                                double gamma) {
        unsigned least = 0;
        for (unsigned d = 1; d < SpatialDim; ++d)
            if (std::abs(waveNumber[d]) < std::abs(waveNumber[least]))
                least = d;
        ippl::Vector<double, SpatialDim> axis(0.0);
        axis[least]                                 = 1;
        const auto raw                              = crossProduct(waveNumber, axis);
        const double norm                           = std::sqrt(dotProduct(waveNumber, waveNumber));
        const ippl::Vector<double, SpatialDim> e1   = raw / std::sqrt(dotProduct(raw, raw));
        const ippl::Vector<double, SpatialDim> e2   = crossProduct(waveNumber, e1) / norm;
        const ippl::Vector<double, SpatialDim> imag = 0.5 * e2 + gamma * (waveNumber / norm);
        return {e1, imag};
    }

    /** @brief Physical source factor of a driven wave at frequency \f$\Omega\f$:
     * \f$\rho = \varepsilon(\lvert\mathbf k\rvert^2 - \varepsilon\Omega^2)\phi\f$,
     * \f$\mathbf J = (\lvert\mathbf k\rvert^2 - \varepsilon\Omega^2)\mathbf A\f$.
     * @details With \f$\hat\phi = \mathbf k\cdot\hat{\mathbf a}/(\varepsilon\Omega)\f$ the wave
     * obeys the medium gauge; inserting it in
     * \f$\partial_t^2\psi - \varepsilon^{-1}\nabla^2\psi = s_\varepsilon\f$ gives
     * \f$s_\varepsilon = (\lvert\mathbf k\rvert^2/\varepsilon - \Omega^2)\psi\f$, and
     * physicalSourceFactor converts it (Deduction). Gauss's law holds:
     * \f$\varepsilon\,\nabla\cdot\mathbf E = \varepsilon(\lvert\mathbf k\rvert^2\hat\phi -
     * \Omega\,\mathbf k\cdot\hat{\mathbf a}) = \hat\rho\f$. Serves the driven-wave convergence
     * check; at \f$\varepsilon = 1\f$ it is the vacuum factor \f$\lvert\mathbf k\rvert^2 -
     * \Omega^2\f$ for all components.
     * @param component 0 for the charge, 1 to 3 for the current components.
     */
    inline double drivenSourceFactor(unsigned component,
                                     ippl::Vector<double, SpatialDim> waveNumber, double eps,
                                     double omega) {
        return (component == 0 ? eps : 1.0)
               * (dotProduct(waveNumber, waveNumber) - eps * omega * omega);
    }

    /** @brief Factor \f$\varepsilon\f$ in the medium forms of the two constraint residuals: Gauss
     * \f$\varepsilon\,D\cdot\mathbf E - \rho\f$ and gauge \f$D\cdot\mathbf A +
     * \varepsilon\,\delta_t\phi\f$.
     * @details They discretize \f$\nabla\cdot(\varepsilon\mathbf E) = \rho\f$
     * [chew2014generalizedgauge, section 2, Eq. (4), p. 2, with \f$\mathbf D = \varepsilon\mathbf
     * E\f$] and \f$\nabla\cdot\mathbf A + \varepsilon\,\partial_t\phi = 0\f$ [ibid., Eq. (11), p.
     * 3,
     * \f$\mu = 1\f$]; the factor multiplies the divergence of E and the time difference of
     * \f$\phi\f$, the time weights of the reconstruction are unchanged. Serves the Gauss and
     * gauge order assertions of the convergence checks.
     */
    inline double constraintFactor(double eps) {
        return eps;
    }

    /** @brief Step of the dispersion record: \f$S\,\Delta t_{\max}(\varepsilon, \mathbf h)\f$, or
     * \f$S\,\Delta t_{\max}(1, \mathbf h)\f$ when a map that holds vacuum sets the step.
     * @details The vacuum-set step is the fraction \f$S/n\f$ of the medium bound; a smaller
     * fraction raises the phase error of this scheme [gedney2011IntroductionFdtd, section 3.6, p.
     * 56].
     */
    inline double recordTimeStep(double eps, ippl::Vector<double, SpatialDim> h, double safety,
                                 bool vacuumSet) {
        return safety * mediumMaxStableTimeStep(vacuumSet ? 1.0 : eps, h);
    }

    /** @brief Relative phase-velocity error \f$\delta v = (\vartheta/\Delta t)/(\lvert\mathbf
     * k\rvert/n) - 1\f$ of a mode; negative when the scheme is slower than the medium.
     * @details Leading term on cubic cells
     * \f$-(\lvert\mathbf k\rvert h)^2(\sum_d\hat k_d^4 - S_{\mathrm{med}}^2/3)/24\f$ (Deduction,
     * series of mediumPhasePerStep). Recorded against the cells per dielectric wavelength; no pass
     * rule.
     */
    inline double phaseVelocityError(ippl::Vector<double, SpatialDim> waveNumber,
                                     ippl::Vector<double, SpatialDim> h, double dt, double eps) {
        return mediumPhasePerStep(waveNumber, h, dt, eps) / dt / mediumFrequency(waveNumber, eps)
               - 1;
    }
    /// @brief Phase error per dielectric wavelength travelled, \f$\varphi_\lambda = 2\pi\,\delta
    /// v\f$ (rad).
    inline double phaseErrorPerWavelength(double velocityError) {
        return 2 * Kokkos::numbers::pi_v<double> * velocityError;
    }
}  // namespace chdr::solver
