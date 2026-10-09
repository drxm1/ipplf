/** @file Operators.h
 * @brief Discrete spatial operators and host reference calculations for the
 * potential-based solver.
 *
 * @details This file provides three groups of helpers:
 *
 * - Cell-level spatial operations. Component adapts one field component to
 *   scalar sampling; centeredDifference and curl are used by magnetic-field
 *   reconstruction and residual calculations.
 *
 * - Wave-operator evaluation. secondDifference and laplacian reproduce the
 *   Standard or NonStandard spatial stencil for Gauss-identity diagnostics,
 *   manufactured-source forcing, and Fourier stencil tests.
 *
 * - Host-side reference and setup calculations. The symbol helpers and
 *   discreteFrequency provide coefficients and frequencies for analytic test
 *   waves and dispersion checks. maxStableTimeStep supplies the bound used
 *   by ChdrSolver's timestep setter and safety-factor calculation.
 *
 * ChdrSolver::solve advances the potentials through the inherited
 * StandardFDTDSolver::step or NonStandardFDTDSolver::step. Those methods
 * implement their own update kernels. The laplacian helper here supports
 * diagnostics and tests; it is not called by those potential-update kernels.
 * Field reconstruction after the update does use the curl helper.
 *
 * Device-callable stencil functions evaluate individual cells. They borrow
 * samples or field views and do not allocate fields, launch kernels, or
 * exchange halos. Their callers prepare halos and perform parallel execution
 * and reductions. Symbol and timestep helpers run on the host.
 *
 * The individual function comments document their formulas, assumptions,
 * and literature sources.
 */
#pragma once

#include "Ippl.h"

#include <cmath>

#include "SpatialDimension.h"

namespace chdr::solver {
    /// Select the inherited periodic vacuum wave operator and its matching Fourier symbol.
    enum class Stencil {
        Standard,
        // The z-averaged one, see mithra2020...
        NonStandard
    };

    /** @brief Borrow one component of a device view as a scalar stencil accessor.
     * @details The wrapper copies a view handle and component index, never field ownership.
     * It lets the same scalar stencil operate on potential components and intermediate
     * monitor fields without embedding a solver or layout object in a device capture.
     */
    template <class View>
    struct Component {
        View view_m;
        unsigned component_m;
        KOKKOS_INLINE_FUNCTION double operator()(int i, int j, int k) const {
            return view_m(i, j, k)[component_m];
        }
    };

    template <class View>
    KOKKOS_INLINE_FUNCTION Component<View> component(View view, unsigned index) {
        return {view, index};
    }

    /** @brief Centred derivative; absolute mode replaces subtraction by absolute addition.
     * @details Along the selected axis, the stencil is
     * \f[
     * (D_cf)_j=\frac{f_{j+1}-f_{j-1}}{2h}.
     * \f]
     * To connect this stencil to its Fourier multiplier, substitute the sampled
     * mode \f$f_j=e^{\mathrm{i}\xi x_j}\f$ with \f$x_j=jh\f$:
     * \f[
     * \begin{aligned}
     * (D_cf)_j
     * &= \frac{e^{\mathrm{i}\xi(x_j+h)}-e^{\mathrm{i}\xi(x_j-h)}}{2h} \\
     * &= \frac{e^{\mathrm{i}\xi h}-e^{-\mathrm{i}\xi h}}{2h}\,f_j
     *  = \frac{\mathrm{i}}{h}\sin(\xi h)\,f_j.
     * \end{aligned}
     * \f]
     * This substitution is a deduction. The multiplier of the unchanged mode
     * is the expression given in [trefethen1996FiniteDifferenceSpectral,
     * Eq. (4.1.7), p. 151]; that equation states the Fourier multiplier,
     * rather than the stencil in terms of neighbouring samples.
     * @pre The requested neighbours are valid owned or synchronized halo cells.
     * @details The absolute mode supplies the local arithmetic scale used by the
     * residual checks [jeannerod2013ImprovedErrorBounds, Eq. (1.6)]. This is a
     * diagnostic scale, not a theorem-bound for FMA-contracted expressions.
     * A nested expression must also propagate the absolute evaluation of its inner terms.
     */
    template <class Sample>
    KOKKOS_INLINE_FUNCTION double centeredDifference(Sample f, int i, int j, int k, int axis,
                                                     double h, bool absolute = false) {
        // Axis determines along which axis the difference is made.
        const int di = axis == 0, dj = axis == 1, dk = axis == 2;
        const double plus = f(i + di, j + dj, k + dk), minus = f(i - di, j - dj, k - dk);
        return (absolute ? Kokkos::abs(plus) + Kokkos::abs(minus) : plus - minus) * (0.5 / h);
    }

    /** @brief Apply the three-point second derivative, or its absolute evaluation.
     * @tparam Absolute Select signed differentiation or absolute evaluation at compile time.
     * @details Along the selected axis, the stencil is
     * \f[
     * (D_{xx}f)_j=\frac{f_{j+1}-2f_j+f_{j-1}}{h^2}.
     * \f]
     * This is the spatial second-derivative formula written directly in
     * [fallahi2020mithra20fullwavesimulation, section 3.1.2,
     * Eqs. (3.11)-(3.13), p. 13], with h the spacing along that axis.
     *
     * For the sampled mode \f$f_j=e^{\mathrm{i}\xi jh}\f$, substitution gives
     * the Fourier multiplier (deduction):
     * \f[
     * (D_{xx}f)_j
     * =\frac{e^{\mathrm{i}\xi h}-2+e^{-\mathrm{i}\xi h}}{h^2}\,f_j
     * =-\frac{4}{h^2}\sin^2\!\left(\frac{\xi h}{2}\right)f_j.
     * \f]
     * [trefethen1996FiniteDifferenceSpectral, Eq. (5.1.9), p. 195] gives the
     * resulting dispersion relation for \f$u_{tt}=D_{xx}u\f$. Its positive
     * \f$\omega^2\f$ follows because the temporal second derivative also
     * contributes a minus sign for a time-harmonic mode.
     *
     * The centre's factor two is retained in the absolute branch; it is part of the
     * evaluated expression, not a field-dependent normalization chosen afterward.
     */
    template <bool Absolute = false, class Sample>
    KOKKOS_INLINE_FUNCTION double secondDifference(Sample f, int i, int j, int k, int axis,
                                                   double h) {
        const int di = axis == 0, dj = axis == 1, dk = axis == 2;
        const double plus = f(i + di, j + dj, k + dk), middle = f(i, j, k),
                     minus = f(i - di, j - dj, k - dk);
        // Absolute ? (|f_{j+1}| + 2|f_j| + |f_{j-1}|)/h^2
        //          : (f_{j+1} - 2f_j + f_{j-1})/h^2.
        // The absolute sum is the stencil's arithmetic scale; see centeredDifference().
        if constexpr (Absolute)
            return (Kokkos::abs(plus) + 2 * Kokkos::abs(middle) + Kokkos::abs(minus)) / (h * h);
        else
            return (plus - 2 * middle + minus) / (h * h);
    }

    /** @brief Standard discrete Laplacian or z-oriented NonStandard spatial wave operator.
     * @tparam Scheme Spatial stencil, fixed for each kernel specialization.
     * @tparam Absolute Select the absolute arithmetic scale instead of the signed operator.
     * @details Call sites are the Gauss-identity diagnostic (GaussIdentity via
     * monitorDetail::laplace), manufactured-source forcing, and Fourier stencil
     * tests. This helper does not advance the potential histories: ChdrSolver::solve()
     * calls the inherited StandardFDTDSolver::step() or NonStandardFDTDSolver::step(),
     * which implement their own update kernels. Diagnostics may run after a solve
     * and evaluate this operator over the cells, but are separate from that update.
     * Scheme and Absolute are template parameters so these choices are fixed
     * before entering the cell kernel; secondDifference receives the same Absolute.
     *
     * For Absolute=false, let \f$D_{dd}\f$ denote secondDifference()
     * along axis d, and define \f$Q_\ell=(D_{xx}f+D_{yy}f)_{i,j,\ell}\f$.
     * Here \f$\ell\f$ denotes the code's spatial z index k, not a wave number.
     * The real-space operator is
     * \f[
     * (L_hf)_{i,j,\ell}
     * = (D_{zz}f)_{i,j,\ell}+(1-2a)Q_\ell+a(Q_{\ell-1}+Q_{\ell+1}).
     * \f]
     * Standard sets \f$a=0\f$, giving the sum of the three second differences.
     * NonStandard uses weight = \f$\mathcal A\f$ from
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.25)]:
     * \f[
     * a=\mathcal A=\frac14\left(1+
     * \frac{0.02}{(h_z/h_x)^2+(h_z/h_y)^2}\right).
     * \f]
     * This is the stencil of MITHRA Eqs. (3.21)-(3.23): transverse x/y
     * differences are averaged across neighbouring z planes, while the z
     * second difference remains unaveraged. The paper averages before taking
     * transverse differences; this code averages afterward. Constant weights
     * and linearity make the two forms equal (deduction).
     * The distinguished direction is code z (axis 2), chosen as the beam axis
     * in MITHRA section 3.1.3. Fields may vary along all three axes; this is
     * not restricted to waves travelling along z.
     *
     * To connect the stencil to dispersion, substitute a spatial Fourier mode
     * with wave-number components \f$q_d\f$. Its second-difference and
     * z-averaging factors give (deduction from the stencil above)
     * \f[
     * \begin{aligned}
     * K_d &= \frac{2}{h_d}\sin\!\left(\frac{q_dh_d}{2}\right), \\
     * S_z &= (1-2a)+a(e^{-\mathrm{i}q_zh_z}+e^{\mathrm{i}q_zh_z})
     *      = 1-4a\sin^2\!\left(\frac{q_zh_z}{2}\right), \\
     * \lambda_h &= -\left[S_z(K_x^2+K_y^2)+K_z^2\right].
     * \end{aligned}
     * \f]
     * laplacianSymbol() returns this multiplier \f$\lambda_h\f$.
     * Combining it with the centred source-free time update
     * \f$D_{tt}f=c^2L_hf\f$ gives
     * \f[
     * \frac{4}{\Delta t^2}\sin^2\!\left(\frac{\omega\Delta t}{2}\right)
     * =c^2\left[S_z(K_x^2+K_y^2)+K_z^2\right].
     * \f]
     * Dividing by \f$4c^2\f$ gives MITHRA Eq. (3.24); setting \f$a=0\f$
     * gives Eq. (3.17). These are dispersion relations for the full time/space
     * scheme; this function supplies only its spatial operator. The solver
     * uses normalized \f$c=1\f$. Sources: [fallahi2020mithra20fullwavesimulation,
     * Eqs. (3.11)-(3.17), (3.21)-(3.25)]; IPPL NonStandardFDTDSolver::step
     * uses the same coefficient convention.
     *
     * Absolute evaluation instead follows the nested expression with absolute
     * samples and weights to supply the residual arithmetic scale described
     * by centeredDifference(). It is not generally the absolute value of L_h f;
     * the Fourier and dispersion formulas above apply only to Absolute=false.
     * @pre The caller has prepared the one-cell faces and, for NonStandard, the
     * mixed transverse/longitudinal neighbours used by the smoothed second differences.
     * This is a read-only device operation: ownership and halo exchanges remain outside it.
     */
    template <Stencil Scheme, bool Absolute = false, class Sample>
    KOKKOS_INLINE_FUNCTION double laplacian(Sample f, int i, int j, int k,
                                            ippl::Vector<double, SpatialDim> h) {
        static_assert(SpatialDim == 3, "laplacian requires three spatial dimensions (x, y, z)");
        double result = secondDifference<Absolute>(f, i, j, k, 2, h[2]);
        double weight = 0.0;
        if constexpr (Scheme == Stencil::NonStandard) {
            const double ratio = h[2] * h[2] / (h[0] * h[0]) + h[2] * h[2] / (h[1] * h[1]);
            weight             = 0.25 * (1 + 0.02 / ratio);
        }
        // Axis 2 is z; the loop below covers the transverse x and y axes.
        for (int d = 0; d < 2; ++d) {
            const double centre = secondDifference<Absolute>(f, i, j, k, d, h[d]);
            result += (Absolute ? Kokkos::abs(1 - 2 * weight) : 1 - 2 * weight) * centre;

            // For NonStandard solver, z is treated differently.
            // The x- and y-difference contributions are averaged across neighbouring z-planes, and
            // the second difference in z remains unaveraged.
            if constexpr (Scheme == Stencil::NonStandard) {
                result += weight
                          * (secondDifference<Absolute>(f, i, j, k - 1, d, h[d])
                             + secondDifference<Absolute>(f, i, j, k + 1, d, h[d]));
            }
        }
        return result;
    }

    /** @brief Real factor of the centered-derivative symbol, omitting the imaginary unit.
     * @details Source: [trefethen1996FiniteDifferenceSpectral, Eq. (4.1.7),
     * p. 151; Eq. (5.1.8), p. 195].
     */
    inline double centeredSymbol(double waveNumber, double h) {
        // std::sin instead of Kokkos::sin because the call sites are host-only reference
        // calculations
        return std::sin(waveNumber * h) / h;
    }
    /// Signed spatial factor in the Laplacian symbol [fallahi2020mithra20fullwavesimulation, Eq.
    /// (3.17)].
    inline double waveSymbol(double waveNumber, double h) {
        // std::sin instead of Kokkos::sin because the call sites are host-only reference
        return 2 * std::sin(waveNumber * h / 2) / h;
    }
    /// One-step centred-time symbol [fallahi2020mithra20fullwavesimulation, Eqs. (3.17),
    /// (3.53)-(3.55)].
    inline double timeSymbol(double omega, double dt) {
        // std::sin instead of Kokkos::sin because the call sites are host-only reference
        return 2 * std::sin(omega * dt / 2) / dt;
    }
    /// Match real-space smoothing weights [fallahi2020mithra20fullwavesimulation, Eqs.
    /// (3.21)-(3.25)].
    inline double laplacianSymbol(ippl::Vector<double, SpatialDim> waveNumber /*k*/,
                                  ippl::Vector<double, SpatialDim> h, Stencil stencil) {
        const double kx = waveSymbol(waveNumber[0], h[0]), ky = waveSymbol(waveNumber[1], h[1]);
        const double kz     = waveSymbol(waveNumber[2], h[2]);
        const double ratio  = h[2] * h[2] / (h[0] * h[0]) + h[2] * h[2] / (h[1] * h[1]);
        const double weight = stencil == Stencil::Standard ? 0.0 : 0.25 * (1 + 0.02 / ratio);
        const double sine   = std::sin(waveNumber[2] * h[2] / 2);
        return -((1 - 4 * weight * sine * sine) * (kx * kx + ky * ky) + kz * kz);
    }
    /** @brief Return the positive leapfrog frequency on its principal stable branch.
     * @pre The supplied mode and step lie inside the real dispersion branch.
     * @details This is the dispersion relation [fallahi2020mithra20fullwavesimulation,
     * Eqs. (3.17), (3.24)],
     * solved for frequency. It is a host reference calculation, not a field update.
     */
    inline double discreteFrequency(ippl::Vector<double, SpatialDim> waveNumber,
                                    ippl::Vector<double, SpatialDim> h, double dt,
                                    Stencil stencil) {
        return 2 * std::asin(dt * std::sqrt(-laplacianSymbol(waveNumber, h, stencil)) / 2) / dt;
    }

    /** @brief CFL upper bound; equality is excluded for the checked Standard setter.
     * @details Standard: fallahi2020mithra20fullwavesimulation, Eqs. (3.19)-(3.20).
     * NonStandard: deduction from Eq. (3.24), under the base constructor's aspect-ratio guard.
     */
    inline double maxStableTimeStep(ippl::Vector<double, SpatialDim> h, Stencil stencil) {
        return stencil == Stencil::NonStandard
                   ? h[2]
                   : 1 / std::sqrt(1 / (h[0] * h[0]) + 1 / (h[1] * h[1]) + 1 / (h[2] * h[2]));
    }
}  // namespace chdr::solver
