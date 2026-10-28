/** @file VacuumEnergy.h
 * @brief Field energies of the vacuum wave: the raw-field sums and the conserved
 * cross-level energy of the three-level update.
 * @details Diagnostic only: no acceptance check uses these values. In the normalized
 * units of this demo (\f$c=\varepsilon_0=\mu_0=1\f$, so \f$\mathbf H=\mathbf B\f$) the
 * continuum energies \f$W_e=\tfrac12\int\varepsilon\lvert\mathbf E\rvert^2\,dv\f$ and
 * \f$W_m=\tfrac12\int\mu\lvert\mathbf H\rvert^2\,dv\f$
 * [ellingson2020ElectromagneticsVol2, section 3.1, Eqs. (3.20), (3.22)] are approximated by
 * the owned-cell sums
 * \f[W_E=\frac{\Delta V}{2}\sum_{\mathrm{owned}}\lvert\mathbf E\rvert^2,\qquad
 *    W_B=\frac{\Delta V}{2}\sum_{\mathrm{owned}}\lvert\mathbf B\rvert^2,\qquad
 *    \Delta V=h_xh_yh_z.\f]
 * All components are the collocated cell-centred raw arrays of
 * FDTDSolverBase::evaluate_EB. The raw E belongs to \f$t_n-\Delta t/2\f$ and the raw B to
 * \f$t_n\f$ (README, "Reference, storage and time labels"), so \f$W_E+W_B\f$ adds energies
 * of two different times and is not conserved exactly.
 *
 * The Standard update \f$v^{n+1}-2v^n+v^{n-1}=\Delta t^2L_hv^n\f$,
 * \f$L_h=\sum_dD_d^-D_d^+\f$, applied to each potential component \f$v=\mathcal A_c\f$
 * without source, conserves the cross-level energy
 * \f$L(t_{n+1},L_h)=\lVert v^{n+1}-v^n\rVert^2-\Delta t^2\langle v^{n+1},L_hv^n\rangle\f$
 * of the explicit scheme (2.8) [kreiss2002SecondOrderWave, section 2, \f$L(t_n,A)\f$ and
 * Eq. (2.10); not \f$L_1\f$, which belongs to the semi-implicit scheme (2.11) and reduces
 * to \f$L\f$ for \f$A_1=0\f$]. Summation by parts on the
 * periodic grid and the factor \f$\Delta V/(2\Delta t^2)\f$ give, for the pair
 * \f$(\mathcal A^n,\mathcal A^{n-1})\f$ held after \f$n\f$ solves (time label
 * \f$t_n-\Delta t/2\f$),
 * \f[W_\times=\frac{\Delta V}{2}\sum_{\mathrm{owned}}\sum_{c=0}^{3}\Bigl[
 *    \Bigl(\frac{\mathcal A_c^n-\mathcal A_c^{n-1}}{\Delta t}\Bigr)^2
 *    +\sum_{d}\bigl(D_d^+\mathcal A_c^n\bigr)\bigl(D_d^+\mathcal A_c^{n-1}\bigr)\Bigr],
 *    \qquad D_d^+u=\frac{u_{\mathbf i+\mathbf e_d}-u_{\mathbf i}}{h_d}.\f]
 * Kreiss et al. show positive definiteness for \f$A=A^*<0\f$ when
 * \f$(\max_j\lvert\lambda_j\rvert+\min_j\lvert\lambda_j\rvert)\Delta t^2\le4\f$ (ibid.); \f$L_h\f$
 * has the eigenvalue zero (constant mode), so here the energy is positive semidefinite
 * under \f$\max_j\lvert\lambda_j\rvert\Delta t^2\le4\f$. It is an invariant of the update for each
 * component, the scalar potential included; it is not the physical field energy in general.
 * The forward differences read one ghost layer in \f$+d\f$ of both histories. These ghosts
 * are current: FDTDSolverBase::finishPotentialStep fills the halo of \f$\mathcal A^{n+1}\f$
 * before timeShift copies whole views, and initializePotentials fills both initial halos.
 * No field is changed, copied to the host, or halo-exchanged.
 */
#pragma once

#include <Kokkos_Core.hpp>

#include <mpi.h>
#include <type_traits>

namespace chdr::apps::vacuumWave {
    /// @brief Global energies of one measured state; the parts refer to different times.
    struct EnergySample {
        double electric_m  = 0.0;  ///< \f$W_E\f$ at the raw-E time \f$t_n-\Delta t/2\f$.
        double magnetic_m  = 0.0;  ///< \f$W_B\f$ at the raw-B time \f$t_n\f$.
        double total_m     = 0.0;  ///< \f$W_E+W_B\f$, a sum over two times.
        double invariant_m = 0.0;  ///< \f$W_\times\f$ of \f$(\mathcal A^n,\mathcal A^{n-1})\f$.
    };

    namespace energy_detail {
        /// @brief Owned-cell sums before the factor \f$\Delta V/2\f$.
        struct Sums {
            double electric_m  = 0.0;
            double magnetic_m  = 0.0;
            double invariant_m = 0.0;
        };

        /// @brief One fused owned-cell reduction of \f$\lvert\mathbf E\rvert^2\f$,
        /// \f$\lvert\mathbf B\rvert^2\f$ and the cross-level terms of \f$W_\times\f$.
        template <class FieldView, class PotentialView>
        struct EnergyKernel {
            using value_type = Sums;
            FieldView electric_m;
            FieldView magnetic_m;
            PotentialView current_m;   ///< \f$\mathcal A^n\f$ with current ghosts.
            PotentialView previous_m;  ///< \f$\mathcal A^{n-1}\f$ with current ghosts.
            double inverseDt_m;
            Kokkos::Array<double, 3> inverseSpacing_m;

            /// @brief Identity of the sum reduction.
            KOKKOS_INLINE_FUNCTION void init(value_type& value) const { value = Sums{}; }

            /// @brief Combine two partial sums.
            KOKKOS_INLINE_FUNCTION void join(value_type& out, const value_type& in) const {
                out.electric_m += in.electric_m;
                out.magnetic_m += in.magnetic_m;
                out.invariant_m += in.invariant_m;
            }

            /// @brief \f$(D_d^+\mathcal A_c^n)(D_d^+\mathcal A_c^{n-1})\f$ in one cell.
            KOKKOS_INLINE_FUNCTION double gradientProduct(int i, int j, int k, unsigned d,
                                                          int c) const {
                const int di = d == 0, dj = d == 1, dk = d == 2;
                const double a = current_m(i + di, j + dj, k + dk)[c] - current_m(i, j, k)[c];
                const double b = previous_m(i + di, j + dj, k + dk)[c] - previous_m(i, j, k)[c];
                return a * b * inverseSpacing_m[d] * inverseSpacing_m[d];
            }

            /// @brief Bracket of \f$W_\times\f$ summed over the four potential components.
            KOKKOS_INLINE_FUNCTION double crossLevelTerm(int i, int j, int k) const {
                double term = 0.0;
                for (int c = 0; c < 4; ++c) {
                    const double rate =
                        (current_m(i, j, k)[c] - previous_m(i, j, k)[c]) * inverseDt_m;
                    term += rate * rate;
                    for (unsigned d = 0; d < 3; ++d)
                        term += gradientProduct(i, j, k, d, c);
                }
                return term;
            }

            /// @brief Add the squared E and B components and the cross-level terms of one cell.
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, value_type& out) const {
                for (int c = 0; c < 3; ++c) {
                    const double e = electric_m(i, j, k)[c];
                    const double b = magnetic_m(i, j, k)[c];
                    out.electric_m += e * e;
                    out.magnetic_m += b * b;
                }
                out.invariant_m += crossLevelTerm(i, j, k);
            }
        };

        /// @brief Cell volume \f$h_xh_yh_z\f$ of the field mesh; the spacing need not be cubic.
        template <class Field>
        double cellVolume(const Field& field) {
            const auto spacing = field.get_mesh().getMeshSpacing();
            return spacing[0] * spacing[1] * spacing[2];
        }

        /// @brief Kernel over the raw fields and the two newest potential histories.
        template <class Solver, class Fields>
        auto makeKernel(Solver& solver, Fields& fields) {
            using FieldView     = std::decay_t<decltype(fields.getE().getView())>;
            using PotentialView = std::decay_t<decltype(solver.A_n.getView())>;
            const auto spacing  = solver.A_n.get_mesh().getMeshSpacing();
            return EnergyKernel<FieldView, PotentialView>{
                fields.getE().getView(), fields.getB().getView(),
                solver.A_n.getView(),    solver.A_nm1.getView(),
                1.0 / solver.getDt(),    {1.0 / spacing[0], 1.0 / spacing[1], 1.0 / spacing[2]}};
        }
    }  // namespace energy_detail

    /** @brief Measure \f$W_E\f$, \f$W_B\f$, their sum and \f$W_\times\f$ over all owned cells.
     * @param solver Periodic Standard solver after initialization or a completed solve().
     * @param fields Allocated E/B after evaluate_EB(); not modified.
     * @return The same global energies on every rank.
     * @pre All layout ranks call in matching order; E, B and the potentials share mesh and
     * layout; source-free run, as the cross-level invariant requires.
     * @details Each rank sums its owned cells in one Kokkos reduction; one MPI sum
     * combines the three partial sums; the factor \f$\Delta V/2\f$ is applied last.
     */
    template <class Solver, class Fields>
    EnergySample measureEnergy(Solver& solver, Fields& fields) {
        energy_detail::Sums local;
        Kokkos::parallel_reduce("Vacuum field energy", fields.getE().getFieldRangePolicy(),
                                energy_detail::makeKernel(solver, fields), local);
        Kokkos::fence();
        double sums[3] = {local.electric_m, local.magnetic_m, local.invariant_m};
        MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_DOUBLE, MPI_SUM, fields.getFL().comm);
        const double factor = 0.5 * energy_detail::cellVolume(fields.getE());
        return {factor * sums[0], factor * sums[1], factor * (sums[0] + sums[1]), factor * sums[2]};
    }
}  // namespace chdr::apps::vacuumWave
