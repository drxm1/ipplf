/** @file AlignedFieldEnergy.h
 * @brief Owned-cell electromagnetic energy from time-aligned staggered R2 fields.
 */
#pragma once
#include <stdexcept>

#include "core/base/SpatialDimension.h"
#include "core/material/O1VSCoefficients.h"
#include "core/solver/ChdrSolver.h"

static_assert(
    chdr::SpatialDim == 3,
    "AlignedFieldEnergy.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::diagnostics {
    using solver::ChdrSolver;

    namespace alignedEnergyDetail {
        /// Sum field squares at their staggered locations on the field execution space.
        struct Reduction {
            EMFieldsType::view_type::const_type electric_m, magnetic_m;
            Kokkos::Array<material::ScalarFieldConstView, 3> epsilon_m;
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& sum) const {
                const auto e = electric_m(i, j, k), b = magnetic_m(i, j, k);
                for (unsigned d = 0; d < 3; ++d)
                    sum += epsilon_m[d](i, j, k) * e[d] * e[d] + b[d] * b[d];
            }
        };
    }  // namespace alignedEnergyDetail

    /** @brief Estimate total electromagnetic energy from the reconstructed staggered fields.
     * @details In normalized units with unit permeability, this approximates the continuum
     * energy [ellingson2020ElectromagneticsVol2, section 3.1, Eqs. (3.20), (3.22)] by
     * \f[W_h=\frac{\Delta V}{2}\sum_{\mathrm{owned},d}\varepsilon_d E_d^2
     *       +\frac{\Delta V}{2}\sum_{\mathrm{owned},d}\overline B_d^{\,2},
     *       \qquad\Delta V=h_xh_yh_z.\f]
     * Electric components occupy positive cell faces; magnetic components occupy edges
     * (LevelContract). Each square is sampled at its own component position, with harmonic
     * face permittivity at the electric position. These are separate volume sums combined
     * in one loop, not a local energy density at a common point. Component-wise staggered
     * sums follow the discrete-inner-product construction in
     * [sakkaplangkul2021YeeDispersiveMedia, section 3, Eq. (30)]. No spatial interpolation occurs.
     *
     * "Aligned" refers only to time: R2 has already formed
     * \f$\overline{\mathbf B}^{n-1/2}=(\mathbf B^n+\mathbf B^{n-1})/2\f$
     * at the electric-field time; see reconstructStaggered and
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.50)-(3.52)].
     * All owned samples receive the full volume weight, including Mur boundary samples;
     * no boundary-specific quadrature weights or boundary fluxes are included. This
     * diagnostic is distinct from the cross-level invariant in FieldEnergy.h and is not
     * asserted to be exactly conserved. MPI reduction counts each owned sample once.
     * @param state A ChdrSolver in the staggered reading with R2, its fields reconstructed for
     * the current state. Every rank of its layout calls this function.
     * @param owner The solver's MaterialOwner; the function does not check this identity. Its
     * harmonic face permittivities are derived on the first call (collective) and cached.
     * @return \f$W_h\f$, identical on every rank.
     */
    template <class Base>
    double alignedFieldEnergy(ChdrSolver<Base>& state, const material::MaterialOwner& owner) {
        if (state.spatialDiscretization() != SpatialDiscretization::Staggered
            || state.ebReconstructionConvention() != EBReconstructionConvention::R2)
            throw std::invalid_argument("Aligned field energy requires staggered R2 fields");
        const alignedEnergyDetail::Reduction reduction{
            state.electricField().getView(),
            state.magneticField().getView(),
            {owner.derive(material::HarmonicFace<0>{}), owner.derive(material::HarmonicFace<1>{}),
             owner.derive(material::HarmonicFace<2>{})}};
        double local = 0, total = 0;
        Kokkos::parallel_reduce("Time-aligned field energy", state.A_n.getFieldRangePolicy(),
                                reduction, local);
        MPI_Allreduce(&local, &total, 1, MPI_DOUBLE, MPI_SUM, state.A_n.getLayout().comm);
        const auto h = state.A_n.get_mesh().getMeshSpacing();
        return 0.5 * h[0] * h[1] * h[2] * total;
    }
}  // namespace chdr::diagnostics
