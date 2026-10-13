/** @file HomogeneousUpdate.h
 * @brief Uniform dielectric update using IPPL's existing vacuum kernel.
 */
#pragma once

#include <cmath>

#include "MaterialStepContext.h"
#include "Operators.h"

namespace chdr::solver {
    /** @brief Evolve homogeneous Lorenz-gauge potentials with physical sources.
     * @details With unit permeability, Chew's equations give
     * \f$\partial_t^2\psi=\varepsilon^{-1}\nabla^2\psi+
     * (\rho/\varepsilon^2,\mathbf J/\varepsilon)\f$
     * [chew2014generalizedgauge, section 2, Eqs. (11)-(13)]. Substitution of
     * \f$\tau=\Delta t/\sqrt\varepsilon\f$ makes the discrete step the vacuum
     * step with source \f$(\rho/\varepsilon,\mathbf J)\f$ (deduction).
     * Reusing that kernel preserves its operation order, including contraction;
     * at unit permittivity every added division is exact. The source copy uses
     * four doubles per allocated cell and no additional halo exchange.
     * Here epsilon is relative permittivity and psi contains the scalar potential
     * followed by the three vector-potential components. The caller continues to
     * supply physical charge/current and the physical step; only the inherited
     * kernel sees the transformed source and step.
     */
    class HomogeneousUpdate final : public MaterialUpdate {
    public:
        /** @brief Allocate one private source copy; retain only the scalar permittivity.
         * @pre The factory has established a finite homogeneous permittivity at least one.
         * @details The cell material field is never read during evolution: one scalar
         * describes the already validated uniform map. The source copy owns its data
         * but borrows the prototype's mesh/layout; those objects must outlive it.
         * Keeping this field as a member lets successive steps reuse the allocation
         * and keeps its storage alive after a device kernel has been dispatched.
         */
        HomogeneousUpdate(double epsilon, AJFourFieldsType& prototype)
            : epsilon_m(epsilon)
            , refractiveIndex_m(std::sqrt(epsilon)) {
            scaledSource_m.initialize(prototype.get_mesh(), prototype.getLayout());
        }

        /** @brief Evolve next() and its halos with the homogeneous-medium coefficients.
         * @details The internal step contributes a factor
         * \f$\tau^2=\Delta t^2/\varepsilon\f$ to both the Laplacian and source.
         * The extra charge division in copySource therefore gives the required
         * \f$\Delta t^2\rho/\varepsilon^2\f$; unchanged current receives
         * \f$\Delta t^2\mathbf J/\varepsilon\f$ (deduction from
         * [chew2014generalizedgauge, Eqs. (12)-(13)]). Do not pre-scale the caller's
         * source or apply another current division. baseStep completes halos and
         * boundaries; history rotation and E/B reconstruction occur outside this call.
         */
        void advance(MaterialStepContext& context) override {
            copySource(context.source());
            context.baseStep(context.dt() / refractiveIndex_m, scaledSource_m);
        }

        /** @brief Return the Standard bound multiplied by the refractive index.
         * @details \f$\Delta t_{\max}=\sqrt\varepsilon
         * (\sum_d h_d^{-2})^{-1/2}\f$; the endpoint is excluded by the shell
         * [fallahi2020mithra20fullwavesimulation, Eq. (3.20), with speed
         * \f$c/\sqrt\varepsilon\f$; kreiss2002SecondOrderWave, section 2].
         * spacing is the physical cell spacing. The return value bounds the caller's
         * physical step, before the internal time rescaling used by the vacuum kernel.
         */
        double maxStableTimeStep(const ippl::Vector<double, SpatialDim>& spacing) const override {
            return refractiveIndex_m * chdr::solver::maxStableTimeStep(spacing, Stencil::Standard);
        }

        /** @brief Copy owned physical sources and divide charge by permittivity.
         * @details The inherited step reads source values only on owned cells, so
         * source ghosts are unnecessary. Views and epsilon are captured by value;
         * this kernel and the inherited update use the field's execution space.
         * Components one through three remain physical current. Only this private
         * copy is modified, so a source reused by the caller is never scaled twice.
         * The owned-cell policy excludes ghosts; StandardFDTDSolver::step reads its
         * source at the destination cell, without neighbour-source accesses.
         * Public visibility permits CUDA extended-lambda compilation.
         */
        void copySource(const AJFourFieldsType& source) {
            const auto input     = source.getView();
            const auto output    = scaledSource_m.getView();
            const double epsilon = epsilon_m;
            // TODO: This later has to be benchmarked for the computational and memory overhead it
            // might make. This is called in every step so this might be too expensive.
            // For now we'll leave the implementation like this but beware that it might be
            // expensive.
            Kokkos::parallel_for(
                "Homogeneous source copy", source.getFieldRangePolicy(),
                KOKKOS_LAMBDA(int i, int j, int k) {
                    auto value = input(i, j, k);
                    value[0] /= epsilon;
                    output(i, j, k) = value;
                });
        }

    private:
        double epsilon_m, refractiveIndex_m;
        AJFourFieldsType scaledSource_m;
    };
}  // namespace chdr::solver
