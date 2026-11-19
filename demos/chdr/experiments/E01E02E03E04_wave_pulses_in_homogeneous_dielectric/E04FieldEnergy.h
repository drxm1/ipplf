/** @file E04FieldEnergy.h
 * @brief E04.meas01 energy of the staggered and the collocated solver configurations.
 */
#pragma once

#include <Kokkos_Core.hpp>

#include <mpi.h>

#include "core/base/FieldContainer.h"
#include "core/solver/operators/Reconstruction.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03E04Setup.h"
#include "opt-in/diagnostics/AlignedFieldEnergy.h"

namespace chdr::experiments::e04 {
    namespace fieldEnergyDetail {
        using View = EMFieldsType::view_type;

        /// \f$\sum(\varepsilon|E|^2+|\bar B|^2)\f$ over the owned cells of one rank.
        struct CollocatedSquares {
            View::const_type electric_m, magnetic_m;
            double epsilon_m;

            /// Add the squares of one owned cell to the partial sum of this rank.
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& sum) const {
                for (unsigned d = 0; d < 3; ++d)
                    sum += epsilon_m * electric_m(i, j, k)[d] * electric_m(i, j, k)[d]
                           + magnetic_m(i, j, k)[d] * magnetic_m(i, j, k)[d];
            }
        };
        /// \f$\sum|\mathbf F|^2\f$ of one field over the owned cells of one rank.
        struct Squares {
            View::const_type field_m;

            /// Add the squares of one owned cell to the partial sum of this rank.
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& sum) const {
                for (unsigned d = 0; d < 3; ++d)
                    sum += field_m(i, j, k)[d] * field_m(i, j, k)[d];
            }
        };
    }  // namespace fieldEnergyDetail

    /** @brief The parts of the energy of one state, each with the weight \f$\Delta V/2\f$.
     * @details electric_m is \f$\varepsilon\sum|\mathbf E|^2\f$ and magnetic_m is
     * \f$\sum|\bar{\mathbf B}|^2\f$; their sum is \f$W_h\f$ of E04FieldEnergy. rawMagnetic_m is
     * \f$\sum(|\mathbf B^n|^2+|\mathbf B^{n-1}|^2)/2\f$ with the unaveraged curls of the two
     * potential slots. Its difference to magnetic_m is the effect of the average:
     * \f$|(a+b)/2|^2=(|a|^2+|b|^2)/2-|(a-b)/2|^2\f$.
     */
    struct EnergyParts {
        double electric_m, magnetic_m, rawMagnetic_m;
    };

    /** @brief The energy of E04.meas01 and the magnetic field it uses, for every configuration.
     * @details The staggered configurations use R2's averaged B and alignedFieldEnergy,
     * unchanged from E04. The collocated R1 configuration has \f$\mathbf E\f$ at
     * \f$t_{n-1/2}\f$ up to R1's first-order term (LevelContract) and \f$\mathbf B^n\f$ at
     * \f$t_n\f$. The energy aligns them with
     * \f$\bar{\mathbf B}=(\mathbf B^n+\mathbf B^{n-1})/2\f$ at \f$t_{n-1/2}\f$, as R2 does on
     * the staggered layout. The collocated R2Prime configuration has \f$\mathbf E\f$ at
     * \f$t_{n-1/2}\f$ alone and takes the same branch: the averaged B computed here equals
     * the B of its reconstruction. The energy is
     * \f$W_h=\tfrac{\Delta V}{2}\sum(\varepsilon|\mathbf E|^2+|\bar{\mathbf B}|^2)\f$ over owned
     * cell centers [ellingson2020ElectromagneticsVol2, Eqs. (3.20), (3.22)], with unit
     * permeability and a uniform \f$\varepsilon\f$.
     * @pre Both potential slots have current halos and closed open-face ghosts, as after
     * loading and after every solve.
     */
    class E04FieldEnergy {
        using Solver = e03::SphericalSimulation::Solver;

    public:
        /// Allocate the averaged B only for the collocated layout, where the solver lacks it.
        E04FieldEnergy(Solver& state, double epsilon)
            : epsilon_m(epsilon)
            , staggered_m(state.spatialDiscretization() == SpatialDiscretization::Staggered) {
            if (!staggered_m)
                averaged_m.initialize(state.A_n.get_mesh(), state.A_n.getLayout());
        }

        /// Whether the run uses the staggered layout and R2's averaged B.
        bool staggered() const { return staggered_m; }

        /** @brief The magnetic field at the electric comparison time.
         * @details On the collocated layout this is \f$(\mathbf B^n+\mathbf B^{n-1})/2\f$:
         * reconstructMagnetic averages the central curls of both potential slots, the curl that
         * IPPL's StandardFDTDSolver::evaluate_EB applies to the current slot.
         */
        EMFieldsType& magnetic(Solver& state) {
            if (staggered_m)
                return state.magneticField();
            solver::reconstructMagnetic(state.A_n, state.A_nm1, averaged_m, true);
            return averaged_m;
        }

        /// \f$W_h\f$ of the current state, identical on every rank.
        double energy(Solver& state, const material::MaterialOwner& owner) {
            if (staggered_m)
                return diagnostics::alignedFieldEnergy(state, owner);
            double local = 0, total = 0;
            Kokkos::parallel_reduce(
                "E04 collocated field energy", state.A_n.getFieldRangePolicy(),
                fieldEnergyDetail::CollocatedSquares{state.electricField().getView(),
                                                     magnetic(state).getView(), epsilon_m},
                local);
            MPI_Allreduce(&local, &total, 1, MPI_DOUBLE, MPI_SUM, state.A_n.getLayout().comm);
            const auto h = state.A_n.get_mesh().getMeshSpacing();
            return 0.5 * h[0] * h[1] * h[2] * total;
        }

        /** @brief The electric part, the magnetic part and the mean energy of the two raw curls.
         * @details The raw curls use the curl of the placement: rawStaggeredMagnetic on the
         * staggered layout and the central curl of reconstructMagnetic on the collocated one.
         * A uniform \f$\varepsilon\f$ is assumed, as in energy() for the collocated layout.
         */
        EnergyParts parts(Solver& state) {
            if (raw_m.getView().size() == 0)
                raw_m.initialize(state.A_n.get_mesh(), state.A_n.getLayout());
            double raw = 0;
            for (auto* slot : {&state.A_n, &state.A_nm1}) {
                if (staggered_m)
                    solver::rawStaggeredMagnetic(*slot, raw_m);
                else
                    solver::reconstructMagnetic(*slot, *slot, raw_m, false);
                raw += 0.5 * squareSum(state, raw_m);
            }
            return {epsilon_m * squareSum(state, state.electricField()),
                    squareSum(state, magnetic(state)), raw};
        }

    private:
        /// \f$(\Delta V/2)\sum|\mathbf F|^2\f$ over the owned cells of all ranks.
        double squareSum(Solver& state, EMFieldsType& field) const {
            double local = 0, total = 0;
            Kokkos::parallel_reduce("E04 energy part", state.A_n.getFieldRangePolicy(),
                                    fieldEnergyDetail::Squares{field.getView()}, local);
            MPI_Allreduce(&local, &total, 1, MPI_DOUBLE, MPI_SUM, state.A_n.getLayout().comm);
            const auto h = state.A_n.get_mesh().getMeshSpacing();
            return 0.5 * h[0] * h[1] * h[2] * total;
        }

        double epsilon_m;
        bool staggered_m;
        EMFieldsType averaged_m, raw_m;
    };
}  // namespace chdr::experiments::e04
