/** @file FieldEnergy.h
 * @brief Compensated O1-VS energy sums: periodic invariants or open-box diagnostics.
 */
#pragma once
#include <stdexcept>
#include <vector>

#include "core/base/SpatialDimension.h"
#include "core/formulation/O1VSOperator.h"
#include "core/solver/ChdrSolver.h"

static_assert(chdr::solver::SpatialDim == 3,
              "FieldEnergy.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    /** @brief Whole-domain sums for one discrete potential invariant, not per-cell values.
     * @details vacuum_m and dielectric_m partition the owned cell-centre, face and edge
     * contributions by material region; total() adds both regions. absolute_m is
     * the sum of absolute elementary terms before cancellation, a numerical
     * cancellation scale rather than an additional physical energy.
     */
    struct EnergyParts {
        double vacuum_m, dielectric_m, absolute_m;
        double total() const { return vacuum_m + dielectric_m; }
    };
    /** @brief Whole-domain vector and scalar leapfrog invariants and their time labels.
     * @details potentialEnergy() combines owned contributions across the communicator
     * and returns the same global record on every participating MPI rank. No per-cell
     * energy array is stored here. After solve k, the vector invariant is labelled
     * \f$(k-1/2)\Delta t\f$ and the scalar invariant \f$k\Delta t\f$.
     * These are the discrete cross-level invariants defined in energyDetail::Reduction,
     * not generally the instantaneous electromagnetic field energy.
     * The accumulation proceeds as follows:
     * 1. Compute contributions from each owned cell: its centre and the associated staggered
     * faces and edges.
     * 2. Sum locally on each rank, excluding halo cells from the ownership count.
     * 3. Combine all ranks and apply the factor \f$\Delta V/2\f$.
     */
    struct EnergyRecord {
        EnergyParts vector_m, scalar_m;
        double vectorTime_m, scalarTime_m;
    };
    namespace energyDetail {
        /** @brief Six compensated sums, each stored as its leading value and correction.
         * @details TwoSum retains the rounding residual of each addition
         * [ogita2005AccurateSumDot, Algorithm 3.1]. The join adds both parts;
         * compensation reduces, but does not remove, partition-order dependence.
         */
        struct Sums {
            double values_m[12]{};
            KOKKOS_INLINE_FUNCTION void add(unsigned c, double value) {
                const double old = values_m[2 * c], sum = old + value, back = sum - old;
                values_m[2 * c + 1] += (old - (sum - back)) + (value - back);
                values_m[2 * c] = sum;
            }
            KOKKOS_INLINE_FUNCTION void merge(const Sums& other) {
                for (unsigned c = 0; c < 6; ++c) {
                    add(c, other.values_m[2 * c]);
                    add(c, other.values_m[2 * c + 1]);
                }
            }
            double value(unsigned c) const { return values_m[2 * c] + values_m[2 * c + 1]; }
        };
        /** @brief Accumulate the two cross-level quadratic invariants in one device pass.
         * @details With a and b the current and previous vector histories,
         *      \f$2W_A/V=\sum[\varepsilon_d((a_d-b_d)/dt)^2+(Ca)_d(Cb)_d]+g(a)g(b)/\chi\f$.
         * The scalar counterpart is
         *      \f$2W_\Phi/V=\chi((\phi-\phi')/dt)^2 + \sum\varepsilon_d(D_d^+\phi)(D_d^+\phi')\f$.
         * These are the cross-level wave invariants [kreiss2002SecondOrderWave,
         * section 2, before Eq. (2.10)], specialized to Chew's masses and stiffness
         * [chew2014generalizedgauge, Eqs. (16)--(18)]. Both are conserved without
         * sources on a periodic grid; no initial gauge constraint is assumed.
         * With open faces these same sums are diagnostics, not conserved quantities;
         * they include boundary samples but no boundary-flux term.
         */
        struct Reduction {
            using value_type = Sums;
            O1VSOperator current_m, previous_m;
            material::ScalarFieldConstView epsilon_m;
            double dt_m;
            KOKKOS_INLINE_FUNCTION void init(Sums& sums) const { sums = Sums{}; }
            KOKKOS_INLINE_FUNCTION void join(Sums& a, const Sums& b) const { a.merge(b); }
            /** @brief Split terms by their staggered locations, counting each once.
             * A face belongs to the dielectric when its permittivity differs from one.
             * An edge belongs there when any of its four adjoining cells is dielectric;
             * cell centres use the cell value. Thus both slab interfaces count as dielectric.
             */
            KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, Sums& sums) const {
                const unsigned region = epsilon_m(i, j, k) == 1 ? 0 : 1;
                const auto a = current_m.potential_m, b = previous_m.potential_m;
                const auto h  = current_m.spacing_m;
                const auto ca = staggeredCurl(a, i, j, k, h, 1),
                           cb = staggeredCurl(b, i, j, k, h, 1);
                term(sums, region,
                     current_m.g(i, j, k) * previous_m.g(i, j, k)
                         * current_m.inverseMass(0, i, j, k));
                const double scalarVelocity = (a(i, j, k)[0] - b(i, j, k)[0]) / dt_m;
                term(sums, 2 + region,
                     scalarVelocity * scalarVelocity / current_m.inverseMass(0, i, j, k));
                for (unsigned d = 0; d < 3; ++d) {
                    material::CellIndex p(i, j, k);
                    ++p[d];
                    const double eps          = current_m.face_m[d](i, j, k);
                    const double velocity     = (a(i, j, k)[d + 1] - b(i, j, k)[d + 1]) / dt_m;
                    const unsigned faceRegion = eps != 1;
                    material::CellIndex edge(i, j, k);
                    ++edge[(d + 2) % 3];
                    const unsigned edgeRegion =
                        current_m.face_m[(d + 1) % 3](i, j, k) != 1
                        || current_m.face_m[(d + 2) % 3](i, j, k) != 1
                        || current_m.face_m[(d + 1) % 3](edge[0], edge[1], edge[2]) != 1;
                    term(sums, faceRegion, eps * velocity * velocity);
                    term(sums, edgeRegion, ca[d] * cb[d]);
                    term(sums, 2 + faceRegion,
                         eps * ((a(p[0], p[1], p[2])[0] - a(i, j, k)[0]) / h[d])
                             * ((b(p[0], p[1], p[2])[0] - b(i, j, k)[0]) / h[d]));
                }
            }
            /// Preserve absolute elementary terms before cancellation in the invariant.
            KOKKOS_INLINE_FUNCTION static void term(Sums& sums, unsigned region, double value) {
                sums.add(region, value);
                sums.add(region < 2 ? 4 : 5, Kokkos::abs(value));
            }
        };
        /** @brief Combine rank-local pairs in increasing communicator-rank order everywhere.
         * MPI reductions may choose their operation order [mpi2023Standard41,
         * section 6.9.1]; gathering twelve doubles per rank fixes this final order.
         */
        inline Sums combine(const Sums& local, MPI_Comm comm) {
            int ranks;
            MPI_Comm_size(comm, &ranks);
            std::vector<Sums> gathered(ranks);
            MPI_Allgather(local.values_m, 12, MPI_DOUBLE, gathered.data(), 12, MPI_DOUBLE, comm);
            Sums total;
            for (const auto& rank : gathered)
                total.merge(rank);
            return total;
        }
    }  // namespace energyDetail

    /** @brief Read both leapfrog invariants, their material split, and cancellation scales.
     * @param solver A ChdrSolver with the O1-VS formulation, both potential slots loaded and
     * their halos valid. Every rank of its layout calls this function.
     * @param owner The MaterialOwner that was passed to the solver constructor. The function
     * does not check this identity; another owner gives sums with other coefficients.
     * @return The two invariants with their vacuum/dielectric split, their cancellation scales
     * and their time labels, identical on every rank.
     * @details No fields are copied to host: only twelve compensated reduction
     * scalars cross MPI. Region membership follows the cell-centre, face and edge rule in
     * Reduction. A positive initial total is required before forming relative drift; an
     * identically zero scalar invariant is reported without normalization.
     */
    template <class Base>
    EnergyRecord potentialEnergy(const ChdrSolver<Base>& solver,
                                 const material::MaterialOwner& owner) {
        if (solver.formulationParameters().formulation_m != MaterialFormulation::O1VS)
            throw std::invalid_argument("Potential energy requires the O1-VS formulation");
        const double alpha = solver.formulationParameters().alpha_m, dt = solver.getDt();
        const energyDetail::Reduction reducer{
            {solver.A_n, owner, alpha}, {solver.A_nm1, owner, alpha}, owner.epsilon(), dt};
        energyDetail::Sums local;
        Kokkos::parallel_reduce("O1-VS potential energy", solver.A_n.getFieldRangePolicy(), reducer,
                                local);
        const auto sum      = energyDetail::combine(local, owner.layout().comm);
        const auto h        = owner.mesh().getMeshSpacing();
        const double factor = 0.5 * h[0] * h[1] * h[2], time = solver.completedSolves() * dt;
        return {{factor * sum.value(0), factor * sum.value(1), factor * sum.value(4)},
                {factor * sum.value(2), factor * sum.value(3), factor * sum.value(5)},
                time - dt / 2,
                time};
    }
}  // namespace chdr::solver
