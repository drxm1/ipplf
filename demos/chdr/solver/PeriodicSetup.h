/** @file PeriodicSetup.h
 * @brief Runtime periodic boxes with a single selectable decomposition axis.
 */
#pragma once

#include <array>

#include "ChdrSolver.h"
#include "SpatialDimension.h"

static_assert(chdr::solver::SpatialDim == 3,
              "PeriodicSetup.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver {
    /** @brief Construct mesh/layout ownership; allocate with initializeFieldStorage.
     * @pre All ranks must call it with identical arguments; lengths are finite and positive.
     * @param cells     Global cell counts in x, y and z; each must be at least two.
     * @param lengths   Physical box lengths in x, y and z; finite and positive.
     *                  Size of the whole container, in each axis.
     * @param splitAxis Axis eligible for MPI decomposition: 0=x, 1=y, 2=z.
     *                  The other two axes remain undistributed. Defaults to x.
     * @return Prvalue container; retain its address once its fields are allocated.
     * @details Construction defines mesh and MPI ownership only. Allocation is a separate
     * operation because the field objects retain pointers to this owner's mesh and layout.
     * Do not return or move an already allocated container through a setup helper.
     */
    ChDRFieldContainer makePeriodicFieldContainer(std::array<int, SpatialDim> cells,
                                                  std::array<double, SpatialDim> lengths,
                                                  unsigned splitAxis = 0);

    /// Reusing from the vacuum wave example
    using chdr::vacuum::checkLocalDomain;

    /// Reusing from the vacuum wave example
    using chdr::vacuum::initializeFieldStorage;

    /// Reusing from the vacuum wave example
    using chdr::vacuum::preparePotentialHalo;

    /** @brief Write owned analytic components at their declared positions and times.
     * @tparam AnalyticSampler Device-copyable functor returning one component value from
     * sample(x,y,z,time,component), callable inside a Kokkos kernel.
     * @param field Allocated four-component potential or source destination: component
     * zero stores scalar potential or charge density; components 1 through SpatialDim
     * store vector potential or current density. Mesh/layout owners remain with the caller.
     * @param sample Analytic sampler with the same component ordering as field.
     * @param times Absolute times for component zero followed by each Cartesian vector
     * component; scalar and vector times may differ.
     * @param contract Supplies the spatial offsets for the chosen reading; times are explicit.
     * @details All coordinates come from global integer indices; no rank-local phase origin.
     * The component offsets come from LevelContract and are applied before calling the
     * analytic functor. The kernel captures only views, scalar/vector metadata and the
     * device-callable functor. Its coordinate-to-sample operation is kept together so
     * the arithmetic grouping is visible. Owned cells are written and fenced; no halo
     * communication is performed by this helper.
     */
    template <class AnalyticSampler>
    void writeSamples(AJFourFieldsType& field, AnalyticSampler sample, std::array<double, 4> times,
                      LevelContract contract = {}) {
        const auto view = field.getView();
        const auto h = field.get_mesh().getMeshSpacing(), origin = field.get_mesh().getOrigin();
        const auto local = field.getLayout().getLocalNDIndex();
        const ippl::Vector<int, SpatialDim> first(local[0].first(), local[1].first(),
                                                  local[2].first());
        const ippl::Vector<double, 4> labels(times[0], times[1], times[2], times[3]);
        const int ghost = field.getNghost();
        Kokkos::parallel_for(
            "Labelled analytic samples", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                const double x = origin[0] + (i + first[0] - ghost + 0.5) * h[0];
                const double y = origin[1] + (j + first[1] - ghost + 0.5) * h[1];
                const double z = origin[2] + (k + first[2] - ghost + 0.5) * h[2];
                for (unsigned c = 0; c < 4; ++c) {
                    const auto offset = contract.positionOffset(
                        c == 0 ? Quantity::Phi : Quantity::VectorPotential, c == 0 ? 0 : c - 1);
                    view(i, j, k)[c] = sample(x + offset[0] * h[0], y + offset[1] * h[1],
                                              z + offset[2] * h[2], labels[c], c);
                }
            });
        Kokkos::fence();
    }

    /** @brief Load both potential slots and synchronize their halos, without reconstructing E/B.
     * @tparam PotentialSolverBase StandardSolverPeriodic or NonStandardSolverPeriodic;
     * ChdrSolver enforces this choice and supplies the required history and label interface.
     * @tparam AnalyticSampler Device-copyable scalar-valued analytic potential sampler.
     * @param solver Initialized shell whose final time spacing is already selected.
     * @param sample Device-callable sample(x,y,z,time,component), with scalar potential
     * in component zero and Cartesian vector potential in components 1 through SpatialDim.
     * @details [fallahi2020mithra20fullwavesimulation, section 3.1.4]; component
     * time labels are in LevelContract. The load freezes the solver's time step.
     * @pre The solver has not evolved since initialization and its final time step is set.
     * @details Scalar and vector labels are obtained separately. Each previous slot is
     * sampled one step before its corresponding current slot; the source labels are not
     * inferred from whichever potential component happens to be loaded first.
     * A_n receives the current slot and A_nm1 its predecessor; E/B remain unevaluated.
     * Analytic sampling alone does not enforce a discrete initial gauge condition.
     * Exact sourced staggered Gauss propagation additionally needs consistent
     * initial gauge and continuity histories [christlieb2024GaugeConserving,
     * Theorems 3.3 and 3.5]; continuum-sampled histories may leave a free gauge wave.
     */
    template <class PotentialSolverBase, class AnalyticSampler>
    void loadHistory(ChdrSolver<PotentialSolverBase>& solver, AnalyticSampler sample) {
        solver.markHistoryLoaded();
        const double a   = solver.labelTime(Quantity::VectorPotential, 0);
        const double phi = solver.labelTime(Quantity::Phi, 0), dt = solver.getDt();
        writeSamples(solver.A_n, sample, {phi, a, a, a}, solver.levelContract());
        writeSamples(solver.A_nm1, sample, {phi - dt, a - dt, a - dt, a - dt},
                     solver.levelContract());
        preparePotentialHalo(solver.A_n);
        preparePotentialHalo(solver.A_nm1);
    }

    /** @brief Write the source consumed by the next solve; stepping reads owned source cells only.
     * @tparam PotentialSolverBase StandardSolverPeriodic or NonStandardSolverPeriodic;
     * ChdrSolver enforces this choice and supplies the required source and label interface.
     * @tparam AnalyticSampler Device-copyable scalar-valued charge/current sampler.
     * @param solver Shell whose sourceField is overwritten on owned cells.
     * @param sample Device-callable sample(x,y,z,time,component), returning charge density
     * for component zero and Cartesian current density for components 1 through SpatialDim.
     * @param sourceStep Source index for the next solve, normally solver.completedSolves().
     * @details Charge and current use their separate LevelContract times. Their spatial
     * offsets match the scalar and vector potential locations used by writeSamples.
     * The functor and integer source index remain explicit so a test can deliberately
     * plant a mislabelled source without changing the production level table.
     */
    template <class PotentialSolverBase, class AnalyticSampler>
    void writeSource(ChdrSolver<PotentialSolverBase>& solver, AnalyticSampler sample,
                     long long sourceStep) {
        const double current = solver.labelTime(Quantity::Current, sourceStep);
        const double charge  = solver.labelTime(Quantity::Charge, sourceStep);
        writeSamples(solver.sourceField(), sample, {charge, current, current, current},
                     solver.levelContract());
    }
}  // namespace chdr::solver
