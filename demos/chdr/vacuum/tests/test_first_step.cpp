/** @file test_first_step.cpp
 * @brief Initial-state and one-step regression for the fixed coarse vacuum demo.
 * @details Tests 8x8x64 periodic cells, source-free Standard evolution and the
 * production setup helpers on one or two ranks. Our independent host oracle
 * substitutes the declared wave into Fallahi, MITHRA 2.0, arXiv:2009.13645v1,
 * "Wave Equation" (3.6)-(3.9) and "FDTD for Wave Equation" (3.11)-(3.15),
 * printed pp.12-13. It never calls referencePotential or solver reconstruction.
 * Code contract: StandardFDTDSolver::step (StandardFDTDSolver.hpp:65-75) and
 * FDTDSolverBase::solve/timeShift/evaluate_EB (FDTDSolverBase.hpp:23-26,60-63,109-127).
 */
#include "Ippl.h"

#include <Kokkos_MathematicalConstants.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <memory>

#include "VacuumSetup.h"

namespace {
    using chdr::vacuum::checkLocalDomain;
    using chdr::vacuum::Dim;
    using chdr::vacuum::Fields;
    using chdr::vacuum::initializeFieldStorage;
    using chdr::vacuum::initializePotentials;
    using chdr::vacuum::makePeriodicFieldContainer;
    using chdr::vacuum::VacuumSolver;
    constexpr double WaveNumber = 2.0 * Kokkos::numbers::pi_v<double>;
    using Potential             = ippl::Vector<double, 4>;
    using FieldValue            = ippl::Vector<double, Dim>;
    using PotentialHost         = SourceField_t<double, Dim>::host_mirror_type;
    using FieldHost             = VField_t<double, Dim>::host_mirror_type;

    /// @brief Independent host copies; keep the initial Ax values alive across solve().
    struct Snapshot {
        PotentialHost current_m, previous_m, next_m, source_m;
        FieldHost electric_m, magnetic_m;
    };

    /// @brief One owned cell, with initial Ax history available for the exact-copy check.
    struct Cell {
        Potential current_m, previous_m, next_m, source_m, original_m;
        FieldValue electric_m, magnetic_m;
    };

    /// @brief Discrete expected values and continuum E/B at their distinct sample times.
    struct Reference {
        double current_m, previous_m, next_m, electric_m, magnetic_m;
        double continuumE_m, continuumB_m;
    };

    /// @brief Local maxima and unnormalized sums; globalErrors performs the MPI reductions.
    struct Errors {
        double potential_m = 0.0;
        double field_m     = 0.0;
        double zero_m      = 0.0;
        double source_m    = 0.0;
        double shift_m     = 0.0;
        double duplicate_m = 0.0;
        double errorE_m = 0.0, normE_m = 0.0, errorB_m = 0.0, normB_m = 0.0;
        unsigned long long count_m = 0;
    };

    /** @brief Independent normalized continuum potential, with k=2*pi and c=1.
     * @details Our wave is \f$A_x(z,t)=\sin(k(z-t))/k\f$, with all other components zero.
     * Starting relations: Fallahi, MITHRA 2.0, Eqs. (3.6)-(3.9).
     */
    double initialAx(double z, double time) {
        return std::sin(WaveNumber * (z - time)) / WaveNumber;
    }

    /** @brief Apply the axial stencil to analytic histories without calling the solver.
     * @details Our substitution gives \f$A^1=2A^0-A^{-1}+dt^2\delta_{zz,h}A^0\f$;
     * transverse differences vanish for this wave. Starting stencil: Fallahi,
     * MITHRA 2.0, Eqs. (3.11)-(3.15). dt and hz come from the actual solver and mesh.
     */
    double updatedAx(double z, double dt, double hz) {
        const double centre = initialAx(z, 0.0);
        const double laplacian =
            (initialAx(z + hz, 0.0) - 2.0 * centre + initialAx(z - hz, 0.0)) / (hz * hz);
        return 2.0 * centre - initialAx(z, -dt) + dt * dt * laplacian;
    }

    /// @brief One-step oracle; step must be -1, 0 or 1, never an arbitrary evolution time.
    double potentialAt(double z, double dt, double hz, int step) {
        return step == 1 ? updatedAx(z, dt, hz) : initialAx(z, step * dt);
    }

    /** @brief Independently reconstruct expected fields for step=0 or step=1.
     * @details For this phi=0 case, \f$E_x=-(A_x^n-A_x^{n-1})/dt\f$ and
     * \f$B_y=(A_x^n(z+h_z)-A_x^n(z-h_z))/(2h_z)\f$.
     * These are our discrete substitutions in Fallahi, MITHRA 2.0, (3.8)-(3.9).
     * IPPL's FDTDSolverBase::evaluate_EB (FDTDSolverBase.hpp:109-127) sets raw E/B times
     * to (n-1/2)dt and n*dt;
     * the continuum entries use those respective times. At step=0, next is zero scratch.
     */
    Reference referenceAt(double z, double dt, double hz, int step) {
        const double current  = potentialAt(z, dt, hz, step);
        const double previous = potentialAt(z, dt, hz, step - 1);
        const double electric = -(current - previous) / dt;
        const double magnetic =
            (potentialAt(z + hz, dt, hz, step) - potentialAt(z - hz, dt, hz, step)) / (2.0 * hz);
        return {current,
                previous,
                step == 0 ? 0.0 : current,
                electric,
                magnetic,
                std::cos(WaveNumber * (z - (step - 0.5) * dt)),
                std::cos(WaveNumber * (z - step * dt))};
    }

    /// @brief Deep-copy a field; a host-accessible alias would not preserve the old history.
    template <class Field>
    auto copyToHost(Field& field) {
        auto host = field.getHostMirror();
        Kokkos::deep_copy(host, field.getView());
        return host;
    }

    /// @brief Copy all solver histories and E/B/J into independent host storage.
    Snapshot snapshot(VacuumSolver& solver, Fields& fields) {
        return {copyToHost(solver.A_n),    copyToHost(solver.A_nm1),  copyToHost(solver.A_np1),
                copyToHost(fields.getJ()), copyToHost(fields.getE()), copyToHost(fields.getB())};
    }

    /// @brief Read one local view index, including halos in the indexing offset.
    Cell cellAt(const Snapshot& state, const PotentialHost& initial, int i, int j, int k) {
        return {state.current_m(i, j, k), state.previous_m(i, j, k), state.next_m(i, j, k),
                state.source_m(i, j, k),  initial(i, j, k),          state.electric_m(i, j, k),
                state.magnetic_m(i, j, k)};
    }

    /** @brief Maximum scaled component error, or infinity for any nonfinite value.
     * @details active identifies the sole nonzero reference component; -1 means all zero.
     * Potential errors use scale=k (amplitude 1/k); E/B use their unit amplitude.
     */
    template <unsigned Components>
    double componentError(const ippl::Vector<double, Components>& value, int active,
                          double expected, double scale) {
        double error = 0.0;
        for (unsigned c = 0; c < Components; ++c) {
            if (!std::isfinite(value[c]))
                return std::numeric_limits<double>::infinity();
            const double reference = static_cast<int>(c) == active ? expected : 0.0;
            error                  = std::max(error, scale * std::abs(value[c] - reference));
        }
        return error;
    }

    /// @brief Ignore a valid active component's magnitude while still checking finiteness.
    template <unsigned Components>
    double inactiveError(const ippl::Vector<double, Components>& value, int active, double scale) {
        return componentError(value, active, value[active], scale);
    }

    /// @brief Maximum scaled component difference; zero means an exact history copy.
    double pairError(const Potential& first, const Potential& second) {
        const Potential difference = first - second;
        return componentError(difference, -1, 0.0, WaveNumber);
    }

    /// @brief Maximum inactive-component error across all histories and reconstructed fields.
    double zeroError(const Cell& cell) {
        return std::max({inactiveError(cell.current_m, 1, WaveNumber),
                         inactiveError(cell.previous_m, 1, WaveNumber),
                         inactiveError(cell.next_m, 1, WaveNumber),
                         inactiveError(cell.electric_m, 0, 1.0),
                         inactiveError(cell.magnetic_m, 1, 1.0)});
    }

    /** @brief Accumulate one owned cell's continuum-error and reference-norm squares.
     * @details Read Ex/By from cell and their own-time references from ref; update error's
     * rank-local sums with \f$(F-F_{\rm ref})^2\f$ and \f$F_{\rm ref}^2\f$.
     * Our relative-error metric is \f$\sqrt{\sum(F-F_{\rm ref})^2/\sum F_{\rm ref}^2}\f$;
     * the caller sums globally and checks a positive denominator before normalization.
     * The continuum reference follows Fallahi, MITHRA 2.0, Eqs. (3.8)-(3.9).
     */
    void checkContinuum(const Cell& cell, const Reference& ref, Errors& error) {
        const double electric = cell.electric_m[0] - ref.continuumE_m;
        const double magnetic = cell.magnetic_m[1] - ref.continuumB_m;
        error.errorE_m += electric * electric;
        error.normE_m += ref.continuumE_m * ref.continuumE_m;
        error.errorB_m += magnetic * magnetic;
        error.normB_m += ref.continuumB_m * ref.continuumB_m;
    }

    /** @brief Add this cell's observations to error; step must be 0 or 1.
     * @details Check fields, sources and inactive components; after one solve also
     * compare the shifted history to its deep snapshot and the duplicated newest arrays.
     * History-copy contract: FDTDSolverBase::timeShift (FDTDSolverBase.hpp:60-63).
     * This helper makes no MPI calls.
     */
    void checkCell(const Cell& cell, const Reference& ref, int step, Errors& error) {
        error.potential_m = std::max(
            {error.potential_m, componentError(cell.current_m, 1, ref.current_m, WaveNumber),
             componentError(cell.previous_m, 1, ref.previous_m, WaveNumber),
             componentError(cell.next_m, 1, ref.next_m, WaveNumber)});
        error.field_m =
            std::max({error.field_m, componentError(cell.electric_m, 0, ref.electric_m, 1.0),
                      componentError(cell.magnetic_m, 1, ref.magnetic_m, 1.0)});
        error.zero_m   = std::max(error.zero_m, zeroError(cell));
        error.source_m = std::max(error.source_m, componentError(cell.source_m, -1, 0.0, 1.0));
        checkContinuum(cell, ref, error);
        if (step == 1) {
            error.shift_m = std::max(error.shift_m, pairError(cell.previous_m, cell.original_m));
            error.duplicate_m = std::max(error.duplicate_m, pairError(cell.current_m, cell.next_m));
        } else {
            error.duplicate_m =
                std::max(error.duplicate_m, componentError(cell.next_m, -1, 0.0, 1.0));
        }
        ++error.count_m;
    }

    /** @brief Compare owned cells only against the independent step=0/1 reference.
     * @details Global cell indices determine positions; subtract the owned first index
     * and add the halo width to address each local view. Halos do not enter the norms.
     */
    Errors compareOwned(VacuumSolver& solver, Fields& fields, const Snapshot& state,
                        const PotentialHost& initial, int step) {
        const auto& local = fields.getFL().getLocalNDIndex();
        const auto& mesh  = fields.getMesh();
        const int ghost   = solver.A_n.getNghost();
        Errors error;
        for (int i = local[0].first(); i <= local[0].last(); ++i)
            for (int j = local[1].first(); j <= local[1].last(); ++j)
                for (int k = local[2].first(); k <= local[2].last(); ++k) {
                    const double z = mesh.getOrigin()[2] + (k + 0.5) * mesh.getMeshSpacing()[2];
                    const auto ref = referenceAt(z, solver.getDt(), mesh.getMeshSpacing()[2], step);
                    const auto cell =
                        cellAt(state, initial, i - local[0].first() + ghost,
                               j - local[1].first() + ghost, k - local[2].first() + ghost);
                    checkCell(cell, ref, step, error);
                }
        return error;
    }

    /** @brief Combine observations before any rank evaluates GoogleTest assertions.
     * @pre All layout ranks participate in the same order.
     * @return Identical global maxima, squared-error/reference sums and cell count on all ranks.
     */
    Errors globalErrors(const Errors& local, Fields& fields) {
        double maxima[6]         = {local.potential_m, local.field_m, local.zero_m,
                                    local.source_m,    local.shift_m, local.duplicate_m};
        double sums[4]           = {local.errorE_m, local.normE_m, local.errorB_m, local.normB_m};
        unsigned long long count = local.count_m;
        MPI_Allreduce(MPI_IN_PLACE, maxima, 6, MPI_DOUBLE, MPI_MAX, fields.getFL().comm);
        MPI_Allreduce(MPI_IN_PLACE, sums, 4, MPI_DOUBLE, MPI_SUM, fields.getFL().comm);
        MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                      fields.getFL().comm);
        return {maxima[0], maxima[1], maxima[2], maxima[3], maxima[4], maxima[5],
                sums[0],   sums[1],   sums[2],   sums[3],   count};
    }

    /** @brief Apply the fixed coarse-grid benchmark thresholds to global observations.
     * @details Expected count is 8*8*64; source and history copies must be exact.
     * Potential errors are normalized by 1/k; E/B have unit amplitude. Relative L2
     * uses the global squared-error/reference sums described in checkContinuum.
     * The cutoffs are our chosen regression requirements, explained in this demo's
     * README.md, "Tests and source reading"; they are not thresholds stated by Fallahi.
     * This tests initialization and one step only.
     */
    void expectErrors(const Errors& error) {
        EXPECT_LE(error.potential_m, 1e-10);
        EXPECT_LE(error.field_m, 1e-10);
        EXPECT_LE(error.zero_m, 1e-12);
        EXPECT_EQ(error.source_m, 0.0);
        EXPECT_EQ(error.shift_m, 0.0);
        EXPECT_EQ(error.duplicate_m, 0.0);
        EXPECT_EQ(error.count_m, 4096ULL);
        ASSERT_GT(error.normE_m, 0.0);
        ASSERT_GT(error.normB_m, 0.0);
        EXPECT_LE(std::sqrt(error.errorE_m / error.normE_m), 0.01);
        EXPECT_LE(std::sqrt(error.errorB_m / error.normB_m), 0.01);
    }

    /** @brief Fresh production setup for each test; solver is destroyed before fields.
     * @details NaN poisoning makes missed initialization visible. E/B reference times
     * for phi=0 are -dt/2 and 0 initially, dt/2 and dt after solve, as implemented in
     * FDTDSolverBase::evaluate_EB (FDTDSolverBase.hpp:109-127), using Fallahi's (3.8)-(3.9).
     * Current-history face halos feed reconstruction
     * and update; all owned values are checked. This is not a long-time convergence test.
     */
    class VacuumFirstStepTest : public ::testing::Test {
    protected:
        void SetUp() override {
            checkLocalDomain(fields_m.getFL());
            initializeFieldStorage(fields_m);
            solver_m =
                std::make_unique<VacuumSolver>(fields_m.getJ(), fields_m.getE(), fields_m.getB());
            solver_m->A_n   = std::numeric_limits<double>::quiet_NaN();
            solver_m->A_nm1 = std::numeric_limits<double>::quiet_NaN();
            initializePotentials(*solver_m);
            fields_m.getE() = std::numeric_limits<double>::quiet_NaN();
            fields_m.getB() = std::numeric_limits<double>::quiet_NaN();
            solver_m->evaluate_EB();
        }

        Fields fields_m = makePeriodicFieldContainer();
        std::unique_ptr<VacuumSolver> solver_m;
    };

    /// @brief Verify both initial histories, zero scratch/source and initial raw E/B.
    TEST_F(VacuumFirstStepTest, InitializesHistoriesAndFields) {
        const auto initial = copyToHost(solver_m->A_n);
        const auto state   = snapshot(*solver_m, fields_m);
        const auto local   = compareOwned(*solver_m, fields_m, state, initial, 0);
        expectErrors(globalErrors(local, fields_m));
    }

    /// @brief Verify one real solve, including exact history shift and E/B reconstruction.
    TEST_F(VacuumFirstStepTest, AdvancesOneStep) {
        const auto initial = copyToHost(solver_m->A_n);
        solver_m->solve();
        const auto state = snapshot(*solver_m, fields_m);
        const auto local = compareOwned(*solver_m, fields_m, state, initial, 1);
        expectErrors(globalErrors(local, fields_m));
    }
}  // namespace
