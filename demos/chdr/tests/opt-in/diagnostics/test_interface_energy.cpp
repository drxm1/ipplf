/** @file test_interface_energy.cpp
 * @brief Conservation of the two discrete potential invariants of O1-VS, and their
 * independence of the MPI decomposition.
 * @details Random potential histories do not impose the generalized-gauge constraint.
 * This checks the two leapfrog invariants defined in energyDetail::Reduction,
 * not physical electromagnetic energy conservation for these initial data.
 * The conserved quantities are cross-level quadratic forms
 * [kreiss2002SecondOrderWave, section 2, before Eq. (2.10)], specialized to the
 * potential equations [chew2014generalizedgauge, Eqs. (16)--(18)].
 *
 * Configuration: a periodic box of 8 by 6 by 10 cells with lengths 8, 7.8 and 8, and the
 * dielectric slab \f$0 \le x \le 4\f$ (cells 0 to 3) with \f$\varepsilon_r = 2.13\f$.
 * ChdrSolver uses R2, the staggered reading, O1-VS and the safety factor 0.9. The
 * histories are random (seeds 23 and 71), the source is zero, and the run has 50 solves.
 * Each test requires a relative drift of at most \f$10^{-12}\f$ for the vector and for
 * the scalar invariant. Rank 0 then repeats the run on a communicator-self layout. The
 * retained histories \f$A^n\f$ and \f$A^{n-1}\f$ must be bitwise equal to the distributed
 * run, and each part of each energy record (total, vacuum, dielectric, absolute sum) must
 * agree within \f$10^{-14}\f$ times the initial one-rank total.
 *
 * CTest: chdr.diagnostics.energy.mpi2 runs NormalSplit and TangentialSplit on two ranks.
 * Each test includes a one-rank replay.
 */
#include <bit>
#include <gtest/gtest.h>

#include "core/material/PlanarSlab.h"
#include "core/solver/BoxSetup.h"
#include "opt-in/diagnostics/FieldEnergy.h"
#include "tests/core/solver/support/RandomFields.h"
#include "tests/core/solver/support/SerialReference.h"

namespace {
    using namespace chdr;
    using namespace chdr::solver;
    using namespace chdr::solver::test;
    using namespace chdr::diagnostics;
    using Solver = ChdrSolver<StandardSolverPeriodic>;
    struct State {
        std::vector<EnergyRecord> energy_m;
        std::vector<double> current_m, previous_m;
    };
    /// Identical global-index samples preserve bits across decompositions without a mean reduction.
    State evolve(AJFourFieldsType& source, EMFieldsType& electric, EMFieldsType& magnetic) {
        material::MaterialOwner owner(
            source.get_mesh(), source.getLayout(),
            material::TwoMaterialMap{material::PlanarSlab{0, 0, 4}, 2.13, 1.0});
        Solver solver(source, electric, magnetic, EBReconstructionConvention::R2,
                      SpatialDiscretization::Staggered, owner, {MaterialFormulation::O1VS, 1});
        solver.setSafetyFactor(0.9);
        randomize(solver.A_n, 23, false);
        randomize(solver.A_nm1, 71, false);
        solver.prepareHistories();
        source = 0;
        State result;
        result.energy_m.push_back(potentialEnergy(solver, owner));
        for (int step = 0; step < 50; ++step) {
            solver.solve();
            result.energy_m.push_back(potentialEnergy(solver, owner));
        }
        result.current_m  = solver::test::gatherOwned(solver.A_n);
        result.previous_m = solver::test::gatherOwned(solver.A_nm1);
        return result;
    }
    /** @brief Bound each potential invariant's relative drift by 1e-12 over the run.
     * @details Compare each whole-domain total with its initial value, separately
     * for the vector and scalar invariants. Both initial totals must be positive;
     * all sampled totals must be finite. Print cancellation scales with the drift.
     */
    void checkEnergy(const State& state) {
        for (auto member : {&EnergyRecord::vector_m, &EnergyRecord::scalar_m}) {
            const auto initial = state.energy_m.front().*member;
            ASSERT_GT(initial.total(), 0);
            double drift = 0;
            for (const auto& record : state.energy_m) {
                ASSERT_TRUE(std::isfinite((record.*member).total()));
                drift = std::max(
                    drift, std::abs((record.*member).total() - initial.total()) / initial.total());
            }
            std::cout << "energy=" << initial.total()
                      << " S_abs/W0=" << initial.absolute_m / initial.total()
                      << " max_relative_drift=" << drift << '\n';
            EXPECT_LE(drift, 1e-12);
        }
    }
    /// Check both retained histories and every energy record against communicator-self evolution.
    void compare(const State& distributed, const State& serial) {
        for (auto member : {&State::current_m, &State::previous_m}) {
            const auto& a = distributed.*member;
            const auto& b = serial.*member;
            ASSERT_EQ(a.size(), b.size());
            std::size_t differences = 0;
            for (std::size_t i = 0; i < a.size(); ++i)
                differences +=
                    std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i]);
            EXPECT_EQ(differences, 0);
        }
        for (auto member : {&EnergyRecord::vector_m, &EnergyRecord::scalar_m}) {
            const double scale = (serial.energy_m.front().*member).total();
            double error       = 0;
            for (std::size_t k = 0; k < serial.energy_m.size(); ++k) {
                const auto a = distributed.energy_m[k].*member, b = serial.energy_m[k].*member;
                error = std::max({error, std::abs(a.total() - b.total()),
                                  std::abs(a.vacuum_m - b.vacuum_m),
                                  std::abs(a.dielectric_m - b.dielectric_m),
                                  std::abs(a.absolute_m - b.absolute_m)});
            }
            std::cout << "rank_energy_relative=" << error / scale << '\n';
            EXPECT_LE(error / scale, 1e-14);
        }
    }
    /// x splits place the seam at x=4; y splits cross the same half-space tangentially.
    void checkSplit(unsigned axis) {
        auto fields = makeBoxFieldContainer({8, 6, 10}, {8, 7.8, 8}, axis);
        initializeFieldStorage(fields);
        const auto distributed = evolve(fields.getJ(), fields.getE(), fields.getB());
        checkEnergy(distributed);
        if (ippl::Comm->rank() == 0) {
            solver::test::SerialFields fieldsSerial(fields.getJ());
            const auto serial =
                evolve(fieldsSerial.source_m, fieldsSerial.electric_m, fieldsSerial.magnetic_m);
            checkEnergy(serial);
            compare(distributed, serial);
        }
    }
    /// Split along x: the rank seam lies on the slab face at x = 4. Drift at most 1e-12.
    TEST(InterfaceEnergy, NormalSplit) {
        checkSplit(0);
    }
    /// Split along y: the rank seam crosses the slab and the vacuum. Drift at most 1e-12.
    TEST(InterfaceEnergy, TangentialSplit) {
        checkSplit(1);
    }
}  // namespace
