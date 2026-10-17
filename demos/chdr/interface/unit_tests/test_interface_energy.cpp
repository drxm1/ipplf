/** @file test_interface_energy.cpp
 * @brief Conservation of discrete potential invariants and rank decomposition checks.
 * @details Random potential histories do not impose the generalized-gauge constraint.
 * This checks the two leapfrog invariants defined in energyDetail::Reduction,
 * not physical electromagnetic energy conservation for these initial data.
 * The conserved quantities are cross-level quadratic forms
 * [kreiss2002SecondOrderWave, section 2, before Eq. (2.10)], specialized to the
 * potential equations [chew2014generalizedgauge, Eqs. (16)--(18)].
 */
#include <gtest/gtest.h>

#include "FieldEnergy.h"
#include "PlanarSlab.h"
#include "RandomFields.h"

namespace {
    using namespace chdr;
    using namespace chdr::interfaceTest;
    using Solver = ChdrSolver<StandardSolverPeriodic>;
    struct State {
        std::vector<EnergyRecord> energy_m;
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
    TEST(InterfaceEnergy, Serial) {
        auto fields = makeBoxFieldContainer({8, 6, 10}, {8, 7.8, 8});
        initializeFieldStorage(fields);
        checkEnergy(evolve(fields.getJ(), fields.getE(), fields.getB()));
    }
}  // namespace
