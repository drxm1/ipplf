/** @file test_simulation.cpp
 * @brief Shared-manager configuration and initial/completed-state observation checks.
 */
#include <gtest/gtest.h>
#include <optional>
#include <vector>

#include "../../simulation/ChDRSimulationManager.h"

namespace {
    using namespace chdr;
    using Simulation = simulation::ChDRSimulationManager<>;

    /// Supply a uniform material independently of box and solver choices.
    struct UniformMaterial {
        double epsilon_m = 1.0;
        KOKKOS_INLINE_FUNCTION double operator()(material::CellIndex, material::Position,
                                                 bool& marked) const {
            marked = false;
            return epsilon_m;
        }
    };

    /// Exercise lifecycle transitions with zero initial potentials in every component.
    struct ZeroPotential {
        KOKKOS_INLINE_FUNCTION double operator()(double, double, double, double, unsigned) const {
            return 0.0;
        }
    };

    /// Check initial observation, exact solve counts, observer restoration and empty observers.
    void checkObservations(Simulation& run) {
        std::vector<long long> states;
        const auto record = [&](long long k, auto&, const auto&) {
            states.push_back(k);
        };
        run.setObserver([&](long long k, auto&, const auto&) {
            states.push_back(100 + k);
        });
        run.advance(0, record);
        run.advance(2, record);
        run.run(1);
        run.setObserver({});
        run.run(1);
        EXPECT_EQ(states, (std::vector<long long>{0, 0, 1, 2, 103}));
        EXPECT_EQ(run.solver().completedSolves(), 4);
        EXPECT_THROW(run.advance(-1, record), std::invalid_argument);
    }

    TEST(Simulation, OpenBoxObserver) {
        const simulation::RunBox box{{8, 2, 8}, {1, 1, 1}, 0, {false, true, false}, 2};
        const simulation::SolverSettings settings{
            .safetyFactor_m   = 0.7,
            .reconstruction_m = solver::EBReconstructionConvention::R1,
            .spatial_m        = solver::SpatialDiscretization::Collocated,
            .formulation_m =
                solver::FormulationParameters{solver::MaterialFormulation::Homogeneous, 1}};
        Simulation run(box, UniformMaterial{}, settings);
        run.load(ZeroPotential{});
        const auto& layout = run.solver().A_n.getLayout();
        EXPECT_FALSE(layout.isPeriodic(0));
        EXPECT_TRUE(layout.isPeriodic(1));
        EXPECT_FALSE(layout.isPeriodic(2));
        EXPECT_EQ(layout.getLocalNDIndex()[0].length(), 8);
        EXPECT_EQ(layout.getLocalNDIndex()[1].length(), 2);
        EXPECT_NEAR(run.solver().safetyFactor(), 0.7, 1e-15);
        checkObservations(run);
    }

    TEST(Simulation, VacuumSelection) {
        const simulation::RunBox box{{8, 2, 8}, {1, 1, 1}};
        const simulation::SolverSettings settings{.formulation_m = std::nullopt};
        Simulation run(box, UniformMaterial{}, settings);
        EXPECT_FALSE(run.solver().hasMaterial());
        run.load(ZeroPotential{});
        run.run(1);
        EXPECT_EQ(run.solver().completedSolves(), 1);
        EXPECT_THROW((Simulation{box, UniformMaterial{2.13}, settings}), std::invalid_argument);
    }
}  // namespace
