/** @file test_simulation.cpp
 * @brief Call contract of ChDRSimulationManager: observer calls, solve counts, solver
 * selection, time step, observation schedules and refusals.
 * @details Most tests use an 8 by 2 by 8 box of unit spacing with zero potentials. For
 * advance(0), advance(2) and two run(1) calls (observer set for the first, cleared for the
 * second), the manager must record the states 0, 0, 1, 2, 103 and 4 completed solves. It must
 * keep the open-periodic-open layout and the safety factor 0.7 (within \f$10^{-15}\f$), select
 * the vacuum update without a material owner, use an explicit step 0.5, and keep two
 * observation schedules independent (frames 0, 2, 4, 5 and samples 1, 2, 3, 4 over 5 steps).
 * It must throw for advance(-1), for \f$\varepsilon = 2.13\f$ on the vacuum path, for a step
 * equal to the stability bound or equal to \f$-1\f$, and for a 4 by 2 by 2 box whose seam at
 * cell 1 or 3 leaves one rank one cell. No field value is compared.
 *
 * Reference: self-consistency (call contract, refusal). CTest entries:
 * chdr.simulation.unit.serial and chdr.simulation.unit.mpi2. RejectsThinSeamCollectively
 * requires two ranks and skips on one rank.
 */
#include <gtest/gtest.h>
#include <optional>
#include <vector>

#include "core/simulation/ChDRSimulationManager.h"
#include "core/simulation/ObservationSchedule.h"

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

    /// Box open in x and z, periodic in y, split along z: layout, safety factor 0.7 and the
    /// observer call sequence of checkObservations().
    TEST(Simulation, OpenBoxObserver) {
        const simulation::ChdrCartesianGridSettings box{
            {8, 2, 8}, {1, 1, 1}, 0, {false, true, false}, 2};
        const simulation::ChdrSolverSettings settings{
            .safetyFactor_m   = 0.7,
            .reconstruction_m = EBReconstructionConvention::R1,
            .spatial_m        = SpatialDiscretization::Collocated,
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

    /// Without a formulation the solver has no material owner and completes 1 solve for run(1);
    /// a permittivity of 2.13 on this vacuum path throws std::invalid_argument.
    TEST(Simulation, VacuumSelection) {
        const simulation::ChdrCartesianGridSettings box{{8, 2, 8}, {1, 1, 1}};
        const simulation::ChdrSolverSettings settings{.formulation_m = std::nullopt};
        Simulation run(box, UniformMaterial{}, settings);
        EXPECT_FALSE(run.solver().hasMaterial());
        run.load(ZeroPotential{});
        run.run(1);
        EXPECT_EQ(run.solver().completedSolves(), 1);
        EXPECT_THROW((Simulation{box, UniformMaterial{2.13}, settings}), std::invalid_argument);
    }

    /** @brief Reject collectively when either rank owns a one-cell partition.
     * @details The automatic two-by-two split is valid. Moving its seam makes only one
     * rank invalid, so validation must follow repartitioning and communicate the failure.
     * The expected runtime_error identifies the pre-allocation field-container guard;
     * the later MaterialOwner validation uses invalid_argument.
     */
    TEST(Simulation, RejectsThinSeamCollectively) {
        if (ippl::Comm->size() != 2)
            GTEST_SKIP() << "The asymmetric seam check requires two ranks";
        for (int seam : {1, 3}) {
            const simulation::ChdrCartesianGridSettings box{{4, 2, 2}, {1, 1, 1}, seam};
            EXPECT_THROW((Simulation{box, UniformMaterial{}}), std::runtime_error);
        }
    }

    /// The solver uses the explicit step 0.5, and the manager's E field is the solver's E field;
    /// a step equal to the stability bound or equal to -1 throws std::invalid_argument.
    TEST(Simulation, ExplicitStep) {
        const simulation::ChdrCartesianGridSettings box{{8, 2, 8}, {1, 1, 1}};
        simulation::ChdrSolverSettings settings{.formulation_m = std::nullopt, .timeStep_m = 0.5};
        Simulation run(box, UniformMaterial{}, settings);
        EXPECT_DOUBLE_EQ(run.solver().getDt(), 0.5);
        EXPECT_EQ(&run.fields().getE(), &run.solver().electricField());
        settings.timeStep_m = run.solver().maxStableTimeStep();
        EXPECT_THROW((Simulation{box, UniformMaterial{}, settings}), std::invalid_argument);
        settings.timeStep_m = -1;
        EXPECT_THROW((Simulation{box, UniformMaterial{}, settings}), std::invalid_argument);
    }

    /// Over 5 steps, a field schedule (every 2, final included) selects 0, 2, 4, 5 and a
    /// diagnostics schedule (first 1, last 4) selects 1, 2, 3, 4.
    TEST(Simulation, IndependentObservationSchedules) {
        Simulation run({{8, 2, 8}, {1, 1, 1}}, UniformMaterial{});
        run.load(ZeroPotential{});
        const simulation::ObservationSchedule fields{.every_m = 2, .includeFinal_m = true};
        const simulation::ObservationSchedule diagnostics{.first_m = 1, .last_m = 4};
        std::vector<long long> frames, samples;
        run.advance(5, [&](long long k, auto&, const auto&) {
            if (fields.selects(k, 5))
                frames.push_back(k);
            if (diagnostics.selects(k, 5))
                samples.push_back(k);
        });
        EXPECT_EQ(frames, (std::vector<long long>{0, 2, 4, 5}));
        EXPECT_EQ(samples, (std::vector<long long>{1, 2, 3, 4}));
    }
}  // namespace
