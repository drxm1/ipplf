/** @file ChDRSimulationManager.h
 * @brief Shared field-simulation setup and IPPL manager lifecycle for experiment drivers.
 */
#pragma once

#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>

#include "BoxSetup.h"
#include "Manager/BaseManager.h"
#include "SimulationConfig.h"

namespace chdr::simulation {
    /** @brief Construct and advance a configured field problem and invoke experiment observers.
     * @tparam FieldSolver ChDR solver type, including any experiment-specific boundary hook.
     * @details The experiment supplies box and solver settings, a material map, initial
     * potential samples and observations. This manager owns their mesh, fields, material
     * and solver. Construction prepares storage; load() supplies initial histories.
     * BaseManager::run advances without an initial observation; advance(steps, observer)
     * additionally observes the starting state and temporarily installs that observer.
     * All ranks supply identical settings and enter collective observations together.
     * The object cannot move because its fields and solver refer to its mesh and layout.
     */
    template <class FieldSolver = solver::ChdrSolver<solver::StandardSolverPeriodic>>
    class ChDRSimulationManager : public ippl::BaseManager {
    public:
        using Solver = FieldSolver;
        /// Host callback; captured resources must remain alive whenever it is invoked.
        using Observer = std::function<void(long long, Solver&, const material::MaterialOwner&)>;

        /** @brief Allocate the problem and select its timestep before loading initial histories.
         * @details Material and field storage share one mesh and layout. The solver constructor
         * binds these fields and selects its update; initialize() clears the potential histories
         * and solve count. The explicit safety factor then sets the step used by load().
         */
        template <class Map>
        ChDRSimulationManager(RunBox box, Map map, SolverSettings settings = {})
            : fields_m(solver::makeBoxFieldContainer(box.cells_m, box.lengths(), box.splitAxis_m,
                                                     box.periodic_m)) {
            applyPartitionSeam(box);
            solver::initializeFieldStorage(fields_m);
            owner_m  = std::make_unique<material::MaterialOwner>(fields_m.getMesh(),
                                                                 fields_m.getFL(), map);
            solver_m = makeSolver(settings);
            solver_m->initialize();
            solver_m->setSafetyFactor(settings.safetyFactor_m);
        }

        ChDRSimulationManager(const ChDRSimulationManager&)            = delete;
        ChDRSimulationManager& operator=(const ChDRSimulationManager&) = delete;
        ChDRSimulationManager(ChDRSimulationManager&&)                 = delete;
        ChDRSimulationManager& operator=(ChDRSimulationManager&&)      = delete;

        /** @brief Sample the initial potential histories, prepare halos and reconstruct fields.
         * @details The device-callable sampler follows solver::loadHistory. Observations are
         * separate so the caller can create writers and probes after the fields are initialized.
         */
        template <class Sampler>
        void load(Sampler sample) {
            solver::loadHistory(*solver_m, sample);
            solver_m->evaluateInitialFields();
        }

        /// Advance the configured field solver by one timestep.
        void advance() override { solver_m->solve(); }

        /// Observe the fields reconstructed by the completed solve.
        void post_step() override { observeCurrentState(); }

        /// Install a persistent observer for inherited run(); an empty function disables
        /// observation.
        void setObserver(Observer observer) { observer_m = std::move(observer); }

        /// Observe the current state without advancing or changing the completed-solve count.
        void observeCurrentState() {
            if (observer_m)
                observer_m(solver_m->completedSolves(), *solver_m, *owner_m);
        }

        /** @brief Observe the starting state, then advance and observe exactly steps solves.
         * @details Restore the previous observer on return or exception, so caller-local captures
         * need only remain alive during this synchronous call. Restoration affects the callback,
         * not numerical state. Step labels continue across consecutive calls.
         */
        template <class Callback>
        void advance(int steps, Callback sample) {
            if (steps < 0)
                throw std::invalid_argument("The number of simulation steps must be nonnegative");
            auto previous = std::exchange(observer_m, Observer(std::move(sample)));
            try {
                observeCurrentState();
                ippl::BaseManager::run(steps);
            } catch (...) {
                observer_m = std::move(previous);
                throw;
            }
            observer_m = std::move(previous);
        }

        /// Access the solver for experiment-specific fields, sources and diagnostics.
        Solver& solver() { return *solver_m; }
        /// Access the immutable material distribution and its cached coefficient fields.
        material::MaterialOwner& material() { return *owner_m; }

    private:
        /** @brief Optionally place an explicit boundary between two MPI ranks before allocation.
         * @details MPI_Comm_size queries the field layout's communicator. The chosen axis is
         * already validated by makeBoxFieldContainer. All field allocations subsequently use
         * this same partition, including the material distribution.
         */
        void applyPartitionSeam(const RunBox& box) {
            if (!box.seam_m)
                return;
            int size = 1;
            MPI_Comm_size(fields_m.getFL().comm, &size);
            const auto axis = box.splitAxis_m;
            if (size != 2 || box.seam_m < 1 || box.seam_m >= box.cells_m[axis])
                throw std::invalid_argument(
                    "A seam requires two ranks and an interior cell boundary");
            auto left = fields_m.getFL().getDomain(), right = left;
            left[axis]  = ippl::Index(0, box.seam_m - 1);
            right[axis] = ippl::Index(box.seam_m, box.cells_m[axis] - 1);
            fields_m.getFL().updateLayout({left, right});
        }

        /// Construct the selected material update or the explicitly requested unit-permittivity
        /// path.
        std::unique_ptr<Solver> makeSolver(const SolverSettings& settings) {
            if (settings.formulation_m)
                return std::make_unique<Solver>(
                    fields_m.getJ(), fields_m.getE(), fields_m.getB(), settings.reconstruction_m,
                    settings.spatial_m, *owner_m, *settings.formulation_m, settings.murSpeed_m);
            if (!owner_m->statistics().homogeneous() || owner_m->homogeneousValue() != 1.0)
                throw std::invalid_argument("The vacuum update requires uniform unit permittivity");
            return std::make_unique<Solver>(fields_m.getJ(), fields_m.getE(), fields_m.getB(),
                                            settings.reconstruction_m, settings.spatial_m,
                                            settings.murSpeed_m);
        }

        solver::ChDRFieldContainer fields_m;
        std::unique_ptr<material::MaterialOwner> owner_m;
        std::unique_ptr<Solver> solver_m;
        Observer observer_m;
    };
}  // namespace chdr::simulation
