/** @file InterfaceRun.h
 * @brief Own a periodic material run and dispatch completed states to a callback.
 */
#pragma once
#include "BoxSetup.h"
#include "FieldEnergy.h"
#include "Manager/BaseManager.h"
#include "PlanarSlab.h"

namespace chdr::interface {

    using InterfaceSolver = solver::ChdrSolver<solver::StandardSolverPeriodic>;

    /** @brief Mesh and clock parameters in normalized units, with an optional two-rank seam.
     * The seam is a global x cell boundary. Zero selects the normal balanced partition.
     */
    struct RunBox {
        std::array<int, 3> cells_m;
        std::array<double, 3> spacing_m;
        double courant_m = 0.9;
        int seam_m       = 0;
        std::array<double, 3> lengths() const {
            return {cells_m[0] * spacing_m[0], cells_m[1] * spacing_m[1],
                    cells_m[2] * spacing_m[2]};
        }
    };

    /** @brief Keep mesh, immutable material, source and solver addresses stable for a whole run.
     * @details All ranks supply identical inputs. Samplers are device-copyable and use
     * the positions and time levels supplied by loadHistory. Each callback sees the initial
     * state and every completed solve; it chooses its own diagnostic/output cadence.
     */
    class InterfaceRun : public ippl::BaseManager {
    public:
        template <class Map>
        InterfaceRun(RunBox box, Map map)
            : fields_m(solver::makeBoxFieldContainer(box.cells_m, box.lengths())) {
            // TODO: What does this do?
            int size = 1;
            MPI_Comm_size(fields_m.getFL().comm, &size);

            // TODO: Some MPI checks it seems
            if (box.seam_m && (size != 2 || box.seam_m < 1 || box.seam_m >= box.cells_m[0]))
                throw std::invalid_argument("A seam requires two ranks and an interior x boundary");
            if (box.seam_m) {
                auto left = fields_m.getFL().getDomain(), right = left;
                left[0]  = ippl::Index(0, box.seam_m - 1);
                right[0] = ippl::Index(box.seam_m, box.cells_m[0] - 1);
                fields_m.getFL().updateLayout({left, right});
            }

            // Initialize fields A, J, E, B by initializing a FELFieldContainer.
            solver::initializeFieldStorage(fields_m);

            // Create a Material with mesh and layout taken from the field container, and a Map
            // providing geometry.
            owner_m = std::make_unique<material::MaterialOwner>(fields_m.getMesh(),
                                                                fields_m.getFL(), map);

            // Create a ChDRSolver and choose its parameters:
            //  - EBReconstructionConvention
            //  - SpatialDiscretization
            //  - MaterialFormulationParameters
            solver_m = std::make_unique<InterfaceSolver>(
                fields_m.getJ(), fields_m.getE(), fields_m.getB(),
                solver::EBReconstructionConvention::R2, solver::SpatialDiscretization::Staggered,
                *owner_m, solver::FormulationParameters{solver::MaterialFormulation::O1VS, 1});

            // Initialize the solver (TODO: Why do we need a separate initialize step that is
            // separate from the solvers constructor?)
            solver_m->initialize();

            // Explicitly choose a safety factor (TODO: Why is that not covered by solver
            // initialization?)
            solver_m->setSafetyFactor(box.courant_m);
        }

        template <class Sampler>
        void load(Sampler sample) {
            solver::loadHistory(*solver_m, sample);
            solver_m->evaluateInitialFields();
        }

        void advance() override { this->solver_m->solve(); }

        template <class Callback>
        void advance(int steps, Callback sample) {
            for (int k = 0; k <= steps; ++k) {
                sample(k, *solver_m, *owner_m);
                if (k < steps)
                    this->advance();
            }
        }
        InterfaceSolver& solver() { return *solver_m; }
        material::MaterialOwner& material() { return *owner_m; }

    private:
        solver::ChDRFieldContainer fields_m;
        std::unique_ptr<material::MaterialOwner> owner_m;
        std::unique_ptr<InterfaceSolver> solver_m;
    };
}  // namespace chdr::interface
