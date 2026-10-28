/** @file test_serial_reference.cpp
 * @brief Sourced-wave test of the communicator-self comparison on both spatial readings.
 * @details The separate serial mesh, layout and field storage must reproduce the
 * distributed evolution while each communicator performs its own halo exchanges.
 */
#include "LorenzWave.h"
#include "SerialReference.h"
#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"

static_assert(
    chdr::SpatialDim == 3,
    "test_serial_reference.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr;
    using namespace chdr::solver;
    using namespace chdr::solver::test;

    void sourcedWitness(SpatialDiscretization spatialDiscretization) {
        auto fields = makeBoxFieldContainer({16, 16, 16}, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  EBReconstructionConvention::R2,
                                                  spatialDiscretization);
        const auto wave = makeSourcedWave();
        loadHistory(solver, wave);
        solver.evaluateInitialFields();
        constexpr int STEPS = 8;
        for (int step = 0; step < STEPS; ++step) {
            writeSource(solver, wave.sources(), step);
            solver.solve();
        }
        verifySerialReference(solver, wave, STEPS, spatialDiscretization);
    }
    TEST(SerialReference, SourcedN16Witness) {
        for (auto spatialDiscretization :
             {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered})
            sourcedWitness(spatialDiscretization);
    }
}  // namespace
