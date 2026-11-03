/** @file test_serial_reference.cpp
 * @brief Sourced-wave test of the communicator-self comparison on both spatial readings.
 * @details The separate serial mesh, layout and field storage must reproduce the
 * distributed evolution while each communicator performs its own halo exchanges.
 *
 * The test makes 8 solves of the sourced Lorenz wave of makeSourcedWave() on Standard
 * 16x16x16 cells (unit box split along x, R2, default step), first collocated, then
 * staggered. verifySerialReference() repeats the run on MPI_COMM_SELF at global rank zero
 * with the same step. It requires A_n, A_nm1, E and B to agree per cell to \f$10^{-13}\f$
 * of the serial field's absolute maximum, and the maximum, the absolute sum and the square
 * sum to agree to relative \f$10^{-12}\f$. On one rank the replay repeats an identical
 * one-rank run; only the two-rank entry compares two decompositions.
 *
 * CTest entries: chdr.solver.serial_reference.{serial,mpi2} (label unit, timeout 60 s).
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
    /// Eight sourced R2 solves on 16x16x16 cells, collocated and staggered: cells agree with
    /// the rank-zero replay to \f$10^{-13}\f$ of the field maximum, statistics to \f$10^{-12}\f$.
    TEST(SerialReference, SourcedN16Witness) {
        for (auto spatialDiscretization :
             {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered})
            sourcedWitness(spatialDiscretization);
    }
}  // namespace
