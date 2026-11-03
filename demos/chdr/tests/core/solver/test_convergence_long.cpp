/** @file test_convergence_long.cpp
 * @brief Full refinement ladders and the explicitly disabled independent anisotropic case.
 * @details LongConvergence.Standard and LongConvergence.NonStandard each run the six
 * collocated ladders of test_convergence.cpp (mixed and sourced wave; R1, R2, R2Prime) on
 * \f$N=16,32,64,128\f$, with \f$N\times N\times 3N/2\f$ cells for NonStandard. They repeat
 * the assertions on the pair (32, 64) and add the pair (64, 128) with the same rules. On
 * (64, 128) only, R1 must also give E and Gauss orders between 0.8 and 1.2. The reference is
 * the continuum wave of LorenzWave.h; the final identities must hold at \f$10^{-12}\f$ on
 * every level. CTest runs them as chdr.solver.long.Convergence.Standard and
 * chdr.solver.long.Convergence.NonStandard, each serial and mpi2, with the label long. Their
 * records go to the directories long-serial and long-mpi2, which no rank comparison reads.
 * DISABLED_Holdout.Collocated runs the independent anisotropic case of Convergence.h. It was
 * the held-out final acceptance check of earlier work; it is disabled and has no CTest entry.
 */
#include "Convergence.h"
#include "core/base/SpatialDimension.h"

static_assert(
    chdr::SpatialDim == 3,
    "test_convergence_long.cpp requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    /// Six collocated Standard ladders on \f$N^3\f$ cells, \f$N=16\f$ to 128, checked on the
    /// pairs (32, 64) and (64, 128).
    TEST(LongConvergence, Standard) {
        for (bool sourced : {false, true})
            for (auto ebReconstructionConvention :
                 {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
                  EBReconstructionConvention::R2Prime})
                waveLadder<StandardSolverPeriodic>(sourced, ebReconstructionConvention, true);
    }
    /// As LongConvergence.Standard on NonStandard \f$N\times N\times 3N/2\f$ cells.
    TEST(LongConvergence, NonStandard) {
        for (bool sourced : {false, true})
            for (auto ebReconstructionConvention :
                 {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
                  EBReconstructionConvention::R2Prime})
                waveLadder<NonStandardSolverPeriodic>(sourced, ebReconstructionConvention, true);
    }
    /// Disabled: the sourced R2 ladder on the independent anisotropic mesh family, collocated.
    TEST(DISABLED_Holdout, Collocated) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R2, true,
                                           SpatialDiscretization::Collocated, true);
    }
}  // namespace chdr::solver::test
