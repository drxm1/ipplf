/** @file test_convergence_long.cpp
 * @brief Full refinement ladders and the explicitly disabled independent anisotropic case.
 *
 * TODO: REVIEW.
 */
#include "Convergence.h"
#include "SpatialDimension.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_convergence_long.cpp requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    TEST(LongConvergence, Standard) {
        for (bool sourced : {false, true})
            for (auto ebReconstructionConvention :
                 {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
                  EBReconstructionConvention::R2Prime})
                waveLadder<StandardSolverPeriodic>(sourced, ebReconstructionConvention, true);
    }
    TEST(LongConvergence, NonStandard) {
        for (bool sourced : {false, true})
            for (auto ebReconstructionConvention :
                 {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
                  EBReconstructionConvention::R2Prime})
                waveLadder<NonStandardSolverPeriodic>(sourced, ebReconstructionConvention, true);
    }
    TEST(DISABLED_Holdout, Collocated) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R2, true,
                                           SpatialDiscretization::Collocated, true);
    }
}  // namespace chdr::solver::test
