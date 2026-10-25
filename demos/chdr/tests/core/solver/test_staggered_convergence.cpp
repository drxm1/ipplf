/** @file test_staggered_convergence.cpp
 * @brief Continuum convergence of the Standard stencil with staggered reconstruction.
 * @details Scalar potentials stay at cell centres, vector potentials and electric
 * components occupy positive faces, and magnetic components occupy edges.
 * StaggeredReconstruction documents the matching forward/backward differences
 * [ryu2016potentialfdtd, App. A, Eqs. (19)-(24)]. The analytic waves and error
 * norms are shared with the collocated tests through Convergence.h. The disabled
 * holdout runs the same harness in record-only mode for this spatial reading.
 */
#include "Convergence.h"
#include "core/base/SpatialDimension.h"

namespace chdr::solver::test {
    TEST(StaggeredConvergence, Quick) {
        for (bool sourced : {false, true})
            waveLadder<StandardSolverPeriodic>(sourced, EBReconstructionConvention::R2, false,
                                               SpatialDiscretization::Staggered);
    }
    TEST(LongStaggeredConvergence, Standard) {
        for (bool sourced : {false, true})
            waveLadder<StandardSolverPeriodic>(sourced, EBReconstructionConvention::R2, true,
                                               SpatialDiscretization::Staggered);
    }
    TEST(DISABLED_Holdout, StaggeredRecord) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R2, true,
                                           SpatialDiscretization::Staggered, true, true);
    }
}  // namespace chdr::solver::test
