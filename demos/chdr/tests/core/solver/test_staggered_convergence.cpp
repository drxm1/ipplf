/** @file test_staggered_convergence.cpp
 * @brief Continuum convergence of the Standard stencil with staggered reconstruction.
 * @details Scalar potentials stay at cell centres, vector potentials and electric
 * components occupy positive faces, and magnetic components occupy edges.
 * StaggeredReconstruction documents the matching forward/backward differences
 * [ryu2016potentialfdtd, App. A, Eqs. (19)-(24)]. The analytic waves and error
 * norms are shared with the collocated tests through Convergence.h. The independent anisotropic
 * case runs the same harness in record-only mode for this spatial reading.
 *
 * StaggeredConvergence.Quick runs the mixed and the sourced wave with R2 on \f$N^3\f$ cells
 * of the unit box, \f$N=16,32,64\f$ (split along x, default step, one period). The E, B and
 * Gauss orders must be at least 1.8 on the pair (32, 64), the gauge order of the sourced
 * wave too when its fine error exceeds \f$10^{-12}\f$, and the final identities must hold at
 * \f$10^{-12}\f$. The reference is the continuum wave of LorenzWave.h. The test writes the
 * states at \f$N=16\f$ and 32 as rank records for chdr.solver.ranks.
 * LongStaggeredConvergence.Standard adds \f$N=128\f$ and the pair (64, 128). CTest runs them
 * as chdr.solver.convergence.staggered and chdr.solver.long.staggered, each serial and mpi2.
 * AnisotropicConvergence.StaggeredRecord runs as
 * chdr.solver.long.anisotropic.StaggeredRecord.mpi2.
 */
#include "Convergence.h"

namespace chdr::solver::test {
    /// Mixed and sourced wave, staggered R2, \f$N=16,32,64\f$: orders at least 1.8 on
    /// (32, 64).
    TEST(StaggeredConvergence, Quick) {
        for (bool sourced : {false, true})
            waveLadder<StandardSolverPeriodic>(sourced, EBReconstructionConvention::R2, false,
                                               SpatialDiscretization::Staggered);
    }
    /// The same two ladders up to \f$N=128\f$, checked on (32, 64) and (64, 128).
    TEST(LongStaggeredConvergence, Standard) {
        for (bool sourced : {false, true})
            waveLadder<StandardSolverPeriodic>(sourced, EBReconstructionConvention::R2, true,
                                               SpatialDiscretization::Staggered);
    }
    /// The sourced ladder on the anisotropic mesh family, staggered, record only.
    TEST(AnisotropicConvergence, StaggeredRecord) {
        waveLadder<StandardSolverPeriodic>(true, EBReconstructionConvention::R2, true,
                                           SpatialDiscretization::Staggered, true, true);
    }
}  // namespace chdr::solver::test
