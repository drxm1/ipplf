/** @file Convergence.h
 * @brief Execute source-free and sourced Lorenz-wave refinement ladders.
 * @details WaveNorms.h defines continuum comparisons. Monitors verifies the
 * independent discrete identities on the last evolved step. RankRecords.h
 * records owned fields for decomposition comparisons; approximation errors
 * remain separate JSONL records.
 * Helpers parameterized by PotentialSolverBase select the inherited Standard or
 * NonStandard potential update. SolverType instead denotes the complete ChdrSolver
 * interface, including its potential histories, source field and E/B outputs.
 *
 * TODO: REVIEW.
 */
#pragma once
#include <gtest/gtest.h>
#include <iostream>

#include "ChdrSolver.h"
#include "Monitors.h"
#include "PeriodicSetup.h"
#include "RankRecords.h"
#include "SerialReference.h"
#include "SpatialDimension.h"
#include "WaveNorms.h"
#include "Waves.h"

static_assert(chdr::solver::SpatialDim == 3,
              "Convergence.h requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    struct WaveOptions {
        bool sourced_m;
        EBReconstructionConvention ebReconstructionConvention_m;
        SpatialDiscretization spatialDiscretization_m;
        bool holdout_m, recordOnly_m;
    };

}  // namespace chdr::solver::test
