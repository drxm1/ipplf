/** @file test_medium_guard.cpp
 * @brief Check cached validity decisions and collective refusal before a material step.
 */
#include <gtest/gtest.h>

#include "ChdrSolver.h"
#include "MediumCases.h"
#include "MediumOracles.h"
#include "Monitors.h"
#include "PeriodicSetup.h"

namespace chdr::solver::test {
    /** @brief Check map-wide validity and the bound without advancing a solver.
     * @details Refused rows exercise the cached-statistics guard directly. Accepted
     * rows also construct the update and compare its physical-step bound with the
     * independently evaluated mediumMaxStableTimeStep oracle, including vacuum.
     */
    TEST(MediumValidity, CachedTableAndBound) {
        auto fields = makePeriodicFieldContainer({8, 8, 8}, {1, 1, 1}, 0);
        initializeFieldStorage(fields);
        for (const auto row : MediumCases) {
            const material::MaterialOwner owner(fields.getMesh(), fields.getFL(), GuardMap{row});
            const FormulationParameters parameters{MaterialFormulation::Homogeneous, row.alpha_m};
            const SolverDiscretizationSelection selection{
                row.nonstandard_m ? Stencil::NonStandard : Stencil::Standard,
                SpatialDiscretization::Collocated, EBReconstructionConvention::R2};
            if (!row.accepted_m) {
                EXPECT_THROW(validateMaterialUpdate(parameters, owner.statistics(), selection),
                             IpplException);
                continue;
            }
            auto update  = makeMaterialUpdate(parameters, owner, selection, fields.getJ());
            const auto h = fields.getHr();
            EXPECT_LE(
                std::abs(update->maxStableTimeStep(h) / mediumMaxStableTimeStep(row.epsilon_m, h)
                         - 1),
                1e-15);
        }
    }

}  // namespace chdr::solver::test
