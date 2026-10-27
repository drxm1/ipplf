/** @file test_medium_guard.cpp
 * @brief Check cached validity decisions and collective refusal before a material step.
 */
#include <gtest/gtest.h>

#include "MediumCases.h"
#include "MediumOracles.h"
#include "core/solver/BoxSetup.h"
#include "core/solver/ChdrSolver.h"
#include "opt-in/diagnostics/Monitors.h"

namespace chdr::solver::test {
    using namespace chdr::diagnostics;
    /** @brief Check map-wide validity and the bound without advancing a solver.
     * @details Refused rows exercise the cached-statistics guard directly. Accepted
     * rows also construct the update and compare its physical-step bound with the
     * independently evaluated mediumMaxStableTimeStep oracle, including vacuum.
     */
    TEST(MediumValidity, CachedTableAndBound) {
        auto fields = makeBoxFieldContainer({8, 8, 8}, {1, 1, 1}, 0);
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

    /** @brief Require the same constructor verdict on every rank, then check monitor eligibility.
     * @details An IpplException counts as a constructor refusal. Accepted unit material
     * must also enter the vacuum monitor; accepted non-unit material must refuse that
     * monitor with std::logic_error. These are separate contracts: rejecting the
     * diagnostic does not count as rejecting an otherwise valid material solver.
     * The final reduction requires either zero constructor refusals or one per rank.
     * It observes the result of the cached guard; it is not part of that guard.
     */
    template <class Base, class Fields>
    void constructorRow(const MediumCase& row, SpatialDiscretization reading, Fields& fields,
                        const material::MaterialOwner& owner) {
        int refused = 0;
        try {
            ChdrSolver<Base> solver(fields.getJ(), fields.getE(), fields.getB(),
                                    EBReconstructionConvention::R2, reading, owner,
                                    {MaterialFormulation::Homogeneous, row.alpha_m});
            Monitors monitors(solver.A_n);
            if (row.epsilon_m == 1) {
                EXPECT_NO_THROW(monitors.beforeSolve(solver));
            } else {
                EXPECT_THROW(monitors.beforeSolve(solver), std::logic_error);
            }
        } catch (const IpplException&) {
            refused = 1;
        }
        int count = 0, size = 0;
        MPI_Comm_size(fields.getFL().comm, &size);
        MPI_Allreduce(&refused, &count, 1, MPI_INT, MPI_SUM, fields.getFL().comm);
        EXPECT_EQ(count, row.accepted_m ? 0 : size);
    }

    /** @brief Exercise full construction for each map defect and both spatial placements.
     * @details Each rank constructs the same owner and selection. A one-cell defect
     * is owned by only one rank, so this also checks that global material statistics
     * make every rank refuse before entering the inherited solver's communication.
     */
    TEST(MediumGuard, ConstructorTableAndMonitor) {
        auto fields = makeBoxFieldContainer({8, 8, 12}, {1, 1, 1}, 0);
        initializeFieldStorage(fields);
        for (const auto row : MediumCases) {
            const material::MaterialOwner owner(fields.getMesh(), fields.getFL(), GuardMap{row});
            for (auto reading :
                 {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered}) {
                if (row.nonstandard_m)
                    constructorRow<NonStandardSolverPeriodic>(row, reading, fields, owner);
                else
                    constructorRow<StandardSolverPeriodic>(row, reading, fields, owner);
            }
        }
    }
}  // namespace chdr::solver::test
