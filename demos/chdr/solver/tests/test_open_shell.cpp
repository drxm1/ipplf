/** @file test_open_shell.cpp
 * @brief Check open-box construction, physical Mur speeds and supported observations.
 */
#include <gtest/gtest.h>

#include "BoxSetup.h"
#include "CoreTestSupport.h"
#include "opt-in/diagnostics/FourierProbe.h"
#include "opt-in/diagnostics/Monitors.h"

namespace chdr::solver::test {
    /** @brief Require the boundary description and its speed to survive a storage reset.
     * Each topology is constructed identically on every rank. Unsupported
     * NonStandard open boxes must throw before entering a halo collective.
     */
    TEST(OpenShell, TopologyAndReset) {
        for (auto axes : {std::array{true, true, true}, std::array{false, false, false},
                          std::array{false, true, true}, std::array{false, false, true}}) {
            auto fields = makeBoxFieldContainer({8, 6, 10}, {8, 6, 10}, 2, axes);
            initializeFieldStorage(fields);
            ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                      EBReconstructionConvention::R2);
            const auto boundary = solver.boundary();
            solver.initialize();
            EXPECT_EQ(solver.boundary().kind, boundary.kind);
            EXPECT_EQ(solver.boundary().speed, boundary.speed);
            for (unsigned d = 0; d < 3; ++d) {
                EXPECT_EQ(fields.getFL().isPeriodic(d), axes[d]);
                EXPECT_EQ(solver.boundary().kind[2 * d],
                          axes[d] ? ippl::FaceKind::periodic : ippl::FaceKind::mur);
            }
            if (!fields.getFL().isAllPeriodic_m) {
                EXPECT_THROW((ChdrSolver<NonStandardSolverPeriodic>(
                                 fields.getJ(), fields.getE(), fields.getB(),
                                 EBReconstructionConvention::R1)),
                             IpplException);
            }
        }
    }

    /** @brief Match the medium's internal step without applying its refractive index twice.
     * The expected kernel factor follows directly from the physical/internal
     * timestep ratio; Mur weights depend on speed times timestep
     * [mur1981AbsorbingBoundaryFdtd, section IV, equation (16)].
     */
    TEST(OpenShell, PhysicalSpeedConversion) {
        auto fields = makeBoxFieldContainer({8, 6, 10}, {8, 6, 10}, 2, {false, true, true});
        initializeFieldStorage(fields);
        const material::MaterialOwner owner(
            fields.getMesh(), fields.getFL(),
            KOKKOS_LAMBDA(material::CellIndex, material::Position, bool& marked) {
                marked = false;
                return 2.13;
            });
        ChdrSolver<StandardSolverPeriodic> matched(fields.getJ(), fields.getE(), fields.getB(),
                                                   EBReconstructionConvention::R2,
                                                   SpatialDiscretization::Collocated, owner);
        EXPECT_EQ(matched.boundary().speed[0], 1.0);
        ChdrSolver<StandardSolverPeriodic> physical(
            fields.getJ(), fields.getE(), fields.getB(), EBReconstructionConvention::R2,
            SpatialDiscretization::Collocated, owner, {}, MurSpeed::physical(0.7));
        EXPECT_EQ(physical.boundary().speed[0], 0.7 * std::sqrt(2.13));
        physical.initialize();
        EXPECT_EQ(physical.boundary().speed[0], 0.7 * std::sqrt(2.13));
        EXPECT_THROW(MurSpeed::physical(0), std::invalid_argument);
        EXPECT_THROW(MurSpeed::physical(std::numeric_limits<double>::quiet_NaN()),
                     std::invalid_argument);
    }

}  // namespace chdr::solver::test
