/** @file test_mur_faces.cpp
 * @brief Check the Mur face description of ChdrSolver and its physical speed conversion.
 * @details Box: 8 x 6 x 10 cells, unit spacing, split along z on two ranks.
 * DescriptionSurvivesStorageReset uses four periodicity patterns: all periodic, all open,
 * x open, and x and y open. PhysicalSpeedConversion opens x and sets a uniform relative
 * permittivity 2.13. The checks are exact (EXPECT_EQ) or require an exception; there is no
 * numerical tolerance. The references are the layout and the closed form of the speed slot:
 * the kernel Mur speed is the physical speed times the refractive index
 * \f$n = \sqrt{\varepsilon_r}\f$, and the matched speed is exactly 1.
 * CTest: chdr.solver.mur.faces.serial and chdr.solver.mur.faces.mpi2 (filter MurFaces.*).
 */
#include <gtest/gtest.h>

#include "core/solver/BoxSetup.h"

namespace chdr::solver::test {
    /** @brief Require the boundary description and its speed to survive a storage reset.
     * For each layout, the lower face of each axis is periodic or Mur as the layout says,
     * and initialize() leaves the face kinds and speeds unchanged.
     * Each topology is constructed identically on every rank. Unsupported
     * NonStandard open boxes (the three layouts with a Mur axis) must throw
     * IpplException before entering a halo collective.
     */
    TEST(MurFaces, DescriptionSurvivesStorageReset) {
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
     * At \f$\varepsilon_r = 2.13\f$ the matched speed slot is exactly 1.0, and
     * MurSpeed::physical(0.7) stores exactly \f$0.7\sqrt{2.13}\f$, also after initialize().
     * Speeds 0 and NaN throw std::invalid_argument.
     */
    TEST(MurFaces, PhysicalSpeedConversion) {
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
