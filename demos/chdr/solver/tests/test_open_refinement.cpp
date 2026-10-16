/** @file test_open_refinement.cpp
 * @brief Compare the centre-plane discrete reflection with its continuum limit.
 */
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <string_view>

#include "OpenBoxCases.h"

namespace chdr::solver::test::reflection {
    namespace {
        struct RefinementRow {
            const char* name_m;
            double epsilon_m, omega_m;
            int period_m;  // Zero denotes normal incidence, with no tangential mode.
        };

        /// Coarse integer periods are fixed; the two finer periods are exact multiples.
        constexpr std::array<RefinementRow, 11> refinementRows{
            {{"v0", 1., .16755160819145562, 0},
             {"v30", 1., .16755160819145565, 75},
             {"v45", 1., .1676559599305044, 53},
             {"v60", 1., .16872552225434587, 43},
             {"v1.2", 1., .16890283083816093, 31},
             {"v1.5", 1., .16755160819145565, 25},
             {"m0", 2.13, .16755160819145562, 0},
             {"m30", 2.13, .16883010330328174, 51},
             {"m43.252", 2.13, .16981093559895621, 37},
             {"m60", 2.13, .16570597950208546, 30},
             {"m1.2", 2.13, .17083998548546367, 21}}};

        /** @brief Absolute centre-plane coefficient error at an exact integer refinement.
         * @details Cell-unit spacing and step stay fixed while frequency and wavevector divide
         * by the refinement factor. This represents physical mesh refinement at fixed mode
         * and Courant number. Both coefficients come from the independent oracle; using
         * centre instead of layerZero removes the reference-plane phase described by
         * discreteReflection() [Mur (1981), section IV, Eq. (16)].
         */
        double refinementError(const RefinementRow& row, int refinement) {
            const double speed = 1 / std::sqrt(row.epsilon_m), omega = row.omega_m / refinement;
            const double kz =
                row.period_m == 0 ? 0. : 2 * std::numbers::pi / (row.period_m * refinement);
            const auto lattice = discreteReflection(
                {1., 1., 1.}, .5, speed, speed, OracleStencil::Standard, 0, {0., 0., kz}, omega);
            return std::abs(lattice.centre - continuumReflection(speed, speed, kz, omega));
        }
    }  // namespace

    /** Finest-pair errors must give second order without changing the coarse period by rounding.
     * The interval allows finite-mesh terms in the centred Mur discretization; the discrete
     * and continuum coefficients are independently defined in ReflectionOracles.h.
     */
    TEST(OpenReflectionRefinement, CentrePlaneOrder) {
        for (const auto& row : refinementRows) {
            SCOPED_TRACE(row.name_m);
            const double coarse = refinementError(row, 2), fine = refinementError(row, 4);
            ASSERT_GT(coarse, 0.);
            ASSERT_GT(fine, 0.);
            const double order = observedOrder(coarse, fine, 2.);
            EXPECT_GE(order, 1.9);
            EXPECT_LE(order, 2.1);
            std::cout << row.name_m << " centre-plane finest-pair order=" << order << '\n';
        }
    }

}  // namespace chdr::solver::test::reflection
