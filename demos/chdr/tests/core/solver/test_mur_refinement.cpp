/** @file test_mur_refinement.cpp
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
    TEST(MurReflectionRefinement, CentrePlaneOrder) {
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

    /// Three-mesh angular and evanescent tables; evolved marks only an independently checked coarse
    /// row.
    TEST(MurReflectionRefinement, CodeConfirmedRows) {
        const auto precision = std::cout.precision(17);
        const auto record    = [](const Case& row, bool normal, bool evolved, int refinement) {
            auto kt = normal ? std::array<double, 3>{} : row.tangent();
            for (auto& value : kt)
                value /= refinement;
            const double omega = row.omega_m / refinement, speed = row.speed();
            const auto discrete = discreteReflection(row.h_m, row.dt_m, speed, speed,
                                                        OracleStencil::Standard, 0, kt, omega);
            const auto continuum =
                continuumReflection(speed, speed, std::hypot(kt[1], kt[2]), omega);
            // A unit physical Mur speed in a medium is a closed-form record
            // only; it is outside the supported matched-intersection policy.
            const auto unitSpeed = discreteReflection(row.h_m, row.dt_m, speed, 1.,
                                                         OracleStencil::Standard, 0, kt, omega);
            std::cout << "reflection-record " << row.name_m << " refinement=" << refinement
                      << " evolved=" << (evolved && refinement == 1) << " omega=" << omega
                      << " dt=" << row.dt_m << " nominal_ppw=" << 37.5 * refinement
                      << " tangential_ratio=" << speed * std::hypot(kt[1], kt[2]) / omega
                      << " medium_ppw_z=" << 2 * std::numbers::pi * speed / (omega * row.h_m[2])
                      << " h=" << row.h_m[0] << ',' << row.h_m[1] << ',' << row.h_m[2]
                      << " layer-zero=" << discrete.layerZero << " centre=" << discrete.centre
                      << " continuum=" << continuum
                      << " unit-physical-speed-oracle=" << unitSpeed.layerZero << '\n';
        };
        const std::array<RefinementRow, 24> rows{
            {{"v0", 1, .16755160819145562, 0},          {"v15", 1, .16742319888454554, 145},
             {"v30", 1, .16755160819145565, 75},        {"v45", 1, .1676559599305044, 53},
             {"v60", 1, .16872552225434587, 43},        {"v75", 1, .16679056629809322, 39},
             {"v85", 1, .16597858037651997, 38},        {"vbeta", 1, .16980972104961076, 37},
             {"v1.01109", 1, .16795313592753203, 37},   {"v1.05", 1, .16622183352326947, 36},
             {"v1.1325", 1, .1681232271638126, 33},     {"v1.2", 1, .16890283083816093, 31},
             {"v1.5", 1, .16755160819145565, 25},       {"m0", 2.13, .16755160819145562, 0},
             {"m30", 2.13, .16883010330328174, 51},     {"m43.252", 2.13, .16981093559895621, 37},
             {"m46.748", 2.13, .16888211355418753, 35}, {"m50.894", 2.13, .16812232917332198, 33},
             {"m52.923", 2.13, .16862885799775393, 32}, {"m60", 2.13, .16570597950208546, 30},
             {"m70", 2.13, .16968384085682203, 27},     {"m75", 2.13, .16507546312572452, 27},
             {"m1.05", 2.13, .17083998548546364, 24},   {"m1.2", 2.13, .17083998548546367, 21}}};
        for (const auto& row : rows)
            for (double dt : {.5, .98 / std::sqrt(3.)}) {
                if (row.epsilon_m != 1 && dt != .5)
                    continue;
                bool evolved = false;
                for (const auto& checked : cases)
                    evolved |= std::string_view(checked.name_m) == row.name_m && checked.dt_m == dt;
                const Case mode{row.name_m,    2,           std::max(4, row.period_m),
                                row.epsilon_m, row.omega_m, dt};
                for (int refinement : {1, 2, 4})
                    record(mode, row.period_m == 0, evolved, refinement);
            }
        for (const auto& row : cases)
            if (row.ny_m != 2)
                for (int refinement : {1, 2, 4})
                    record(row, false, true, refinement);
        std::cout.precision(precision);
    }
}  // namespace chdr::solver::test::reflection
