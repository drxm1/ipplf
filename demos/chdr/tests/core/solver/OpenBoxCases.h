/** @file OpenBoxCases.h
 * @brief Fixed integer-period boxes for the thirty-step reflection checks.
 * @details The 13 rows of cases are 6 vacuum rows (names v...), 5 rows at
 * \f$\varepsilon_r = 2.13\f$ with matched Mur speed (m...) and 2 vacuum rows with
 * h = (1.25, 0.8, 1) (a...). The step is 0.5; the two anisotropic rows use
 * 0.5029188806044014. A row name carries the incidence angle in degrees,
 * \f$\arcsin(v\lvert\mathbf k_t\rvert/\omega)\f$ (for example v30, m43.252, a40), or, for an
 * evanescent row, the ratio \f$v\lvert\mathbf k_t\rvert/\omega > 1\f$ (for example v1.05,
 * m1.2, a1.2). Users: test_mur_reflection.cpp, test_mur_field_reflection.cpp,
 * test_open_gauge_records.cpp and test_mur_refinement.cpp.
 */
#pragma once
#include "ReflectionOracles.h"

namespace chdr::solver::test::reflection {
    /** @brief One source-free mode, with physical time step and permittivity.
     * @details The normal extent is 40 cells. Each nonzero tangential component
     * spans exactly one periodic box length; two y cells denote zero y mode.
     */
    struct Case {
        const char* name_m;
        int ny_m, nz_m;
        double epsilon_m, omega_m, dt_m = .5;
        std::array<double, 3> h_m{1, 1, 1};
        double speed() const { return 1 / std::sqrt(epsilon_m); }
        std::array<int, 3> cells() const { return {40, ny_m, nz_m}; }
        std::array<double, 3> lengths() const {
            return {40 * h_m[0], ny_m * h_m[1], nz_m * h_m[2]};
        }
        std::array<double, 3> tangent() const {
            return {0, ny_m == 2 ? 0 : 2 * std::numbers::pi / (ny_m * h_m[1]),
                    2 * std::numbers::pi / (nz_m * h_m[2])};
        }
        DiscreteReflection oracle() const {
            return discreteReflection(h_m, dt_m, speed(), speed(), OracleStencil::Standard, 0,
                                      tangent(), omega_m);
        }
    };

    /// Vacuum, medium and anisotropic rows; frequencies and integer periods are fixed.
    inline constexpr std::array<Case, 13> cases{
        {{"v30", 2, 75, 1, .16755160819145565},
         {"v60", 2, 43, 1, .16872552225434587},
         {"v85", 2, 38, 1, .16597858037651997},
         {"vbeta", 2, 37, 1, .16980972104961076},
         {"v1.05", 2, 36, 1, .16622183352326947},
         {"v1.5", 2, 25, 1, .16755160819145565},
         {"m43.252", 2, 37, 2.13, .16981093559895621},
         {"m46.748", 2, 35, 2.13, .16888211355418753},
         {"m60", 2, 30, 2.13, .16570597950208546},
         {"m1.05", 2, 24, 2.13, .17083998548546364},
         {"m1.2", 2, 21, 2.13, .17083998548546367},
         {"a40", 120, 72, 1, .16970314394025246, .5029188806044014, {1.25, .8, 1}},
         {"a1.2", 65, 39, 1, .16782012038407013, .5029188806044014, {1.25, .8, 1}}}};
}  // namespace chdr::solver::test::reflection
