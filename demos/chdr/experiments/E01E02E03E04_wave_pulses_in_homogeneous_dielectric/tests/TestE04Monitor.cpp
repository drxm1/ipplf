/** @file TestE04Monitor.cpp
 * @brief Check E04 monitor colocation, face ownership and normalization with affine fields.
 */
#include <array>
#include <gtest/gtest.h>
#include <limits>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03E04Setup.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04Monitor.h"

namespace {
    namespace spherical = chdr::experiments::e03;
    using Solver        = spherical::SphericalSimulation::Solver;

    /** @brief Set owned Ez=x+z and By=-1 at their native component positions.
     * @details E's z offset makes its monitor colocation read across a z-rank seam.
     * Clear halos first, so a missing exchange cannot supply valid stale samples.
     * Component positions follow LevelContract [ryu2016potentialfdtd, Appendix A].
     */
    void loadAffineFields(Solver& state, bool nonfinite = false) {
        auto& e             = state.electricField();
        auto& b             = state.magneticField();
        e                   = 0;
        b                   = 0;
        const auto electric = e.getView(), magnetic = b.getView();
        const auto first = e.getLayout().getLocalNDIndex().first();
        const auto h = e.get_mesh().getMeshSpacing(), origin = e.get_mesh().getOrigin();
        const auto offset    = state.levelContract().positionOffset(chdr::Quantity::Electric, 2);
        const int ghost      = e.getNghost();
        const double invalid = std::numeric_limits<double>::quiet_NaN();
        Kokkos::parallel_for(
            "E04 affine monitor input", e.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                const int gx = i + first[0] - ghost, gy = j + first[1] - ghost;
                const int gz         = k + first[2] - ghost;
                const double x       = origin[0] + (gx + 0.5 + offset[0]) * h[0];
                const double z       = origin[2] + (gz + 0.5 + offset[2]) * h[2];
                electric(i, j, k)[2] = nonfinite && gx == 3 && gy == 4 && gz == 7 ? invalid : x + z;
                magnetic(i, j, k)[1] = -1;
            });
        Kokkos::fence();
    }

    /** @brief Integrate an affine Poynting flux whose midpoint quadrature is exact.
     * @details E cross B=(x+z,0,0) [ellingson2020ElectromagneticsVol2, section 3.2].
     * On x=3 and x=9 the mean fluxes are 9 and 15 over area 36, giving inward
     * power 324 and outward power 540. Tangential face powers are zero.
     * The same expected values apply on one rank and two z slabs.
     */
    TEST(E04Monitor, AffineFluxAcrossNativeOffsetsAndRankSeam) {
        auto run    = spherical::makeSphericalSimulation(spherical::sphericalGrid(16, true),
                                                         spherical::sphericalReference("TEr"));
        auto& state = run.solver();
        loadAffineFields(state);
        const auto powers = chdr::experiments::e04::measureE04Monitor(state, 1, 1);
        std::array<double, 12> expected{};
        expected[0] = 324;
        expected[3] = 540;
        for (unsigned i = 0; i < expected.size(); ++i)
            EXPECT_NEAR(powers[i], expected[i], 1e-12) << "power index " << i;
        loadAffineFields(state, true);
        EXPECT_THROW(chdr::experiments::e04::measureE04Monitor(state, 1, 1), std::runtime_error);
    }
}  // namespace
