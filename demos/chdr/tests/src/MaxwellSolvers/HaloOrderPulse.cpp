/** @file HaloOrderPulse.cpp
 * @brief Pulse checks of the executable HaloOrder: one rank against two ranks for a Gaussian
 * pulse that crosses the rank seam near an absorbing face.
 * @details The all-absorbing boxes have 16 by 16 by 48 cells and \f$\Delta t = 0.5\f$, with
 * spacing (1, 1, 1) for Standard and (1, 1, 0.5) for NonStandard. The world layout splits
 * along \f$z\f$. Over 256 steps, the maximum difference of the owned potential samples from a
 * separate MPI_COMM_SELF solver, divided by the maximum of the one-rank solution, must be at
 * most \f$10^{-13}\f$. The old-order control of ShippedOrderSolver.h must differ by at least
 * \f$10^{-3}\f$. The test prints whether the two runs agree bit for bit; it does not assert it.
 *
 * Both tests require exactly two ranks (ASSERT_EQ on the communicator size). They run only in
 * chdr.src.HaloOrder.mpi2; chdr.src.HaloOrder.serial compiles this file but filters out its
 * tests. HaloOrderMixed.cpp also uses difference(), which this file defines.
 */
#include <algorithm>
#include <bit>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>

#include "HaloOrderTest.h"
#include "ShippedOrderSolver.h"

namespace fdtd_test {
    /** @brief Load a Gaussian intersecting an absorbing face near a rank seam.
     * @details Both histories have the same values: sigma=2.5 cells, center
     * x=4 and z=24 in cell coordinates, independent of y. All four components
     * share the scalar pulse. Initial halo exchange is part of load().
     */
    template <class Solver>
    void pulse(Solver& solver) {
        for (auto* field : {&solver.A_nm1, &solver.A_n})
            load(*field, [](int i, int, int k, unsigned) {
                const double x = (i + .5 - 4) / 2.5, z = (k + .5 - 24) / 2.5;
                return std::exp(-.5 * (x * x + z * z));
            });
    }

    /** @brief Accumulate per-cell differences from a separate serial solver.
     * @details Numerator and denominator are maxima across all owned samples
     * and steps, not a ratio formed after the pulse has left the domain.
     */
    void difference(Potential& distributed, Potential& serial, std::array<double, 2>& maximum,
                    int& identical) {
        const auto reference =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), serial.getView());
        const int g = serial.getNghost();
        owned(distributed, [&](int i, int j, int k, const auto& value) {
            for (unsigned q = 0; q < 4; ++q) {
                const double exact = reference(i + g, j + g, k + g)[q];
                identical &=
                    std::bit_cast<std::uint64_t>(value[q]) == std::bit_cast<std::uint64_t>(exact);
                maximum[0] = std::isfinite(value[q]) && std::isfinite(exact)
                                 ? std::max(maximum[0], std::abs(value[q] - exact))
                                 : std::numeric_limits<double>::infinity();
                maximum[1] = std::max(maximum[1], std::abs(exact));
            }
        });
    }

    /** @brief Compare the same 256-step all-absorbing run on one and two ranks.
     * @details The independent serial layout uses MPI_COMM_SELF. The global
     * domain has 16x16x48 cells with dt=0.5. Standard uses spacing (1,1,1);
     * NonStandard uses (1,1,0.5). The world layout splits z, so its seam cuts
     * the side faces updated by Mur. Both use the same pulse and source data.
     */
    template <class Solver>
    std::pair<double, bool> rankDifference(ippl::Vector<double, 3> spacing) {
        Fields fields({16, 16, 48}, spacing, false);
        Fields reference({16, 16, 48}, spacing, false, MPI_COMM_SELF);
        Solver distributed(fields.source, fields.electric, fields.magnetic);
        Solver serial(reference.source, reference.electric, reference.magnetic);
        EXPECT_EQ(distributed.getDt(), .5);
        EXPECT_EQ(serial.getDt(), .5);
        pulse(distributed);
        pulse(serial);
        std::array<double, 2> local{}, global{};
        int identical = 1, globallyIdentical = 0;
        for (int step = 0; step < 256; ++step) {
            distributed.step();
            serial.step();
            difference(distributed.A_np1, serial.A_np1, local, identical);
            distributed.timeShift();
            serial.timeShift();
        }
        MPI_Allreduce(local.data(), global.data(), 2, MPI_DOUBLE, MPI_MAX, fields.layout.comm);
        MPI_Allreduce(&identical, &globallyIdentical, 1, MPI_INT, MPI_MIN, fields.layout.comm);
        EXPECT_GT(global[1], 0);
        return {global[0] / global[1], globallyIdentical != 0};
    }

    /// Apply the same rank-agreement threshold and stale-halo control to either stencil.
    template <class Solver>
    void checkPulse(ippl::Vector<double, 3> spacing) {
        ASSERT_EQ(ippl::Comm->size(), 2);
        const auto repaired = rankDifference<Solver>(spacing);
        const auto control  = rankDifference<ShippedOrderSolver<Solver>>(spacing);
        EXPECT_LE(repaired.first, 1e-13);
        EXPECT_GE(control.first, 1e-3);
        std::cout << "pulse repaired=" << repaired.first << " control=" << control.first
                  << " bitwise=" << repaired.second << '\n';
    }

    /// Standard stencil, spacing (1, 1, 1): rank difference at most 1e-13, control at least 1e-3.
    TEST(HaloOrder, StandardPulseAcrossSeam) {
        checkPulse<Standard<ippl::absorbing>>({1, 1, 1});
    }

    /// NonStandard stencil, spacing (1, 1, 0.5): rank difference at most 1e-13, control at least
    /// 1e-3.
    TEST(HaloOrder, NonStandardPulseAcrossSeam) {
        checkPulse<NonStandard<ippl::absorbing>>({1, 1, .5});
    }
}  // namespace fdtd_test
