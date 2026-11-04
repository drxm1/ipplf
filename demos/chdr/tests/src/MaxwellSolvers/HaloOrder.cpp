/** @file HaloOrder.cpp
 * @brief Drift checks of the finishing order of the potential step in StandardFDTDSolver and
 * NonStandardFDTDSolver with Mur faces: Mur update first, halo exchange second.
 * @details The potentials start at \f$\psi^{-1} = -1\f$ and \f$\psi^0 = 0\f$ with zero source,
 * so the exact lattice solution is \f$\psi^k = k\f$ in every component (see drift()). The test
 * requires \f$\max\lvert\psi^k - k\rvert/64 \le 10^{-12}\f$ over 64 steps in two all-absorbing
 * boxes: Standard with 16 by 16 by 16 cells of spacing 1, and NonStandard with 16 by 16 by 48
 * cells of spacing (1, 1, 0.5). The layouts split along \f$z\f$. Each box runs once on
 * MPI_COMM_SELF and once on MPI_COMM_WORLD. ShippedOrderSolver.h restores the original order
 * (halo exchange first, Mur update second) as a negative control. On two ranks the control
 * must miss by at least \f$10^{-2}\f$. On one rank it has no stale ghost copies and also
 * meets \f$10^{-12}\f$. The repaired order is the owner-to-ghost order of
 * [oskooi2010meep, section 2.3].
 *
 * This file holds main() of the executable HaloOrder. HaloOrderMixed.cpp adds the box that is
 * open in \f$x\f$ and periodic in \f$y\f$ and \f$z\f$; HaloOrderPulse.cpp adds the Gaussian
 * pulses across the rank seam. Reference: exact lattice value. The CTest entry
 * chdr.src.HaloOrder.mpi2 runs all six tests on two ranks, including the one-rank drift replays.
 */
#include <gtest/gtest.h>
#include <limits>

#include "HaloOrderTest.h"
#include "ShippedOrderSolver.h"

namespace fdtd_test {
    /** @brief Measure the maximum error of the affine-in-step solution.
     * @details Substitution of a spatial constant linear in time into IPPL's
     * potential leapfrog update and Mur's recurrence gives \f$k\f$ at step \f$k\f$ (deduction
     * from [mur1981AbsorbingBoundaryFdtd, section IV, Eq. (16), p. 380]). Both histories are
     * prepared before the first step. The maximum over all owned components and all 64 steps
     * is divided by 64.
     */
    template <class Solver>
    double drift(Fields& fields) {
        Solver solver(fields.source, fields.electric, fields.magnetic);
        load(solver.A_nm1, [](int, int, int, unsigned) {
            return -1.;
        });
        double error = 0;
        for (int step = 1; step <= 64; ++step) {
            solver.step();
            owned(solver.A_np1, [&](int, int, int, const auto& value) {
                for (unsigned q = 0; q < 4; ++q)
                    error = std::isfinite(value[q]) ? std::max(error, std::abs(value[q] - step))
                                                    : std::numeric_limits<double>::infinity();
            });
            solver.timeShift();
        }
        double maximum = 0;
        MPI_Allreduce(&error, &maximum, 1, MPI_DOUBLE, MPI_MAX, fields.layout.comm);
        return maximum / 64;
    }

    /** @brief Check drift on serial and distributed all-absorbing boxes.
     * @details The old order has no stale copies in serial. With two ranks its
     * error must be at least 0.01, which shows that the repaired tolerance of
     * 1e-12 distinguishes the ordering defect from roundoff. On one world rank the
     * distributed case repeats the serial case, and the control bound is not checked.
     */
    template <class Solver>
    void checkDrift(std::array<int, 3> n, ippl::Vector<double, 3> h) {
        Fields serial(n, h, false, MPI_COMM_SELF), distributed(n, h, false);
        EXPECT_LE(drift<Solver>(serial), 1e-12);
        EXPECT_LE(drift<ShippedOrderSolver<Solver>>(serial), 1e-12);
        const double repaired = drift<Solver>(distributed);
        const double control  = drift<ShippedOrderSolver<Solver>>(distributed);
        EXPECT_LE(repaired, 1e-12);
        if (distributed.layout.comm.size() == 2) {
            EXPECT_GE(control, 1e-2);
        }
        std::cout << "drift repaired=" << repaired << " control=" << control << '\n';
    }

    /// All-absorbing Standard box, 16 by 16 by 16 cells of spacing 1: drift at most 1e-12.
    TEST(HaloOrder, StandardDrift) {
        checkDrift<Standard<ippl::absorbing>>({16, 16, 16}, {1, 1, 1});
    }
    /// All-absorbing NonStandard box, 16 by 16 by 48 cells of spacing (1, 1, 0.5): drift at
    /// most 1e-12.
    TEST(HaloOrder, NonStandardDrift) {
        checkDrift<NonStandard<ippl::absorbing>>({16, 16, 48}, {1, 1, .5});
    }
}  // namespace fdtd_test

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
