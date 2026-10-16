#include "HaloOrderTest.h"

namespace fdtd_test {
    /** @brief Measure the maximum error of the affine-in-step solution.
     * @details Substitution of a spatial constant linear in time into IPPL's
     * potential leapfrog update and Mur's recurrence gives \f$k\f$ at step \f$k\f$ (deduction
     * from Mur, 1981, "Absorbing boundary conditions for the finite-difference
     * approximation of the time-domain electromagnetic-field equations",
     * section IV, Eq. (16), p. 380). Both histories are prepared before the first step. The
     * maximum over all owned components and all 64 steps is divided by 64.
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
     * error must exceed 0.01, demonstrating that the repaired tolerance of
     * 1e-12 distinguishes the ordering defect from roundoff.
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

    TEST(HaloOrder, StandardDrift) {
        checkDrift<Standard<ippl::absorbing>>({16, 16, 16}, {1, 1, 1});
    }
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
