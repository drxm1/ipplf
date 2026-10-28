#include <cmath>
#include <gtest/gtest.h>

#include "HaloOrderTest.h"
#include "ShippedOrderSolver.h"

namespace fdtd_test {
    /** @brief Set the two prescribed Standard steps without changing its kernel. */
    struct MixedStandard : Standard<ippl::configured> {
        using Base = Standard<ippl::configured>;
        using Base::Base;
        void setStep(double step) { this->dt = step; }
    };

    /** @brief Reuse the owned-sample comparison from HaloOrderPulse.cpp.
     * @details It records maximum discrepancy, reference magnitude and exact
     * binary64 equality; neither physical exterior ghosts nor rank ghosts enter.
     */
    void difference(Potential&, Potential&, std::array<double, 2>&, int&);

    /** @brief Load the affine-time or Gaussian mixed-boundary fixture.
     * @details The pulse has width 2.5 cells, center x=4,z=16, and no y variation.
     * Both pulse histories are equal. Drift uses constant levels -1 and 0.
     * load() exchanges both histories before the first stencil evaluation.
     */
    template <class Solver>
    void mixedHistory(Solver& solver, bool drift) {
        for (int level = -1; level <= 0; ++level)
            load(level == -1 ? solver.A_nm1 : solver.A_n, [=](int i, int, int k, unsigned) {
                const double x = (i + .5 - 4) / 2.5, z = (k + .5 - 16) / 2.5;
                return drift ? double(level) : std::exp(-.5 * (x * x + z * z));
            });
    }

    /** @brief Measure mixed-box drift or decomposition error over the full run.
     * @details The 32x4x32 box absorbs on x and wraps on y,z. Drift compares 64
     * steps with the exact affine-time solution; see HaloOrder.cpp for its
     * derivation from Mur's recurrence. Pulse compares 256 steps with a separate
     * serial solver. The normalization is the maximum reference amplitude over
     * the run: exactly 64 for drift. No post-step fill repairs the old-order
     * control's stale periodic wraps or seam copies.
     */
    template <class Solver>
    std::pair<double, bool> mixedRun(double dt, bool drift, MPI_Comm comm) {
        const std::array<bool, 3> periodic{false, true, true};
        Fields fields({32, 4, 32}, {1, 1, 1}, periodic, comm);
        Fields reference({32, 4, 32}, {1, 1, 1}, periodic, MPI_COMM_SELF);
        Solver solver(fields.source, fields.electric, fields.magnetic,
                      ippl::FDTDBoundary::fromLayout(fields.layout));
        Solver serial(reference.source, reference.electric, reference.magnetic,
                      ippl::FDTDBoundary::fromLayout(reference.layout));
        solver.setStep(dt);
        serial.setStep(dt);
        mixedHistory(solver, drift);
        mixedHistory(serial, drift);
        std::array<double, 2> local{}, global{};
        int identical = 1, allIdentical = 0;
        for (int step = 1; step <= (drift ? 64 : 256); ++step) {
            solver.step();
            if (drift)
                serial.A_np1 = double(step);
            else
                serial.step();
            difference(solver.A_np1, serial.A_np1, local, identical);
            solver.timeShift();
            serial.timeShift();
        }
        MPI_Allreduce(local.data(), global.data(), 2, MPI_DOUBLE, MPI_MAX, comm);
        MPI_Allreduce(&identical, &allIdentical, 1, MPI_INT, MPI_MIN, comm);
        EXPECT_GT(global[1], 0);
        return {global[0] / global[1], allIdentical != 0};
    }

    TEST(HaloOrder, MixedStandardDrift) {
        for (MPI_Comm comm : {MPI_COMM_SELF, MPI_COMM_WORLD}) {
            const auto repaired = mixedRun<MixedStandard>(.5, true, comm);
            const auto control  = mixedRun<ShippedOrderSolver<MixedStandard>>(.5, true, comm);
            EXPECT_LE(repaired.first, 1e-12);
            EXPECT_GE(control.first, 1e-2);
            std::cout << "mixed drift repaired=" << repaired.first << " control=" << control.first
                      << '\n';
        }
    }

    TEST(HaloOrder, MixedStandardPulseAcrossSeam) {
        ASSERT_EQ(ippl::Comm->size(), 2);
        for (double dt : {.5, .98 / std::sqrt(3.)}) {
            const auto repaired = mixedRun<MixedStandard>(dt, false, MPI_COMM_WORLD);
            const auto control =
                mixedRun<ShippedOrderSolver<MixedStandard>>(dt, false, MPI_COMM_WORLD);
            EXPECT_LE(repaired.first, 1e-13);
            EXPECT_GE(control.first, 1e-3);
            std::cout << "mixed pulse dt=" << dt << " repaired=" << repaired.first
                      << " control=" << control.first << " bitwise=" << repaired.second << '\n';
        }
    }
}  // namespace fdtd_test
