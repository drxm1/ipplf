#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <gtest/gtest.h>

#include "HaloOrderTest.h"
#include "ShippedOrderSolver.h"

namespace fdtd_test {
    /** @brief Assemble owned samples in canonical global order on every rank.
     * @details Exclude the six outer layers when requested. Bitwise OR on the
     * single owner's binary64 representation preserves signed zeros and avoids
     * floating-point reduction arithmetic. No ghost enters the comparison.
     */
    template <class F>
    std::vector<double> samples(F& field, bool interior) {
        constexpr unsigned components = F::value_type::dim;
        const auto n                  = field.getLayout().getDomain().length();
        const unsigned o              = interior ? 1 : 0;
        std::vector<std::uint64_t> bits((n[0] - 2 * o) * (n[1] - 2 * o) * (n[2] - 2 * o)
                                        * components);
        owned(field, [&](std::size_t i, std::size_t j, std::size_t k, const auto& value) {
            if (interior
                && (i == 0 || j == 0 || k == 0 || i == n[0] - 1 || j == n[1] - 1 || k == n[2] - 1))
                return;
            const std::size_t index = ((i - o) * (n[1] - 2 * o) + j - o) * (n[2] - 2 * o) + k - o;
            for (unsigned q = 0; q < components; ++q)
                bits[index * components + q] = std::bit_cast<std::uint64_t>(value[q]);
        });
        if (field.getLayout().comm.size() > 1)
            MPI_Allreduce(MPI_IN_PLACE, bits.data(), bits.size(), MPI_UINT64_T, MPI_BOR,
                          field.getLayout().comm);
        std::vector<double> result(bits.size());
        std::transform(bits.begin(), bits.end(), result.begin(), [](std::uint64_t value) {
            return std::bit_cast<double>(value);
        });
        return result;
    }

    /** @brief Make deterministic histories without an external data generator.
     * @details Unsigned 32-bit arithmetic supplies a fixed sequence; scaling by
     * a power of two gives exactly representable binary64 samples in \f$[-1,1)\f$.
     * Both solvers receive the same histories in global cell/component order.
     */
    std::vector<double> portableHistories(std::size_t cells) {
        std::vector<double> result(8 * cells);
        std::uint32_t state = 271828;
        for (auto& value : result) {
            state = 1664525u * state + 1013904223u;
            value = (double(state) - 2147483648.) / 2147483648.;
        }
        return result;
    }

    /** @brief Load both independent histories and their initial ghost copies. */
    template <class Solver>
    void loadHistories(Solver& solver, const std::vector<double>& data, std::array<int, 3> n) {
        const std::size_t cells = n[0] * n[1] * n[2];
        int level               = 0;
        for (auto* field : {&solver.A_nm1, &solver.A_n}) {
            load(*field, [&](int i, int j, int k, unsigned q) {
                return data[(4 * level + q) * cells + (i * n[1] + j) * n[2] + k];
            });
            ++level;
        }
        solver.evaluate_EB();
    }

    /** @brief Assert exact in-build equality and finite samples.
     * @details memcmp also distinguishes signed zero. Finiteness is checked
     * separately, so a pair of identical invalid samples cannot pass.
     */
    void compare(const std::vector<double>& actual, const std::vector<double>& control) {
        ASSERT_EQ(actual.size(), control.size());
        EXPECT_EQ(std::memcmp(actual.data(), control.data(), actual.size() * sizeof(double)), 0);
        EXPECT_TRUE(std::all_of(actual.begin(), actual.end(), [](double x) {
            return std::isfinite(x);
        }));
    }

    /** @brief Check all histories and off-boundary electric and magnetic fields.
     * @details Raw B is the output of evaluate_EB(); averaged B is the arithmetic
     * mean of its current and previous curls. This test asserts preservation of
     * those samples, not the accuracy of either reconstruction. The reference
     * solver uses identical kernels with the original finishing order.
     */
    template <class Solver>
    void regression(std::array<int, 3> n, ippl::Vector<double, 3> h, bool periodic,
                    const std::vector<double>& data) {
        const MPI_Comm comm = periodic ? MPI_COMM_WORLD : MPI_COMM_SELF;
        Fields fields(n, h, periodic, comm), reference(n, h, periodic, comm);
        Solver solver(fields.source, fields.electric, fields.magnetic);
        ShippedOrderSolver<Solver> control(reference.source, reference.electric,
                                           reference.magnetic);
        loadHistories(solver, data, n);
        loadHistories(control, data, n);
        for (int step = 0; step < 256; ++step) {
            auto oldB        = samples(fields.magnetic, !periodic);
            auto oldControlB = samples(reference.magnetic, !periodic);
            solver.step();
            control.step();
            compare(samples(solver.A_nm1, false), samples(control.A_nm1, false));
            compare(samples(solver.A_n, false), samples(control.A_n, false));
            compare(samples(solver.A_np1, false), samples(control.A_np1, false));
            solver.timeShift();
            control.timeShift();
            solver.evaluate_EB();
            control.evaluate_EB();
            compare(samples(fields.electric, !periodic), samples(reference.electric, !periodic));
            auto rawB     = samples(fields.magnetic, !periodic);
            auto controlB = samples(reference.magnetic, !periodic);
            compare(rawB, controlB);
            for (std::size_t i = 0; i < rawB.size(); ++i) {
                rawB[i]     = .5 * (rawB[i] + oldB[i]);
                controlB[i] = .5 * (controlB[i] + oldControlB[i]);
            }
            compare(rawB, controlB);
        }
    }

    /// Exercise both stencils and boundary paths with the same two input histories.
    void runRegression(const std::vector<double>& standard,
                       const std::vector<double>& nonstandard) {
        regression<Standard<ippl::periodic>>({24, 24, 24}, {1, 1, 1}, true, standard);
        regression<Standard<ippl::absorbing>>({24, 24, 24}, {1, 1, 1}, false, standard);
        regression<NonStandard<ippl::periodic>>({16, 16, 48}, {1, 1, .5}, true, nonstandard);
        regression<NonStandard<ippl::absorbing>>({16, 16, 48}, {1, 1, .5}, false, nonstandard);
    }

    TEST(OpenBoundaryRegression, OriginalOrder) {
        runRegression(portableHistories(24 * 24 * 24), portableHistories(16 * 16 * 48));
    }

}  // namespace fdtd_test

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
