/** @file test_mur_stability.cpp
 * @brief Check long-run potential envelopes from fixed inner-cell random histories.
 */
#include "OpenRunSupport.h"
#include "RandomPotentialHistoryFixture.h"
#include "RankWitness.h"

namespace chdr::solver::test {
    namespace {
        struct GrowthWindow {
            double first_m = 0, last_m = 0, maximum_m = 0;
        };

        /** @brief Observe owned potential maxima after steps 1 through count.
         * @details Compare the maximum in the final tenth with the first tenth. This is a
         * bounded-growth diagnostic, not a discrete energy or a proof of stability. The passing
         * envelope suite zeros two layers at open faces; the separate all-cell record retains
         * unrestricted random boundary values. This observer accepts either initialization.
         */
        GrowthWindow evolveWindow(OpenStandardSolver& solver, int count) {
            GrowthWindow window;
            for (int step = 1; step <= count; ++step) {
                advancePotentials(solver);
                const double maximum = potentialMaximum(solver.A_n);
                window.maximum_m     = std::max(window.maximum_m, maximum);
                if (step <= count / 10)
                    window.first_m = std::max(window.first_m, maximum);
                if (step > count - count / 10)
                    window.last_m = std::max(window.last_m, maximum);
            }
            return window;
        }

        /// Replay the four distributed fixtures on a rank-zero self communicator with identical
        /// data.
        void checkStabilityRanks(OpenStandardSolver& distributed, double epsilon, const char* name,
                                 std::array<int, 3> size, const GrowthWindow& observed) {
            const auto current  = gatherOwned(distributed.A_n),
                       previous = gatherOwned(distributed.A_nm1);
            RankAgreement report{};
            if (ippl::Comm->rank() == 0) {
                SerialFields fields(distributed.A_n);
                auto serial = makeOpenTestSolver(fields.source_m, fields.electric_m,
                                                 fields.magnetic_m, epsilon);
                serial->setTimeStepUnchecked(distributed.getDt());
                loadRandomPotentialHistoryFixture(*serial, name, size, true);
                const auto reference = evolveWindow(*serial, 1000);
                const auto a = gatherOwned(serial->A_n), b = gatherOwned(serial->A_nm1);
                const auto now    = rankDifference(current, a, reference.maximum_m);
                const auto before = rankDifference(previous, b, reference.maximum_m);
                report            = {std::max(now[0], before[0]), now[1] * before[1]};
                EXPECT_LE(std::abs(observed.first_m - reference.first_m),
                          1.e-12 * reference.maximum_m);
                EXPECT_LE(std::abs(observed.last_m - reference.last_m),
                          1.e-12 * reference.maximum_m);
            }
            checkRankAgreement(report, MPI_COMM_WORLD, "stability");
        }

        /// Fixed box, physical timestep, exact NumPy histories, matched speed unless deliberately
        /// overridden.
        void stabilityCase(bool periodicY, std::array<double, 3> h, double epsilon, bool nearCfl,
                           bool mismatch = false) {
            const std::array<int, 3> size{16, periodicY ? 8 : 16, 16};
            const char* name = periodicY ? "stability-mixed" : "stability-open";
            auto fields = makeBoxFieldContainer(size, {16 * h[0], size[1] * h[1], 16 * h[2]}, 2,
                                                {false, periodicY, false});
            initializeFieldStorage(fields);
            auto solver = makeOpenTestSolver(fields.getJ(), fields.getE(), fields.getB(), epsilon);
            if (nearCfl)
                solver->setSafetyFactor(.98);
            if (mismatch)
                solver->setKernelMurSpeedForTest(std::sqrt(epsilon));
            loadRandomPotentialHistoryFixture(*solver, name, size, true);
            const auto window = evolveWindow(*solver, ippl::Comm->size() == 1 ? 10000 : 1000);
            ASSERT_GT(window.first_m, 0.);
            ASSERT_TRUE(std::isfinite(window.maximum_m));
            ASSERT_TRUE(std::isfinite(window.last_m));
            if (mismatch)
                EXPECT_GE(window.last_m / window.first_m, 1.e6);
            else
                // A factor of two detects growth while allowing the initially excited envelope.
                EXPECT_LE(window.last_m / window.first_m, 2.);
            if (ippl::Comm->rank() == 0)
                std::cout << "stability periodic-y=" << periodicY << " hx=" << h[0]
                          << " epsilon=" << epsilon << " dt=" << solver->getDt()
                          << " mismatch=" << mismatch << " ratio=" << window.last_m / window.first_m
                          << '\n';
            if (ippl::Comm->size() == 2)
                checkStabilityRanks(*solver, epsilon, name, size, window);
        }
    }  // namespace

    /** Sixteen serial fixtures cover two boundary patterns, two aspects, two media and two
     * physical steps. The MPI subset uses cubic cells and the near-CFL step: four fixtures,
     * each split along z.
     * Exact NumPy data use seed 81413, independent histories and component-major global order.
     */
    TEST(MurStability, InnerDataEnvelope) {
        for (bool periodicY : {false, true})
            for (const auto h : {std::array{1., 1., 1.}, std::array{1.25, .8, 1.}})
                for (double epsilon : {1., 2.13})
                    for (bool nearCfl : {false, true}) {
                        if (ippl::Comm->size() == 2 && (h[0] != 1. || !nearCfl))
                            continue;
                        SCOPED_TRACE(::testing::Message() << periodicY << ' ' << h[0] << ' '
                                                          << epsilon << ' ' << nearCfl);
                        stabilityCase(periodicY, h, epsilon, nearCfl);
                    }
    }

    /** Physical Mur speed one in the dielectric gives kernel speed \f$\sqrt{\varepsilon_r}\f$.
     * This violates the supported boundary-speed interval. The 1e6 floor is far below this
     * fixture's measured growth of 1.9e22 and separates it from the bounded matched runs.
     */
    TEST(MurStability, MismatchedSpeedControl) {
        if (ippl::Comm->size() != 1)
            GTEST_SKIP() << "The growing control is a serial fixture";
        stabilityCase(false, {1., 1., 1.}, 2.13, true, true);
    }
    /** @brief Positive vacuum potential diagnostic, summed over all four components.
     * @details Define \f$Q=\Delta V\sum_q[\sum_i(\delta_t\psi_q)^2+
     * \sum_{d,\mathrm{links}}(D_d^+\bar\psi_q)^2]\f$ with
     * \f$\bar\psi=(\psi^k+\psi^{k-1})/2\f$. Links include periodic closing links and omit
     * links beyond physical faces. This chosen diagnostic uses the difference operators of
     * Kreiss, Petersson and Ystrom (2002), Eq. (1.4); it is not their conserved energy.
     */
    double positivePotentialNorm(OpenStandardSolver& solver) {
        const auto a = solver.A_n.getView(), b = solver.A_nm1.getView();
        const auto& layout = solver.A_n.getLayout();
        const auto n       = layout.getDomain().length();
        const auto first   = layout.getLocalNDIndex().first();
        const auto h       = solver.A_n.get_mesh().getMeshSpacing();
        const double dt    = solver.getDt();
        const Kokkos::Array<bool, 3> periodic{layout.isPeriodic(0), layout.isPeriodic(1),
                                              layout.isPeriodic(2)};
        double result = 0;
        Kokkos::parallel_reduce(
            "Positive potential diagnostic", solver.A_n.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k, double& sum) {
                const int cell[3]{i, j, k};
                for (unsigned q = 0; q < 4; ++q) {
                    const double time = (a(i, j, k)[q] - b(i, j, k)[q]) / dt;
                    sum += time * time;
                    for (unsigned d = 0; d < 3; ++d) {
                        if (!periodic[d] && cell[d] - 1 + first[d] == int(n[d]) - 1)
                            continue;
                        int next[3]{i, j, k};
                        ++next[d];
                        const double gradient =
                            (a(next[0], next[1], next[2])[q] + b(next[0], next[1], next[2])[q]
                             - a(i, j, k)[q] - b(i, j, k)[q])
                            / (2 * h[d]);
                        sum += gradient * gradient;
                    }
                }
            },
            result);
        MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, layout.comm);
        return result * h[0] * h[1] * h[2];
    }

    /// Unrestricted initial boundary data: report affine drift and the independently defined Q
    /// ratio.
    TEST(MurStabilityRecords, AllCellDrift) {
        if (ippl::Comm->size() != 1)
            GTEST_SKIP() << "This record uses one rank";
        auto fields = makeBoxFieldContainer({16, 16, 16}, {16, 16, 16}, 2, {false, false, false});
        initializeFieldStorage(fields);
        auto solver = makeOpenTestSolver(fields.getJ(), fields.getE(), fields.getB());
        solver->setSafetyFactor(.98);
        loadRandomPotentialHistoryFixture(*solver, "stability-open", {16, 16, 16}, false);
        const double initial = positivePotentialNorm(*solver);
        const auto window    = evolveWindow(*solver, 10000);
        const auto current = gatherOwned(solver->A_n), previous = gatherOwned(solver->A_nm1);
        std::array<double, 4> drift{};
        for (std::size_t i = 0; i < current.size(); ++i)
            drift[i % 4] += (current[i] - previous[i]) / 4096;
        std::cout << "all-cell-drift steps=10000 seed=81413 envelope_ratio="
                  << window.last_m / window.first_m
                  << " Q_N/Q_0=" << positivePotentialNorm(*solver) / initial
                  << " mean_step_increment=" << drift[0] << ',' << drift[1] << ',' << drift[2]
                  << ',' << drift[3] << '\n';
    }
}  // namespace chdr::solver::test
