/** @file test_mur_families.cpp
 * @brief Check local plane-wave compatibility of Mur faces, edges and corners.
 */
#include <cstring>

#include "OpenRunSupport.h"
#include "RankWitness.h"
#include "ReflectionOracles.h"
#include "tests/src/MaxwellSolvers/OuterLayerTest.h"

namespace chdr::solver::test {
    namespace {
        /** @brief Four or eight exact lattice waves near one selected edge or corner.
         * @details Each wave amplitude and the frequency come from the independent oracle.
         * Sum their phases in inward layer coordinates, with the outermost owned sample at
         * zero. The local family satisfies the selected faces and their intersections by
         * linear superposition [Bamberger, Joly and Roberts (1987), section II, p. 05], as
         * explained in familyAmplitude(). Opposite faces are excluded from the causal region.
         */
        struct Family {
            std::array<double, 3> h_m, k_m;
            std::array<bool, 3> open_m;
            std::array<OracleComplex, 8> amplitude_m{};
            double dt_m, omega_m;
            unsigned corner_m;
            Family(std::array<double, 3> spacing, double step, double epsilon, unsigned selected)
                : h_m(spacing)
                , k_m{.6 * std::numbers::pi / 4, .48 * std::numbers::pi / 4,
                      .64 * std::numbers::pi / 4}
                , open_m{true, selected != 8, true}
                , dt_m(step)
                , corner_m(selected) {
                if (!open_m[1]) {
                    k_m[1] = 2 * std::numbers::pi / (16 * h_m[1]);
                    const double transverse =
                        std::sqrt(std::pow(std::numbers::pi / 4, 2) - k_m[1] * k_m[1]);
                    k_m[0] = .8 * transverse;
                    k_m[2] = .6 * transverse;
                }
                const double speed = 1 / std::sqrt(epsilon);
                omega_m            = latticeFrequency(h_m, dt_m, speed, k_m);
                for (unsigned mask = 0; mask < 8; ++mask)
                    if (open_m[1] || !(mask & 2))
                        amplitude_m[mask] = familyAmplitude(h_m, dt_m, speed, k_m, open_m, mask);
            }

            /// Convert global cell labels to distance from the selected physical sample plane.
            int inward(int index, unsigned axis) const {
                return corner_m < 8 && ((corner_m >> axis) & 1u) ? 15 - index : index;
            }

            /// Fixed three-step comparison region: inward layer indices zero through eight.
            bool causal(std::array<int, 3> index) const {
                for (unsigned d = 0; d < 3; ++d)
                    if (open_m[d] && inward(index[d], d) >= 9)
                        return false;
                return true;
            }

            /// Real part of the independently weighted family at potential level -1, 0 or 3.
            double value(std::array<int, 3> index, int level) const {
                OracleComplex spatial = 0;
                for (unsigned mask = 0; mask < 8; ++mask) {
                    if (!open_m[1] && (mask & 2))
                        continue;
                    double phase = 0;
                    for (unsigned d = 0; d < 3; ++d)
                        phase += (open_m[d] && !(mask & (1u << d)) ? -1 : 1) * k_m[d] * h_m[d]
                                 * inward(index[d], d);
                    spatial += amplitude_m[mask] * std::polar(1., phase);
                }
                return (spatial * std::polar(1., -omega_m * level * dt_m)).real();
            }
        };

        /// Load identical exact scalar families in all components, then prepare both histories.
        void loadFamily(OpenStandardSolver& solver, const Family& family) {
            int level = -1;
            for (auto* field : {&solver.A_nm1, &solver.A_n}) {
                auto host = coreTest::copyToHost(*field);
                outer_layer_test::visit(*field, false, [&](int i, int j, int k, auto index) {
                    host(i, j, k) = family.value(index, level);
                });
                Kokkos::deep_copy(field->getView(), host);
                ++level;
            }
            solver.prepareHistories();
        }

        /** @brief Compare globally ordered samples only in the fixed three-step causal region.
         * @details The positive maximum of the exact real family at level three normalizes
         * every component. A separate serial witness isolates decomposition error; the
         * oracle comparison supplies the independent test of the boundary equations.
         */
        double familyError(const std::vector<double>& actual, const Family& family) {
            double error = 0, scale = 0;
            for (int i = 0; i < 16; ++i)
                for (int j = 0; j < 16; ++j)
                    for (int k = 0; k < 16; ++k) {
                        if (!family.causal({i, j, k}))
                            continue;
                        const double expected = family.value({i, j, k}, 3);
                        scale                 = std::max(scale, std::abs(expected));
                        for (unsigned q = 0; q < 4; ++q) {
                            const unsigned index    = 4 * ((i * 16 + j) * 16 + k) + q;
                            const double difference = std::abs(actual[index] - expected);
                            error                   = std::max(error, std::isfinite(difference)
                                                                          ? difference
                                                                          : std::numeric_limits<double>::infinity());
                        }
                    }
            EXPECT_GT(scale, 0.);
            return error / scale;
        }

        /// Execute one family using the production finish, or the explicit stale-halo control.
        void familyCase(std::array<double, 3> h, double epsilon, bool nearCfl, unsigned corner,
                        bool shipped = false) {
            auto fields = makeBoxFieldContainer({16, 16, 16}, {16 * h[0], 16 * h[1], 16 * h[2]}, 1,
                                                {false, corner == 8, false});
            initializeFieldStorage(fields);
            auto solver = makeOpenTestSolver(fields.getJ(), fields.getE(), fields.getB(), epsilon);
            if (nearCfl)
                solver->setSafetyFactor(.98);
            solver->useShippedOrder(shipped);
            const Family family(h, solver->getDt(), epsilon, corner);
            loadFamily(*solver, family);
            for (int step = 0; step < 3; ++step)
                advancePotentials(*solver);
            const auto actual = gatherOwned(solver->A_n);
            RankAgreement report{};
            if (ippl::Comm->rank() == 0) {
                const double error = familyError(actual, family);
                if (shipped && (corner == 8 || ippl::Comm->size() == 2))
                    EXPECT_GE(error, 1.e-3);
                else
                    EXPECT_LE(error, 1.e-13);
                if (!shipped && ippl::Comm->size() == 2) {
                    SerialFields serialFields(solver->A_n);
                    auto serial = makeOpenTestSolver(serialFields.source_m, serialFields.electric_m,
                                                     serialFields.magnetic_m, epsilon);
                    serial->setTimeStepUnchecked(solver->getDt());
                    loadFamily(*serial, family);
                    for (int step = 0; step < 3; ++step)
                        advancePotentials(*serial);
                    const auto reference = gatherOwned(serial->A_n);
                    report               = rankDifference(actual, reference);
                    std::cout << "family corner=" << corner << " hx=" << h[0]
                              << " epsilon=" << epsilon << " near-cfl=" << nearCfl << '\n';
                }
            }
            if (!shipped && ippl::Comm->size() == 2)
                checkRankAgreement(report, MPI_COMM_WORLD, "family");
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }  // namespace

    /// Two aspects, two media, two physical steps and eight corners plus one mixed edge: 72
    /// fixtures.
    TEST(MurEdgeCornerFamilies, LocalCompatibility) {
        for (const auto h : {std::array{1., 1., 1.}, std::array{1.25, .8, 1.}})
            for (double epsilon : {1., 2.13})
                for (bool nearCfl : {false, true})
                    for (unsigned corner = 0; corner < 9; ++corner) {
                        SCOPED_TRACE(::testing::Message()
                                     << h[0] << ' ' << epsilon << ' ' << nearCfl << ' ' << corner);
                        familyCase(h, epsilon, nearCfl, corner);
                    }
    }

    /// Serial all-open fields have no stale wrap; only distributed all-open and mixed cases must
    /// fail.
    TEST(MurEdgeCornerFamilies, ShippedOrderControls) {
        familyCase({1., 1., 1.}, 1., false, 0, true);
        familyCase({1., 1., 1.}, 1., false, 8, true);
    }
}  // namespace chdr::solver::test
