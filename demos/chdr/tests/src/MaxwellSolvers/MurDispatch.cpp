#include <gtest/gtest.h>

#include "BoundaryOracles.h"
#include "HaloOrderTest.h"

namespace {
    /** @brief Enumerate the primitive associated with each paired face mask.
     * This table does not inspect coordinates. The independent cell oracle chooses
     * its row, while MurSpeed checks the primitive's weights separately.
     */
    template <class Visit>
    void faces(Visit visit, ippl::Vector<double, 3> h, double dt) {
        visit(1u, second_order_abc_face<double, 0, 1, 2>(h, dt, 1));
        visit(2u, second_order_abc_face<double, 0, 1, 2>(h, dt, -1));
        visit(4u, second_order_abc_face<double, 1, 0, 2>(h, dt, 1));
        visit(8u, second_order_abc_face<double, 1, 0, 2>(h, dt, -1));
        visit(16u, second_order_abc_face<double, 2, 0, 1>(h, dt, 1));
        visit(32u, second_order_abc_face<double, 2, 0, 1>(h, dt, -1));
    }

    /// Edge masks pair the two inward normals; the remaining axis runs along the edge.
    template <class Visit>
    void edges(Visit visit, ippl::Vector<double, 3> h, double dt) {
        visit(4u | 16u, second_order_abc_edge<double, 0, 1, 2, true, true>(h, dt));
        visit(4u | 32u, second_order_abc_edge<double, 0, 1, 2, true, false>(h, dt));
        visit(8u | 16u, second_order_abc_edge<double, 0, 1, 2, false, true>(h, dt));
        visit(8u | 32u, second_order_abc_edge<double, 0, 1, 2, false, false>(h, dt));
        visit(16u | 1u, second_order_abc_edge<double, 1, 2, 0, true, true>(h, dt));
        visit(16u | 2u, second_order_abc_edge<double, 1, 2, 0, true, false>(h, dt));
        visit(32u | 1u, second_order_abc_edge<double, 1, 2, 0, false, true>(h, dt));
        visit(32u | 2u, second_order_abc_edge<double, 1, 2, 0, false, false>(h, dt));
        visit(1u | 4u, second_order_abc_edge<double, 2, 0, 1, true, true>(h, dt));
        visit(1u | 8u, second_order_abc_edge<double, 2, 0, 1, true, false>(h, dt));
        visit(2u | 4u, second_order_abc_edge<double, 2, 0, 1, false, true>(h, dt));
        visit(2u | 8u, second_order_abc_edge<double, 2, 0, 1, false, false>(h, dt));
    }

    /// Each corner row selects the inward sign of all three axes.
    template <class Visit>
    void corners(Visit visit, ippl::Vector<double, 3> h, double dt) {
        visit(1u | 4u | 16u, second_order_abc_corner<double, true, true, true>(h, dt));
        visit(2u | 4u | 16u, second_order_abc_corner<double, false, true, true>(h, dt));
        visit(1u | 8u | 16u, second_order_abc_corner<double, true, false, true>(h, dt));
        visit(2u | 8u | 16u, second_order_abc_corner<double, false, false, true>(h, dt));
        visit(1u | 4u | 32u, second_order_abc_corner<double, true, true, false>(h, dt));
        visit(2u | 4u | 32u, second_order_abc_corner<double, false, true, false>(h, dt));
        visit(1u | 8u | 32u, second_order_abc_corner<double, true, false, false>(h, dt));
        visit(2u | 8u | 32u, second_order_abc_corner<double, false, false, false>(h, dt));
    }

    struct Solver : fdtd_test::Standard<ippl::configured> {
        using Base = fdtd_test::Standard<ippl::configured>;
        using Base::applyBCs;
        using Base::Base;
    };

    /** @brief Check all owned samples against their independent paired face mask.
     * The negative control executes the dispatcher with all faces open. Mixed
     * layouts must expose both changed periodic samples and wrong open-face values.
     */
    void checkDispatch(unsigned pattern, bool raw = false) {
        const std::array<bool, 3> periodic{bool(pattern & 1), bool(pattern & 2), bool(pattern & 4)};
        fdtd_test::Fields fields({6, 5, 4}, {1, 1.25, .5}, periodic);
        Solver solver(fields.source, fields.electric, fields.magnetic,
                      ippl::FDTDBoundary::fromLayout(fields.layout));
        for (auto* field : {&solver.A_n, &solver.A_nm1, &solver.A_np1}) {
            fdtd_test::load(*field, [](int i, int j, int k, unsigned q) {
                return 1 + i * i + 3 * j * j + 5 * k * k + q * (i + j + k);
            });
        }
        solver.prepareHistories();
        if (raw) {
            second_order_mur_boundary_conditions{}.apply(solver.A_n, solver.A_nm1, solver.A_np1,
                                                         solver.getDt(), {8, 7, 6},
                                                         fields.layout.getLocalNDIndex(), 63);
        } else {
            solver.applyBCs();
        }
        const auto now =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), solver.A_n.getView());
        const auto old =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), solver.A_nm1.getView());
        const auto next =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), solver.A_np1.getView());
        const auto local = fields.layout.getLocalNDIndex();
        unsigned errors[2]{};
        fdtd_test::owned(solver.A_np1, [&](int i, int j, int k, auto actual) {
            const auto mask = open_boundary_oracle::murFaceMask(
                {i + 1, j + 1, k + 1}, {6, 5, 4}, {periodic[0], periodic[1], periodic[2]});
            const ippl::Vector<uint32_t, 3> p{unsigned(i + 1 - local[0].first()),
                                              unsigned(j + 1 - local[1].first()),
                                              unsigned(k + 1 - local[2].first())};
            if (mask == 0) {
                for (unsigned q = 0; q < 4; ++q) {
                    errors[0] += actual[q] != 1 + i * i + 3 * j * j + 5 * k * k + q * (i + j + k);
                }
            }
            auto check = [&](unsigned required, auto operation) {
                if (mask == required) {
                    const auto expected = operation(now, old, next, p);
                    for (unsigned q = 0; q < 4; ++q) {
                        errors[1] += actual[q] != expected[q];
                    }
                }
            };
            faces(check, fields.mesh.getMeshSpacing(), solver.getDt());
            edges(check, fields.mesh.getMeshSpacing(), solver.getDt());
            corners(check, fields.mesh.getMeshSpacing(), solver.getDt());
        });
        unsigned globalErrors[2]{};
        MPI_Allreduce(errors, globalErrors, 2, MPI_UNSIGNED, MPI_SUM,
                      fields.layout.comm.getCommunicator());
        EXPECT_EQ(globalErrors[0] > 0, raw);
        EXPECT_EQ(globalErrors[1] > 0, raw);
        if (raw && fields.layout.comm.rank() == 0) {
            std::cout << "raw dispatch pattern=" << pattern
                      << " changed periodic components=" << globalErrors[0]
                      << " wrong open components=" << globalErrors[1] << '\n';
        }
    }
}  // namespace

TEST(MurDispatch, EveryOpenFaceEdgeAndCorner) {
    for (unsigned pattern = 0; pattern < 8; ++pattern) {
        SCOPED_TRACE(pattern);
        checkDispatch(pattern);
        if (pattern != 0 && pattern != 7) {
            checkDispatch(pattern, true);
        }
    }
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
