/** @file MurDispatch.cpp
 * @brief Dispatch of the Mur face, edge and corner updates to the owned boundary cells of a
 * StandardFDTDSolver with a configured FDTDBoundary.
 * @details The box has 6 by 5 by 4 cells of spacing (1, 1.25, 0.5) and splits along \f$z\f$.
 * For all 8 axis patterns, the test loads distinct potential histories with an integer
 * polynomial of the global cell index and calls applyBCs() once. Every owned cell with oracle
 * mask 0 (murFaceMask() in BoundaryOracles.h) must keep its loaded value. Every face, edge and
 * corner cell must equal exactly the output of the Mur struct that its mask selects. In the
 * negative control, the dispatcher treats all 6 faces as open; on the 6 mixed patterns it must
 * change periodic cells and write wrong values into open cells. MurSpeed.cpp checks the
 * weights of the structs. In one more run per pattern, setMurCellSpeeds() gives every cell its
 * own speed between 0.55 and 0.95, and each face, edge and corner row must use the speed of
 * the cell that it writes. CellSpeedRefusals checks that the setter refuses a view with other
 * extents and a zero speed.
 *
 * Reference: self-consistency (routing against an independent mask). CTest entries:
 * chdr.src.MurDispatch.serial and chdr.src.MurDispatch.mpi2.
 */
#include <gtest/gtest.h>
#include <optional>

#include "BoundaryOracles.h"
#include "HaloOrderTest.h"

namespace {
    /** @brief Enumerate the primitive associated with each paired face mask.
     * This table does not inspect coordinates. The independent cell oracle chooses
     * its row, while MurSpeed checks the primitive's weights separately.
     */
    template <class Visit>
    void faces(Visit visit, ippl::Vector<double, 3> h, double dt,
               const std::array<double, 2>& endSpeeds) {
        visit(1u, second_order_abc_face<double, 0, 1, 2>(h, dt, 1, endSpeeds[0]));
        visit(2u, second_order_abc_face<double, 0, 1, 2>(h, dt, -1, endSpeeds[1]));
        visit(4u, second_order_abc_face<double, 1, 0, 2>(h, dt, 1, endSpeeds[0]));
        visit(8u, second_order_abc_face<double, 1, 0, 2>(h, dt, -1, endSpeeds[1]));
        visit(16u, second_order_abc_face<double, 2, 0, 1>(h, dt, 1, endSpeeds[0]));
        visit(32u, second_order_abc_face<double, 2, 0, 1>(h, dt, -1, endSpeeds[1]));
    }

    /// Edge masks pair the two inward normals; the remaining axis runs along the edge.
    template <class Visit>
    void edges(Visit visit, ippl::Vector<double, 3> h, double dt, double speed) {
        visit(4u | 16u, second_order_abc_edge<double, 0, 1, 2, true, true>(h, dt, speed));
        visit(4u | 32u, second_order_abc_edge<double, 0, 1, 2, true, false>(h, dt, speed));
        visit(8u | 16u, second_order_abc_edge<double, 0, 1, 2, false, true>(h, dt, speed));
        visit(8u | 32u, second_order_abc_edge<double, 0, 1, 2, false, false>(h, dt, speed));
        visit(16u | 1u, second_order_abc_edge<double, 1, 2, 0, true, true>(h, dt, speed));
        visit(16u | 2u, second_order_abc_edge<double, 1, 2, 0, true, false>(h, dt, speed));
        visit(32u | 1u, second_order_abc_edge<double, 1, 2, 0, false, true>(h, dt, speed));
        visit(32u | 2u, second_order_abc_edge<double, 1, 2, 0, false, false>(h, dt, speed));
        visit(1u | 4u, second_order_abc_edge<double, 2, 0, 1, true, true>(h, dt, speed));
        visit(1u | 8u, second_order_abc_edge<double, 2, 0, 1, true, false>(h, dt, speed));
        visit(2u | 4u, second_order_abc_edge<double, 2, 0, 1, false, true>(h, dt, speed));
        visit(2u | 8u, second_order_abc_edge<double, 2, 0, 1, false, false>(h, dt, speed));
    }

    /// Each corner row selects the inward sign of all three axes.
    template <class Visit>
    void corners(Visit visit, ippl::Vector<double, 3> h, double dt, double speed) {
        visit(1u | 4u | 16u, second_order_abc_corner<double, true, true, true>(h, dt, speed));
        visit(2u | 4u | 16u, second_order_abc_corner<double, false, true, true>(h, dt, speed));
        visit(1u | 8u | 16u, second_order_abc_corner<double, true, false, true>(h, dt, speed));
        visit(2u | 8u | 16u, second_order_abc_corner<double, false, false, true>(h, dt, speed));
        visit(1u | 4u | 32u, second_order_abc_corner<double, true, true, false>(h, dt, speed));
        visit(2u | 4u | 32u, second_order_abc_corner<double, false, true, false>(h, dt, speed));
        visit(1u | 8u | 32u, second_order_abc_corner<double, true, false, false>(h, dt, speed));
        visit(2u | 8u | 32u, second_order_abc_corner<double, false, false, false>(h, dt, speed));
    }

    struct Solver : fdtd_test::Standard<ippl::configured> {
        using Base = fdtd_test::Standard<ippl::configured>;
        using Base::applyBCs;
        using Base::Base;
    };

    /// Distinguish spatial index, time level and component in each owned sample.
    double label(int i, int j, int k, unsigned component, int level) {
        return 1 + i * i + 3 * j * j + 5 * k * k + component * (i + j + k + 1)
               + level * (16.0 + i + 2 * j + k + component);
    }

    /// Distinct admissible Mur speeds below one, indexed by the global owned cell.
    double cellSpeed(int i, int j, int k) {
        return .55 + .04 * ((i + 2 * j + 3 * k) % 11);
    }

    /// Load cellSpeed() into a view with the extents of the local A_n view.
    Solver::MurCellSpeedView::non_const_type cellSpeedView(const fdtd_test::Potential& field) {
        const auto first = field.getLayout().getLocalNDIndex().first();
        Solver::MurCellSpeedView::non_const_type speeds(
            "Mur cell speeds", field.getView().extent(0), field.getView().extent(1),
            field.getView().extent(2));
        auto host = Kokkos::create_mirror_view(speeds);
        for (int a = 0; a < int(host.extent(0)); ++a)
            for (int b = 0; b < int(host.extent(1)); ++b)
                for (int c = 0; c < int(host.extent(2)); ++c)
                    host(a, b, c) = cellSpeed(a - 1 + first[0], b - 1 + first[1], c - 1 + first[2]);
        Kokkos::deep_copy(speeds, host);
        return speeds;
    }

    void loadLevels(Solver& solver) {
        int level = 0;
        for (auto* field : {&solver.A_nm1, &solver.A_n, &solver.A_np1}) {
            fdtd_test::load(*field, [level](int i, int j, int k, unsigned q) {
                return label(i, j, k, q, level);
            });
            ++level;
        }
        solver.prepareHistories();
    }

    /** @brief Check all owned samples against their independent paired face mask.
     * The negative control executes the dispatcher with all faces open. Mixed
     * layouts must expose both changed periodic samples and wrong open-face values.
     * With perCell, every face, edge and corner row must use cellSpeed() of its own cell.
     */
    void checkDispatch(unsigned pattern, bool raw = false, std::array<double, 2> endSpeeds = {1, 1},
                       std::optional<unsigned> selected = std::nullopt, bool perCell = false) {
        const std::array<bool, 3> periodic{bool(pattern & 1), bool(pattern & 2), bool(pattern & 4)};
        fdtd_test::Fields fields({6, 5, 4}, {1, 1.25, .5}, periodic);
        auto boundary = ippl::FDTDBoundary::fromLayout(fields.layout);
        for (unsigned face = 0; face < 6; ++face)
            boundary.speed[face] = periodic[face / 2] ? -1 : endSpeeds[face % 2];
        Solver solver(fields.source, fields.electric, fields.magnetic, boundary);
        loadLevels(solver);
        if (perCell)
            solver.setMurCellSpeeds(cellSpeedView(solver.A_n));
        if (raw) {
            second_order_mur_boundary_conditions{}.apply(solver.A_n, solver.A_nm1, solver.A_np1,
                                                         solver.getDt(), {8, 7, 6},
                                                         fields.layout.getLocalNDIndex(), 63);
        } else if (selected) {
            second_order_mur_boundary_conditions{}.apply(
                solver.A_n, solver.A_nm1, solver.A_np1, solver.getDt(), {8, 7, 6},
                fields.layout.getLocalNDIndex(), *selected, boundary.kernelSpeeds<double>());
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
            const auto mask =
                selected.value_or(63u)
                & open_boundary_oracle::murFaceMask({i + 1, j + 1, k + 1}, {6, 5, 4},
                                                    {periodic[0], periodic[1], periodic[2]});
            const ippl::Vector<uint32_t, 3> p{unsigned(i + 1 - local[0].first()),
                                              unsigned(j + 1 - local[1].first()),
                                              unsigned(k + 1 - local[2].first())};
            if (mask == 0) {
                for (unsigned q = 0; q < 4; ++q) {
                    errors[0] += actual[q] != label(i, j, k, q, 2);
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
            const double speed = perCell ? cellSpeed(i, j, k) : endSpeeds[0];
            faces(check, fields.mesh.getMeshSpacing(), solver.getDt(),
                  perCell ? std::array{speed, speed} : endSpeeds);
            edges(check, fields.mesh.getMeshSpacing(), solver.getDt(), speed);
            corners(check, fields.mesh.getMeshSpacing(), solver.getDt(), speed);
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

/// For all 8 axis patterns, interior cells keep their values and boundary cells equal the
/// output of their mask's struct; on the 6 mixed patterns the all-open control must fail.
TEST(MurDispatch, EveryOpenFaceEdgeAndCorner) {
    for (unsigned pattern = 0; pattern < 8; ++pattern) {
        SCOPED_TRACE(pattern);
        checkDispatch(pattern);
        if (pattern == 0 || pattern == 1 || pattern == 2 || pattern == 4)
            checkDispatch(pattern, false, {.8, .8});
        if (pattern == 3 || pattern == 5 || pattern == 6)
            checkDispatch(pattern, false, {.8, 1.2});
        checkDispatch(pattern, false, {1, 1}, std::nullopt, true);
        if (pattern != 0 && pattern != 7) {
            checkDispatch(pattern, true);
        }
    }
}

/// Check direct calls with paired, one-sided, upper-only and empty selections.
TEST(MurDispatch, DirectSelectedFaces) {
    for (unsigned selected : {63u, 1u | 8u | 32u, 32u, 0u}) {
        SCOPED_TRACE(selected);
        checkDispatch(0, false, {.8, .8}, selected);
    }
}

/// setMurCellSpeeds refuses a view with other extents and a zero owned speed on all ranks.
TEST(MurDispatch, CellSpeedRefusals) {
    fdtd_test::Fields fields({6, 5, 4}, {1, 1.25, .5}, {false, false, false});
    Solver solver(fields.source, fields.electric, fields.magnetic,
                  ippl::FDTDBoundary::fromLayout(fields.layout));
    EXPECT_THROW(solver.setMurCellSpeeds(Solver::MurCellSpeedView::non_const_type("s", 2, 2, 2)),
                 IpplException);
    auto speeds = cellSpeedView(solver.A_n);
    Kokkos::deep_copy(speeds, 0.0);
    EXPECT_THROW(solver.setMurCellSpeeds(speeds), IpplException);
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
