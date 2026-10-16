#include "Ippl.h"

#include "MaxwellSolvers/NonStandardFDTDSolver.h"
#include "OuterLayerTest.h"
#include "gtest/gtest.h"

using namespace outer_layer_test;

namespace {
    /// Reconstruct the same loaded histories on rank zero's self communicator and compare fields.
    void checkPolynomialRanks(const Box& box, bool periodicY) {
        if (box.layout.comm.size() == 1)
            return;
        const auto electric = ownedFieldBits(box.electric), magnetic = ownedFieldBits(box.magnetic);
        if (box.layout.comm.rank() == 0) {
            Box serialBox(box.size, spacing, periodicY, 2, MPI_COMM_SELF);
            Solver serial(serialBox.source, serialBox.electric, serialBox.magnetic);
            serial.usePolynomialTimeStep();
            load(serial.A_n, 0.);
            load(serial.A_nm1, -serial.getDt());
            serial.prepareHistories();
            serial.evaluate_EB();
            std::cout << "outer-layer R1 periodic-y=" << periodicY << '\n';
            expectSameField(electric, serialBox.electric, "E");
            expectSameField(magnetic, serialBox.magnetic, "B");
        }
        MPI_Barrier(box.layout.comm);
    }
}  // namespace

/** The forty-coefficient oracle exercises every component at all six faces, edges and corners.
 * Zero ghosts must produce a relative error of at least 0.1 in each field, while prepared
 * histories must agree within 1e-13 of each field's positive analytic global maximum.
 */
TEST(OpenFaceClosure, PolynomialFieldsAndZeroGhostControl) {
    for (bool periodicY : {false, true}) {
        Box box({8, 6, 10}, spacing, periodicY);
        Solver solver(box.source, box.electric, box.magnetic);
        solver.usePolynomialTimeStep();
        load(solver.A_n, 0.);
        load(solver.A_nm1, -solver.getDt());
        solver.A_n.fillHalo();
        solver.A_nm1.fillHalo();
        solver.evaluate_EB();
        EXPECT_GE(fieldError(box.electric, box, false), .1);
        EXPECT_GE(fieldError(box.magnetic, box, true), .1);
        solver.prepareHistories();
        solver.evaluate_EB();
        EXPECT_LE(fieldError(box.electric, box, false), 1.e-13);
        EXPECT_LE(fieldError(box.magnetic, box, true), 1.e-13);
        checkPolynomialRanks(box, periodicY);
    }
}

/** With two ranks, the third normal sample lies in an MPI halo. The tensor polynomial also
 * checks seam/face intersections and full ghost corners, at both loaded time levels.
 */
TEST(OpenFaceClosure, TwoOwnedNormalSamplesAndCorners) {
    for (bool periodicY : {false, true}) {
        Box box({4, 6, 10}, spacing, periodicY, 0);
        Solver solver(box.source, box.electric, box.magnetic);
        solver.usePolynomialTimeStep();
        load(solver.A_n, 0.);
        load(solver.A_nm1, -solver.getDt());
        solver.prepareHistories();
        EXPECT_LE(historyError(solver.A_n, 0.), 1.e-13);
        EXPECT_LE(historyError(solver.A_nm1, -solver.getDt()), 1.e-13);
    }
}

/** Both production step tails must close the newly computed level. A constant polynomial
 * makes the expected face ghost equal to its owned samples; omitted closure leaves zero.
 */
TEST(OpenFaceClosure, BothStepsCloseTheNewLevel) {
    Box box({8, 6, 10}, {1., .8, .5});
    ippl::StandardFDTDSolver<Fields, Potentials, ippl::absorbing> standard(box.source, box.electric,
                                                                           box.magnetic);
    ippl::NonStandardFDTDSolver<Fields, Potentials, ippl::absorbing> nonstandard(
        box.source, box.electric, box.magnetic);
    using Base = ippl::FDTDSolverBase<Fields, Potentials, ippl::absorbing>;
    for (auto* solver : std::array<Base*, 2>{&standard, &nonstandard}) {
        solver->A_n   = 2.;
        solver->A_nm1 = 2.;
        solver->step();
        auto host = solver->A_np1.getHostMirror();
        Kokkos::deep_copy(host, solver->A_np1.getView());
        for (unsigned q = 0; q < 4; ++q)
            EXPECT_NEAR(host(1, 0, 1)[q], 2., 1.e-13);
    }
}

/** A periodic template on a nonperiodic layout needs PeriodicFace conditions in addition to
 * fillHalo(). Check both loaded levels at every face centre, including MPI-exchanged faces.
 */
TEST(OpenFaceClosure, PreparePeriodicHistoriesOnNonperiodicLayout) {
    Box box;
    ippl::StandardFDTDSolver<Fields, Potentials, ippl::periodic> solver(box.source, box.electric,
                                                                        box.magnetic);
    load(solver.A_n, 0.);
    load(solver.A_nm1, -solver.getDt());
    solver.prepareHistories();
    const auto c = open_boundary_oracle::polynomialFixture(false);
    for (auto* field : {&solver.A_n, &solver.A_nm1}) {
        const double time = field == &solver.A_n ? 0. : -solver.getDt();
        auto host         = field->getHostMirror();
        Kokkos::deep_copy(host, field->getView());
        visit(*field, true, [&](int i, int j, int k, auto index) {
            const auto v     = field->getView();
            const int ghosts = (i == 0 || i == int(v.extent(0)) - 1)
                               + (j == 0 || j == int(v.extent(1)) - 1)
                               + (k == 0 || k == int(v.extent(2)) - 1);
            if (ghosts != 1)
                return;
            for (unsigned d = 0; d < 3; ++d)
                index[d] = (index[d] + box.size[d]) % box.size[d];
            const auto x = open_boundary_oracle::fieldPosition(index, spacing, false, false, 0);
            for (unsigned q = 0; q < 4; ++q)
                EXPECT_DOUBLE_EQ(host(i, j, k)[q], open_boundary_oracle::polynomialPotential(
                                                       c, q, x, time, -1, false));
        });
    }
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    ippl::finalize();
    return result;
}
