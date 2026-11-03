/** @file test_monitor_records.cpp
 * @brief Gauss formulas and retained-source safety on random data with Nyquist content.
 * @details GaussIdentity explains the discrete formula derived from
 * [christlieb2024GaugeConserving, Theorem 3.5]. Distinct histories exercise each
 * reconstruction's scalar slot and the delayed R2Prime electric output.
 *
 * Configuration: periodic boxes of unit length per axis, Standard with 12 by 10 by 8 cells
 * and NonStandard with 8 by 8 by 12 cells, for R1, R2 and R2Prime (collocated reading).
 * The data add a hash of the global cell index, the checkerboard \f$(-1)^{x+y+z}\f$ and
 * the x-Nyquist mode \f$(-1)^x\f$, so the input bits do not depend on the rank count.
 * Over five monitor steps, each with a new source, the normalized mismatch between the Gauss
 * defect and its discrete formula (cell mismatch divided by its uncancelled absolute scale,
 * maximum over cells) is at most \f$10^{-12}\f$. After the live source is overwritten with
 * 12345.0, a second evaluation returns the same normalized formula mismatch and the same L2
 * norm of the Gauss defect (EXPECT_DOUBLE_EQ, 4 units in the last place). Reference: exact
 * lattice value (the
 * discrete formula at round-off).
 * The records are the in-memory reports of MonitorRecords.h; the test writes no file.
 *
 * CTest: chdr.diagnostics.monitor-records.serial and chdr.diagnostics.monitor-records.mpi2.
 */
#include <gtest/gtest.h>
#include <iostream>

#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "core/solver/ChdrSolver.h"
#include "opt-in/diagnostics/Monitors.h"

static_assert(
    chdr::SpatialDim == 3,
    "test_monitor_records.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr;
    using namespace chdr::solver;
    using namespace chdr::diagnostics;

    void fillRecordData(AJFourFieldsType& field, double phase) {
        const auto view  = field.getView();
        const auto local = field.getLayout().getLocalNDIndex();
        const int firstX = local[0].first() - 1, firstY = local[1].first() - 1,
                  firstZ = local[2].first() - 1;
        Kokkos::parallel_for(
            "Gauss random Nyquist fixture", field.getFieldRangePolicy(),
            KOKKOS_LAMBDA(int i, int j, int k) {
                const int x = i + firstX, y = j + firstY, z = k + firstZ;
                for (unsigned c = 0; c < 4; ++c) {
                    const unsigned hash = (unsigned(x) * 1103515245u + unsigned(y) * 12345u
                                           + unsigned(z) * 2654435761u + c * 374761393u)
                                          ^ 0x85ebca6bu;
                    view(i, j, k)[c] = double(hash % 1021) / 1021.0 + phase * (c + 1)
                                       + ((x + y + z) & 1 ? -0.2 : 0.2) + (x & 1 ? -0.3 : 0.3);
                }
            });
        field.fillHalo();
    }

    void checkFormulaReport(const MonitorReport& report, double& maximum) {
        EXPECT_LE(report.gauss_m.cell_m, 1e-12);
        maximum = std::max(maximum, report.gauss_m.cell_m);
    }

    /// The second evaluation must use the source consumed by the step, despite this overwrite.
    template <class SolverType>
    void checkLateSource(SolverType& solver, Monitors& monitors, const MonitorReport& original) {
        solver.sourceField() = 12345.0;
        const auto again     = monitors.evaluate(solver);
        EXPECT_DOUBLE_EQ(original.gauss_m.cell_m, again.gauss_m.cell_m);
        EXPECT_DOUBLE_EQ(original.gaussDefect_m.l2_m, again.gaussDefect_m.l2_m);
    }

    template <class PotentialSolverBase>
    void checkGaussReconstruction(std::array<int, SpatialDim> cells,
                                  EBReconstructionConvention ebReconstructionConvention) {
        auto fields = makeBoxFieldContainer(cells, {1.0, 1.0, 1.0}, 0);
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               ebReconstructionConvention);
        fillRecordData(solver.A_n, 0.13);
        fillRecordData(solver.A_nm1, 0.31);
        Monitors monitors(fields.getJ());
        double maximum = 0.0;
        for (unsigned step = 0; step < 5; ++step) {
            fillRecordData(fields.getJ(), 0.17 * step);
            const auto report = monitors.advance(solver);
            checkFormulaReport(report, maximum);
            checkLateSource(solver, monitors, report);
        }
        // The printed line reports the Gauss-defect versus discrete-formula comparison.
        if (ippl::Comm->rank() == 0)
            std::cout << "gauss stencil=" << int(solver.stencil())
                      << " reconstruction=" << int(ebReconstructionConvention)
                      << " Gauss=" << maximum << " threshold=1e-12\n";
    }

    template <class PotentialSolverBase>
    void checkGauss(std::array<int, SpatialDim> cells) {
        for (auto r : {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
                       EBReconstructionConvention::R2Prime})
            checkGaussReconstruction<PotentialSolverBase>(cells, r);
    }
    /// Standard stencil, 12 by 10 by 8 cells: Gauss formula check of R1, R2 and R2Prime.
    TEST(ChdrMonitorRecords, GaussStandard) {
        checkGauss<StandardSolverPeriodic>({12, 10, 8});
    }
    /// NonStandard stencil, 8 by 8 by 12 cells: the same checks.
    TEST(ChdrMonitorRecords, GaussNonStandard) {
        checkGauss<NonStandardSolverPeriodic>({8, 8, 12});
    }

}  // namespace
