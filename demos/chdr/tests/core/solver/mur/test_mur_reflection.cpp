/** @file test_mur_reflection.cpp
 * @brief Exact two-wave continuation, sensitive controls and MPI witnesses for Mur faces.
 * @details Each row of OpenBoxCases.h starts from an exact incident-plus-reflected lattice
 * state (TwoWaveState.h) and advances 30 steps. Box: 40 normal cells, Mur faces on x,
 * periodic y and z, split along z on two ranks, R2 in both readings. The 13 rows are 6 vacuum
 * rows, 5 rows at \f$\varepsilon_r = 2.13\f$ with matched Mur speed and 2 vacuum rows with
 * h = (1.25, 0.8, 1). The step is 0.5; the two anisotropic rows use 0.5029188806044014.
 * The reference is the exact lattice value of the discrete reflection coefficient of Mur's
 * condition, discreteReflection() in ReflectionOracles.h. Pass: the absolute coefficient error
 * and the cell error on layers 0 to 9 (relative to the exact maximum on these layers) are each
 * at most 1e-11. Controls: a kernel Mur speed scaled by \f$1 + 10^{-6}\f$ and, on row v60, the
 * old halo order must exceed 1e-11 in at least one criterion. verifyRanks() replays each run
 * without a control on MPI_COMM_SELF: cells within 1e-13, projections within 1e-12. In the
 * serial entry this replay compares a run with an identical run.
 *
 * CTest: chdr.solver.mur.reflection.serial and chdr.solver.mur.reflection.mpi2
 * (MurReflection.ExactTwoWaveRows). chdr.solver.mur.anisotropic.mpi2 selects
 * MurReflection.AnisotropicMedium; see that test.
 */
#include <gtest/gtest.h>

#include "ReflectionRun.h"
#include "core/solver/BoxSetup.h"
#include "tests/core/solver/support/RankWitness.h"
#include "tests/core/solver/support/SerialReference.h"

namespace chdr::solver::test::reflection {
    /** @brief Compare every sample in the first ten layers with exact continuation.
     * @details With forty normal cells and thirty updates, the far face cannot
     * reach these samples [trefethen1996FiniteDifferenceSpectral, section 4.3,
     * p. 161]. The scale is the exact maximum on these same ten layers, not
     * the possibly much larger evanescent amplitude near the opposite face.
     */
    std::array<double, 2> potentialErrors(const Case& row, SpatialDiscretization reading,
                                          const Result& result) {
        const TwoWaveState real(row, reading, false, 0), imaginary(row, reading, false, 1);
        const auto offset = LevelContract{EBReconstructionConvention::R2, reading}.positionOffset(
            Quantity::VectorPotential, 1);
        double error = 0, scale = 0;
        for (int i = 0; i < 10; ++i)
            for (int j = 0; j < row.ny_m; ++j)
                for (int k = 0; k < row.nz_m; ++k) {
                    const double x = (i + .5) * row.h_m[0], y = (j + .5 + offset[1]) * row.h_m[1],
                                 z = (k + .5 + offset[2]) * row.h_m[2], t = 30 * row.dt_m;
                    const Complex exact(real(x, y, z, t, 2), imaginary(x, y, z, t, 2));
                    const double difference =
                        std::abs(result.potential_m[(i * row.ny_m + j) * row.nz_m + k] - exact);
                    error = std::isfinite(difference) ? std::max(error, difference) : INFINITY;
                    scale = std::max(scale, std::abs(exact));
                }
        const auto measured = fit(project(result.potential_m, row, reading, false, false, 0, 10),
                                  row.oracle().reflectedRoot);
        return {std::abs(measured - row.oracle().layerZero), error / scale};
    }

    /** @brief Compare projected modes against a non-cancelling serial input scale.
     * @details Each plane's scale is the transverse mean of absolute input values;
     * use its largest value over the same ten causal planes as the projection.
     * This measures the reduction separately from per-cell agreement.
     */
    double projectionError(const std::vector<Complex>& a, const std::vector<Complex>& b,
                           const Case& row, SpatialDiscretization reading, bool magnetic,
                           bool field) {
        const auto projected = project(a, row, reading, magnetic, field, 0, 10);
        const auto reference = project(b, row, reading, magnetic, field, 0, 10);
        double error = 0, scale = 0;
        for (int i = 0; i < 10; ++i) {
            error              = std::max(error, std::abs(projected[i] - reference[i]));
            double absoluteSum = 0;
            for (int j = 0; j < row.ny_m * row.nz_m; ++j)
                absoluteSum += std::abs(b[i * row.ny_m * row.nz_m + j]);
            scale = std::max(scale, absoluteSum / (row.ny_m * row.nz_m));
        }
        return error / scale;
    }

    /** @brief Replay on global rank zero with MPI_COMM_SELF, without modifying global
     * communicators.
     * @details Cells must agree within 1e-13 and projections within 1e-12. On one rank the
     * replay repeats an identical run, so the check then tests only determinism.
     */
    void verifyRanks(const Case& row, SpatialDiscretization reading, bool magnetic,
                     AJFourFieldsType& source, const Result& distributed) {
        RankAgreement report{};
        double projection = 0;
        if (source.getLayout().comm.rank() == 0) {
            SerialFields fields(source);
            const auto serial = run(row, reading, magnetic, 0, fields.source_m, fields.electric_m,
                                    fields.magnetic_m);
            const auto potential = rankDifference(distributed.potential_m, serial.potential_m);
            const auto field     = rankDifference(distributed.field_m, serial.field_m);
            report               = {std::max(potential[0], field[0]), potential[1] * field[1]};
            projection           = std::max(
                projectionError(distributed.potential_m, serial.potential_m, row, reading, magnetic,
                                          false),
                projectionError(distributed.field_m, serial.field_m, row, reading, magnetic, true));
            std::cout << row.name_m << " reading=" << int(reading)
                      << " projection_error=" << projection << '\n';
        }
        checkRankAgreement(report, source.getLayout().comm, row.name_m);
        MPI_Bcast(&projection, 1, MPI_DOUBLE, 0, source.getLayout().comm);
        EXPECT_LE(projection, 1e-12);
    }

    /** @brief Check one prescribed row and both deliberately defective finishes.
     * @details The two passing criteria are absolute reflection-coefficient error
     * and normalized cell error, each at most 1e-11: above the measured 4.2e-14 cell
     * round-off floor and below the 1.3e-6 speed-control error. Each control (speed scaled by
     * \f$1 + 10^{-6}\f$, old halo order) must fail at least one criterion. The stale-wrap
     * control uses the vacuum 60-degree row; its one-rank run retains the serial
     * periodic-wrap defect.
     */
    void checkPotentialCase(const Case& row, SpatialDiscretization reading, bool controls) {
        auto fields = makeBoxFieldContainer(row.cells(), row.lengths(), 2, {false, true, true});
        initializeFieldStorage(fields);
        for (int control = 0; control <= (controls ? 2 : 0); ++control) {
            if (control == 2 && std::string_view(row.name_m) != "v60")
                continue;
            const auto result =
                run(row, reading, false, control, fields.getJ(), fields.getE(), fields.getB());
            if (control == 0)
                verifyRanks(row, reading, false, fields.getJ(), result);
            if (ippl::Comm->rank() == 0) {
                const auto error = potentialErrors(row, reading, result);
                SCOPED_TRACE(std::string(row.name_m) + " reading=" + std::to_string(int(reading))
                             + " control=" + std::to_string(control));
                if (control == 0) {
                    EXPECT_LE(error[0], 1e-11);
                    EXPECT_LE(error[1], 1e-11);
                } else {
                    EXPECT_GT(std::max(error[0], error[1]), 1e-11);
                }
                std::cout << row.name_m << " reading=" << int(reading) << " control=" << control
                          << " reflection_error=" << error[0] << " cell_error=" << error[1] << '\n';
            }
        }
    }

    /** @brief Check the 13 rows of OpenBoxCases.h in both readings, with both controls.
     * Pass: coefficient and cell errors at most 1e-11 without a control, above 1e-11 with
     * each control, and the MPI_COMM_SELF replay within 1e-13 (cells) and 1e-12 (projections).
     */
    TEST(MurReflection, ExactTwoWaveRows) {
        for (const auto& row : cases)
            for (auto reading :
                 {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered})
                checkPotentialCase(row, reading, true);
    }

    /** @brief Independent anisotropic medium rows with a one-rank replay.
     * @details The two modes use a separate permittivity, aspect ratio and integer
     * periods. The propagating row also checks both field polarizations.
     * Rows: \f$\varepsilon_r = 3.1\f$, h = (1.15, 1, 0.85),
     * \f$\Delta t = 0.9\sqrt{\varepsilon_r/\sum_d h_d^{-2}}\f$, staggered R2; one propagating row
     * (36 z cells) and one evanescent row (19 z cells). Criteria: those of checkPotentialCase()
     * without controls (1e-11), the field rows of checkFieldReflection() (1e-10 and 1e-9) and
     * the rank replay (1e-13). CTest runs these rows on two ranks as
     * chdr.solver.mur.anisotropic.mpi2.
     */
    TEST(MurReflection, AnisotropicMedium) {
        const std::array<double, 3> h{1.15, 1, .85};
        const double dt =
            .9 * std::sqrt(3.1) / std::sqrt(1 / (h[0] * h[0]) + 1 + 1 / (h[2] * h[2]));
        for (const Case row :
             {Case{"anisotropic-propagating", 2, 36, 3.1, .17776026053686098, dt, h},
              Case{"anisotropic-evanescent", 2, 19, 3.1, .17677322362748854, dt, h}}) {
            checkPotentialCase(row, SpatialDiscretization::Staggered, false);
            if (row.nz_m != 36)
                continue;
            auto fields = makeBoxFieldContainer(row.cells(), row.lengths(), 2, {false, true, true});
            initializeFieldStorage(fields);
            for (bool magnetic : {false, true}) {
                const auto result = run(row, SpatialDiscretization::Staggered, magnetic, 0,
                                        fields.getJ(), fields.getE(), fields.getB());
                verifyRanks(row, SpatialDiscretization::Staggered, magnetic, fields.getJ(), result);
                if (ippl::Comm->rank() == 0)
                    checkFieldReflection(row, SpatialDiscretization::Staggered, magnetic, result);
            }
        }
    }
}  // namespace chdr::solver::test::reflection
