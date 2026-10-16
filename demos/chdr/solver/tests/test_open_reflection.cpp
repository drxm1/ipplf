/** @file test_open_reflection.cpp
 * @brief Exact two-wave continuation, sensitive controls and MPI witnesses for Mur faces.
 */
#include "RankWitness.h"
#include "ReflectionRun.h"

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

    /// Replay on global rank zero with MPI_COMM_SELF, without modifying global communicators.
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
     * round-off floor and below the 1.3e-6 speed-control error. Each speed perturbation must
     * fail at least one criterion. The stale-wrap control uses the vacuum 60-degree
     * row; its one-rank run retains the serial periodic-wrap defect.
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

    TEST(OpenReflection, ExactTwoWaveRows) {
        for (const auto& row : {cases[0], cases[1], cases[6], cases[7], cases[8], cases[11]})
            for (auto reading :
                 {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered})
                checkPotentialCase(row, reading, true);
    }

}  // namespace chdr::solver::test::reflection
