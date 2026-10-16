/** @file test_open_field_reflection.cpp
 * @brief Check signed electric and magnetic reflection at their declared sample positions.
 */
#include "ReflectionRun.h"

namespace chdr::solver::test::reflection {
    namespace {
        /** @brief Evaluate the incident oracle wave under its documented field sampling rule.
         * @details For TE, \f$E_y=-(A_y^{30}-A_y^{29})/\Delta t\f$. For TM0,
         * \f$B_y=-D_x(A_z^{30}+A_z^{29})/2\f$, with centred or forward difference for
         * collocated or staggered samples. These evaluate the incident root explicitly;
         * no production reconstruction supplies an expected value. fieldRowRatio() documents
         * the corresponding reflected-wave signs and sample offsets, from the potential-field
         * relations [Fallahi (2020), MITHRA 2.0, section 3.1.1, Eqs. (3.6)-(3.8)].
         */
        Complex incidentField(const Case& row, Complex theta, bool magnetic, bool staggered,
                              int j) {
            const Complex imaginary(0, 1);
            const auto normal = [&](int index) {
                return std::exp(-imaginary * theta * double(index));
            };
            const Complex current  = std::polar(1., -30 * row.omega_m * row.dt_m),
                          previous = std::polar(1., -29 * row.omega_m * row.dt_m);
            if (!magnetic)
                return -normal(j) * (current - previous) / row.dt_m;
            const Complex derivative = staggered
                                           ? (normal(j + 1) - normal(j)) / row.h_m[0]
                                           : (normal(j + 1) - normal(j - 1)) / (2 * row.h_m[0]);
            return -derivative * (current + previous) / 2.;
        }
    }  // namespace

    /** @brief Compare all transverse cells on layers one through eight, then fit their signed
     * ratio.
     * @details The field error uses the positive incident-component amplitude at layer one.
     * Both waves have time label \f$29.5\Delta t\f$. The ratio expectation comes directly
     * from fieldRowRatio(); fit() refers its answer to the first supplied plane, also layer one.
     * The transverse projection respects the selected field component's position offsets.
     * This routine only examines an already evolved result and never starts a trajectory.
     */
    void checkFieldReflection(const Case& row, SpatialDiscretization reading, bool magnetic,
                              const Result& result) {
        if (ippl::Comm->rank() != 0)
            return;
        ASSERT_EQ(result.field_m.size(), std::size_t(40 * row.ny_m * row.nz_m));
        const auto reflection = row.oracle();
        const bool staggered  = reading == SpatialDiscretization::Staggered;
        const auto offset = LevelContract{EBReconstructionConvention::R2, reading}.positionOffset(
            magnetic ? Quantity::Magnetic : Quantity::Electric, 1);
        const auto kt = row.tangent();
        const double scale =
            std::abs(incidentField(row, reflection.reflectedRoot, magnetic, staggered, 1));
        ASSERT_GT(scale, 0.);
        double error = 0;
        for (int i = 1; i <= 8; ++i) {
            const Complex normal =
                incidentField(row, reflection.reflectedRoot, magnetic, staggered, i)
                * (1. + fieldRowRatio(reflection, magnetic, staggered, true, i));
            for (int j = 0; j < row.ny_m; ++j)
                for (int k = 0; k < row.nz_m; ++k) {
                    const double phase = kt[1] * (j + .5 + offset[1]) * row.h_m[1]
                                         + kt[2] * (k + .5 + offset[2]) * row.h_m[2];
                    const double difference =
                        std::abs(result.field_m[(i * row.ny_m + j) * row.nz_m + k]
                                 - normal * std::polar(1., phase));
                    error = std::isfinite(difference) ? std::max(error, difference)
                                                      : std::numeric_limits<double>::infinity();
                }
        }
        const auto measured = fit(project(result.field_m, row, reading, magnetic, true, 1, 8),
                                  reflection.reflectedRoot);
        const double coefficientError =
            std::abs(measured - fieldRowRatio(reflection, magnetic, staggered, true, 1));
        EXPECT_LE(error / scale, 1.e-10);
        EXPECT_LE(coefficientError, 1.e-9);
        std::cout << row.name_m << " reading=" << int(reading) << " magnetic=" << magnetic
                  << " field_error=" << error / scale
                  << " field_coefficient_error=" << coefficientError << '\n';
    }

    /// Two propagating rows, two readings and two polarizations, with identical MPI/serial
    /// fixtures.
    TEST(OpenFieldReflection, SignedElectricAndMagnetic) {
        for (const auto& row : {cases[1], cases[7]})
            for (auto reading :
                 {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered})
                for (bool magnetic : {false, true}) {
                    SCOPED_TRACE(::testing::Message()
                                 << row.name_m << ' ' << int(reading) << ' ' << magnetic);
                    auto fields =
                        makeBoxFieldContainer(row.cells(), row.lengths(), 2, {false, true, true});
                    initializeFieldStorage(fields);
                    const auto result =
                        run(row, reading, magnetic, 0, fields.getJ(), fields.getE(), fields.getB());
                    checkFieldReflection(row, reading, magnetic, result);
                    verifyRanks(row, reading, magnetic, fields.getJ(), result);
                }
    }
}  // namespace chdr::solver::test::reflection
