/** @file test_outer_layer.cpp
 * @brief Check reconstructed outer-layer fields against independent polynomial values.
 */
#include <complex>
#include <gtest/gtest.h>

#include "OuterLayerTest.h"
#include "RankWitness.h"
#include "SerialReference.h"
#include "core/solver/BoxSetup.h"

namespace chdr::solver::test {
    namespace {
        using Shell      = ChdrSolver<StandardSolverPeriodic>;
        namespace oracle = open_boundary_oracle;
        constexpr std::array<int, 3> Cells{8, 6, 10};
        constexpr std::array<double, 3> Spacing{1, .8, 1.25};

        /** @brief Load owned samples at each potential component's spatial and temporal labels.
         * Ghosts start at zero so the unchanged reconstruction can first demonstrate
         * the missing-extension error. Expected samples come only from the oracle.
         */
        void loadPolynomial(Shell& solver, bool periodicY) {
            const auto coefficients = oracle::polynomialFixture(periodicY);
            for (auto* field : {&solver.A_n, &solver.A_nm1}) {
                *field    = 0.;
                auto host = field->getHostMirror();
                Kokkos::deep_copy(host, field->getView());
                outer_layer_test::visit(*field, false, [&](int i, int j, int k, auto index) {
                    for (unsigned q = 0; q < 4; ++q) {
                        const auto quantity = q == 0 ? Quantity::Phi : Quantity::VectorPotential;
                        const auto offset =
                            solver.levelContract().positionOffset(quantity, q == 0 ? 0 : q - 1);
                        std::array<double, 3> x;
                        for (unsigned d = 0; d < 3; ++d)
                            x[d] = (index[d] + .5 + offset[d]) * Spacing[d];
                        const double time = solver.labelTime(quantity, 0)
                                            - (field == &solver.A_nm1 ? solver.getDt() : 0.);
                        host(i, j, k)[q] =
                            oracle::polynomialPotential(coefficients, q, x, time, -1, false);
                    }
                });
                Kokkos::deep_copy(field->getView(), host);
                field->fillHalo();
            }
        }

        /** @brief Normalize the maximum layer error by the positive analytic full-box maximum.
         * Include only the outermost owned layers of open axes, including their
         * edges and corners. Each component uses its own sample location.
         */
        double layerError(EMFieldsType& field, Shell& solver, bool periodicY, Quantity quantity) {
            const auto c = oracle::polynomialFixture(periodicY);
            const bool staggered =
                solver.spatialDiscretization() == SpatialDiscretization::Staggered;
            const bool magnetic = quantity != Quantity::Electric;
            const double time   = quantity == Quantity::RawMagnetic ? 0. : -solver.getDt() / 2;
            auto host           = field.getHostMirror();
            Kokkos::deep_copy(host, field.getView());
            double error = 0;
            outer_layer_test::visit(field, false, [&](int i, int j, int k, auto index) {
                bool layer = false;
                for (unsigned d = 0; d < 3; ++d)
                    layer |= !(d == 1 && periodicY) && (index[d] == 0 || index[d] == Cells[d] - 1);
                if (!layer)
                    return;
                for (unsigned q = 0; q < 3; ++q) {
                    const auto x = oracle::fieldPosition(index, Spacing, staggered, magnetic, q);
                    const double delta = std::abs(
                        host(i, j, k)[q] - oracle::polynomialField(c, magnetic, q, x, time));
                    error = std::max(error, std::isfinite(delta)
                                                ? delta
                                                : std::numeric_limits<double>::infinity());
                }
            });
            return outer_layer_test::maximum(error)
                   / oracle::polynomialFieldMaximum(c, Cells, Spacing, staggered, magnetic, time);
        }
        /** @brief Compare all owned E, averaged B and raw B samples with an independent-rank
         * replay.
         * @details Retain the physical step, reading and periodic axes. Only global rank zero
         * reconstructs the self-communicator fields; the shared comparison uses each field's
         * nonzero serial maximum and records binary equality separately from its tolerance.
         */
        void checkFieldRanks(Shell& solver, const EMFieldsType& raw, bool periodicY) {
            if (solver.A_n.getLayout().comm.size() == 1)
                return;
            const auto electric    = gatherOwned(solver.electricField());
            const auto magnetic    = gatherOwned(solver.magneticField());
            const auto rawMagnetic = gatherOwned(raw);
            RankAgreement electricRank{}, magneticRank{}, rawRank{};
            if (solver.A_n.getLayout().comm.rank() == 0) {
                SerialFields fields(solver.A_n);
                Shell serial(fields.source_m, fields.electric_m, fields.magnetic_m,
                             EBReconstructionConvention::R2, solver.spatialDiscretization());
                serial.setTimeStepUnchecked(solver.getDt());
                EMFieldsType serialRaw(fields.mesh_m, fields.layout_m);
                loadPolynomial(serial, periodicY);
                serial.prepareHistories();
                serial.evaluateInitialFields();
                serial.computeRawMagneticField(serialRaw);
                std::cout << "outer-layer reading=" << int(solver.spatialDiscretization())
                          << " periodic-y=" << periodicY << '\n';
                electricRank = rankDifference(electric, gatherOwned(fields.electric_m));
                magneticRank = rankDifference(magnetic, gatherOwned(fields.magnetic_m));
                rawRank      = rankDifference(rawMagnetic, gatherOwned(serialRaw));
            }
            checkRankAgreement(electricRank, solver.A_n.getLayout().comm, "outer-layer E");
            checkRankAgreement(magneticRank, solver.A_n.getLayout().comm, "outer-layer averaged B");
            checkRankAgreement(rawRank, solver.A_n.getLayout().comm, "outer-layer raw B");
        }
    }  // namespace

    /** @brief Require quadratic closure accuracy for collocated and staggered half-step fields.
     * The zero-ghost control must exceed 0.1 of the analytic maximum for E,
     * averaged B and raw B. Prepared fields must be within 1e-13 on both
     * all-open and mixed boxes. CTest runs the same fixtures on one and two ranks.
     */
    TEST(OpenOuterLayer, PolynomialFieldsAndZeroGhostControl) {
        for (bool periodicY : {false, true})
            for (auto reading :
                 {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered}) {
                auto fields =
                    makeBoxFieldContainer(Cells, {8, 4.8, 12.5}, 2, {false, periodicY, false});
                initializeFieldStorage(fields);
                Shell solver(fields.getJ(), fields.getE(), fields.getB(),
                             EBReconstructionConvention::R2, reading);
                solver.setSafetyFactor(.9);
                EMFieldsType raw(fields.getMesh(), fields.getFL());
                loadPolynomial(solver, periodicY);
                for (bool closed : {false, true}) {
                    if (closed)
                        solver.prepareHistories();
                    solver.evaluateInitialFields();
                    solver.computeRawMagneticField(raw);
                    const auto errors =
                        std::array{layerError(fields.getE(), solver, periodicY, Quantity::Electric),
                                   layerError(fields.getB(), solver, periodicY, Quantity::Magnetic),
                                   layerError(raw, solver, periodicY, Quantity::RawMagnetic)};
                    for (double error : errors) {
                        if (closed)
                            EXPECT_LE(error, 1e-13);
                        else
                            EXPECT_GE(error, .1);
                    }
                }
                checkFieldRanks(solver, raw, periodicY);
            }
    }
    /** @brief Record the missing-ghost magnetic error on unit-amplitude normal TE snapshots.
     * @details Compare raw \f$B_z=D_x A_y\f$ with the exact lattice derivative of the
     * continued plane wave. Normalize by \f$k|A_y|\f$, where \f$k=2\pi/\mathrm{PPW}\f$.
     * The missing term has magnitude \f$1/(2kh)\f$ for the centred lower face and
     * \f$1/(kh)\f$ for the forward upper face (deduction from the difference rows in
     * Fornberg (1988), Table 1, and Ryu et al. (2016), Appendix A). No evolution is used.
     */
    TEST(OpenOuterLayerRecords, ZeroGhostMagnitude) {
        if (ippl::Comm->size() != 1)
            GTEST_SKIP() << "This record uses one rank";
        for (auto reading : {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered})
            for (int ppw : {10, 20, 40, 80}) {
                auto fields = makeBoxFieldContainer({24, 2, 4}, {24, 2, 4}, 2, {false, true, true});
                initializeFieldStorage(fields);
                Shell solver(fields.getJ(), fields.getE(), fields.getB(),
                             EBReconstructionConvention::R2, reading);
                EMFieldsType raw(fields.getMesh(), fields.getFL());
                const bool staggered = reading == SpatialDiscretization::Staggered;
                const double k       = 2 * Kokkos::numbers::pi_v<double> / ppw;
                std::complex<double> measured{};
                for (int part = 0; part < 2; ++part) {
                    solver.A_n = 0.;
                    writeSamples(
                        solver.A_n,
                        KOKKOS_LAMBDA(double x, double, double, double, unsigned q) {
                            return q == 2 ? (part == 0 ? Kokkos::cos(k * x) : Kokkos::sin(k * x))
                                          : 0.;
                        },
                        {0., 0., 0., 0.}, solver.levelContract());
                    solver.A_n.fillHalo();
                    solver.computeRawMagneticField(raw);
                    const auto host    = coreTest::copyToHost(raw);
                    const double value = host(staggered ? 24 : 1, 1, 1)[2];
                    measured +=
                        std::complex<double>(part == 0 ? value : 0., part == 1 ? value : 0.);
                }
                const auto expected =
                    std::complex<double>(0, staggered ? 2 * std::sin(k / 2) : std::sin(k))
                    * std::polar(1., k * (staggered ? 24. : .5));
                std::cout << "zero-ghost ppw=" << ppw << " field=raw-Bz reading=" << int(reading)
                          << " face=" << (staggered ? "upper" : "lower")
                          << " normalization=k*abs(Ay)=" << k
                          << " relative_error=" << std::abs(measured - expected) / k << '\n';
            }
    }
}  // namespace chdr::solver::test
