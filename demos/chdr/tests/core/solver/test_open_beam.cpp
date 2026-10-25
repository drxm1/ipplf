/** @file test_open_beam.cpp
 * @brief Compare the central plane-wave reflection of a finite beam and a periodic witness.
 * The spatial projection and two-wave fit distinguish the selected tangential
 * mode from the beam's other angles [jones1988MultiPointImpedance, equations
 * (1)-(3), p. 3; oskooi2010meep, equation (6)].
 */
#include <Kokkos_Complex.hpp>
#include <gtest/gtest.h>

#include "CoreTestSupport.h"
#include "RankWitness.h"
#include "SerialReference.h"
#include "TwoWaveState.h"
#include "core/solver/BoxSetup.h"

namespace chdr::solver::test {
    namespace beam {
        using Complex       = std::complex<double>;
        using Projection    = std::array<Complex, 20>;
        constexpr double Dt = .5196152422706632, Omega = .6346975625940523, Kz = .4487989505128276;
        constexpr double Width = 31.51107106549871, WaistY = 19.79898987322333,
                         WaistZ = 39.59797974644666;
        constexpr int Steps = 1436, Center = 364;

        /** @brief Write one real quadrature of the Gaussian sheet source at global x index 20.
         * The temporal phase uses the consumed source time; the field projection
         * below uses the completed potential time. Other components remain zero.
         */
        void drive(AJFourFieldsType& source, int step, bool finiteBeam, int quadrature) {
            const auto view  = source.getView();
            const auto local = source.getLayout().getLocalNDIndex();
            const int y0 = local[1].first(), z0 = local[2].first();
            const double pulse = std::exp(-.5 * std::pow((step - Center) * Dt / Width, 2));
            Kokkos::parallel_for(
                "Beam source sheet",
                Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
                    {1, 1}, {int(local[1].length()) + 1, int(local[2].length()) + 1}),
                KOKKOS_LAMBDA(int j, int k) {
                    const double y = y0 + j - .5, z = z0 + k - .5,
                                 phase = Kz * z - Omega * Dt * step;
                    const double envelope =
                        finiteBeam ? Kokkos::exp(-.5 * Kokkos::pow((y - 80) / WaistY, 2)
                                                 - .5 * Kokkos::pow((z - 160) / WaistZ, 2))
                                   : 1.;
                    view(21, j, k)[1] =
                        pulse * envelope * (quadrature ? Kokkos::sin(phase) : Kokkos::cos(phase));
                });
        }

        /** @brief Add the spatial and temporal Fourier sum on the twenty source-free layers.
         * Each rank accumulates its owned transverse cells. One MPI sum at the end
         * completes the projection; ghosts never enter the sum. The complex source
         * is represented by two independent real trajectories.
         */
        void project(AJFourFieldsType& field, int step, Projection& output) {
            const auto view  = field.getView();
            const auto local = field.getLayout().getLocalNDIndex();
            const int z0     = local[2].first();
            const auto time  = std::polar(1., Omega * Dt * (step + 1));
            for (int x = 0; x < 20; ++x) {
                Kokkos::complex<double> sum;
                Kokkos::parallel_reduce(
                    "Beam transverse projection",
                    Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
                        {1, 1}, {int(local[1].length()) + 1, int(local[2].length()) + 1}),
                    KOKKOS_LAMBDA(int j, int k, Kokkos::complex<double>& value) {
                        const double phase = Kz * (z0 + k - .5), sample = view(x + 1, j, k)[1];
                        value += Kokkos::complex<double>(sample * Kokkos::cos(phase),
                                                         -sample * Kokkos::sin(phase));
                    },
                    sum);
                output[x] += Complex(sum.real(), sum.imag()) * time;
            }
        }

        struct Result {
            Projection projection_m{};
            std::vector<double> cells_m;
            double maximum_m = 0;
        };

        /** @brief Run the prescribed source from zero histories, optionally retaining rank
         * evidence. The non-cancelling reduction scale is the number of projection terms times the
         * largest owned potential during this real trajectory.
         */
        Result run(AJFourFieldsType& source, EMFieldsType& electric, EMFieldsType& magnetic,
                   bool finiteBeam, int steps, int quadrature, bool retainCells) {
            ChdrSolver<StandardSolverPeriodic> solver(source, electric, magnetic,
                                                      EBReconstructionConvention::R2);
            solver.setTimeStepUnchecked(Dt);
            source = 0.;
            Result result;
            for (int step = 0; step < steps; ++step) {
                drive(source, step, finiteBeam, quadrature);
                solver.step();
                solver.timeShift();
                project(solver.A_n, step, result.projection_m);
                if (retainCells) {
                    const auto view  = solver.A_n.getView();
                    result.maximum_m = std::max(result.maximum_m,
                                                coreTest::localMaximum(
                                                    solver.A_n, KOKKOS_LAMBDA(int i, int j, int k) {
                                                        return Kokkos::abs(view(i, j, k)[1]);
                                                    }));
                }
            }
            MPI_Allreduce(MPI_IN_PLACE, result.projection_m.data(), 40, MPI_DOUBLE, MPI_SUM,
                          source.getLayout().comm);
            if (retainCells) {
                result.maximum_m = coreTest::globalMaximum(solver.A_n, result.maximum_m);
                result.cells_m   = gatherOwned(solver.A_n);
            }
            return result;
        }

        /** @brief Fit the two propagating roots, both unit-modulus, to the projected twenty layers.
         * Solving the two-column normal equations gives the reflected/incident
         * coefficient at layer zero [jones1988MultiPointImpedance, equations (1)-(3)].
         * The normal phase comes only from the independent dispersion oracle.
         */
        Complex reflection(const Projection& values) {
            const auto mode = discreteReflection({1, 1, 1}, Dt, 1, 1, OracleStencil::Standard, 0,
                                                 {0, 0, Kz}, Omega);
            return chdr::solver::test::reflection::fit(
                std::vector<Complex>(values.begin(), values.end()), mode.reflectedRoot);
        }

        /// Combine the two real source quadratures before fitting a complex coefficient.
        Complex coefficient(std::array<int, 3> cells, bool finiteBeam) {
            auto fields =
                makeBoxFieldContainer(cells, {double(cells[0]), double(cells[1]), double(cells[2])},
                                      2, {false, !finiteBeam, !finiteBeam});
            initializeFieldStorage(fields);
            Projection combined{};
            for (int quadrature = 0; quadrature < 2; ++quadrature) {
                const auto result = run(fields.getJ(), fields.getE(), fields.getB(), finiteBeam,
                                        Steps, quadrature, false);
                for (int x = 0; x < 20; ++x)
                    combined[x] +=
                        (quadrature ? Complex(0, 1) : Complex(1, 0)) * result.projection_m[x];
            }
            return reflection(combined);
        }
    }  // namespace beam

    /** @brief Require a finite beam and periodic witness to give the same central-mode reflection.
     * Both use the same carrier and time gate. Bounds of 1e-4 absolute and one percent of
     * the witness amplitude resolve a meaningful beam-width bias while admitting the measured
     * 4e-6 finite-window difference; they are not a round-off equality criterion.
     */
    TEST(OpenBeam, FullBoxAgainstPeriodicWitness) {
        const auto witness      = beam::coefficient({40, 2, 14}, false);
        const auto finite       = beam::coefficient({40, 160, 400}, true);
        const double difference = std::abs(finite - witness);
        std::cout << "beam cells=40x160x400 witness=40x2x14 steps=" << beam::Steps
                  << " dt=" << beam::Dt << " witness=" << witness << " finite=" << finite
                  << " difference=" << difference << '\n';
        EXPECT_LE(difference, 1e-4);
        EXPECT_LE(difference, .01 * std::abs(witness));
    }

    /** @brief Compare the reduced beam trajectory on the active communicator and on one rank.
     * Compare every final owned potential with a run-maximum scale. Compare Fourier
     * sums using their non-cancelling input scale, including when a sum is near zero.
     * Record binary64 equality, including signed zeros, separately from the tolerances.
     */
    TEST(OpenBeam, ReducedRankWitness) {
        auto fields = makeBoxFieldContainer({40, 80, 200}, {40, 80, 200}, 2, {false, false, false});
        initializeFieldStorage(fields);
        for (int quadrature = 0; quadrature < 2; ++quadrature) {
            const auto distributed =
                beam::run(fields.getJ(), fields.getE(), fields.getB(), true, 300, quadrature, true);
            RankAgreement cells{}, projection{};
            if (ippl::Comm->rank() == 0) {
                SerialFields self(fields.getJ());
                const auto serial = beam::run(self.source_m, self.electric_m, self.magnetic_m, true,
                                              300, quadrature, true);
                cells      = rankDifference(distributed.cells_m, serial.cells_m, serial.maximum_m);
                projection = rankDifference(distributed.projection_m, serial.projection_m,
                                            300. * 80 * 200 * serial.maximum_m);
                std::cout << "reduced beam quadrature=" << quadrature << '\n';
            }
            checkRankAgreement(cells, MPI_COMM_WORLD, "reduced beam cells");
            checkRankAgreement(projection, MPI_COMM_WORLD, "reduced beam projection", 1e-12);
        }
    }
}  // namespace chdr::solver::test
