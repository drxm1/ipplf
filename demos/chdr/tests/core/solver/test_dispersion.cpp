/** @file test_dispersion.cpp
 * @brief Measure dispersion by projecting the evolved potential arrays.
 * @details Single-root histories and the independent dispersion prediction use
 * [fallahi2020mithra20fullwavesimulation, Eqs. (3.17), (3.24)-(3.25)].
 * T5 identifies these dispersion records. Projection uses the spatial kernel
 * \f$\exp(-i k\cdot x)\f$ for modes \f$\exp(i(k\cdot x-\omega t))\f$.
 * All potential components share the same scalar update, so the test projects
 * A_x in every direction; field reconstruction is checked by the Convergence suites.
 */
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <iostream>

#include "CoreTestSupport.h"
#include "MediumOracles.h"
#include "MonitorTestSupport.h"
#include "RankRecords.h"
#include "SerialReference.h"
#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "core/solver/ChdrSolver.h"

static_assert(
    chdr::SpatialDim == 3,
    "test_dispersion.cpp requires three spatial dimensions for its Cartesian implementation");

namespace chdr::solver::test {
    using namespace chdr::diagnostics;
    struct TravellingMode {
        ippl::Vector<double, SpatialDim> waveNumber_m;
        double omega_m;
        KOKKOS_FUNCTION double operator()(double x, double y, double z, double t,
                                          unsigned c) const {
            return c == 1 ? Kokkos::cos(waveNumber_m[0] * x + waveNumber_m[1] * y
                                        + waveNumber_m[2] * z - omega_m * t)
                          : 0;
        }
    };

    /** @brief Project owned potential samples using separate real and imaginary reductions.
     * @tparam PotentialView Three-dimensional view of four-component potential values;
     * component 1 stores the x component of the vector potential used by this test.
     */
    template <class PotentialView>
    struct Projection {
        PotentialView view_m;
        ippl::Vector<int, SpatialDim> first_m;
        ippl::Vector<double, SpatialDim> spacing_m, waveNumber_m;
        int part_m;
        KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, double& result) const {
            const double phase = waveNumber_m[0] * (i + first_m[0] - 0.5) * spacing_m[0]
                                 + waveNumber_m[1] * (j + first_m[1] - 0.5) * spacing_m[1]
                                 + waveNumber_m[2] * (k + first_m[2] - 0.5) * spacing_m[2];
            result += view_m(i, j, k)[1] * (part_m == 0 ? Kokkos::cos(phase) : -Kokkos::sin(phase));
        }
    };

    /** @brief Measure the current x-potential mode without changing its samples.
     * @tparam FourPotentialA IPPL field with scalar potential at component 0 and
     * vector potential at components 1 through SpatialDim.
     * @param currentA Current potential history, sampled at cell centres in this test.
     * @param waveNumber Physical wave vector of the Fourier mode to project.
     */
    template <class FourPotentialA>
    std::complex<double> projectMode(const FourPotentialA& currentA,
                                     ippl::Vector<double, SpatialDim> waveNumber) {
        const auto view    = currentA.getView();
        const auto first   = currentA.getLayout().getLocalNDIndex().first();
        const auto spacing = currentA.get_mesh().getMeshSpacing();
        double sums[2]     = {};
        for (int part = 0; part < 2; ++part)
            Kokkos::parallel_reduce("spatial Fourier mode", currentA.getFieldRangePolicy(),
                                    Projection{view, first, spacing, waveNumber, part}, sums[part]);
        MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, currentA.getLayout().comm);
        return {sums[0], sums[1]};
    }

    /** @brief Geometry and measurement policy for one periodic potential mode.
     * @details epsilon_m equal to zero selects the owner-free vacuum solver;
     * a positive value constructs a uniform material. recordOnly_m suppresses
     * the symbol-agreement assertion, but still measures and records the phase.
     * splitAxis_m overrides the usual decomposition along a varying direction.
     */
    struct DispersionInput {
        std::array<int, SpatialDim> cells_m, mode_m;
        std::array<double, SpatialDim> lengths_m;
        double safety_m;
        std::string name_m;
        bool rankRecord_m;
        double epsilon_m  = 0;
        bool recordOnly_m = false;
        int splitAxis_m   = -1;
        /// Split an axis where the wave varies, so the MPI seam is exercised.
        unsigned splitAxis() const {
            if (splitAxis_m >= 0)
                return unsigned(splitAxis_m);
            unsigned split = 0;
            while (split < 2 && mode_m[split] == 0)
                ++split;
            return split;
        }
    };
    struct DispersionPrediction {
        TravellingMode wave_m;
        double omegaStep_m, waveNormSquared_m;
    };

    /** @brief Independent host substitution into the selected spatial symbol.
     * @details The frequency is selected before evolution and used only for
     * initialization and comparison. The measurement below projects the actual
     * updated arrays. Formula source: [fallahi2020mithra20fullwavesimulation,
     * Eqs. (3.17), (3.24)-(3.25)].
     */
    template <class SolverType>
    DispersionPrediction prediction(SolverType& solver, const DispersionInput& input) {
        const auto h    = solver.A_n.get_mesh().getMeshSpacing();
        const double dt = solver.getDt(), pi = Kokkos::numbers::pi_v<double>;
        ippl::Vector<double, SpatialDim> k;
        double compact2[SpatialDim], kn2 = 0;
        for (int d = 0; d < static_cast<int>(SpatialDim); ++d) {
            k[d]        = 2 * pi * input.mode_m[d] / input.lengths_m[d];
            compact2[d] = std::pow(2 * std::sin(k[d] * h[d] / 2) / h[d], 2);
            kn2 += k[d] * k[d];
        }
        if (input.epsilon_m > 0) {
            const double theta = mediumPhasePerStep(k, h, dt, input.epsilon_m);
            return {{k, theta / dt}, theta, kn2};
        }
        const double alpha =
            0.25 * (1 + 0.02 / (std::pow(h[2] / h[0], 2) + std::pow(h[2] / h[1], 2)));
        const double sz = SolverType::stencil() == Stencil::NonStandard
                              ? 1 - 4 * alpha * std::pow(std::sin(k[2] * h[2] / 2), 2)
                              : 1;
        const double omegaStep =
            2 * std::asin(dt * std::sqrt(sz * (compact2[0] + compact2[1]) + compact2[2]) / 2);
        return {{k, omegaStep / dt}, omegaStep, kn2};
    }

    /** @brief Cumulative phase of the evolved Fourier coefficient, in radians.
     * @details The travelling mode has decreasing Fourier phase, so frequency()
     * negates the accumulated increments. Despite its name it returns radians
     * per step, not per unit time; divide by the actual step for angular frequency.
     * expected is used only for the increment-error diagnostic, not the measurement.
     */
    struct PhaseMeasurement {
        double phase_m = 0, largestIncrementError_m = 0;
        int steps_m = 0;
        /// Accumulate principal increments; every frozen per-step phase lies below pi.
        void observe(std::complex<double> current, std::complex<double> previous, double expected) {
            const double increment = std::arg(current * std::conj(previous));
            phase_m += increment;
            largestIncrementError_m =
                std::max(largestIncrementError_m, std::abs(increment + expected));
        }
        double frequency() const { return -phase_m / steps_m; }
    };

    /** @brief Select one travelling root before measuring the evolved phase.
     * @details Sampling both initial slots with the discrete frequency selects
     * the same root of the two-level recurrence (deduction from
     * [fallahi2020mithra20fullwavesimulation, Eqs. (3.17), (3.24)]).
     * A continuum-frequency history would mix the two
     * discrete roots and introduce a transient into the finite-window phase.
     */
    template <class SolverType>
    void measurePhase(SolverType& solver, const DispersionPrediction& expected,
                      PhaseMeasurement& measured) {
        loadHistory(solver, expected.wave_m);
        solver.evaluateInitialFields();
        measured.steps_m = std::ceil(2 * Kokkos::numbers::pi_v<double> / expected.omegaStep_m);
        auto previous    = projectMode(solver.A_n, expected.wave_m.waveNumber_m);
        // Monitors implements vacuum residuals; medium phase needs only the evolved potentials.
        std::unique_ptr<Monitors> monitors;
        if (!solver.hasMaterial())
            monitors = std::make_unique<Monitors>(solver.A_n);
        for (int n = 0; n < measured.steps_m; ++n) {
            if (monitors)
                advanceWithFinalMonitor(solver, *monitors, n == measured.steps_m - 1);
            else
                solver.solve();
            const auto current = projectMode(solver.A_n, expected.wave_m.waveNumber_m);
            ASSERT_GT(std::abs(current), 0);
            measured.observe(current, previous, expected.omegaStep_m);
            previous = current;
        }
    }

    template <class SolverType>
    void requirePhase(SolverType& solver, const DispersionInput& input,
                      const DispersionPrediction& expected, const PhaseMeasurement& measured) {
        EXPECT_TRUE(std::isfinite(measured.frequency()));
        EXPECT_LE(std::abs(measured.frequency() / expected.omegaStep_m - 1), 1e-12);
        if (SolverType::stencil() == Stencil::NonStandard && input.mode_m[0] == 0
            && input.mode_m[1] == 0 && input.safety_m == 1) {
            const double velocity =
                measured.frequency() / (solver.getDt() * std::sqrt(expected.waveNormSquared_m));
            EXPECT_LE(std::abs(velocity - 1), 1e-12);
        }
    }

    /** @brief Compare measured and discrete-symbol phase with the continuum medium wave.
     * @details frequency is the continuum angular frequency from mediumFrequency;
     * measured.frequency() and omegaStep_m are respectively measured and predicted
     * phase advances per step. Dividing either by denominator gives phase velocity
     * relative to the continuum velocity. phaseErrorPerWavelength converts its
     * relative error to radians per dielectric wavelength; see that oracle's
     * derivation from [fallahi2020mithra20fullwavesimulation, Eq. (3.17)].
     */
    template <class SolverType>
    nlohmann::json dispersionRecord(SolverType& solver, const DispersionInput& input,
                                    const DispersionPrediction& expected,
                                    const PhaseMeasurement& measured) {
        const double frequency =
            input.epsilon_m > 0 ? mediumFrequency(expected.wave_m.waveNumber_m, input.epsilon_m)
                                : std::sqrt(expected.waveNormSquared_m);
        const double denominator = solver.getDt() * frequency;
        return {{"name", input.name_m},
                {"cells", input.cells_m},
                {"lengths", input.lengths_m},
                {"mode", input.mode_m},
                {"dt", solver.getDt()},
                {"safety", solver.safetyFactor()},
                {"epsilon", input.epsilon_m > 0 ? input.epsilon_m : 1},
                {"formulation", solver.formulationName()},
                {"alpha", solver.formulationParameters().alpha_m},
                {"phase_per_wavelength",
                 phaseErrorPerWavelength(measured.frequency() / denominator - 1)},
                {"phase_per_wavelength_prediction",
                 phaseErrorPerWavelength(expected.omegaStep_m / denominator - 1)},
                {"steps", measured.steps_m},
                {"omega", measured.frequency()},
                {"omega_prediction", expected.omegaStep_m},
                {"phase_velocity_error", measured.frequency() / denominator - 1},
                {"phase_velocity_prediction", expected.omegaStep_m / denominator - 1},
                {"max_increment_error", measured.largestIncrementError_m}};
    }

    /** @brief Keep velocity errors and individual increments beside the mean phase advance.
     * @details The JSON omega value is the mean angular phase advance per step;
     * dividing it by the time step gives angular frequency. The decomposition
     * comparison uses that cumulative measurement, not individual increments.
     */
    inline void writeDispersion(const DispersionInput& input, const nlohmann::json& record) {
        std::cout << "T5 " << input.name_m << " omega=" << record["omega"].get<double>()
                  << " predicted=" << record["omega_prediction"].get<double>()
                  << " v_minus_1=" << record["phase_velocity_error"].get<double>()
                  << " predicted_v_minus_1=" << record["phase_velocity_prediction"].get<double>()
                  << " max_increment_error=" << record["max_increment_error"].get<double>() << '\n';
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank != 0)
            return;
        std::filesystem::create_directories(recordDirectory());
        std::ofstream out(recordDirectory() / ("T5-" + input.name_m + ".jsonl"));
        out << record.dump() << '\n';
        if (input.rankRecord_m) {
            std::ofstream comparison(recordDirectory() / ("T5-" + input.name_m + ".json"));
            comparison << nlohmann::json{{"omega", record["omega"]}}.dump() << '\n';
        }
    }

    /** @brief Evolve one discrete mode at the requested step and record its phase.
     * @details The production path applies safety to the material step bound.
     * vacuumSet instead selects the smaller vacuum-bound step through the unchecked
     * test entry. This records the phase error from retaining a vacuum-bound step;
     * the material remains uniform with permittivity eps.
     * Both initialization and prediction use getDt() after this selection.
     * recordOnly cases record dispersion without imposing the development threshold.
     */
    template <class PotentialSolverBase>
    void dispersionCase(std::array<int, SpatialDim> cells, std::array<double, SpatialDim> lengths,
                        std::array<int, SpatialDim> mode, double safety, const std::string& name,
                        bool rankRecord, double eps = 0, bool recordOnly = false,
                        bool vacuumSet = false, int splitAxis = -1) {
        const DispersionInput input{cells,      mode, lengths,    safety,   name,
                                    rankRecord, eps,  recordOnly, splitAxis};
        auto fields = makeBoxFieldContainer(cells, lengths, input.splitAxis());
        initializeFieldStorage(fields);
        auto handle = coreTest::makeTestSolver<PotentialSolverBase>(
            fields.getJ(), fields.getE(), fields.getB(), EBReconstructionConvention::R1,
            SpatialDiscretization::Collocated, eps);
        auto& solver = *handle;
        if (safety > 0)
            solver.setSafetyFactor(safety);
        if (vacuumSet)
            solver.setTimeStepUnchecked(recordTimeStep(eps, fields.getHr(), safety, true));
        const auto expected = prediction(solver, input);
        PhaseMeasurement measured;
        ASSERT_NO_FATAL_FAILURE(measurePhase(solver, expected, measured));
        if (!recordOnly)
            requirePhase(solver, input, expected, measured);
        writeDispersion(input, dispersionRecord(solver, input, expected, measured));
    }

    void standardDirection(int n, int direction, std::array<int, SpatialDim> mode, double safety,
                           bool ranks, double eps = 0, bool recordOnly = false,
                           bool vacuumSet = false) {
        std::array<int, SpatialDim> cells      = {n, n, n};
        std::array<double, SpatialDim> lengths = {1, 1, 1};
        if (direction < static_cast<int>(SpatialDim))
            for (int d = 0; d < static_cast<int>(SpatialDim); ++d)
                if (d != direction) {
                    cells[d]   = 8;
                    lengths[d] = 8.0 / n;
                }
        dispersionCase<StandardSolverPeriodic>(
            cells, lengths, mode, safety,
            (eps > 0 ? "Medium-" + std::to_string(eps) + "-" : "") + "Std-" + std::to_string(n)
                + "-" + std::to_string(direction) + "-" + std::to_string(safety)
                + (vacuumSet ? "-vacuum" : ""),
            ranks, eps, recordOnly, vacuumSet);
    }
    void standardDispersion(int n, bool ranks) {
        const std::array<std::array<int, SpatialDim>, 5> modes = {
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}, {1, 1, 1}}};
        for (int direction = 0; direction < 5; ++direction)
            for (double safety : {0.0, 0.98})
                standardDirection(n, direction, modes[direction], safety, ranks);
    }
    void nonstandardDispersion(int nz, bool ranks) {
        const double h = 1.0 / nz;
        for (double safety : {1.0, 0.99, 0.98})
            dispersionCase<NonStandardSolverPeriodic>(
                {8, 8, nz}, {12 * h, 12 * h, 1}, {0, 0, 1}, safety,
                "NS-axis-" + std::to_string(nz) + "-" + std::to_string(safety), ranks);
        for (int length : {2, 4, 8})
            for (double safety : {1.0, 0.98})
                dispersionCase<NonStandardSolverPeriodic>(
                    {length * nz * 2 / 3, 4, nz}, {double(length), 6 * h, 1}, {1, 0, 1}, safety,
                    "NS-near-" + std::to_string(nz) + "-" + std::to_string(length) + "-"
                        + std::to_string(safety),
                    ranks);
    }
    TEST(Dispersion, StandardRankCases) {
        standardDispersion(32, true);
    }
    TEST(Dispersion, NonStandardRankCases) {
        nonstandardDispersion(24, true);
    }
    TEST(RecordDispersion, Standard) {
        standardDispersion(64, false);
    }
    TEST(RecordDispersion, NonStandard) {
        for (int n : {24, 48, 96})
            nonstandardDispersion(n, false);
    }

    /** @brief Check the discrete phase for axial, face-diagonal and body-diagonal medium modes.
     * @details standardDirection keeps eight cells in directions where an axial
     * wave is constant, while preserving the same spacing as the full cubic cases.
     * ranks requests a separate cumulative-phase record for rank-count comparison.
     */
    void mediumDispersion(int n, bool ranks) {
        const std::array<std::array<int, 3>, 5> modes = {
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}, {1, 1, 1}}};
        for (double eps : {2.13, 4.0})
            for (int direction = 0; direction < 5; ++direction)
                for (double safety : {std::sqrt(3.0) / 2, 0.98})
                    standardDirection(n, direction, modes[direction], safety, ranks, eps);
    }
    TEST(MediumDispersion, RankCases) {
        mediumDispersion(32, true);
    }
    TEST(MediumDispersionLong, ResolvedModes) {
        mediumDispersion(64, false);
    }
    // Pair the medium-bound and vacuum-bound steps at each resolution and direction.
    // Record-only cases omit the measured-versus-symbol agreement assertion.
    TEST(MediumDispersionRecord, WavelengthError) {
        const std::array<std::array<int, 3>, 3> modes = {{{1, 0, 0}, {1, 1, 0}, {1, 1, 1}}};
        for (int n : {16, 32, 64, 128})
            for (int d = 0; d < 3; ++d)
                for (double safety : {0.9, std::sqrt(3.0) / 2})
                    for (bool vacuumSet : {false, true})
                        standardDirection(n, d == 0 ? 0 : d + 2, modes[d], safety, false, 2.13,
                                          true, vacuumSet);
    }

    /** @brief Independently held anisotropic single-root phase and rank comparison.
     * @details The serial reference repeats the same measured cumulative phase on
     * MPI_COMM_SELF. It is deliberately absent from ordinary test registration.
     * The mesh is the finest held continuum-wave mesh, but initialization uses
     * the discrete prediction to isolate one travelling root. First require
     * agreement with that independent symbol, then compare the two measured
     * phases across rank counts. The root-only replay uses exactly the distributed
     * step; serialFrequency also means radians per step, as in PhaseMeasurement.
     */
    TEST(MediumHoldout, DISABLED_Phase) {
        const DispersionInput input{{120, 144, 160}, {-1, 2, 3}, {1, 1, 1}, 0.9, "medium-held",
                                    false,           3.7,        false,     1};
        auto fields = makeBoxFieldContainer(input.cells_m, input.lengths_m, input.splitAxis());
        initializeFieldStorage(fields);
        auto handle = coreTest::makeTestSolver<StandardSolverPeriodic>(
            fields.getJ(), fields.getE(), fields.getB(), EBReconstructionConvention::R2,
            SpatialDiscretization::Collocated, input.epsilon_m);
        auto& solver = *handle;
        solver.setSafetyFactor(input.safety_m);
        const auto expected = prediction(solver, input);
        PhaseMeasurement measured;
        measurePhase(solver, expected, measured);
        requirePhase(solver, input, expected, measured);
        double serialFrequency = 0;
        if (ippl::Comm->rank() == 0) {
            SerialFields serial(solver.A_n);
            auto reference = coreTest::makeTestSolver<StandardSolverPeriodic>(
                serial.source_m, serial.electric_m, serial.magnetic_m,
                EBReconstructionConvention::R2, SpatialDiscretization::Collocated, input.epsilon_m);
            reference->setTimeStepUnchecked(solver.getDt());
            PhaseMeasurement result;
            measurePhase(*reference, expected, result);
            serialFrequency = result.frequency();
        }
        MPI_Bcast(&serialFrequency, 1, MPI_DOUBLE, 0, fields.getFL().comm);
        EXPECT_LE(std::abs(measured.frequency() / serialFrequency - 1), 1e-12);
        auto record                   = dispersionRecord(solver, input, expected, measured);
        record["serial_omega"]        = serialFrequency;
        record["rank_phase_relative"] = std::abs(measured.frequency() / serialFrequency - 1);
        writeDispersion(input, record);
    }
}  // namespace chdr::solver::test
