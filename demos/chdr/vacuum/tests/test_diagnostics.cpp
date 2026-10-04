/** @file test_diagnostics.cpp
 * @brief Regression and fault tests for owned-cell vacuum diagnostics.
 * @details CPU checks for the default coarse mode on one/two ranks. The independent
 * host recurrence and derivative references follow our substitution in Fallahi,
 * MITHRA 2.0, arXiv:2009.13645v1, Eqs. (3.6)-(3.17), printed pp.12-14.
 * Production referencePotential, evaluate_EB and the closed-form helper are not
 * used to create expected arrays. This does not replace the full propagation test.
 */
#include "Ippl.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <memory>

#include "VacuumDiagnostics.h"
#include "VacuumSetup.h"

namespace {
    using namespace chdr::vacuum;
    using Potential             = ippl::Vector<double, 4>;
    using FieldValue            = ippl::Vector<double, 3>;
    using Solver                = VacuumSolver;
    constexpr double WaveNumber = 2.0 * Kokkos::numbers::pi_v<double>;
    using Complex               = std::complex<double>;

    /// Test-only independent scalar recurrence, from Fallahi (3.11)-(3.17).
    std::pair<Complex, Complex> scalarState(double dt, double h, int step) {
        Complex previous         = std::polar(1.0, WaveNumber * dt);
        Complex current          = 1.0;
        const double ratio       = dt / h;
        const double sine        = std::sin(WaveNumber * h / 2.0);
        const double coefficient = 2.0 - 4.0 * ratio * ratio * sine * sine;
        for (int n = 0; n < step; ++n) {
            const Complex next = coefficient * current - previous;
            previous           = current;
            current            = next;
        }
        return {current, previous};
    }

    /// Host fixture fill deliberately leaves every halo NaN; production must ignore it.
    template <class Field, class Function>
    void fillOwned(Field& field, Function value) {
        auto host = field.getHostMirror();
        Kokkos::deep_copy(host,
                          typename Field::value_type(std::numeric_limits<double>::quiet_NaN()));
        const auto& local = field.getLayout().getLocalNDIndex();
        const int ghost   = field.getNghost();
        const auto& mesh  = field.get_mesh();
        for (int i = 0; i < static_cast<int>(local[0].length()); ++i)
            for (int j = 0; j < static_cast<int>(local[1].length()); ++j)
                for (int k = 0; k < static_cast<int>(local[2].length()); ++k) {
                    const double z = mesh.getOrigin()[2]
                                     + (local[2].first() + k + 0.5) * mesh.getMeshSpacing()[2];
                    host(i + ghost, j + ghost, k + ghost) = value(z);
                }
        Kokkos::deep_copy(field.getView(), host);
    }

    /// Test potentials come from the independently advanced scalar recurrence.
    Potential potentialValue(double z, Complex mode) {
        return Potential(0.0, std::imag(std::polar(1.0, WaveNumber * z) * mode) / WaveNumber, 0.0,
                         0.0);
    }

    /// Exact stencil reconstruction of the synthetic mode; no solver calls or oracle reuse.
    void seedFields(Fields& fields, Complex current, Complex previous, double dt) {
        const double h = fields.getMesh().getMeshSpacing()[2];
        fillOwned(fields.getJ(), [](double) {
            return Potential(0.0);
        });
        fillOwned(fields.getE(), [=](double z) {
            const double value = -std::imag(std::polar(1.0, WaveNumber * z) * (current - previous));
            return FieldValue(value / (WaveNumber * dt), 0.0, 0.0);
        });
        fillOwned(fields.getB(), [=](double z) {
            const double value = std::real(std::polar(1.0, WaveNumber * z) * current);
            return FieldValue(0.0, value * std::sin(WaveNumber * h) / (WaveNumber * h), 0.0);
        });
    }

    void seed(Solver& solver, Fields& fields, int step) {
        const double dt = solver.getDt();
        const auto [current, previous] =
            scalarState(dt, fields.getMesh().getMeshSpacing()[2], step);
        fillOwned(solver.A_n, [=](double z) {
            return potentialValue(z, current);
        });
        fillOwned(solver.A_nm1, [=](double z) {
            return potentialValue(z, previous);
        });
        fillOwned(solver.A_np1, [=](double z) {
            return step == 0 ? Potential(0.0) : potentialValue(z, current);
        });
        seedFields(fields, current, previous, dt);
    }

    template <class Field>
    void defectOnLastRank(Field& field, int component, double value) {
        const auto& comm = field.getLayout().comm;
        if (comm.rank() != comm.size() - 1)
            return;
        auto host = field.getHostMirror();
        Kokkos::deep_copy(host, field.getView());
        const int g              = field.getNghost();
        host(g, g, g)[component] = value;
        Kokkos::deep_copy(field.getView(), host);
    }

    void checkInitial(Solver& solver, Fields& fields) {
        seed(solver, fields, 0);
        const auto result = chdr::vacuum::measure(solver, fields, 0);
        const double x    = WaveNumber * solver.getDt() / 2.0;
        const double q    = WaveNumber * fields.getMesh().getMeshSpacing()[2];
        EXPECT_TRUE(result.finite_m && result.count_m == 4096) << "initial validity/count";
        EXPECT_TRUE(result.discreteError_m < 1e-12 && result.zeroError_m == 0.0)
            << "initial oracle/zero";
        EXPECT_TRUE(std::abs(result.relativeE_m - (1.0 - std::sin(x) / x)) < 1e-13)
            << "initial E norm";
        EXPECT_TRUE(std::abs(result.relativeB_m - (1.0 - std::sin(q) / q)) < 1e-13)
            << "initial B norm";
        EXPECT_TRUE(std::abs(result.phaseReal_m) < 1e-13) << "initial Fourier real";
        EXPECT_TRUE(std::abs(result.phaseImag_m + 1.0) < 1e-13) << "initial Fourier imaginary";
        EXPECT_TRUE(result.timeE_m == -solver.getDt() / 2.0 && result.timeB_m == 0.0)
            << "raw times";
    }

    void checkSourceDefect(Solver& solver, Fields& fields) {
        seed(solver, fields, 0);
        defectOnLastRank(fields.getJ(), 2, 1e-5);
        const auto result = chdr::vacuum::measure(solver, fields, 0);
        EXPECT_TRUE(result.finite_m && result.sourceError_m == 1e-5 && result.zeroError_m == 1e-5)
            << "source fault";
    }

    template <class Field>
    void checkFiniteDefect(Solver& solver, Fields& fields, Field& target, int component,
                           bool inactive) {
        seed(solver, fields, 0);
        defectOnLastRank(target, component, 0.4);
        const auto result = chdr::vacuum::measure(solver, fields, 0);
        EXPECT_TRUE(result.finite_m && result.discreteError_m > 1e-10)
            << "finite field/history fault";
        EXPECT_TRUE(inactive ? result.zeroError_m > 1e-12 : result.zeroError_m == 0.0)
            << "zero-channel fault";
    }

    void checkFiniteDefects(Solver& solver, Fields& fields) {
        checkFiniteDefect(solver, fields, solver.A_n, 1, false);
        checkFiniteDefect(solver, fields, solver.A_nm1, 1, false);
        checkFiniteDefect(solver, fields, solver.A_np1, 1, false);
        checkFiniteDefect(solver, fields, fields.getE(), 0, false);
        checkFiniteDefect(solver, fields, fields.getB(), 1, false);
        checkFiniteDefect(solver, fields, solver.A_n, 0, true);
        checkFiniteDefect(solver, fields, fields.getE(), 1, true);
        checkFiniteDefect(solver, fields, fields.getE(), 2, true);
        checkFiniteDefect(solver, fields, fields.getB(), 0, true);
        checkFiniteDefect(solver, fields, fields.getB(), 2, true);
        checkSourceDefect(solver, fields);
    }

    template <class Field>
    void checkNonfinite(Solver& solver, Fields& fields, Field& target, int component) {
        seed(solver, fields, 0);
        defectOnLastRank(target, component, std::numeric_limits<double>::quiet_NaN());
        EXPECT_TRUE(!chdr::vacuum::measure(solver, fields, 0).finite_m)
            << "non-root NaN propagation";
    }

    void checkNonfiniteFields(Solver& solver, Fields& fields) {
        checkNonfinite(solver, fields, solver.A_n, 1);
        checkNonfinite(solver, fields, solver.A_nm1, 1);
        checkNonfinite(solver, fields, solver.A_np1, 1);
        checkNonfinite(solver, fields, fields.getE(), 1);
        checkNonfinite(solver, fields, fields.getB(), 1);
        checkNonfinite(solver, fields, fields.getJ(), 0);
    }

    void checkStorageMismatch(Solver& solver, Fields& fields) {
        seed(solver, fields, 0);
        const auto original = solver.A_np1.getView();
        const auto& comm    = fields.getFL().comm;
        if (comm.rank() == comm.size() - 1)
            solver.A_np1.resize(original.extent(0) + 1, original.extent(1), original.extent(2));
        const auto result = chdr::vacuum::measure(solver, fields, 0);
        EXPECT_TRUE(!result.finite_m && result.count_m == 0)
            << "collective storage mismatch preflight";
        solver.A_np1.getView() = original;
    }

    void checkMeasuredPhase(const chdr::vacuum::Diagnostics& result,
                            const chdr::vacuum::PhaseSample& sample, double& accumulated) {
        const auto [current, previous] = scalarState(result.dt_m, 1.0 / 64.0, result.step_m);
        const Complex expected         = Complex(0.0, -1.0) * current;
        EXPECT_TRUE(std::abs(result.phaseReal_m - expected.real()) < 1e-12)
            << "Fourier real recurrence";
        EXPECT_TRUE(std::abs(result.phaseImag_m - expected.imag()) < 1e-12)
            << "Fourier imaginary recurrence";
        EXPECT_TRUE(sample.valid_m && sample.hasSpeed_m == (result.step_m > 0))
            << "phase structure";
        if (result.step_m == 0)
            return;
        const double increment = std::arg(current / previous);
        accumulated += increment;
        const double speed = -accumulated / (WaveNumber * result.timeA_m);
        EXPECT_TRUE(std::abs(sample.increment_m - increment) < 1e-12)
            << "independent phase increment";
        EXPECT_TRUE(std::abs(sample.speed_m - speed) < 1e-12)
            << "independent cumulative phase speed";
    }

    void checkSteps(Solver& solver, Fields& fields) {
        chdr::vacuum::PhaseTracker phase;
        double expectedPhase = 0.0;
        for (int step = 0; step <= 128; ++step) {
            seed(solver, fields, step);
            const auto result = chdr::vacuum::measure(solver, fields, step);
            EXPECT_TRUE(result.finite_m && result.discreteError_m <= 1e-10)
                << "independent recurrence oracle";
            const auto sample = phase.update(result);
            checkMeasuredPhase(result, sample, expectedPhase);
        }
    }

    chdr::vacuum::Diagnostics phaseInput(chdr::vacuum::Diagnostics sample, int step, double speed,
                                         double amplitude) {
        sample.step_m      = step;
        sample.timeA_m     = step * sample.dt_m;
        const double angle = WaveNumber * speed * sample.timeA_m;
        sample.phaseReal_m = -amplitude * std::sin(angle);
        sample.phaseImag_m = -amplitude * std::cos(angle);
        return sample;
    }

    chdr::vacuum::PhaseSample trackOne(const chdr::vacuum::Diagnostics& initial, double speed,
                                       double amplitude) {
        chdr::vacuum::PhaseTracker phase;
        const auto first = phase.update(phaseInput(initial, 0, speed, amplitude));
        EXPECT_TRUE(first.valid_m && !first.hasSpeed_m && first.speed_m == 0.0)
            << "initial speed skipped";
        return phase.update(phaseInput(initial, 1, speed, amplitude));
    }

    void checkPhasePolicySeparation(const chdr::vacuum::Diagnostics& initial) {
        const auto slow = trackOne(initial, 0.9, 1.0);
        EXPECT_TRUE(slow.valid_m && slow.increment_m < 0.0) << "valid direction at wrong speed";
        EXPECT_TRUE(std::abs(slow.speed_m - 0.9) < 1e-12) << "wrong speed measured accurately";
        EXPECT_TRUE(std::abs(slow.speed_m - 1.0) > 0.002) << "speed policy belongs to caller";
        const auto backwards = trackOne(initial, -1.0, 1.0);
        EXPECT_TRUE(backwards.valid_m && backwards.increment_m > 0.0)
            << "direction policy belongs to caller";
        const auto weak = trackOne(initial, 1.0, 0.5);
        EXPECT_TRUE(weak.valid_m && std::abs(weak.amplitude_m - 0.5) < 1e-15)
            << "amplitude policy belongs to caller";
        const auto tiny = trackOne(initial, 1.0, 1e-250);
        EXPECT_TRUE(tiny.valid_m && std::abs(tiny.speed_m - 1.0) < 1e-12)
            << "tiny amplitude phase measurement";
    }

    void checkPhaseStructure(const chdr::vacuum::Diagnostics& initial) {
        chdr::vacuum::PhaseTracker phase;
        EXPECT_TRUE(phase.update(initial).valid_m) << "initial phase accepted";
        EXPECT_TRUE(!phase.update(initial).valid_m) << "repeated step rejected";
        EXPECT_TRUE(!phase.update(phaseInput(initial, 2, 1.0, 1.0)).valid_m)
            << "skipped step rejected";
        auto next = phaseInput(initial, 1, 1.0, 1.0);
        next.timeA_m += initial.dt_m / 2.0;
        EXPECT_TRUE(!phase.update(next).valid_m) << "inconsistent time rejected";
        next = phaseInput(initial, 1, 1.0, 1.0);
        next.dt_m *= 2.0;
        EXPECT_TRUE(!phase.update(next).valid_m) << "changed dt rejected";
        EXPECT_TRUE(!phase.update(phaseInput(initial, 1, 1.0, 0.0)).valid_m)
            << "zero amplitude undefined";
        next             = phaseInput(initial, 1, 1.0, 1.0);
        next.phaseReal_m = std::numeric_limits<double>::quiet_NaN();
        EXPECT_TRUE(!phase.update(next).valid_m) << "nonfinite phase rejected";
    }

    void checkPhaseFailures(Solver& solver, Fields& fields) {
        seed(solver, fields, 0);
        const auto result = chdr::vacuum::measure(solver, fields, 0);
        checkPhasePolicySeparation(result);
        checkPhaseStructure(result);
        EXPECT_TRUE(!chdr::vacuum::measure(solver, fields, -1).finite_m) << "negative step";
    }

    /** @brief Synthetic owned data isolate diagnostic behavior from solver evolution.
     * @details Histories use an independent host recurrence from Fallahi (3.11)-(3.17).
     * Halos remain NaN. All mutations are on the last rank; measured errors are
     * globally reduced before assertions, so a two-rank run checks non-root failures.
     */
    class VacuumDiagnosticsTest : public ::testing::Test {
    protected:
        void SetUp() override {
            checkLocalDomain(fields_m.getFL());
            initializeFieldStorage(fields_m);
            solver_m = std::make_unique<Solver>(fields_m.getJ(), fields_m.getE(), fields_m.getB());
        }
        Fields fields_m = makePeriodicFieldContainer();
        std::unique_ptr<Solver> solver_m;
    };

    TEST_F(VacuumDiagnosticsTest, MatchesInitialReferenceAndOwnership) {
        checkInitial(*solver_m, fields_m);
    }

    TEST_F(VacuumDiagnosticsTest, DetectsFiniteFieldAndHistoryDefects) {
        checkFiniteDefects(*solver_m, fields_m);
    }

    TEST_F(VacuumDiagnosticsTest, RejectsNonfiniteOwnedValuesOnAnyRank) {
        checkNonfiniteFields(*solver_m, fields_m);
    }

    TEST_F(VacuumDiagnosticsTest, RejectsStorageMismatchCollectively) {
        checkStorageMismatch(*solver_m, fields_m);
    }

    TEST_F(VacuumDiagnosticsTest, RejectsOverflowedReferenceFromFiniteMetadata) {
        seed(*solver_m, fields_m, 0);
        ippl::Vector<double, Dim> origin(0.0);
        origin[2] = std::numeric_limits<double>::max();
        fields_m.getMesh().setOrigin(origin);
        EXPECT_FALSE(measure(*solver_m, fields_m, 0).finite_m);
    }

    TEST_F(VacuumDiagnosticsTest, MatchesIndependentModesAndUnwrappedPhase) {
        checkSteps(*solver_m, fields_m);
    }

    TEST_F(VacuumDiagnosticsTest, SeparatesPhaseStructureFromAcceptancePolicy) {
        checkPhaseFailures(*solver_m, fields_m);
    }
}  // namespace
