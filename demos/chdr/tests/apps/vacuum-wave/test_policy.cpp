/** @file test_policy.cpp
 * @brief Test the acceptance policy of chdr-vacuum-wave (VacuumChecks.h) with controlled
 * measurements.
 * @details These tests check which tolerance validateStep() applies (coarse or fine,
 * at a checkpoint or at every step) and which inputs make validateStep() throw
 * std::runtime_error or acceptedPhase() return false. They do not run the solver.
 * Thresholds are our fixed benchmark contract in apps/vacuum-wave/README.md,
 * based on the references from Fallahi, MITHRA 2.0, Eqs. (3.6)-(3.20)
 * [fallahi2020mithra20fullwavesimulation]: discrete error 1e-10, zero error 1e-12,
 * phase amplitude above 0.5, phase increment strictly between \f$-\pi\f$ and 0, speed
 * within 0.002 of 1, and at checkpoints relative L2 errors 0.01 (coarse) and 0.003 (fine).
 * No expected decision calls a duplicated local implementation of the policy.
 *
 * CTest: chdr.vacuum.unit.serial and chdr.vacuum.unit.mpi2.
 */
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

#include "VacuumChecks.h"

namespace {
    using namespace chdr::apps::vacuumWave;

    Diagnostics validDiagnostics(int step = 1) {
        Diagnostics value;
        value.step_m   = step;
        value.dt_m     = 1.0 / 128.0;
        value.timeA_m  = step * value.dt_m;
        value.timeE_m  = value.timeA_m - value.dt_m / 2.0;
        value.timeB_m  = value.timeA_m;
        value.count_m  = 4096;
        value.finite_m = true;
        return value;
    }

    PhaseSample validPhase(bool hasSpeed = true) {
        PhaseSample value;
        value.amplitude_m = 1.0;
        value.increment_m = hasSpeed ? -0.05 : 0.0;
        value.speed_m     = hasSpeed ? 1.0 : 0.0;
        value.hasSpeed_m  = hasSpeed;
        value.valid_m     = true;
        return value;
    }

    /// A valid step-0 sample without speed and a valid step-1 sample with speed pass: coarse at
    /// a checkpoint, coarse between checkpoints and fine at a checkpoint.
    TEST(VacuumPolicy, AcceptsConsistentInitialAndEvolvedSamples) {
        EXPECT_NO_THROW(validateStep(validDiagnostics(0), validPhase(false), true, false));
        EXPECT_NO_THROW(validateStep(validDiagnostics(), validPhase(), false, false));
        EXPECT_NO_THROW(validateStep(validDiagnostics(), validPhase(), true, true));
    }

    /// A relative E or B error of 0.005 passes the coarse checkpoint (0.01), fails the fine
    /// checkpoint (0.003) and passes between checkpoints.
    TEST(VacuumPolicy, SelectsCoarseAndFineToleranceForEachField) {
        for (double Diagnostics::* member :
             {&Diagnostics::relativeE_m, &Diagnostics::relativeB_m}) {
            auto value    = validDiagnostics();
            value.*member = 0.005;
            EXPECT_NO_THROW(validateStep(value, validPhase(), true, false));
            EXPECT_THROW(validateStep(value, validPhase(), true, true), std::runtime_error);
            EXPECT_NO_THROW(validateStep(value, validPhase(), false, true));
        }
    }

    /// A relative E error equal to the tolerance passes; the next larger double fails at a
    /// checkpoint and passes between checkpoints, for coarse and fine.
    TEST(VacuumPolicy, EnforcesContinuumLimitsOnlyAtCheckpoints) {
        for (bool fine : {false, true}) {
            auto value        = validDiagnostics();
            value.relativeE_m = fine ? FineFieldTolerance : CoarseFieldTolerance;
            EXPECT_NO_THROW(validateStep(value, validPhase(), true, fine));
            value.relativeE_m =
                std::nextafter(value.relativeE_m, std::numeric_limits<double>::infinity());
            EXPECT_THROW(validateStep(value, validPhase(), true, fine), std::runtime_error);
            EXPECT_NO_THROW(validateStep(value, validPhase(), false, fine));
        }
    }

    /// Discrete error 1e-10 and zero error 1e-12 pass; twice either value fails, also
    /// between checkpoints.
    TEST(VacuumPolicy, RejectsDiscreteAndZeroChannelErrorsEveryStep) {
        auto value            = validDiagnostics();
        value.discreteError_m = DiscreteTolerance;
        value.zeroError_m     = ZeroTolerance;
        EXPECT_NO_THROW(validateStep(value, validPhase(), false, false));
        value.discreteError_m = 2.0 * DiscreteTolerance;
        EXPECT_THROW(validateStep(value, validPhase(), false, false), std::runtime_error);
        value.discreteError_m = 0.0;
        value.zeroError_m     = 2.0 * ZeroTolerance;
        EXPECT_THROW(validateStep(value, validPhase(), false, true), std::runtime_error);
    }

    /// A measurement with finite_m false or a phase sample with valid_m false fails.
    TEST(VacuumPolicy, RejectsInvalidMeasurementAndPhaseFlags) {
        auto value     = validDiagnostics();
        value.finite_m = false;
        EXPECT_THROW(validateStep(value, validPhase(), false, false), std::runtime_error);
        auto phase    = validPhase();
        phase.valid_m = false;
        EXPECT_THROW(validateStep(validDiagnostics(), phase, false, false), std::runtime_error);
    }

    /// acceptedPhase() requires no speed at step 0 and a speed after step 0.
    TEST(VacuumPolicy, RequiresCorrectInitialSpeedAvailability) {
        EXPECT_FALSE(acceptedPhase(validDiagnostics(0), validPhase(true)));
        EXPECT_FALSE(acceptedPhase(validDiagnostics(1), validPhase(false)));
        EXPECT_TRUE(acceptedPhase(validDiagnostics(0), validPhase(false)));
    }

    /// The amplitude 0.5 (the minimum itself) fails; the phase increments 0, 0.1 and
    /// \f$-\pi\f$ fail.
    TEST(VacuumPolicy, RejectsBoundaryAmplitudeAndWrongDirection) {
        auto phase        = validPhase();
        phase.amplitude_m = MinimumPhaseAmplitude;
        EXPECT_THROW(validateStep(validDiagnostics(), phase, false, false), std::runtime_error);
        phase = validPhase();
        for (double increment : {0.0, 0.1, -Kokkos::numbers::pi_v<double>}) {
            phase.increment_m = increment;
            EXPECT_THROW(validateStep(validDiagnostics(), phase, false, false), std::runtime_error);
        }
    }

    /// Speed 0.9 fails although the increment is negative; a NaN speed fails acceptedPhase().
    TEST(VacuumPolicy, RejectsWrongSpeedEvenWithValidDirection) {
        auto phase    = validPhase();
        phase.speed_m = 0.9;
        EXPECT_LT(phase.increment_m, 0.0);
        EXPECT_THROW(validateStep(validDiagnostics(), phase, false, false), std::runtime_error);
        phase.speed_m = std::numeric_limits<double>::quiet_NaN();
        EXPECT_FALSE(acceptedPhase(validDiagnostics(), phase));
    }
}  // namespace
