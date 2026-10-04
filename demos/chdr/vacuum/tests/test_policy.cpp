/** @file test_policy.cpp
 * @brief Test the real driver acceptance policy using controlled measurements.
 * @details These tests validate gate selection and failure behavior, not the
 * physics solver. Thresholds are our fixed benchmark contract in README.md,
 * based on the references from Fallahi, MITHRA 2.0, Eqs. (3.6)-(3.20).
 * No expected decision calls a duplicated local implementation of the policy.
 */
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

#include "VacuumChecks.h"

namespace {
    using namespace chdr::vacuum;

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

    TEST(VacuumPolicy, AcceptsConsistentInitialAndEvolvedSamples) {
        EXPECT_NO_THROW(validateStep(validDiagnostics(0), validPhase(false), true, false));
        EXPECT_NO_THROW(validateStep(validDiagnostics(), validPhase(), false, false));
        EXPECT_NO_THROW(validateStep(validDiagnostics(), validPhase(), true, true));
    }

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

    TEST(VacuumPolicy, RejectsInvalidMeasurementAndPhaseFlags) {
        auto value     = validDiagnostics();
        value.finite_m = false;
        EXPECT_THROW(validateStep(value, validPhase(), false, false), std::runtime_error);
        auto phase    = validPhase();
        phase.valid_m = false;
        EXPECT_THROW(validateStep(validDiagnostics(), phase, false, false), std::runtime_error);
    }

    TEST(VacuumPolicy, RequiresCorrectInitialSpeedAvailability) {
        EXPECT_FALSE(acceptedPhase(validDiagnostics(0), validPhase(true)));
        EXPECT_FALSE(acceptedPhase(validDiagnostics(1), validPhase(false)));
        EXPECT_TRUE(acceptedPhase(validDiagnostics(0), validPhase(false)));
    }

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

    TEST(VacuumPolicy, RejectsWrongSpeedEvenWithValidDirection) {
        auto phase    = validPhase();
        phase.speed_m = 0.9;
        EXPECT_LT(phase.increment_m, 0.0);
        EXPECT_THROW(validateStep(validDiagnostics(), phase, false, false), std::runtime_error);
        phase.speed_m = std::numeric_limits<double>::quiet_NaN();
        EXPECT_FALSE(acceptedPhase(validDiagnostics(), phase));
    }
}  // namespace
