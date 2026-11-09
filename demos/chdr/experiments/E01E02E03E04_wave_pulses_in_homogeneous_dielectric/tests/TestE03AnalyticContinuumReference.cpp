/** @file TestE03AnalyticContinuumReference.cpp
 * @brief Check E03 fields by differentiating each potential pair independently.
 */
#include <array>
#include <cmath>
#include <gtest/gtest.h>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03AnalyticContinuumReference.h"

namespace {
    namespace analytic = chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric::
        E03_analytic_continuum_reference;
    using Reference          = analytic::E03AnalyticContinuumReferenceSampler;
    using Point              = std::array<double, 4>;
    constexpr double Epsilon = 2.13;

    /// Evaluate phi or one Cartesian A component at the supplied space-time point.
    double potential(const Reference& reference, const Point& point, unsigned component) {
        return reference(point[0], point[1], point[2], point[3], component);
    }

    /** @brief Numerically differentiate a potential component along x, y, z or time.
     * @details The centred five-point derivative avoids using the implementation's
     * radial derivatives or TE/TM duality identities. Its coefficients follow by
     * Taylor expansion through fourth order [fornberg1988FiniteDifferenceFormulas,
     * section 2].
     */
    double derivative(const Reference& reference, Point point, unsigned component,
                      unsigned direction) {
        constexpr double Step = 1e-4;
        constexpr std::array shifts{-2, -1, 1, 2};
        constexpr std::array weights{1, -8, 8, -1};
        double result = 0;
        for (unsigned sample = 0; sample < shifts.size(); ++sample) {
            Point shifted = point;
            shifted[direction] += shifts[sample] * Step;
            result += weights[sample] * potential(reference, shifted, component);
        }
        return result / (12 * Step);
    }

    /** @brief Apply E=-dA/dt-grad(phi), B=curl(A), and the homogeneous Lorenz gauge.
     * @details These identities are [chew2014generalizedgauge, Eqs. (5), (6), (11)].
     * Both potential pairs are differentiated separately, so a shared duality
     * implementation cannot make the test pass merely by repeating itself.
     */
    void checkFieldsAndGauge(const Reference& reference, const Point& point) {
        double gauge = Epsilon * derivative(reference, point, 0, 3);
        for (unsigned c = 0; c < 3; ++c) {
            const unsigned a = (c + 1) % 3, b = (c + 2) % 3;
            const double electric =
                -derivative(reference, point, c + 1, 3) - derivative(reference, point, 0, c);
            const double magnetic =
                derivative(reference, point, b + 1, a) - derivative(reference, point, a + 1, b);
            EXPECT_NEAR(reference.referenceE(point[0], point[1], point[2], point[3], c), electric,
                        2e-8 * (1 + std::abs(electric)));
            EXPECT_NEAR(reference.referenceB(point[0], point[1], point[2], point[3], c), magnetic,
                        2e-8 * (1 + std::abs(magnetic)));
            gauge += derivative(reference, point, c + 1, c);
        }
        EXPECT_NEAR(gauge, 0, 2e-8);
    }

    /// Require the exact zero continuation for all potentials and both reconstructed fields.
    void checkZero(const Reference& reference, double radius, double time) {
        for (unsigned c = 0; c < 4; ++c)
            EXPECT_DOUBLE_EQ(reference(radius, 0, 0, time, c), 0);
        for (unsigned c = 0; c < 3; ++c) {
            EXPECT_DOUBLE_EQ(reference.referenceE(radius, 0, 0, time, c), 0);
            EXPECT_DOUBLE_EQ(reference.referenceB(radius, 0, 0, time, c), 0);
        }
    }

    /// Compare translated samples without using the radial helper's private state.
    void checkTranslation(const Reference& reference, const Point& point) {
        Reference shifted = reference;
        shifted.centre_m  = 6;
        const double x = point[0], y = point[1], z = point[2], time = point[3];
        for (unsigned c = 0; c < 4; ++c)
            EXPECT_NEAR(shifted(x + 6, y + 6, z + 6, time, c), reference(x, y, z, time, c), 5e-14);
        for (unsigned c = 0; c < 3; ++c) {
            EXPECT_NEAR(shifted.referenceE(x + 6, y + 6, z + 6, time, c),
                        reference.referenceE(x, y, z, time, c), 5e-14);
            EXPECT_NEAR(shifted.referenceB(x + 6, y + 6, z + 6, time, c),
                        reference.referenceB(x, y, z, time, c), 5e-14);
        }
    }

    TEST(E03AnalyticReference, FieldsAndGaugeFromIndependentPotentialDerivatives) {
        const double speed = 1 / std::sqrt(Epsilon);
        const std::array<Point, 3> points{
            {{1.1, 0.9, -1.2, -0.04}, {-1.8, 0.8, 1.1, 0.3}, {1.4, -1.6, 1.9, 1.1}}};
        for (const auto polarization : {analytic::Polarization::TEr, analytic::Polarization::TMr})
            for (const auto& point : points) {
                SCOPED_TRACE(static_cast<int>(polarization));
                SCOPED_TRACE(point[3]);
                checkFieldsAndGauge(Reference{speed, polarization}, point);
            }
    }

    TEST(E03AnalyticReference, CompactSupportAndTranslation) {
        const double speed = 1 / std::sqrt(Epsilon);
        for (const auto polarization : {analytic::Polarization::TEr, analytic::Polarization::TMr}) {
            const Reference reference{speed, polarization};
            for (const double radius : {0.0, 0.5, 1.0, 3.0, 3.5})
                checkZero(reference, radius, 0);
            checkZero(reference, 0, -0.0625);
            checkZero(reference, 0.5, 1);
            checkZero(reference, 5, 1);
            checkTranslation(reference, Point{1.1, 0.9, -1.2, 0.3});
        }
    }
}  // namespace
