/** @file test_probe.cpp
 * @brief Minimal closed-form and interface checks for the time-labelled Fourier probe.
 *
 * TODO: REVIEW.
 */
#include <gtest/gtest.h>
#include <limits>

#include "FourierProbe.h"
#include "SpatialDimension.h"

static_assert(chdr::solver::SpatialDim == 3,
              "test_probe.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;
    constexpr double PI = Kokkos::numbers::pi_v<double>;

    /** @brief Independent integer remainder for an exactly represented dyadic frequency.
     * @details Binary64 frequencies are integer multiples of powers of two. Repeated
     * modular doubling forms the exact product remainder without a floating product;
     * this is an integer-arithmetic reference for FourierProbe::phase's FMA correction.
     */
    double exactPhase(double frequency, std::uint64_t halfStep) {
        int exponent;
        const auto mantissa    = std::uint64_t(std::ldexp(std::frexp(frequency, &exponent), 53));
        const auto denominator = std::uint64_t{1} << (54 - exponent);
        std::uint64_t result = 0, term = mantissa;
        for (; halfStep; halfStep >>= 1, term = (2 * term) % denominator)
            if (halfStep & 1)
                result = (result + term) % denominator;
        const double cycles = double(result) / double(denominator);
        return 2 * PI * (cycles - std::nearbyint(cycles));
    }

    TEST(FourierProbe, CorrectsLargeProductRounding) {
        for (double frequency : {0.1, 0.23456789, 0.499999999999})
            for (std::uint64_t step : {(std::uint64_t{1} << 51) - 7, (std::uint64_t{1} << 52) - 5})
                EXPECT_NEAR(FourierProbe::phase(frequency, step), exactPhase(frequency, step),
                            2e-15);
    }

}  // namespace
