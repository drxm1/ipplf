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

    chdr::vacuum::Fields probeGrid() {
        ippl::Vector<double, SpatialDim> lower(0.0), upper(1.0, 1.0, 2.0), spacing(0.125);
        ippl::NDIndex<SpatialDim> domain(ippl::Index(8), ippl::Index(8), ippl::Index(16));
        return {spacing, lower, upper, {false, false, true}, domain, lower, true};
    }

    constexpr double DT = 0.0625;
    const ProbeBox POINT{{0, 0, 0}, {0, 0, 0}};
    const LevelContract R2{EBReconstructionConvention::R2}, R1{EBReconstructionConvention::R1};
    const ProbeChannel ELECTRIC{Quantity::Electric, 0};

    /// Keep the allocated field owner at its original address throughout each interface test.
    struct ProbeFixture {
        chdr::vacuum::Fields fields_m = probeGrid();
        ProbeFixture() { chdr::vacuum::initializeFieldStorage(fields_m); }
        FourierProbe make(LevelContract levels, ProbeBox box, ProbeChannel channel,
                          double frequency, std::int64_t first, std::int64_t count,
                          bool phiIsZero = false) {
            return FourierProbe(fields_m.getE(), levels, DT, box, {channel}, {frequency}, first,
                                count, ProbeWindow::Rectangle, phiIsZero);
        }
    };
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

    /** @brief Check Hann metadata independently of the recorded window choice.
     * @details The sine-square reference is a half-angle rewrite of
     * [heinzel2002SpectrumSpectralDensity, Eq. (17)]; sums follow Eqs. (19)-(20).
     */
    TEST(FourierProbe, HannWeightsAndSums) {
        ProbeFixture fixture;
        for (int count : {1, 2, 3, 5, 4096}) {
            FourierProbe probe(fixture.fields_m.getE(), R2, DT, POINT, {ELECTRIC}, {0.1}, 0, count,
                               ProbeWindow::Hann);
            EXPECT_EQ(probe.metadata().window_m, ProbeWindow::Hann);
            EXPECT_DOUBLE_EQ(probe.metadata().windowSum_m, count == 1 ? 0 : count / 2.0);
            EXPECT_DOUBLE_EQ(probe.metadata().windowSquareSum_m, count == 1   ? 0
                                                                 : count == 2 ? 1
                                                                              : 3 * count / 8.0);
            for (int sample = 0; sample < count; ++sample)
                EXPECT_NEAR(probe.weight(sample), std::pow(std::sin(PI * sample / count), 2),
                            4e-16);
        }
    }

    TEST(FourierProbe, RejectsNonfiniteReportedFrequency) {
        ProbeFixture fixture;
        EXPECT_THROW(
            FourierProbe(fixture.fields_m.getE(), R2, std::numeric_limits<double>::denorm_min(),
                         POINT, {ELECTRIC}, {0.25}, 0, 1),
            std::invalid_argument);
    }
}  // namespace
