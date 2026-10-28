/** @file test_probe.cpp
 * @brief Minimal closed-form and interface checks for the time-labelled Fourier probe.
 */
#include <gtest/gtest.h>
#include <limits>

#include "core/base/SpatialDimension.h"
#include "core/solver/BoxSetup.h"
#include "opt-in/diagnostics/FourierProbe.h"

static_assert(chdr::SpatialDim == 3,
              "test_probe.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr;
    using namespace chdr::solver;
    using namespace chdr::diagnostics;
    constexpr double PI = Kokkos::numbers::pi_v<double>;

    /// The periodic 8 by 8 by 16 box of unit spacing 1/8, decomposed along z.
    ChDRFieldContainer probeGrid() {
        return makeBoxFieldContainer({8, 8, 16}, {1.0, 1.0, 2.0}, 2);
    }

    constexpr double DT = 0.0625;
    const ProbeBox POINT{{0, 0, 0}, {0, 0, 0}};
    const LevelContract R2{EBReconstructionConvention::R2}, R1{EBReconstructionConvention::R1};
    const ProbeChannel ELECTRIC{Quantity::Electric, 0};

    /// Keep the allocated field owner at its original address throughout each interface test.
    struct ProbeFixture {
        ChDRFieldContainer fields_m = probeGrid();
        ProbeFixture() { initializeFieldStorage(fields_m); }
        FourierProbe make(LevelContract levels, ProbeBox box, ProbeChannel channel,
                          double frequency, std::int64_t first, std::int64_t count,
                          bool phiIsZero = false) {
            return FourierProbe(fields_m.getE(), levels, DT, box, {channel}, {frequency}, first,
                                count, ProbeWindow::Rectangle, phiIsZero);
        }
        void sample(FourierProbe& probe, int step) {
            probe.sample(step, fields_m.getE(), fields_m.getB(), fields_m.getB(), fields_m.getJ());
        }
    };

    void checkInvalidGeometryAndChannels(ProbeFixture& fixture) {
        EXPECT_THROW(fixture.make(R2, {{-1, 0, 0}, {0, 0, 0}}, ELECTRIC, 0.1, 0, 10),
                     std::invalid_argument);
        EXPECT_THROW(fixture.make(R2, {{0, 0, 0}, {8, 0, 0}}, ELECTRIC, 0.1, 0, 10),
                     std::invalid_argument);
        EXPECT_THROW(fixture.make(R2, POINT, {Quantity::Charge, 0}, 0.1, 0, 10),
                     std::invalid_argument);
        EXPECT_THROW(fixture.make(R2, POINT, {Quantity::Phi, 1}, 0.1, 0, 10),
                     std::invalid_argument);
        EXPECT_THROW(fixture.make(R2, POINT, {Quantity::Electric, SpatialDim}, 0.1, 0, 10),
                     std::invalid_argument);
        EXPECT_THROW(fixture.make(R1, POINT, ELECTRIC, 0.1, 0, 10), std::exception);
        EXPECT_NO_THROW(fixture.make(R1, POINT, ELECTRIC, 0.1, 0, 10, true));
    }

    void checkInvalidTimeGate(ProbeFixture& fixture) {
        for (double frequency : {0.0, 0.5, std::numeric_limits<double>::quiet_NaN()})
            EXPECT_THROW(fixture.make(R2, POINT, ELECTRIC, frequency, 0, 10),
                         std::invalid_argument);
        EXPECT_THROW(fixture.make(R2, POINT, ELECTRIC, 0.1, -1, 10), std::invalid_argument);
        EXPECT_THROW(fixture.make(R2, POINT, ELECTRIC, 0.1, 0, 0), std::invalid_argument);
    }

    /// Samples outside the gate are ignored; samples inside it must be complete and ordered.
    void checkGateSampling(ProbeFixture& fixture) {
        auto probe = fixture.make(R2, POINT, ELECTRIC, 0.1, 3, 2);
        EXPECT_THROW(probe.gather(), std::logic_error);
        EXPECT_NO_THROW(fixture.sample(probe, 0));
        EXPECT_THROW(fixture.sample(probe, 4), std::logic_error);
        fixture.sample(probe, 3);
        EXPECT_THROW(fixture.sample(probe, 3), std::logic_error);
        fixture.sample(probe, 4);
        EXPECT_NO_THROW(fixture.sample(probe, 5));
        EXPECT_NO_THROW(probe.gather());
    }

    TEST(FourierProbe, ValidatesConfiguration) {
        ProbeFixture fixture;
        checkInvalidGeometryAndChannels(fixture);
        checkInvalidTimeGate(fixture);
        checkGateSampling(fixture);
    }

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

    /// Different per-rank boxes deliberately make one output index lack an owner.
    TEST(FourierProbe, OwnershipFailureIsCollective) {
        if (ippl::Comm->size() < 2)
            GTEST_SKIP() << "Ownership disagreement requires multiple ranks";
        ProbeFixture fixture;
        const auto& local = fixture.fields_m.getFL().getLocalNDIndex();
        const int outside = local[2].first() == 0 ? 15 : 0;
        auto probe = fixture.make(R2, {{0, 0, outside}, {0, 0, outside}}, ELECTRIC, 0.1, 0, 1);
        fixture.sample(probe, 0);
        EXPECT_THROW(probe.gather(), std::runtime_error);
        int reached = 1, ranks = 0;
        MPI_Allreduce(&reached, &ranks, 1, MPI_INT, MPI_SUM, fixture.fields_m.getFL().comm);
        EXPECT_EQ(ranks, ippl::Comm->size());
    }
}  // namespace
