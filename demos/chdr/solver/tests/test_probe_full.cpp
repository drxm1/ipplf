/** @file test_probe_full.cpp
 * @brief Closed-form checks of point/plane spectra and their rank-count invariance.
 * @details FourierProbeT8 is the extended probe test suite. Its output tag "T8"
 * identifies per-gate error summaries and the spectra consumed by compare_probe_ranks.py.
 *
 * TODO: REVIEW.
 */
#include <cstdlib>
#include <fstream>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "ChdrSolver.h"
#include "FourierProbe.h"
#include "PeriodicSetup.h"
#include "SpatialDimension.h"

static_assert(
    chdr::solver::SpatialDim == 3,
    "test_probe_full.cpp requires three spatial dimensions for its Cartesian implementation");

namespace {
    using namespace chdr::solver;
    using Complex            = std::complex<long double>;
    constexpr long double PI = Kokkos::numbers::pi_v<long double>;
    constexpr double DT      = 0.0625;
    const std::array<double, 2> TONES{40.0 / 4096.0, 97.37 / 4096.0};

    std::vector<ProbeChannel> channels() {
        return {{Quantity::Electric, 0},
                {Quantity::Electric, 1},
                {Quantity::Electric, 2},
                {Quantity::Magnetic, 0},
                {Quantity::Magnetic, 1},
                {Quantity::Magnetic, 2},
                {Quantity::Phi, 0},
                {Quantity::VectorPotential, 0},
                {Quantity::VectorPotential, 1},
                {Quantity::VectorPotential, 2},
                {Quantity::RawMagnetic, 0}};
    }
    chdr::vacuum::Fields makeFields() {
        ippl::Vector<double, SpatialDim> origin(0.0), upper(1.0, 1.0, 2.0), spacing(0.125);
        ippl::NDIndex<SpatialDim> domain(ippl::Index(8), ippl::Index(8), ippl::Index(16));
        return {spacing, origin, upper, {false, false, true}, domain, origin, true};
    }

    Complex cycle(long double frequency, long double index) {
        return std::polar(1.0L, 2.0L * PI * std::remainder(frequency * index, 1.0L));
    }

    double amplitude(std::size_t channel, unsigned tone) {
        return tone == 0 ? 1.0 + 0.1 * channel : 0.37;
    }
    double initialPhase(std::size_t channel, unsigned tone) {
        return tone == 0 ? 0.3 + 0.7 * channel : -1.1 + 0.4 * channel;
    }

    /** @brief Write actual IPPL arrays, poisoning every halo before every sample.
     * @details Temporal factors use FourierProbe::phase; spatial factors are
     * precomputed at each channel's position. Only the selected probe points and
     * plane need finite samples. Writing synthetic fields without a solver step
     * isolates transform accuracy from propagation and reconstruction errors.
     */
    class SignalWriter {
        using SpaceView     = Kokkos::View<double****>;
        using TimeView      = Kokkos::View<double***>;
        using VectorView    = FourierProbe::VectorField::view_type;
        using PotentialView = FourierProbe::PotentialField::view_type;

        /** @brief One device work item poisons a ghost or writes one selected owned cell.
         * @details Ghost poisoning occurs on every step, including when the plane gate
         * is closed. The probe therefore cannot silently use halo copies left by a
         * previous step. Component assignment follows channels(): three E components,
         * three output-B components, phi, three A components and one raw-B component.
         */
        struct SampleKernel {
            VectorView electric_m, magnetic_m, raw_m;
            PotentialView potential_m;
            TimeView temporal_m;
            SpaceView spatial_m;
            ippl::Vector<int, SpatialDim> first_m;
            int ghost_m;
            bool plane_m;
            double nan_m;

            KOKKOS_INLINE_FUNCTION bool owned(int i, int j, int k) const {
                return i >= ghost_m && i < int(electric_m.extent(0)) - ghost_m && j >= ghost_m
                       && j < int(electric_m.extent(1)) - ghost_m && k >= ghost_m
                       && k < int(electric_m.extent(2)) - ghost_m;
            }
            KOKKOS_INLINE_FUNCTION bool selected(int i, int j, int k) const {
                const int x = i - ghost_m + first_m[0], y = j - ghost_m + first_m[1],
                          z = k - ghost_m + first_m[2];
                return (plane_m && x == 3) || (x == 0 && y == 0 && z == 7)
                       || (x == 7 && y == 7 && (z == 8 || z == 15));
            }
            KOKKOS_INLINE_FUNCTION void writeFinite(std::size_t item, int i, int j, int k) const {
                double values[11]{};
                for (unsigned c = 0; c < 11; ++c)
                    for (unsigned tone = 0; tone < 2; ++tone)
                        values[c] += temporal_m(c, tone, 0) * spatial_m(item, c, tone, 0)
                                     - temporal_m(c, tone, 1) * spatial_m(item, c, tone, 1);
                for (unsigned c = 0; c < SpatialDim; ++c) {
                    electric_m(i, j, k)[c]      = values[c];
                    magnetic_m(i, j, k)[c]      = values[c + SpatialDim];
                    potential_m(i, j, k)[c + 1] = values[c + 7];
                }
                potential_m(i, j, k)[0] = values[6];
                raw_m(i, j, k)[0]       = values[10];
            }
            KOKKOS_INLINE_FUNCTION void operator()(const std::size_t item) const {
                const int i = int(item / (electric_m.extent(1) * electric_m.extent(2)));
                const int j = int(item / electric_m.extent(2) % electric_m.extent(1));
                const int k = int(item % electric_m.extent(2));
                if (!owned(i, j, k)) {
                    electric_m(i, j, k)  = nan_m;
                    magnetic_m(i, j, k)  = nan_m;
                    potential_m(i, j, k) = nan_m;
                    raw_m(i, j, k)       = nan_m;
                } else if (selected(i, j, k))
                    writeFinite(item, i, j, k);
            }
        };

    public:
        SignalWriter(chdr::vacuum::Fields& fields, FourierProbe::VectorField& raw,
                     LevelContract levels)
            : fields_m(fields)
            , raw_m(raw)
            , levels_m(levels) {
            const auto e   = fields.getE().getView();
            total_m        = e.extent(0) * e.extent(1) * e.extent(2);
            spatial_m      = SpaceView("T8 spatial factors", total_m, 11, 2, 2);
            temporal_m     = TimeView("T8 temporal factors", 11, 2, 2);
            temporalHost_m = Kokkos::create_mirror_view(temporal_m);
            initializeSpatial();
        }

        /// Separate host phase preparation from the one device write/poison operation.
        void write(std::int64_t step, bool plane) {
            prepareTemporal(step);
            const auto local = fields_m.getFL().getLocalNDIndex();
            const SampleKernel kernel{fields_m.getE().getView(),
                                      fields_m.getB().getView(),
                                      raw_m.getView(),
                                      fields_m.getJ().getView(),
                                      temporal_m,
                                      spatial_m,
                                      {local[0].first(), local[1].first(), local[2].first()},
                                      fields_m.getE().getNghost(),
                                      plane,
                                      std::numeric_limits<double>::quiet_NaN()};
            Kokkos::parallel_for("T8 poison ghosts and write samples", total_m, kernel);
        }

    private:
        /// Evaluate the component's centre, face or edge position through LevelContract.
        std::array<double, SpatialDim> position(std::size_t item, std::size_t channel) const {
            const auto e      = fields_m.getE().getView();
            const auto local  = fields_m.getFL().getLocalNDIndex();
            const int ghost   = fields_m.getE().getNghost();
            const auto spec   = channels()[channel];
            const auto offset = levels_m.positionOffset(spec.quantity_m, spec.component_m);
            const int i       = int(item / (e.extent(1) * e.extent(2)));
            const int j = int(item / e.extent(2) % e.extent(1)), k = int(item % e.extent(2));
            return {(i - ghost + local[0].first() + 0.5 + offset[0]) * 0.125,
                    (j - ghost + local[1].first() + 0.5 + offset[1]) * 0.125,
                    (k - ghost + local[2].first() + 0.5 + offset[2]) * 0.125};
        }
        void initializeSpatialCell(SpaceView::host_mirror_type host, std::size_t item) const {
            for (std::size_t c = 0; c < 11; ++c) {
                const auto x = position(item, c);
                for (unsigned tone = 0; tone < 2; ++tone) {
                    const double phase =
                        2.0 * double(PI) * (x[0] + 2.0 * x[1] + 1.5 * x[2]) + initialPhase(c, tone);
                    host(item, c, tone, 0) = std::cos(phase);
                    host(item, c, tone, 1) = std::sin(phase);
                }
            }
        }
        void initializeSpatial() {
            auto host = Kokkos::create_mirror_view(spatial_m);
            for (std::size_t item = 0; item < total_m; ++item)
                initializeSpatialCell(host, item);
            Kokkos::deep_copy(spatial_m, host);
        }

        /** @brief Generate tones at each channel's true half-step time.
         * @details The temporal factor uses integer step counts, as does the probe,
         * while the closed-form reference independently uses long-double geometric
         * sums. Spatial factors were prepared once and are never advanced in time.
         */
        void prepareTemporal(std::int64_t step) {
            const auto specs = channels();
            for (std::size_t c = 0; c < specs.size(); ++c)
                for (unsigned tone = 0; tone < 2; ++tone) {
                    const double phase = FourierProbe::phase(
                        TONES[tone],
                        2 * step + levels_m.halfStepLabel(specs[c].quantity_m).count());
                    temporalHost_m(c, tone, 0) = amplitude(c, tone) * std::cos(phase);
                    temporalHost_m(c, tone, 1) = amplitude(c, tone) * std::sin(phase);
                }
            Kokkos::deep_copy(temporal_m, temporalHost_m);
        }
        chdr::vacuum::Fields& fields_m;
        FourierProbe::VectorField& raw_m;
        LevelContract levels_m;
        std::size_t total_m = 0;
        SpaceView spatial_m;
        TimeView temporal_m;
        TimeView::host_mirror_type temporalHost_m;
    };

    /// Transverse vector potential with zero scalar potential; exact discrete temporal frequency.
    struct LivePotential {
        double omega_m;
        KOKKOS_INLINE_FUNCTION double operator()(double, double, double z, double time,
                                                 unsigned component) const {
            return component == 1 ? Kokkos::cos(Kokkos::numbers::pi_v<double> * z - omega_m * time)
                                  : 0;
        }
    };

    /** @brief Evaluate the travelling mode's independently differentiated channel value.
     * @details Substitution of a cosine into the differences gives electric factor
     * \f$2\sin(\omega\Delta t/2)/\Delta t\f$, centred curl factor
     * \f$\sin(kh)/h\f$, or staggered derivative factor \f$2\sin(kh/2)/h\f$.
     * These follow from [fallahi2020mithra20fullwavesimulation, Eqs. (3.17),
     * (3.50)-(3.55)] and [ryu2016potentialfdtd, App. A, Eqs. (19)-(24)].
     */
    long double liveSample(unsigned channel, int step,
                           EBReconstructionConvention ebReconstructionConvention,
                           SpatialDiscretization spatialDiscretization, long double omega,
                           long double dt, long double hz) {
        const bool raw =
            channel == 2
            || (channel == 1 && ebReconstructionConvention == EBReconstructionConvention::R1);
        const long double label = (channel < 3 && !raw) ? step - 0.5L : step;
        const long double z     = 7.5L * hz
                              + (spatialDiscretization == SpatialDiscretization::Staggered
                                         && channel > 0 && channel < 3
                                     ? hz / 2
                                     : 0);
        const long double theta = PI * z - omega * label * dt;
        if (channel == 3)
            return std::cos(theta);
        if (channel == 4)
            return 0;
        const long double derivative = spatialDiscretization == SpatialDiscretization::Staggered
                                           ? 2 * std::sin(PI * hz / 2) / hz
                                           : std::sin(PI * hz) / hz;
        const long double factor     = channel == 0 ? 2 * std::sin(omega * dt / 2) / dt
                                                    : derivative * (raw ? 1 : std::cos(omega * dt / 2));
        return -factor * std::sin(theta);
    }

    void checkLiveRecord(const ProbeRecord& record,
                         EBReconstructionConvention ebReconstructionConvention,
                         SpatialDiscretization spatialDiscretization, double omega, double dt,
                         double hz, double frequency, int first) {
        if (ippl::Comm->rank() != 0)
            return;
        for (unsigned channel = 0; channel < 5; ++channel) {
            Complex expected(0, 0);
            const bool half =
                channel == 0
                || (channel == 1 && ebReconstructionConvention != EBReconstructionConvention::R1);
            for (int step = first; step < first + 32; ++step)
                expected += static_cast<long double>(dt)
                            * liveSample(channel, step, ebReconstructionConvention,
                                         spatialDiscretization, omega, dt, hz)
                            * cycle(frequency, step - (half ? 0.5L : 0));
            EXPECT_LE(std::abs(Complex(record.value(0, channel, 0)) - expected), 1e-12 * dt * 32);
        }
    }

    /// Exercise reconstructed E/B, independent raw B, and both kinds of potential slot.
    std::vector<ProbeChannel> liveChannels() {
        return {{Quantity::Electric, 0},
                {Quantity::Magnetic, 1},
                {Quantity::RawMagnetic, 1},
                {Quantity::VectorPotential, 0},
                {Quantity::Phi, 0}};
    }

    /// Preserve the same two-unit axial box with a NonStandard-compatible aspect ratio.
    chdr::vacuum::Fields nonStandardLiveFields() {
        ippl::Vector<double, SpatialDim> origin(0.0), upper(1.0, 1.0, 2.0),
            spacing(0.125, 0.125, 0.0625);
        ippl::NDIndex<SpatialDim> domain(ippl::Index(8), ippl::Index(8), ippl::Index(32));
        return {spacing, origin, upper, {false, false, true}, domain, origin, true};
    }

    /** @brief Sample an evolved mode through the bound probe in either stencil family.
     * @details For an axial mode, both wave symbols reduce to the longitudinal
     * second difference [fallahi2020mithra20fullwavesimulation, Eqs. (3.17),
     * (3.24)]. The NonStandard probe uses the resonant frequency so the raw-B
     * check cannot pass merely because a whole-period off-resonant sum vanishes.
     */
    template <class PotentialSolverBase = StandardSolverPeriodic>
    void liveSolverProbe(EBReconstructionConvention ebReconstructionConvention,
                         SpatialDiscretization spatialDiscretization) {
        constexpr bool NONSTANDARD = std::is_same_v<PotentialSolverBase, NonStandardSolverPeriodic>;
        auto fields                = NONSTANDARD ? nonStandardLiveFields() : makeFields();
        initializeFieldStorage(fields);
        ChdrSolver<PotentialSolverBase> solver(fields.getJ(), fields.getE(), fields.getB(),
                                               ebReconstructionConvention, spatialDiscretization);
        const double dt = solver.getDt(), hz = fields.getMesh().getMeshSpacing()[2];
        const double omega     = 2 * std::asin(dt * std::sin(double(PI) * hz / 2) / hz) / dt;
        const double frequency = NONSTANDARD ? 1.0 / 32.0 : 0.125;
        loadHistory(solver, LivePotential{omega});
        SolverFourierProbe initial(solver, {{0, 0, 7}, {0, 0, 7}}, liveChannels(), {frequency}, 0,
                                   32, ProbeWindow::Rectangle, true);
        SolverFourierProbe delayed(solver, {{0, 0, 7}, {0, 0, 7}}, liveChannels(), {frequency}, 3,
                                   32, ProbeWindow::Rectangle, true);
        for (int step = 0; step < 35; ++step) {
            initial.sample();
            delayed.sample();
            if (step < 34)
                solver.solve();
        }
        checkLiveRecord(initial.gather(), ebReconstructionConvention, spatialDiscretization, omega,
                        dt, hz, frequency, 0);
        checkLiveRecord(delayed.gather(), ebReconstructionConvention, spatialDiscretization, omega,
                        dt, hz, frequency, 3);
    }

    TEST(FourierProbeT8, LiveSolverChannels) {
        for (auto ebReconstructionConvention :
             {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
              EBReconstructionConvention::R2Prime})
            liveSolverProbe(ebReconstructionConvention, SpatialDiscretization::Collocated);
        liveSolverProbe(EBReconstructionConvention::R2, SpatialDiscretization::Staggered);
        liveSolverProbe<NonStandardSolverPeriodic>(EBReconstructionConvention::R2,
                                                   SpatialDiscretization::Collocated);
    }

    /// A separately labelled raw B channel carries the label-shift factor for identical data.
    TEST(FourierProbeT8, IdenticalSamplesGiveTheLabelShiftFactor) {
        auto fields = makeFields();
        chdr::vacuum::initializeFieldStorage(fields);
        FourierProbe probe(fields.getE(), {EBReconstructionConvention::R2}, DT,
                           {{0, 0, 7}, {0, 0, 7}},
                           {{Quantity::Magnetic, 0}, {Quantity::RawMagnetic, 0}}, {0.125}, 9, 64);
        for (int step = 9; step < 73; ++step) {
            fields.getB() = std::cos(FourierProbe::phase(0.125, 2 * step - 1) + 0.3);
            probe.sample(step, fields.getE(), fields.getB(), fields.getB(), fields.getJ());
        }
        const auto record = probe.gather();
        if (ippl::Comm->rank() != 0)
            return;
        EXPECT_LE(std::abs(record.value(0, 1, 0)
                           - std::polar(1.0, double(PI) * 0.125) * record.value(0, 0, 0)),
                  1e-12 * DT * 64);
    }

    /// Independent entries of LevelContract's documented time table, without calling its label
    /// code.
    int expectedLabel(EBReconstructionConvention ebReconstructionConvention, std::size_t channel) {
        int expected = channel < SpatialDim ? -1 : 0;
        if (channel >= SpatialDim && channel < 2 * SpatialDim
            && ebReconstructionConvention != EBReconstructionConvention::R1)
            expected = -1;
        if (channel == 6)
            expected = ebReconstructionConvention == EBReconstructionConvention::R2        ? 1
                       : ebReconstructionConvention == EBReconstructionConvention::R2Prime ? -1
                                                                                           : 0;
        return expected;
    }

    /// Independent centre/face/edge coordinates for the fixed cell, using LevelContract's
    /// documented placement.
    std::array<double, SpatialDim> expectedCoordinate(SpatialDiscretization spatialDiscretization,
                                                      std::size_t channel) {
        std::array<double, SpatialDim> coordinate{0.3125, 0.4375, 0.9375};
        if (spatialDiscretization == SpatialDiscretization::Staggered) {
            if (channel < SpatialDim)
                coordinate[channel] += 0.0625;
            if (channel >= 7 && channel <= 9)
                coordinate[channel - 7] += 0.0625;
            if ((channel >= SpatialDim && channel < 2 * SpatialDim) || channel == 10)
                for (unsigned d = 0; d < SpatialDim; ++d)
                    if (d != (channel == 10 ? 0 : channel - SpatialDim))
                        coordinate[d] += 0.0625;
        }
        return coordinate;
    }

    void checkShellWiring(EBReconstructionConvention ebReconstructionConvention,
                          SpatialDiscretization spatialDiscretization) {
        auto fields = makeFields();
        chdr::vacuum::initializeFieldStorage(fields);
        ChdrSolver<StandardSolverPeriodic> solver(fields.getJ(), fields.getE(), fields.getB(),
                                                  ebReconstructionConvention,
                                                  spatialDiscretization);
        FourierProbe probe(solver.electricField(), solver.levelContract(), solver.getDt(),
                           {{2, 3, 7}, {2, 3, 7}}, channels(), {0.125}, 0, 8,
                           ProbeWindow::Rectangle,
                           ebReconstructionConvention == EBReconstructionConvention::R1);
        const auto& record = probe.metadata();
        EXPECT_EQ(record.levels_m.ebReconstructionConvention_m, ebReconstructionConvention);
        EXPECT_EQ(record.levels_m.spatialDiscretization_m, spatialDiscretization);
        EXPECT_EQ(record.phiIsZero_m, ebReconstructionConvention == EBReconstructionConvention::R1);
        for (std::size_t c = 0; c < channels().size(); ++c) {
            EXPECT_EQ(record.channels_m[c].halfStepLabel_m.count(),
                      expectedLabel(ebReconstructionConvention, c));
            EXPECT_EQ(record.coordinates_m[c], expectedCoordinate(spatialDiscretization, c));
        }
    }

    TEST(FourierProbeT8, ShellLabelsAndChannelCoordinates) {
        for (auto ebReconstructionConvention :
             {EBReconstructionConvention::R1, EBReconstructionConvention::R2,
              EBReconstructionConvention::R2Prime})
            for (auto spatialDiscretization :
                 {SpatialDiscretization::Collocated, SpatialDiscretization::Staggered})
                if (spatialDiscretization == SpatialDiscretization::Collocated
                    || ebReconstructionConvention == EBReconstructionConvention::R2)
                    checkShellWiring(ebReconstructionConvention, spatialDiscretization);
    }

    /// Exercise the same accumulation kernel with the Standard staggered channel positions.
}  // namespace
