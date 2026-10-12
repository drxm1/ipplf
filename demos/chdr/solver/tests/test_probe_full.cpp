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
    chdr::vacuum::Fields makeFields() {
        ippl::Vector<double, SpatialDim> origin(0.0), upper(1.0, 1.0, 2.0), spacing(0.125);
        ippl::NDIndex<SpatialDim> domain(ippl::Index(8), ippl::Index(8), ippl::Index(16));
        return {spacing, origin, upper, {false, false, true}, domain, origin, true};
    }

    Complex cycle(long double frequency, long double index) {
        return std::polar(1.0L, 2.0L * PI * std::remainder(frequency * index, 1.0L));
    }

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

}  // namespace
