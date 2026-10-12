/** @file FourierProbe.h
 * @brief Owned-cell, time-labelled running Fourier transforms.
 *
 * This file will provide frequency-domain measurements during a simulation.
 *
 * One can select a point, plane or region, the field components to observe,
 * and the frequencies of interest. As the solver advances, the fourier probe
 * can accumulate their fourier transforms without storing the full time history.
 * It uses each quantities time label (see the LevelContract) and handles the
 * ownership across MPI ranks.
 *
 * Will be used for: measurements for reflection, transmission and radiation spectra.
 */
#pragma once

#include "Ippl.h"

#include <Kokkos_MathematicalConstants.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "LevelContract.h"
#include "SpatialDimension.h"
#include "VacuumSetup.h"

namespace chdr::solver {

    /// Inclusive global-index box: a point, line, plane or slab uses the same storage.
    struct ProbeBox {
        std::array<int, SpatialDim> first_m, last_m;
    };

    /// Finite time gate with unit or periodic-Hann weights.
    enum class ProbeWindow {
        Rectangle,
        Hann
    };

    /// Named solver quantity and Cartesian component; scalar phi uses component zero.
    struct ProbeChannel {
        Quantity quantity_m;
        unsigned component_m = 0;
    };

    /// Resolved half-step label and spatial offset for one channel.
    struct ProbeChannelRecord {
        ProbeChannel channel_m;
        HalfStepOffset halfStepLabel_m;
        Kokkos::Array<double, SpatialDim> offset_m;
    };

    /** @brief Complete probe record, available on communicator rank zero.
     * @details Values are ordered cell, channel, frequency, real/imaginary.
     * Cell order is x, y, z with z contiguous. Coordinates are stored separately
     * per channel because the staggered reading locates components differently
     * [ryu2016potentialfdtd, App. A, Eqs. (19)-(24)].
     */
    struct ProbeRecord {
        static constexpr const char* CONVENTION =
            "dt sum w f exp(+i omega t); spatial exp(-i k.x); forward factor 1; inverse 1/(2pi) "
            "per variable";
        double dt_m = 0.0, windowSum_m = 0.0, windowSquareSum_m = 0.0;
        std::int64_t firstStep_m = 0, samples_m = 0;
        ProbeWindow window_m = ProbeWindow::Rectangle;
        LevelContract levels_m;
        bool phiIsZero_m = false;
        std::array<double, SpatialDim> origin_m{}, spacing_m{};
        std::vector<double> frequencies_m, angularFrequencies_m;
        std::vector<ProbeChannelRecord> channels_m;
        std::vector<std::array<int, SpatialDim>> indices_m;
        std::vector<std::array<double, SpatialDim>> coordinates_m;
        std::vector<int> ownership_m;
        std::vector<double> values_m;

        std::complex<double> value(std::size_t cell, std::size_t channel,
                                   std::size_t frequency) const {
            const auto index =
                2 * ((cell * channels_m.size() + channel) * frequencies_m.size() + frequency);
            return {values_m.at(index), values_m.at(index + 1)};
        }
    };

    /** @brief Accumulate a time transform independently at each owned probe cell.
     * @details The convention is
     * \f$\widehat f_c=\Delta t\sum_n w_n f_c[k_0+n]
     * e^{2\pi i\nu(k_0+n+m_c/2)}\f$
     * [oskooi2010meep, Sec. 5.1, Eq. (6); oskooi2013sources, Sec. 4.3,
     * Eq. (4.6), with its normalization converted to forward factor one].
     * LevelContract supplies the integer half-step offset \f$m_c\f$:
     * \f$f_c[k]\f$ represents the channel at \f$(k+m_c/2)\Delta t\f$.
     * Here \f$\nu\f$ is cycles per step and \f$\omega=2\pi\nu/\Delta t\f$.
     * The kernel has positive time phase. The factor \f$\Delta t\f$ is the
     * rectangle-rule integration weight; no division by the sample/window count
     * and no spatial area or volume factor is applied. The gate is
     * \f$k_0\leq k<k_0+N\f$, and weight() defines its optional window.
     * Every owned cell has one sequential running sum. Phases are evaluated on
     * the host and copied to the default execution space; fields and accumulators
     * stay in that space. Only gather() copies accumulated data back to the host.
     * MPI operates on pairs of doubles, never Kokkos complex types.
     * For P global probe cells, C channels and F frequencies, each rank stores
     * P indices and PC coordinates; gather additionally allocates 2PCF doubles
     * per rank, with another 2PCF doubles on rank zero. Temporal accumulators
     * contain 2P_local CF doubles. Large planes should budget this global buffer.
     * @pre All supplied fields share the constructor field's layout and halo width.
     * Call sample once for each step in the gate, in order, on every rank.
     */
    class FourierProbe {
    public:
        /// E/B storage with Cartesian components indexed from zero through SpatialDim-1.
        using VectorField = VField_t<double, SpatialDim>;
        /// Potential storage: scalar component zero, Cartesian vector components 1 through
        /// SpatialDim.
        using PotentialField                         = SourceField_t<double, SpatialDim>;
        FourierProbe(const FourierProbe&)            = delete;
        FourierProbe& operator=(const FourierProbe&) = delete;
        FourierProbe(FourierProbe&&)                 = default;
        FourierProbe& operator=(FourierProbe&&)      = default;

        /** @brief Reduce the product before forming an angle using its FMA remainder.
         * @details The time-transform phase is periodic in cycles
         * [oskooi2010meep, Eq. (6)]. For integer \f$n=2k+m\f$, compute the
         * rounded product \f$p=(\nu/2)n\f$ and its fused correction
         * \f$e=\mathrm{fma}(\nu/2,n,-p)\f$. Removing the nearest integer from
         * p before adding e gives a small cycle argument for the trigonometric call.
         * This product/correction construction is a deduction under the
         * normal-range round-to-nearest model [jeannerod2013ImprovedErrorBounds,
         * Sec. 1, Eq. (1.1)]. Each sample is phased directly from its integer step.
         * The binary64 gate bound keeps the integer half-step index exactly representable.
         */
        static double phase(double frequency, std::int64_t halfStep) {
            const double index = double(halfStep), halfFrequency = frequency * 0.5;
            const double product   = halfFrequency * index;
            const double remainder = std::fma(halfFrequency, index, -product);
            return TWO_PI * ((product - std::nearbyint(product)) + remainder);
        }

    private:
        static constexpr double TWO_PI         = 2.0 * Kokkos::numbers::pi_v<double>;
        static constexpr std::int64_t MAX_STEP = (std::int64_t{1} << 51) - 2;
        using CellView                         = Kokkos::View<int* [4]>;
        using ChannelView                      = Kokkos::View<ProbeChannel*>;
        using FactorView                       = Kokkos::View<double***>;
        using ValueView                        = Kokkos::View<double****>;
        using VectorFieldView                  = VectorField::view_type::const_type;
        using FourPotentialView                = PotentialField::view_type::const_type;
        struct OutputBuffers {
            std::vector<double> values_m;
            std::vector<int> counts_m;
        };

        /// Store enough temporal metadata to interpret every transform independently.
        void initializeMetadata(LevelContract levels, double dt, std::int64_t first,
                                std::int64_t count, ProbeWindow window, bool phiIsZero) {
            record_m.dt_m        = dt;
            record_m.levels_m    = levels;
            record_m.phiIsZero_m = phiIsZero;
            record_m.firstStep_m = first;
            record_m.samples_m   = count;
            record_m.window_m    = window;
        }

        /// Keep cycles/step as phase input; report angular frequencies in reciprocal units of dt.
        void initializeFrequencies(std::vector<double> frequencies) {
            if (frequencies.empty())
                throw std::invalid_argument("Probe needs frequencies");
            for (double frequency : frequencies)
                if (!(frequency > 0.0 && frequency < 0.5))
                    throw std::invalid_argument("Probe frequency must be in (0, 0.5) cycles/step");
            record_m.frequencies_m = std::move(frequencies);
            for (double frequency : record_m.frequencies_m) {
                const double angular = TWO_PI * frequency / record_m.dt_m;
                if (!std::isfinite(angular))
                    throw std::invalid_argument("Probe angular frequency is not finite");
                record_m.angularFrequencies_m.push_back(angular);
            }
        }

        /// Resolve every label and position once; R1 E remains an explicit zero-phi declaration.
        void initializeChannels(const std::vector<ProbeChannel>& channels) {
            if (channels.empty())
                throw std::invalid_argument("Probe needs channels");
            for (const auto channel : channels) {
                if ((channel.quantity_m == Quantity::Phi && channel.component_m != 0)
                    || channel.component_m > 2 || channel.quantity_m == Quantity::Charge
                    || channel.quantity_m == Quantity::Current)
                    throw std::invalid_argument("Unsupported probe channel or component");
                record_m.channels_m.push_back(
                    {channel,
                     record_m.levels_m.halfStepLabel(channel.quantity_m, record_m.phiIsZero_m),
                     record_m.levels_m.positionOffset(channel.quantity_m, channel.component_m)});
            }
        }
        ProbeBox box_m;
        MPI_Comm comm_m;
        ProbeRecord record_m;
        std::int64_t samplesSeen_m = 0;
        CellView cells_m;
        ChannelView channels_m;
        ValueView sums_m;
        FactorView factors_m;
        FactorView::host_mirror_type factorsHost_m;
    };

}  // namespace chdr::solver
