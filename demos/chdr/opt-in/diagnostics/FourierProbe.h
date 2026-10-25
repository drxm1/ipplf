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

#include "core/base/FieldContainer.h"
#include "core/base/LevelContract.h"
#include "core/base/SpatialDimension.h"

static_assert(chdr::solver::SpatialDim == 3,
              "FourierProbe.h requires three spatial dimensions for its Cartesian implementation");

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
        FourierProbe(const FourierProbe&)            = delete;
        FourierProbe& operator=(const FourierProbe&) = delete;
        FourierProbe(FourierProbe&&)                 = default;
        FourierProbe& operator=(FourierProbe&&)      = default;

        /** @brief Bind immutable geometry, channels, frequencies and a half-open gate.
         * @details Labels and offsets come from one LevelContract. The caller supplies
         * quantities/components rather than phases, so changing a reconstruction cannot
         * silently retain the previous reconstruction's channel time or position.
         * @param geometry Allocated field supplying mesh, layout and owned-cell indices.
         * @param levels Time labels and component positions of the supplied data.
         * @param dt Finite positive sample interval, fixed throughout the gate.
         * @param box Inclusive global-index sampling box.
         * @param channels Named quantities and components to transform.
         * @param frequencies Frequencies in cycles per step, strictly between zero and one half.
         * @param firstStep First completed-solve index included in the gate.
         * @param samples Number of consecutive samples in the gate.
         * @param window Rectangle or periodic Hann weights on gate-local indices.
         * @param phiIsZero Explicit declaration permitting R1 electric samples.
         */
        FourierProbe(const EMFieldsType& geometry, LevelContract levels, double dt, ProbeBox box,
                     std::vector<ProbeChannel> channels, std::vector<double> frequencies,
                     std::int64_t firstStep, std::int64_t samples,
                     ProbeWindow window = ProbeWindow::Rectangle, bool phiIsZero = false)
            : box_m(box)
            , comm_m(geometry.getLayout().comm) {
            validateTimeGate(dt, firstStep, samples, window);
            initializeMetadata(levels, dt, firstStep, samples, window, phiIsZero);
            initializeFrequencies(std::move(frequencies));
            initializeChannels(channels);
            initializeGeometry(geometry);
            allocateOwnedCells(ownedCells(geometry));
            allocateChannels(channels);
            allocateAccumulators();
            initializeWindowSums();
        }

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

        /// Periodic Hann weights use gate-local indices [heinzel2002SpectrumSpectralDensity, Eq.
        /// (17)].
        double weight(std::int64_t sample) const {
            return record_m.window_m == ProbeWindow::Rectangle
                       ? 1.0
                       : 0.5
                             * (1.0
                                - std::cos(TWO_PI * double(sample) / double(record_m.samples_m)));
        }
        const ProbeRecord& metadata() const { return record_m; }

        /** @brief Sample E/B, separately supplied raw B and the current potential slot.
         * @details Calls outside the gate are harmless. Inside it, duplicate or skipped
         * steps are rejected so the recorded window sums still describe the data.
         * Phases are shared across cells; the device kernel reads only owned indices.
         * @param step Completed-solve count after solve(), or zero after evaluateInitialFields().
         * @param electric Reconstructed electric output at this completed-solve count.
         * @param magnetic Reconstructed magnetic output at this completed-solve count.
         * @param rawMagnetic Separately reconstructed unaveraged magnetic field.
         * @param currentA Current A_n potential slot at this completed-solve count:
         * scalar potential in component zero and Cartesian vector potential in components
         * 1 through SpatialDim. Channel times/positions come from the bound LevelContract.
         * @pre At count zero, reconstruct initial E/B after loading histories.
         * Prefer SolverFourierProbe for evolving solver fields: it supplies these arguments.
         */
        void sample(std::int64_t step, const EMFieldsType& electric, const EMFieldsType& magnetic,
                    const EMFieldsType& rawMagnetic, const AJFourFieldsType& currentA) {
            if (step < record_m.firstStep_m || step >= record_m.firstStep_m + record_m.samples_m)
                return;
            if (step != record_m.firstStep_m + samplesSeen_m)
                throw std::logic_error("Probe samples must cover its gate exactly once in order");
            const double scale = record_m.dt_m * weight(samplesSeen_m++);
            preparePhases(step, scale);
            accumulate(electric, magnetic, rawMagnetic, currentA);
        }

        /** @brief Collect completed transforms and counts on communicator rank zero.
         * @pre Collective after the complete gate, in the same order on every rank.
         * @details Each cell's complete temporal sum is computed by its one owner.
         * Other ranks supply zero in that output position, so MPI does not combine
         * partial temporal sums in a different order. Ownership counts check this
         * condition explicitly. Real and imaginary parts are transferred as plain
         * doubles, avoiding any dependency on a complex MPI datatype mapping.
         */
        ProbeRecord gather() const {
            if (samplesSeen_m != record_m.samples_m)
                throw std::logic_error("Cannot gather an incomplete Fourier probe gate");
            auto result        = record_m;
            const auto buffers = packOutput();
            reduceOutput(buffers, result);
            return result;
        }

    private:
        static constexpr double TWO_PI         = 2.0 * Kokkos::numbers::pi_v<double>;
        static constexpr std::int64_t MAX_STEP = (std::int64_t{1} << 51) - 2;
        using CellView                         = Kokkos::View<int* [4]>;
        using ChannelView                      = Kokkos::View<ProbeChannel*>;
        using FactorView                       = Kokkos::View<double***>;
        using ValueView                        = Kokkos::View<double****>;
        using VectorFieldView                  = EMFieldsType::view_type::const_type;
        using FourPotentialView                = AJFourFieldsType::view_type::const_type;
        struct OutputBuffers {
            std::vector<double> values_m;
            std::vector<int> counts_m;
        };

        /// Validate before allocating device buffers or making any collective call.
        static void validateTimeGate(double dt, std::int64_t first, std::int64_t count,
                                     ProbeWindow window) {
            if (!(dt > 0.0) || !std::isfinite(dt) || first < 0 || count < 1
                || first > MAX_STEP - count)
                throw std::invalid_argument("Invalid Fourier probe time gate");
            if (window != ProbeWindow::Rectangle && window != ProbeWindow::Hann)
                throw std::invalid_argument("Unknown Fourier probe window");
        }

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

        /// Reject out-of-domain boxes and bound the two-double output by MPI's integer count.
        void initializeGeometry(const EMFieldsType& geometry) {
            const auto& domain = geometry.getLayout().getDomain();
            std::size_t count  = 1;
            const auto width   = record_m.channels_m.size() * record_m.frequencies_m.size() * 2;
            for (unsigned d = 0; d < SpatialDim; ++d) {
                if (box_m.first_m[d] < domain[d].first() || box_m.last_m[d] > domain[d].last()
                    || box_m.last_m[d] < box_m.first_m[d])
                    throw std::invalid_argument("Probe box lies outside the global domain");
                // The closed outer layer supports reconstruction, but is not an
                // interior field sample for a radiation measurement. Potential
                // channels remain valid there for boundary-operator checks.
                if (!geometry.getLayout().isPeriodic(d)
                    && (box_m.first_m[d] == domain[d].first()
                        || box_m.last_m[d] == domain[d].last())) {
                    for (const auto& channel : record_m.channels_m)
                        if (channel.channel_m.quantity_m == Quantity::Electric
                            || channel.channel_m.quantity_m == Quantity::Magnetic
                            || channel.channel_m.quantity_m == Quantity::RawMagnetic)
                            throw std::invalid_argument(
                                "Field probes must exclude the outermost open layer");
                }
                const auto extent = std::size_t(box_m.last_m[d] - box_m.first_m[d] + 1);
                if (extent > std::size_t(std::numeric_limits<int>::max()) / width / count)
                    throw std::invalid_argument("Probe output exceeds one MPI count");
                count *= extent;
                record_m.origin_m[d]  = geometry.get_mesh().getOrigin()[d];
                record_m.spacing_m[d] = geometry.get_mesh().getMeshSpacing()[d];
                if (!std::isfinite(record_m.origin_m[d]) || !std::isfinite(record_m.spacing_m[d])
                    || record_m.spacing_m[d] <= 0.0)
                    throw std::invalid_argument(
                        "Probe geometry must be finite with positive spacing");
            }
        }

        /// Global indices are ordered x/y/z, with z contiguous in the output.
        std::vector<std::array<int, 4>> ownedCells(const EMFieldsType& geometry) {
            const auto& local = geometry.getLayout().getLocalNDIndex();
            std::vector<std::array<int, 4>> owned;
            for (int i = box_m.first_m[0]; i <= box_m.last_m[0]; ++i)
                for (int j = box_m.first_m[1]; j <= box_m.last_m[1]; ++j)
                    for (int k = box_m.first_m[2]; k <= box_m.last_m[2]; ++k)
                        appendCell({i, j, k}, local, geometry.getNghost(), owned);
            return owned;
        }

        /// Ownership is tested in global indices; the kernel receives local indices including
        /// halos.
        void appendCell(std::array<int, SpatialDim> global, const ippl::NDIndex<SpatialDim>& local,
                        int ghost, std::vector<std::array<int, 4>>& owned) {
            const int index = int(record_m.indices_m.size());
            record_m.indices_m.push_back(global);
            bool contains = true;
            for (unsigned d = 0; d < SpatialDim; ++d)
                contains =
                    contains && global[d] >= local[d].first() && global[d] <= local[d].last();
            if (contains)
                owned.push_back({index, global[0] - local[0].first() + ghost,
                                 global[1] - local[1].first() + ghost,
                                 global[2] - local[2].first() + ghost});
            appendCoordinates(global);
        }

        /// Staggered channels report their own face/edge positions, not a common cell centre.
        void appendCoordinates(std::array<int, SpatialDim> global) {
            for (const auto& channel : record_m.channels_m) {
                std::array<double, SpatialDim> coordinate;
                for (unsigned d = 0; d < SpatialDim; ++d)
                    coordinate[d] =
                        record_m.origin_m[d]
                        + (global[d] + 0.5 + channel.offset_m[d]) * record_m.spacing_m[d];
                record_m.coordinates_m.push_back(coordinate);
            }
        }

        void allocateOwnedCells(const std::vector<std::array<int, 4>>& owned) {
            cells_m   = CellView("probe owned cells", owned.size());
            auto host = Kokkos::create_mirror_view(cells_m);
            for (std::size_t i = 0; i < owned.size(); ++i)
                for (unsigned d = 0; d < 4; ++d)
                    host(i, d) = owned[i][d];
            Kokkos::deep_copy(cells_m, host);
        }
        void allocateChannels(const std::vector<ProbeChannel>& channels) {
            channels_m = ChannelView("probe channels", channels.size());
            auto host  = Kokkos::create_mirror_view(channels_m);
            for (std::size_t c = 0; c < channels.size(); ++c)
                host(c) = channels[c];
            Kokkos::deep_copy(channels_m, host);
        }
        void allocateAccumulators() {
            const auto channels    = record_m.channels_m.size(),
                       frequencies = record_m.frequencies_m.size();
            sums_m = ValueView("probe accumulators", cells_m.extent(0), channels, frequencies, 2);
            Kokkos::deep_copy(sums_m, 0.0);
            factors_m     = FactorView("probe phases", channels, frequencies, 2);
            factorsHost_m = Kokkos::create_mirror_view(factors_m);
        }

        /** @brief Periodic Hann sums include the cases with one or two samples.
         * @details Deduction from the window and sum definitions in
         * [heinzel2002SpectrumSpectralDensity, Eqs. (17),(19),(20)]: summing
         * the roots of unity gives zero for the first cosine harmonic when
         * \f$N>1\f$ and for the second when \f$N>2\f$. Expanding the squared
         * window therefore gives \f$S_1=N/2\f$ and \f$S_2=3N/8\f$ for
         * \f$N>2\f$. Direct weights are \f$(0)\f$ for one sample and
         * \f$(0,1)\f$ for two, giving the explicit exceptional branches.
         */
        void initializeWindowSums() {
            const auto n               = record_m.samples_m;
            const bool rectangle       = record_m.window_m == ProbeWindow::Rectangle;
            record_m.windowSum_m       = rectangle ? double(n) : (n == 1 ? 0.0 : 0.5 * n);
            record_m.windowSquareSum_m = rectangle ? double(n)
                                                   : (n == 1   ? 0.0
                                                      : n == 2 ? 1.0
                                                               : 0.375 * n);
        }

        /// Reuse the first channel with this label; phases do not depend on cell or component.
        std::size_t phaseSource(std::size_t channel) const {
            std::size_t prior = 0;
            while (prior < channel
                   && record_m.channels_m[prior].halfStepLabel_m.count()
                          != record_m.channels_m[channel].halfStepLabel_m.count())
                ++prior;
            return prior;
        }
        void preparePhases(std::int64_t step, double scale) {
            // Host mirrors may alias device-accessible memory on unified-memory backends.
            Kokkos::fence("Complete probe sampling before updating shared phase storage");
            for (std::size_t c = 0; c < record_m.channels_m.size(); ++c)
                prepareChannelPhase(c, phaseSource(c), step, scale);
            Kokkos::deep_copy(factors_m, factorsHost_m);
        }
        void prepareChannelPhase(std::size_t channel, std::size_t prior, std::int64_t step,
                                 double scale) {
            for (std::size_t f = 0; f < record_m.frequencies_m.size(); ++f) {
                const double angle =
                    prior == channel
                        ? phase(record_m.frequencies_m[f],
                                2 * step + record_m.channels_m[channel].halfStepLabel_m.count())
                        : 0.0;
                factorsHost_m(channel, f, 0) =
                    prior == channel ? scale * std::cos(angle) : factorsHost_m(prior, f, 0);
                factorsHost_m(channel, f, 1) =
                    prior == channel ? scale * std::sin(angle) : factorsHost_m(prior, f, 1);
            }
        }

        /** @brief Read one channel directly from its matching field view.
         * @param channel Quantity and Cartesian component; scalar potential uses component zero.
         * @param i First local spatial index, including the halo offset.
         * @param j Second local spatial index, including the halo offset.
         * @param k Third local spatial index, including the halo offset.
         * @param electricView Reconstructed E components at the bound electric time.
         * @param magneticView Reconstructed B components, including convention-specific averaging.
         * @param rawMagneticView Unaveraged curl of the current vector potential.
         * @param currentAView Current four-component potential slot, ordered scalar then vector.
         * @details No interpolation is performed; LevelContract supplies each channel's
         * physical position and time, and sample checks the completed-solve index.
         */
        KOKKOS_INLINE_FUNCTION static double readSample(ProbeChannel channel, int i, int j, int k,
                                                        VectorFieldView electricView,
                                                        VectorFieldView magneticView,
                                                        VectorFieldView rawMagneticView,
                                                        FourPotentialView currentAView) {
            const unsigned component = channel.component_m;
            switch (channel.quantity_m) {
                case Quantity::Electric:
                    return electricView(i, j, k)[component];
                case Quantity::Magnetic:
                    return magneticView(i, j, k)[component];
                case Quantity::RawMagnetic:
                    return rawMagneticView(i, j, k)[component];
                case Quantity::Phi:
                    return currentAView(i, j, k)[0];
                case Quantity::VectorPotential:
                    return currentAView(i, j, k)[component + 1];
                default:
                    return 0.0;
            }
        }

    public:
        /** @brief Internal launcher: one work item owns a cell/channel's temporal sums.
         * @param electric Reconstructed E; Cartesian components start at zero.
         * @param magnetic Reconstructed B at the bound convention's output time.
         * @param rawMagnetic Separately reconstructed unaveraged B at the current A time.
         * @param currentA Current A_n slot: scalar potential then Cartesian vector potential.
         * @details Frequency and real/imaginary accumulation order stay fixed.
         * CUDA extended host/device lambdas require a public enclosing member;
         * this access follows the same restriction documented by IPPL's FFT wrappers.
         * Call sample() for the public sampling protocol, including gate/phase preparation.
         */
        void accumulate(const EMFieldsType& electric, const EMFieldsType& magnetic,
                        const EMFieldsType& rawMagnetic, const AJFourFieldsType& currentA) {
            const VectorFieldView e = electric.getView(), b = magnetic.getView(),
                                  raw = rawMagnetic.getView();
            const FourPotentialView a = currentA.getView();
            const auto cells          = cells_m;
            const auto channels       = channels_m;
            const auto sums           = sums_m;
            const auto factors        = factors_m;
            const auto channelCount = channels.extent(0), frequencyCount = factors.extent(1);
            Kokkos::parallel_for(
                "Fourier probe sample", cells.extent(0) * channelCount,
                KOKKOS_LAMBDA(const std::size_t item) {
                    const auto cell = item / channelCount, channel = item % channelCount;
                    const double value = readSample(channels(channel), cells(cell, 1),
                                                    cells(cell, 2), cells(cell, 3), e, b, raw, a);
                    for (std::size_t f = 0; f < frequencyCount; ++f)
                        for (unsigned part = 0; part < 2; ++part)
                            sums(cell, channel, f, part) += value * factors(channel, f, part);
                });
        }

    private:
        /// Scatter only owned accumulators into globally ordered zero-padded host buffers.
        OutputBuffers packOutput() const {
            const auto cells = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, cells_m);
            const auto sums  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, sums_m);
            const auto frequencies = record_m.frequencies_m.size();
            const auto width       = record_m.channels_m.size() * frequencies * 2;
            OutputBuffers buffers{std::vector<double>(record_m.indices_m.size() * width, 0.0),
                                  std::vector<int>(record_m.indices_m.size(), 0)};
            for (std::size_t cell = 0; cell < cells.extent(0); ++cell) {
                const auto global        = std::size_t(cells(cell, 0));
                buffers.counts_m[global] = 1;
                for (std::size_t c = 0; c < record_m.channels_m.size(); ++c)
                    for (std::size_t f = 0; f < frequencies; ++f)
                        for (unsigned part = 0; part < 2; ++part)
                            buffers.values_m[global * width + 2 * (c * frequencies + f) + part] =
                                sums(cell, c, f, part);
            }
            return buffers;
        }

        /// Counts accompany the same global cell order, detecting gaps or duplicated ownership.
        void reduceOutput(const OutputBuffers& buffers, ProbeRecord& result) const {
            int rank;
            MPI_Comm_rank(comm_m, &rank);
            if (rank == 0) {
                result.values_m.resize(buffers.values_m.size());
                result.ownership_m.resize(buffers.counts_m.size());
            }
            MPI_Reduce(buffers.values_m.data(), result.values_m.data(),
                       int(buffers.values_m.size()), MPI_DOUBLE, MPI_SUM, 0, comm_m);
            MPI_Reduce(buffers.counts_m.data(), result.ownership_m.data(),
                       int(buffers.counts_m.size()), MPI_INT, MPI_SUM, 0, comm_m);
            int invalid =
                rank == 0
                && std::any_of(result.ownership_m.begin(), result.ownership_m.end(), [](int count) {
                       return count != 1;
                   });
            MPI_Bcast(&invalid, 1, MPI_INT, 0, comm_m);
            if (invalid)
                throw std::runtime_error("Fourier probe cells do not have exactly one owner");
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

    /** @brief Bind a Fourier probe to one solver's fields, labels, step and history generation.
     * @tparam SolverType Solver exposing ChdrSolver's field, history and level interface.
     * @pre Construct after setting the time step and loading the potential histories.
     * @details sample() reads completedSolves() directly and reconstructs initial E/B
     * when the gate includes count zero. A requested raw-B channel gets separate storage
     * and is produced by computeRawMagneticField(); the common-time output is preserved.
     * Reloading or reinitializing histories invalidates the binding. The solver and its
     * caller-owned mesh/layout must outlive this object. Sampling remains opt-in.
     */
    template <class SolverType>
    class SolverFourierProbe {
    public:
        /** @brief Capture one solver's immutable sampling configuration.
         * @param solver Initialized shell with its final time step and loaded histories.
         * @param box Inclusive global-index sampling box.
         * @param channels Quantities and Cartesian components; scalar phi uses component zero.
         * @param frequencies Frequencies in cycles per step.
         * @param firstStep First completed-solve count in the half-open gate.
         * @param samples Number of consecutive samples.
         * @param window Rectangle or periodic-Hann weighting.
         * @param phiIsZero Explicit zero-scalar-potential declaration for R1 electric output.
         */
        SolverFourierProbe(SolverType& solver, ProbeBox box, std::vector<ProbeChannel> channels,
                           std::vector<double> frequencies, std::int64_t firstStep,
                           std::int64_t samples, ProbeWindow window = ProbeWindow::Rectangle,
                           bool phiIsZero = false)
            : solver_m(solver)
            , revision_m(solver.historyRevision())
            , probe_m(solver.electricField(), solver.levelContract(), solver.getDt(), box, channels,
                      std::move(frequencies), firstStep, samples, window, phiIsZero) {
            if (firstStep < solver.completedSolves())
                throw std::invalid_argument("Solver has already passed the probe gate start");
            if (std::any_of(channels.begin(), channels.end(), [](const auto& channel) {
                    return channel.quantity_m == Quantity::RawMagnetic;
                })) {
                raw_m.emplace();
                raw_m->initialize(solver.electricField().get_mesh(),
                                  solver.electricField().getLayout());
            }
        }

        /// Sample the current completed solve; reject changed histories before reading fields.
        void sample() {
            if (revision_m != solver_m.historyRevision()
                || probe_m.metadata().dt_m != solver_m.getDt())
                throw std::logic_error(
                    "Construct a new Fourier probe after changing solver histories");
            const auto step    = solver_m.completedSolves();
            const auto& record = probe_m.metadata();
            if (step < record.firstStep_m || step >= record.firstStep_m + record.samples_m)
                return;
            if (step == 0)
                solver_m.evaluateInitialFields();
            if (raw_m)
                solver_m.computeRawMagneticField(*raw_m);
            probe_m.sample(step, solver_m.electricField(), solver_m.magneticField(),
                           raw_m ? *raw_m : solver_m.magneticField(), solver_m.A_n);
        }
        const ProbeRecord& metadata() const { return probe_m.metadata(); }
        ProbeRecord gather() const { return probe_m.gather(); }

    private:
        SolverType& solver_m;
        unsigned long long revision_m;
        FourierProbe probe_m;
        std::optional<EMFieldsType> raw_m;
    };
}  // namespace chdr::solver
