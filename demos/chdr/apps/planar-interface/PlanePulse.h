/** @file PlanePulse.h
 * @brief Exact lattice plane packets and incident-subtracted Fourier measurements.
 */
#pragma once
#include <complex>
#include <functional>

#include "apps/planar-interface/InterfaceRun.h"
#include "apps/planar-interface/PlaneSpectrum.h"
#include "core/material/PlanarSlab.h"
#include "opt-in/diagnostics/FieldEnergy.h"
#include "opt-in/diagnostics/FourierProbe.h"

namespace chdr::apps::planarInterface {
    using PlaneFourierProbe = diagnostics::SolverFourierProbe<InterfaceSolver>;

    /// Positive normal mode magnitudes, with the propagation sign and exact lattice frequency.
    inline std::vector<std::array<double, 3>> planeModes(const PlaneRow& p) {
        std::vector<std::array<double, 3>> modes;
        for (int m = 1; m < p.cellsX_m / 2; ++m) {
            const double k = 2 * Kokkos::numbers::pi * m / (p.cellsX_m * p.h_m[0]);
            if (std::abs(k - p.kx0_m) >= 9 * p.sigmaK_m)
                continue;
            modes.push_back({p.incidentFromVacuum_m ? -k : k, planeFrequency(p, k),
                             std::exp(-0.5 * std::pow((k - p.kx0_m) / p.sigmaK_m, 2))});
        }
        if (modes.empty())
            throw std::invalid_argument("Plane packet has no periodic modes");
        return modes;
    }

    /** @brief Device sampler of an exact periodic plane packet in the x-z plane.
     * @details Each mode obeys Standard leapfrog dispersion
     * [fallahi2020mithra20fullwavesimulation, section 3.1.3, Eq. (3.17)].
     * TE samples \f$A_y=\psi\f$. TM samples the permittivity-weighted discrete curl
     * of \f$\psi\hat{\mathbf y}\f$: \f$A_x=-D_z\psi\f$, \f$A_z=D_x\psi\f$,
     * multiplied by incident epsilon and divided by harmonic face epsilon.
     * Its weighted backward divergence vanishes by commuting differences
     * [ryu2016potentialfdtd, Appendix A; chew2014generalizedgauge, Eq. (18)].
     * The sign gives a positive fixed-Cartesian \f$H_y\f$ incident amplitude.
     * Mode data are mirrored once from host to the default execution space.
     */
    struct PlanePulse {
        Kokkos::View<double* [3]> modes_m;
        double kz_m, centre_m, dt_m, hx_m, hz_m, epsilon_m, incidentEpsilon_m;
        int cellsX_m, dielectricCells_m;
        bool tm_m;

        explicit PlanePulse(const PlaneRow& p)
            : kz_m(p.kz_m)
            , centre_m(p.pulseCentre_m)
            , dt_m(p.dt_m)
            , hx_m(p.h_m[0])
            , hz_m(p.h_m[2])
            , epsilon_m(p.epsilon_m)
            , incidentEpsilon_m(planeMedia(p)[0])
            , cellsX_m(p.cellsX_m)
            , dielectricCells_m(p.dielectricCells_m)
            , tm_m(p.polarization_m == Polarization::TM) {
            const auto modes = planeModes(p);
            modes_m          = decltype(modes_m)("Plane pulse modes", modes.size());
            auto host        = Kokkos::create_mirror_view(modes_m);
            for (size_t m = 0; m < modes.size(); ++m)
                for (int c = 0; c < 3; ++c)
                    host(m, c) = modes[m][c];
            Kokkos::deep_copy(modes_m, host);
        }

        /// Scalar mode sum at a physical point; y is the invariant out-of-plane coordinate.
        KOKKOS_INLINE_FUNCTION double psi(double x, double z, double t) const {
            double result = 0;
            for (size_t m = 0; m < modes_m.extent(0); ++m)
                result +=
                    modes_m(m, 2)
                    * Kokkos::cos(modes_m(m, 0) * (x - centre_m) + kz_m * z - modes_m(m, 1) * t);
            return result;
        }

        /** @brief Look up a material cell using the resolved inclusive-centre map.
         * @details The curl weights query cell centres. Using their integer cell avoids
         * changing a tie through the rounded component-coordinate plus/minus half-spacing.
         * The resolver counted cells with exactly MaterialOwner's centre expression.
         */
        KOKKOS_INLINE_FUNCTION double epsilon(double x) const {
            const int cell    = static_cast<int>(Kokkos::floor(x / hx_m));
            const int wrapped = (cell % cellsX_m + cellsX_m) % cellsX_m;
            return wrapped < dielectricCells_m ? epsilon_m : 1;
        }

        KOKKOS_INLINE_FUNCTION double operator()(double x, double, double z, double t,
                                                 unsigned c) const {
            if (!tm_m)
                return c == 2 ? psi(x, z, t) : 0;
            if (c != 1 && c != 3)
                return 0;
            const double a    = epsilon(x - (c == 1 ? hx_m / 2 : 0));
            const double b    = epsilon(x + (c == 1 ? hx_m / 2 : 0));
            const double curl = c == 1 ? (psi(x, z - hz_m / 2, t) - psi(x, z + hz_m / 2, t)) / hz_m
                                       : (psi(x + hx_m / 2, z, t) - psi(x - hx_m / 2, z, t)) / hx_m;
            return incidentEpsilon_m * (1 / a + 1 / b) * curl / 2;
        }

        /// The measured channel, including the magnetic two-level average [ryu2016potentialfdtd,
        /// Appendix A].
        template <class HostModes>
        double incidentField(const HostModes& modes, double x, double z, double time) const {
            double field = 0;
            for (size_t m = 0; m < modes.extent(0); ++m) {
                const double omega = modes(m, 1), W = 2 * std::sin(omega * dt_m / 2) / dt_m;
                const double phase = modes(m, 0) * (x - centre_m) + kz_m * z - omega * time;
                field += modes(m, 2)
                         * (tm_m ? incidentEpsilon_m * W * W * std::cos(omega * dt_m / 2)
                                       * std::cos(phase)
                                 : -W * std::sin(phase));
            }
            return field;
        }

        /** @brief Sample the analytic incident channel once for all requested Fourier bins.
         * Modes are copied to host once. Cached real values use the measured half-step
         * labels [ryu2016potentialfdtd, Appendix A]; no time trace is published.
         */
        std::vector<double> incidentSamples(double x, double z, int gate) const {
            const auto modes = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), modes_m);
            std::vector<double> samples;
            samples.reserve(gate);
            for (int k = 0; k < gate; ++k)
                samples.push_back(incidentField(modes, x, z, (k - 0.5) * dt_m));
            return samples;
        }
    };

    /// The Fourier probe stores frequencies in cycles per completed simulation step.
    inline std::vector<double> planeProbeFrequencies(const PlaneRow& p, const PlaneResult& result) {
        std::vector<double> frequencies;
        for (double omega : result.omega_m)
            frequencies.push_back(omega * p.dt_m / (2 * Kokkos::numbers::pi));
        return frequencies;
    }

    /// Both polarizations observe the fixed y component at their own physical staggered location.
    inline PlaneFourierProbe planeProbe(InterfaceSolver& solver, const PlaneRow& p, int index,
                                        int gate, const std::vector<double>& frequencies) {
        const auto channel =
            p.polarization_m == Polarization::TM ? Quantity::Magnetic : Quantity::Electric;
        return PlaneFourierProbe(solver, {{index, 0, 0}, {index, 0, 0}}, {{channel, 1}},
                                 frequencies, 0, gate);
    }

    /// Preserve invariant normalization while keeping energy observation separate from spectra.
    inline void observePlaneEnergy(InterfaceSolver& state, const material::MaterialOwner& owner,
                                   long long k, PlaneResult& result, double& baseline) {
        const auto energy     = diagnostics::potentialEnergy(state, owner);
        result.scalarEnergy_m = energy.scalar_m.total();
        if (k == 0)
            baseline = energy.vector_m.total();
        if (!(baseline > 0) || !std::isfinite(energy.vector_m.total()))
            throw std::runtime_error("Plane pulse needs finite energy and a positive baseline");
        result.drift_m = std::max(result.drift_m, std::abs(energy.vector_m.total() / baseline - 1));
        result.absoluteRatio_m =
            std::max(result.absoluteRatio_m, energy.vector_m.absolute_m / baseline);
    }

    /// Observe every completed field state, independently of the energy diagnostic cadence.
    inline void advancePlaneMeasurement(InterfaceRun& run, const PlaneRow& p, PlaneResult& result,
                                        PlaneFourierProbe& incident, PlaneFourierProbe& other,
                                        int cadence, const InterfaceRun::Observer& observer) {
        double baseline = 0;
        run.advance(p.steps_m, [&](long long k, auto& state, const auto& owner) {
            incident.sample();
            other.sample();
            if (k % cadence == 0 || k == p.steps_m)
                observePlaneEnergy(state, owner, k, result, baseline);
            if (observer)
                observer(k, state, owner);
        });
    }

    /// Gather each probe once, retaining raw data and one cached incident signal on rank zero.
    inline void extractPlaneResult(const PlaneRow& p, const PlanePulse& pulse,
                                   PlaneFourierProbe& incident, PlaneFourierProbe& other,
                                   PlaneResult& result) {
        const auto total = incident.gather(), transmitted = other.gather();
        if (ippl::Comm->rank() != 0)
            return;
        const double z     = (p.polarization_m == Polarization::TM ? 1 : 0.5) * p.h_m[2];
        const auto samples = pulse.incidentSamples(p.incidentProbe_m, z, p.gateIncident_m);
        for (size_t bin = 0; bin < result.omega_m.size(); ++bin) {
            result.incident_m[bin] = planeIncidentTransform(samples, p.dt_m, result.omega_m[bin]);
            result.total_m[bin]    = total.value(0, 0, bin);
            result.transmitted_m[bin] = transmitted.value(0, 0, bin);
        }
        completePlaneCoefficients(p, result);
    }

    /// Construct the sampled material and periodic box from the same resolved row as the pulse.
    inline InterfaceRun planeRun(const PlaneRow& p, int cadence, int seam) {
        validatePlaneFrequencies(p);
        if (cadence <= 0 || p.steps_m < std::max(p.gateIncident_m, p.gateOther_m))
            throw std::invalid_argument(
                "Plane measurement needs positive cadence and complete gates");
        return InterfaceRun({{p.cellsX_m, p.cellsY_m, p.cellsZ_m}, p.h_m, seam},
                            material::TwoMaterialMap{material::PlanarSlab{0, 0, p.facePosition_m},
                                                     p.epsilon_m, 1.0},
                            {.safetyFactor_m = p.courant_m});
    }

    /** @brief Initialize and measure a caller-owned run, optionally observing its same states.
     * @details The caller may own a field writer whose lifetime ends before the run, including
     * on exceptions. Output observes the same evolution after the full-rate Fourier samples;
     * its cadence cannot remove a Fourier sample. Coefficients are available only on rank zero.
     */
    inline PlaneResult measurePlaneRun(
        InterfaceRun& run, const PlaneRow& p, int cadence,
        std::function<void(InterfaceSolver&, const material::MaterialOwner&)> initialCheck = {},
        std::function<void(InterfaceSolver&, const material::MaterialOwner&)> finalCheck   = {},
        InterfaceRun::Observer observer                                                    = {}) {
        PlanePulse pulse(p);
        run.load(pulse);
        if (initialCheck)
            initialCheck(run.solver(), run.material());
        auto result            = planeResult(p, pulse.modes_m.extent(0));
        const auto frequencies = planeProbeFrequencies(p, result);
        auto incident =
            planeProbe(run.solver(), p, p.incidentProbeIndex_m, p.gateIncident_m, frequencies);
        auto other = planeProbe(run.solver(), p, p.otherProbeIndex_m, p.gateOther_m, frequencies);
        advancePlaneMeasurement(run, p, result, incident, other, cadence, observer);
        extractPlaneResult(p, pulse, incident, other, result);
        if (finalCheck)
            finalCheck(run.solver(), run.material());
        return result;
    }

    /// Own a complete measurement when no external observer needs to borrow its fields.
    inline PlaneResult measurePlane(
        const PlaneRow& p, int cadence = 1, int seam = 0,
        std::function<void(InterfaceSolver&, const material::MaterialOwner&)> initialCheck = {},
        std::function<void(InterfaceSolver&, const material::MaterialOwner&)> finalCheck   = {}) {
        auto run = planeRun(p, cadence, seam);
        return measurePlaneRun(run, p, cadence, initialCheck, finalCheck);
    }
}  // namespace chdr::apps::planarInterface
