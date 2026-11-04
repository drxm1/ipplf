/** @file PlanePulse.h
 * @brief Exact lattice plane packets and incident-subtracted Fourier measurements.
 */
#pragma once
#include <complex>
#include <functional>

#include "apps/planar-interface/InterfaceRun.h"
#include "apps/planar-interface/PlaneSetup.h"
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
        double kz_m, centre_m, dt_m, hx_m, hz_m, face_m, length_m, epsilon_m, incidentEpsilon_m;
        bool tm_m;

        explicit PlanePulse(const PlaneRow& p)
            : kz_m(p.kz_m)
            , centre_m(p.pulseCentre_m)
            , dt_m(p.dt_m)
            , hx_m(p.h_m[0])
            , hz_m(p.h_m[2])
            , face_m(p.dielectricCells_m * p.h_m[0])
            , length_m(p.cellsX_m * p.h_m[0])
            , epsilon_m(p.epsilon_m)
            , incidentEpsilon_m(planeMedia(p)[0])
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

        /// Periodic cell material; component samples adjoining the aligned face use cell centres.
        KOKKOS_INLINE_FUNCTION double epsilon(double x) const {
            x -= Kokkos::floor(x / length_m) * length_m;
            return x < face_m ? epsilon_m : 1;
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

        /** @brief Analytic incident spectrum in precisely the measured channel and gate.
         * The labels and differences follow [ryu2016potentialfdtd, Appendix A].
         * Incident subtraction uses [oskooi2013sources, sections 4.2.2-4.3].
         */
        std::complex<double> incident(double x, double z, int gate, double omega) const {
            const auto modes = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), modes_m);
            std::complex<double> result = 0;
            for (int k = 0; k < gate; ++k) {
                const double time = (k - 0.5) * dt_m;
                result += dt_m * incidentField(modes, x, z, time)
                          * std::exp(std::complex<double>(0, omega * time));
            }
            return result;
        }
    };

    /** @brief Complex face amplitudes and sampled invariant diagnostics.
     * The three frequencies are the central value and the lattice frequencies of
     * \f$k_{x0}\mp3\sigma_k\f$. Scalar energy and the vector invariant retain
     * their existing meanings; they are not continuum transmitted/reflected power.
     */
    struct PlaneResult {
        std::array<std::complex<double>, 3> r_m, t_m;
        std::array<double, 3> omega_m;
        double drift_m = 0, absoluteRatio_m = 0, scalarEnergy_m = 0;
        size_t modes_m = 0;
    };

    /// Construct the three existing comparison bins from the same resolved incident dispersion.
    inline PlaneResult planeResult(const PlaneRow& p, const PlanePulse& pulse) {
        PlaneResult result;
        result.omega_m = {p.omega0_m, planeFrequency(p, p.kx0_m - 3 * p.sigmaK_m),
                          planeFrequency(p, p.kx0_m + 3 * p.sigmaK_m)};
        result.modes_m = pulse.modes_m.extent(0);
        return result;
    }

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
                                        int cadence) {
        double baseline = 0;
        run.advance(p.steps_m, [&](long long k, auto& state, const auto& owner) {
            incident.sample();
            other.sample();
            if (k % cadence == 0 || k == p.steps_m)
                observePlaneEnergy(state, owner, k, result, baseline);
        });
    }

    /** @brief Refer one pair of incident-normalised spectra to the physical interface.
     * Phase transfer is the propagation phase difference [schneider2010UnderstandingFdtd,
     * section 5.8.2, Eqs. (5.96)-(5.99)], using the lattice normal roots
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.17)]. The TM denominator is fixed-y H.
     */
    inline std::array<std::complex<double>, 2> planeFaceCoefficients(
        const PlaneRow& p, const PlanePulse& pulse, double omega, std::complex<double> reflected,
        std::complex<double> transmitted) {
        const double z      = (p.polarization_m == Polarization::TM ? 1 : 0.5) * p.h_m[2];
        const auto incident = pulse.incident(p.incidentProbe_m, z, p.gateIncident_m, omega);
        const double side   = p.incidentFromVacuum_m ? -1 : 1;
        const double ki     = side * planeNormalK(p, omega, planeMedia(p)[0]);
        const double kt     = side * planeNormalK(p, omega, planeMedia(p)[1]);
        const double face   = p.dielectricCells_m * p.h_m[0];
        return {(reflected / incident - 1.0)
                    * std::exp(std::complex<double>(0, -2 * ki * (face - p.incidentProbe_m))),
                transmitted / incident
                    * std::exp(std::complex<double>(
                        0, kt * (face - p.otherProbe_m) - ki * (face - p.incidentProbe_m)))};
    }

    /// Gather only once and perform the host coefficient extraction on the root rank.
    inline void extractPlaneResult(const PlaneRow& p, const PlanePulse& pulse,
                                   PlaneFourierProbe& incident, PlaneFourierProbe& other,
                                   PlaneResult& result) {
        const auto reflected = incident.gather(), transmitted = other.gather();
        if (ippl::Comm->rank() != 0)
            return;
        for (size_t f = 0; f < result.omega_m.size(); ++f) {
            const auto coefficient = planeFaceCoefficients(
                p, pulse, result.omega_m[f], reflected.value(0, 0, f), transmitted.value(0, 0, f));
            result.r_m[f] = coefficient[0];
            result.t_m[f] = coefficient[1];
        }
    }

    /// Construct the sampled material and periodic box from the same resolved row as the pulse.
    inline InterfaceRun planeRun(const PlaneRow& p, int cadence, int seam) {
        if (cadence <= 0 || p.steps_m < std::max(p.gateIncident_m, p.gateOther_m))
            throw std::invalid_argument(
                "Plane measurement needs positive cadence and complete gates");
        return InterfaceRun(
            {{p.cellsX_m, p.cellsY_m, p.cellsZ_m}, p.h_m, seam},
            material::TwoMaterialMap{material::PlanarSlab{0, 0, p.dielectricCells_m * p.h_m[0]},
                                     p.epsilon_m, 1.0},
            {.safetyFactor_m = p.courant_m});
    }

    /** @brief Run one resolved row through the common initialization and measurement path.
     * Optional checks run while the solver and material owner remain alive. Complex
     * coefficients are available only on rank zero; no reference values enter the evolution.
     */
    inline PlaneResult measurePlane(
        const PlaneRow& p, int cadence = 1, int seam = 0,
        std::function<void(InterfaceSolver&, const material::MaterialOwner&)> initialCheck = {},
        std::function<void(InterfaceSolver&, const material::MaterialOwner&)> finalCheck   = {}) {
        auto run = planeRun(p, cadence, seam);
        PlanePulse pulse(p);
        run.load(pulse);
        if (initialCheck)
            initialCheck(run.solver(), run.material());
        auto result            = planeResult(p, pulse);
        const auto frequencies = planeProbeFrequencies(p, result);
        auto incident =
            planeProbe(run.solver(), p, p.incidentProbeIndex_m, p.gateIncident_m, frequencies);
        auto other = planeProbe(run.solver(), p, p.otherProbeIndex_m, p.gateOther_m, frequencies);
        advancePlaneMeasurement(run, p, result, incident, other, cadence);
        extractPlaneResult(p, pulse, incident, other, result);
        if (finalCheck)
            finalCheck(run.solver(), run.material());
        return result;
    }
}  // namespace chdr::apps::planarInterface
