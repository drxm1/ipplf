/** @file PlanePulse.h
 * @brief Exact lattice plane packets and incident-subtracted Fourier measurements.
 */
#pragma once
#include <complex>
#include <functional>

#include "InterfaceRun.h"
#include "opt-in/diagnostics/FieldEnergy.h"
#include "opt-in/diagnostics/FourierProbe.h"

namespace chdr::interface {
    /** @brief Configuration of a normal-incidence example in units of the x cell width.
     * The long layers prevent periodic images from reaching the probes during their gates.
     * Each rectangular gate includes completed states zero through gate-1.
     * Integers follow the latest-arrival recipe for every contaminating periodic path at
     * its lattice group velocity; keep them synchronized with normalIncidenceRun in the tests.
     */
    struct NormalPlane {
        int cellsX_m = 812, cellsY_m = 2, cellsZ_m = 2, dielectricCells_m = 211;
        std::array<double, 3> h_m{1, 1, 1};
        double epsilon_m = 2.13, courant_m = 0.9, ky_m = 0, dt_m = 0.9 / std::sqrt(3.0);
        double omega0_m = 2 * Kokkos::numbers::pi / (20 * std::sqrt(2.13));
        double kx0_m = 2 * std::asin(std::sin(omega0_m * dt_m / 2) / dt_m), sigmaK_m = kx0_m / 8;
        double pulseCentre_m = 214 + 8 / sigmaK_m, incidentProbe_m = 214.5, otherProbe_m = 207.5;
        int incidentProbeIndex_m = 214, otherProbeIndex_m = 207, gateIncident_m = 1175,
            gateOther_m = 1178, steps_m = 1180;
        bool incidentFromVacuum_m = true;
    };
    /** @brief Device sampler of one signed set of exact periodic lattice modes.
     * @details Each carrier obeys the Standard leapfrog dispersion relation
     * [fallahi2020mithra20fullwavesimulation, section 3.1.3, Eq. (3.17)]. TE samples \f$A_z\f$;
     * TM samples \f$\varepsilon_{\mathrm{inc}}\f$ times the adjoint curl of an edge scalar, divided
     * by the harmonic face epsilon. Its weighted backward divergence vanishes by commuting
     * differences [ryu2016potentialfdtd, Appendix A; chew2014generalizedgauge, Eq. (18)].
     * Mode storage resides in the default execution space; initialization mirrors it once.
     */
    struct PlanePulse {
        Kokkos::View<double* [3]> modes_m;
        double ky_m, centre_m, dt_m, hx_m, hy_m, face_m, length_m, epsilon_m, incidentEpsilon_m;
        bool tm_m;
        template <class Parameters>
        PlanePulse(const Parameters& p, bool tm)
            : ky_m(p.ky_m)
            , centre_m(p.pulseCentre_m)
            , dt_m(p.dt_m)
            , hx_m(p.h_m[0])
            , hy_m(p.h_m[1])
            , face_m(p.dielectricCells_m * p.h_m[0])
            , length_m(p.cellsX_m * p.h_m[0])
            , epsilon_m(p.epsilon_m)
            , incidentEpsilon_m(p.incidentFromVacuum_m ? 1 : p.epsilon_m)
            , tm_m(tm) {
            std::vector<std::array<double, 3>> modes;
            for (int m = 1; m < p.cellsX_m / 2; ++m) {
                const double k = 2 * Kokkos::numbers::pi * m / (p.cellsX_m * p.h_m[0]);
                if (std::abs(k - p.kx0_m) >= 9 * p.sigmaK_m)
                    continue;
                modes.push_back({p.incidentFromVacuum_m ? -k : k, frequency(k),
                                 std::exp(-0.5 * std::pow((k - p.kx0_m) / p.sigmaK_m, 2))});
            }
            modes_m   = decltype(modes_m)("Plane pulse modes", modes.size());
            auto host = Kokkos::create_mirror_view(modes_m);
            for (size_t m = 0; m < modes.size(); ++m)
                for (int c = 0; c < 3; ++c)
                    host(m, c) = modes[m][c];
            Kokkos::deep_copy(modes_m, host);
        }
        /// Incident-medium lattice frequency [fallahi2020mithra20fullwavesimulation, Eq. (3.17)].
        double frequency(double k) const {
            const double Kx = 2 * std::sin(k * hx_m / 2) / hx_m,
                         Ky = 2 * std::sin(ky_m * hy_m / 2) / hy_m;
            return 2 / dt_m
                   * std::asin(dt_m * std::sqrt((Kx * Kx + Ky * Ky) / incidentEpsilon_m) / 2);
        }
        KOKKOS_INLINE_FUNCTION double psi(double x, double y, double t) const {
            double result = 0;
            for (size_t m = 0; m < modes_m.extent(0); ++m)
                result +=
                    modes_m(m, 2)
                    * Kokkos::cos(modes_m(m, 0) * (x - centre_m) + ky_m * y - modes_m(m, 1) * t);
            return result;
        }
        KOKKOS_INLINE_FUNCTION double operator()(double x, double y, double, double t,
                                                 unsigned c) const {
            if (!tm_m)
                return c == 3 ? psi(x, y, t) : 0;
            if (c != 1 && c != 2)
                return 0;
            const double a    = epsilon(x - (c == 1 ? hx_m / 2 : 0));
            const double b    = epsilon(x + (c == 1 ? hx_m / 2 : 0));
            const double curl = c == 1 ? (psi(x, y + hy_m / 2, t) - psi(x, y - hy_m / 2, t)) / hy_m
                                       : (psi(x - hx_m / 2, y, t) - psi(x + hx_m / 2, y, t)) / hx_m;
            return incidentEpsilon_m * (1 / a + 1 / b) * curl / 2;
        }
        KOKKOS_INLINE_FUNCTION double epsilon(double x) const {
            x -= Kokkos::floor(x / length_m) * length_m;
            return x < face_m ? epsilon_m : 1;
        }
        /** @brief Closed-form incident field in the same channel and gate as the probe.
         * TE uses the time difference of A; TM uses the mean of two raw curls, hence
         * the cosine half-step factor [ryu2016potentialfdtd, Appendix A].
         */
        std::complex<double> incident(double x, double y, int gate, double omega) const {
            auto modes = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), modes_m);
            std::complex<double> result = 0;
            for (int k = 0; k < gate; ++k) {
                const double time = (k - 0.5) * dt_m;
                double field      = 0;
                for (size_t m = 0; m < modes.extent(0); ++m) {
                    const double w = modes(m, 1), W = 2 * std::sin(w * dt_m / 2) / dt_m;
                    const double phase = modes(m, 0) * (x - centre_m) + ky_m * y - w * time;
                    field += modes(m, 2)
                             * (tm_m ? incidentEpsilon_m * W * W * std::cos(w * dt_m / 2)
                                           * std::cos(phase)
                                     : -W * std::sin(phase));
                }
                result += dt_m * field * std::exp(std::complex<double>(0, omega * time));
            }
            return result;
        }
    };
    /** @brief Complex amplitudes at the interface and the largest sampled energy drift.
     * Entry 0 is at \f$\omega_0\f$; entries 1 and 2 are records at the lattice frequencies
     * of \f$k_{x0}\mp3\sigma_k\f$. Each amplitude has its own phase transfer to the face.
     * scalarEnergy_m is the last sampled \f$W_\Phi\f$; modes_m is the selected mode count.
     */
    struct PlaneResult {
        std::array<std::complex<double>, 3> r_m, t_m;
        std::array<double, 3> omega_m;
        double drift_m = 0, absoluteRatio_m = 0, scalarEnergy_m = 0;
        size_t modes_m = 0;
    };
    /** @brief Measure one plane packet; gates exclude all reflected periodic images.
     * The incident is subtracted coherently at its probe, then lattice propagation
     * moves both amplitudes to the face [schneider2010UnderstandingFdtd, section 5.8.2;
     * oskooi2013sources, sections 4.2.2-4.3]. Complex amplitudes exist on rank zero only.
     * The centre and two band records use the same gates; only the centre is asserted.
     *
     * TODO: Too long, needs refactoring.
     *
     * 1. Initialize the wave packed and place a Fourier probe on each side
     * 2. Record selected field channel and a potential invariant drift
     * 3. Subtract the known incident amplitude from the incident-side measurement
     * 4. Normalize by that incident amplitude and transfer both phases from the probes
     * to the interface using lattice wavenumbers
     */
    template <class Parameters>
    PlaneResult measurePlane(
        const Parameters& p, bool tm, int cadence = 1, int seam = 0,
        std::function<void(InterfaceSolver&, const material::MaterialOwner&)> initialCheck = {}) {
        InterfaceRun run(
            {{p.cellsX_m, p.cellsY_m, p.cellsZ_m}, p.h_m, seam},
            material::TwoMaterialMap{material::PlanarSlab{0, 0, p.dielectricCells_m * p.h_m[0]},
                                     p.epsilon_m, 1.0},
            {.safetyFactor_m = p.courant_m});
        PlanePulse pulse(p, tm);
        run.load(pulse);
        auto& s = run.solver();
        if (initialCheck)
            initialCheck(s, run.material());
        const auto channel = tm ? solver::Quantity::Magnetic : solver::Quantity::Electric;
        const std::array<double, 3> kx{p.kx0_m, p.kx0_m - 3 * p.sigmaK_m, p.kx0_m + 3 * p.sigmaK_m};
        PlaneResult result;
        result.omega_m = {p.omega0_m, pulse.frequency(kx[1]), pulse.frequency(kx[2])};
        result.modes_m = pulse.modes_m.extent(0);
        std::vector<double> frequencies;
        for (double omega : result.omega_m)
            frequencies.push_back(omega * p.dt_m / (2 * Kokkos::numbers::pi));
        solver::SolverFourierProbe incident(
            s, {{p.incidentProbeIndex_m, 0, 0}, {p.incidentProbeIndex_m, 0, 0}}, {{channel, 2}},
            frequencies, 0, p.gateIncident_m);
        solver::SolverFourierProbe other(s,
                                         {{p.otherProbeIndex_m, 0, 0}, {p.otherProbeIndex_m, 0, 0}},
                                         {{channel, 2}}, frequencies, 0, p.gateOther_m);
        double baseline = 0;
        run.advance(p.steps_m, [&](int k, auto& state, const auto& owner) {
            incident.sample();
            other.sample();
            if (k % cadence == 0 || k == p.steps_m) {
                const auto energy     = solver::potentialEnergy(state, owner);
                result.scalarEnergy_m = energy.scalar_m.total();
                if (k == 0)
                    baseline = energy.vector_m.total();
                if (!(baseline > 0) || !std::isfinite(energy.vector_m.total()))
                    throw std::runtime_error(
                        "Plane pulse needs finite energy and a positive baseline");
                result.drift_m =
                    std::max(result.drift_m, std::abs(energy.vector_m.total() / baseline - 1));
                result.absoluteRatio_m =
                    std::max(result.absoluteRatio_m, energy.vector_m.absolute_m / baseline);
            }
        });
        const auto reflected = incident.gather(), transmitted = other.gather();
        if (ippl::Comm->rank() == 0) {
            for (size_t f = 0; f < result.omega_m.size(); ++f) {
                const double omega  = result.omega_m[f];
                const auto I        = pulse.incident(p.incidentProbe_m, (tm ? 1 : 0.5) * p.h_m[1],
                                                     p.gateIncident_m, omega);
                const double side   = p.incidentFromVacuum_m ? -1 : 1,
                             epsOut = p.incidentFromVacuum_m ? p.epsilon_m : 1;
                const double W      = 2 * std::sin(omega * p.dt_m / 2) / p.dt_m,
                             Ky     = 2 * std::sin(p.ky_m * p.h_m[1] / 2) / p.h_m[1];
                const double kt     = side * 2 / p.h_m[0]
                                  * std::asin(p.h_m[0] * std::sqrt(epsOut * W * W - Ky * Ky) / 2);
                const double ki = side * kx[f], face = p.dielectricCells_m * p.h_m[0];
                result.r_m[f] =
                    (reflected.value(0, 0, f) / I - 1.0)
                    * std::exp(std::complex<double>(0, -2 * ki * (face - p.incidentProbe_m)));
                result.t_m[f] =
                    transmitted.value(0, 0, f) / I
                    * std::exp(std::complex<double>(
                        0, kt * (face - p.otherProbe_m) - ki * (face - p.incidentProbe_m)));
            }
        }
        return result;
    }
}  // namespace chdr::interface
