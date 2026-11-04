/** @file PlaneSpectrum.h
 * @brief Raw plane-probe spectra, explicit bin validity and face-referred coefficients.
 */
#pragma once
#include <complex>
#include <string>

#include "apps/planar-interface/PlaneSetup.h"

namespace chdr::apps::planarInterface {
    /// Validate the requested observation frequencies before constructing distributed fields.
    inline void validatePlaneFrequencies(const PlaneRow& p) {
        if (!std::isfinite(p.dt_m) || p.dt_m <= 0 || p.omega_m.empty()
            || p.omega_m.front() != p.omega0_m)
            throw std::invalid_argument(
                "Plane frequencies need a valid time step and central bin zero");
        for (double omega : p.omega_m)
            if (!std::isfinite(omega) || omega <= 0 || omega * p.dt_m >= Kokkos::numbers::pi)
                throw std::invalid_argument(
                    "Plane observation frequencies must lie in the principal time branch");
    }

    /** @brief Measured spectra and separately derived amplitude ratios.
     * @details Raw spectra retain their physical probe positions. Ratios refer to the
     * requested interface. Invalid ratios remain NaN and carry a reason; the corresponding
     * raw data are preserved. The central frequency is bin zero. Complex data are available
     * only on rank zero. Invariant diagnostics retain their existing normalization.
     */
    struct PlaneResult {
        std::vector<double> omega_m;
        std::vector<std::complex<double>> r_m, t_m, incident_m, total_m, transmitted_m;
        std::vector<std::string> reason_m;
        double incidentPeak_m = 0, drift_m = 0, absoluteRatio_m = 0, scalarEnergy_m = 0;
        size_t modes_m = 0;

        /// A bin is usable only after every domain and normalization check succeeded.
        bool valid(size_t bin) const { return reason_m.at(bin).empty(); }
    };

    /// Allocate result columns without allowing unevaluated coefficients to look like zero.
    inline PlaneResult planeResult(const PlaneRow& p, size_t modes) {
        PlaneResult result;
        result.omega_m   = p.omega_m;
        result.modes_m   = modes;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (auto* column :
             {&result.r_m, &result.t_m, &result.incident_m, &result.total_m, &result.transmitted_m})
            column->assign(p.omega_m.size(), {nan, nan});
        result.reason_m.assign(p.omega_m.size(), "not_evaluated");
        return result;
    }

    /** @brief Rectangular DTFT of cached half-step-labelled incident samples.
     * This uses the same dt-weighted sum as FourierProbe, without a \f$1/(2\pi)\f$ factor.
     * The convention follows [oskooi2013sources, section 4.3, Eq. (4.6)], with its
     * common normalization omitted consistently from measured and incident spectra.
     */
    inline std::complex<double> planeIncidentTransform(const std::vector<double>& samples,
                                                       double dt, double omega) {
        std::complex<double> result = 0;
        for (size_t k = 0; k < samples.size(); ++k) {
            const double time = (static_cast<double>(k) - 0.5) * dt;
            result += dt * samples[k] * std::exp(std::complex<double>(0, omega * time));
        }
        return result;
    }

    /// Reject non-finite raw data before normalization or coefficient phase transfer.
    inline bool finitePlaneSpectrum(std::complex<double> value) {
        return std::isfinite(value.real()) && std::isfinite(value.imag());
    }

    /// Strict lattice propagation domain [fallahi2020mithra20fullwavesimulation, Eq. (3.17)].
    inline std::string planeLatticePropagationReason(const PlaneRow& p, double omega) {
        const auto media = planeMedia(p);
        const double W   = 2 * std::sin(omega * p.dt_m / 2) / p.dt_m;
        const double Kz  = planeSpatialSymbol(p.kz_m, p.h_m[2]);
        for (size_t side = 0; side < media.size(); ++side) {
            const double sine2 = p.h_m[0] * p.h_m[0] * (media[side] * W * W - Kz * Kz) / 4;
            if (!(sine2 > 0 && sine2 < 1))
                return side == 0 ? "lattice_incident_cutoff" : "lattice_transmitted_cutoff";
        }
        return {};
    }

    /** @brief Per-bin continuum and discrete propagation domains, independently of amplitude.
     * Continuum roots require \f$\varepsilon\omega^2-k_z^2>0\f$
     * [ellingson2020ElectromagneticsVol2, sections 5.6-5.7]. Lattice roots require
     * \f$0<\sin^2(k_xh_x/2)<1\f$ on the principal temporal branch
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.17); deduction]. A lattice-propagating
     * bin can fail the continuum test on a coarse oblique mesh.
     */
    inline std::string planePropagationReason(const PlaneRow& p, double omega) {
        if (!std::isfinite(omega) || omega <= 0 || omega * p.dt_m >= Kokkos::numbers::pi)
            return "outside_time_branch";
        const auto media = planeMedia(p);
        for (size_t side = 0; side < media.size(); ++side)
            if (!(media[side] * omega * omega > p.kz_m * p.kz_m))
                return side == 0 ? "continuum_incident_cutoff" : "continuum_transmitted_cutoff";
        return planeLatticePropagationReason(p, omega);
    }

    /** @brief Require covered timing and a usable incident normalization for each bin.
     * Dividing by a weak incident spectrum amplifies measurement error; apply an explicit
     * fraction of the largest selected incident magnitude [schneider2010UnderstandingFdtd,
     * section 5.8.2, Eqs. (5.96)-(5.99); conditioning deduction from the spectral ratio].
     * This floor is a diagnostic criterion, not a proof of the remaining error.
     */
    inline std::string planeBinReason(const PlaneRow& p, const PlaneResult& result, size_t bin) {
        if (!finitePlaneSpectrum(result.incident_m[bin])
            || !finitePlaneSpectrum(result.total_m[bin])
            || !finitePlaneSpectrum(result.transmitted_m[bin]))
            return "nonfinite_spectrum";
        if (auto reason = planePropagationReason(p, result.omega_m[bin]); !reason.empty())
            return reason;
        const double k      = planeNormalK(p, result.omega_m[bin], planeMedia(p)[0]);
        const double extent = p.tailSigmas_m * p.sigmaK_m;
        if (std::abs(k - p.kx0_m) > extent * (1 + 64 * std::numeric_limits<double>::epsilon()))
            return "outside_qualified_band";
        if (!(result.incidentPeak_m > 0)
            || std::abs(result.incident_m[bin]) < p.incidentFloor_m * result.incidentPeak_m)
            return "weak_incident";
        return {};
    }

    /** @brief Refer one incident-normalized pair to the requested material face.
     * Phase transfer is the propagation phase difference [schneider2010UnderstandingFdtd,
     * section 5.8.2, Eqs. (5.96)-(5.99)], using the lattice normal roots
     * [fallahi2020mithra20fullwavesimulation, Eq. (3.17)]. The TM basis is fixed-y H.
     * The caller has already checked propagation and nonzero normalization.
     */
    inline std::array<std::complex<double>, 2> planeFaceCoefficients(
        const PlaneRow& p, double omega, std::complex<double> incident, std::complex<double> total,
        std::complex<double> transmitted) {
        const double side = p.incidentFromVacuum_m ? -1 : 1;
        const double ki   = side * planeNormalK(p, omega, planeMedia(p)[0]);
        const double kt   = side * planeNormalK(p, omega, planeMedia(p)[1]);
        const double face = p.facePosition_m;
        return {(total / incident - 1.0)
                    * std::exp(std::complex<double>(0, -2 * ki * (face - p.incidentProbe_m))),
                transmitted / incident
                    * std::exp(std::complex<double>(
                        0, kt * (face - p.otherProbe_m) - ki * (face - p.incidentProbe_m)))};
    }

    /// Apply the domain/conditioning checks before exposing one pair of finite coefficients.
    inline void completePlaneBin(const PlaneRow& p, PlaneResult& result, size_t bin) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        result.r_m[bin] = result.t_m[bin] = {nan, nan};
        result.reason_m[bin]              = planeBinReason(p, result, bin);
        if (!result.valid(bin))
            return;
        const auto coefficient =
            planeFaceCoefficients(p, result.omega_m[bin], result.incident_m[bin],
                                  result.total_m[bin], result.transmitted_m[bin]);
        if (!finitePlaneSpectrum(coefficient[0]) || !finitePlaneSpectrum(coefficient[1])) {
            result.reason_m[bin] = "nonfinite_coefficients";
            return;
        }
        result.r_m[bin] = coefficient[0];
        result.t_m[bin] = coefficient[1];
    }

    /// Derive ratios only for qualified bins; raw spectra are never overwritten by references.
    inline void completePlaneCoefficients(const PlaneRow& p, PlaneResult& result) {
        result.incidentPeak_m = 0;
        for (const auto value : result.incident_m)
            if (finitePlaneSpectrum(value))
                result.incidentPeak_m = std::max(result.incidentPeak_m, std::abs(value));
        for (size_t bin = 0; bin < result.omega_m.size(); ++bin)
            completePlaneBin(p, result, bin);
    }
}  // namespace chdr::apps::planarInterface
