/** @file TiltedOutput.h
 * @brief Auditable tilted-plane geometry, raw point spectra and projected face amplitudes.
 */
#pragma once
#include <fstream>
#include <iomanip>

#include "apps/planar-interface/TiltedMeasurement.h"
#include "core/output/RunRecord.h"

namespace chdr::apps::planarInterface {
    /// Common physical selectors retain their existing meanings; kappa is explicitly tangential.
    inline void recordTiltedPhysics(conduit_cpp::Node node, const TiltedRow& r) {
        const auto& p                   = r.plane_m;
        node["schema_version"]          = 2;
        node["geometry"]                = "periodic_tilted_plane";
        node["polarization"]            = p.polarization_m == Polarization::TM ? "TM" : "TE";
        node["component"]               = p.polarization_m == Polarization::TM ? "H_y" : "E_y";
        node["incidence"]               = "vacuum";
        node["epsilon"]                 = p.epsilon_m;
        node["incidence_plane"]         = "xz";
        node["requested_angle_degrees"] = p.requestedAngleDegrees_m;
        node["actual_angle_degrees"]    = p.actualAngleDegrees_m;
        node["requested_cells_per_wavelength"] = p.requestedCellsPerWavelength_m;
        node["cells_per_wavelength"]           = p.cellsPerWavelength_m;
        node["omega0"]                         = p.omega0_m;
        node["kappa"]                          = r.kappa_m;
        node["frequency_units"]                = "normalized";
        node["coordinate_units"]               = "normalized";
        node["angle_definition"] = "continuum phase angle relative to interface normal";
    }

    /// Exact physical tilt and the compatible material/box translations are recorded separately.
    inline void recordTiltedGeometry(conduit_cpp::Node node, const TiltedRow& r) {
        node["interface_tilt_degrees"] = r.tiltDegrees_m;
        node["normal"].set(std::vector<double>{r.cosine_m, 0, r.sine_m});
        node["tangent"].set(std::vector<double>{-r.sine_m, 0, r.cosine_m});
        const auto anchor = tiltedAnchor(r);
        node["anchor"].set(std::vector<double>(anchor.begin(), anchor.end()));
        node["requested_face_position"] = r.facePosition_m;
        node["requested_cell_shape"].set(
            std::vector<double>(r.requestedShape_m.begin(), r.requestedShape_m.end()));
        node["primitive_cells_x"]        = r.primitiveX_m;
        node["primitive_cells_z"]        = r.primitiveZ_m;
        node["period_cells_multiplier"]  = r.multiplier_m;
        node["normal_period"]            = r.normalPeriod_m;
        node["tangent_period"]           = r.tangentPeriod_m;
        node["primitive_tangent_period"] = r.primitivePeriod_m;
        node["normal_dielectric_width"]  = r.dielectricWidth_m;
        node["normal_incident_width"]    = r.incidentWidth_m;
        node["tangent_mode"]             = r.tangentMode_m;
        node["material_sampling"]        = "integer_normal_rows_closed_dielectric_interval";
    }

    /// Initial histories are reconstructible from this source anchor and the saved mode table.
    inline void recordTiltedPulse(conduit_cpp::Node node, const TiltedRow& r) {
        const auto& p = r.plane_m;
        auto source   = tiltedAnchor(r);
        source[0] += r.cosine_m * r.sourceDistance_m;
        source[2] += r.sine_m * r.sourceDistance_m;
        node["source_position"].set(std::vector<double>(source.begin(), source.end()));
        node["source_normal_distance"]          = r.sourceDistance_m;
        node["sigma_k"]                         = p.sigmaK_m;
        node["central_normal_k"]                = p.kx0_m;
        node["grazing_normal_k"]                = p.cutoffK_m;
        node["tail_sigmas"]                     = p.tailSigmas_m;
        node["mode_cutoff_sigmas"]              = p.tailSigmas_m;
        node["cutoff_fraction"]                 = p.cutoffFraction_m;
        node["envelope_half_width"]             = p.envelopeHalfWidth_m;
        node["retained_mode_count"]             = static_cast<conduit_int64>(r.modes_m.size());
        node["discarded_field_weight_fraction"] = r.discardedFieldWeightFraction_m;
        node["mode_file"]                       = "modes.csv";
    }

    /// Each component sample carries its own lift; no scalar represented staircase face is claimed.
    inline void recordTiltedProbes(conduit_cpp::Node node, const TiltedRow& r) {
        node["probe_cells"]            = r.probeCells_m;
        node["probe_half_width_cells"] = r.probeHalfWidth_m;
        const int count = (2 * r.probeHalfWidth_m + 1) * (2 * r.probeHalfWidth_m + 1);
        node["probe_point_count_incident"]    = count;
        node["probe_point_count_transmitted"] = count;
        const auto first = tiltedProbeBox(r, true), second = tiltedProbeBox(r, false);
        node["incident_probe_first"].set(
            std::vector<int>(first.first_m.begin(), first.first_m.end()));
        node["incident_probe_last"].set(std::vector<int>(first.last_m.begin(), first.last_m.end()));
        node["transmitted_probe_first"].set(
            std::vector<int>(second.first_m.begin(), second.first_m.end()));
        node["transmitted_probe_last"].set(
            std::vector<int>(second.last_m.begin(), second.last_m.end()));
        node["projection"]               = "unweighted_complex_plane_wave_fit";
        node["projection_basis"]         = "exp(i*(normal_k*u+kappa*v))";
        node["projection_normalization"] = "mean_over_points";
        node["fit_residual_limit"]       = r.fitResidualLimit_m;
        node["raw_point_file"]           = "probe-spectra.csv";
    }

    /// Actual gates and branch/order qualifications accompany every declared usable bin.
    inline void recordTiltedQualification(conduit_cpp::Node node, const TiltedRow& r) {
        const auto& p       = r.plane_m;
        node["central_bin"] = 0;
        node["angular_frequencies"].set(p.omega_m);
        node["frequency_grid"]          = "uniform_normal_wavenumber";
        node["band_sigmas"]             = p.bandSigmas_m;
        node["incident_floor"]          = p.incidentFloor_m;
        node["floquet_margin"]          = r.floquetMargin_m;
        node["retained_floquet_margin"] = r.retainedMargin_m;
        node["timing_samples"]          = p.timingSamples_m;
        node["separation_factor"]       = p.separationFactor_m;
        node["gate_first"]              = 0;
        node["gate_incident"]           = p.gateIncident_m;
        node["gate_other"]              = p.gateOther_m;
        node["field_time_offset_steps"] = -0.5;
        node["window"]                  = "rectangular";
        node["timing_model"]            = "four_direction_group_velocity_envelope_estimate";
    }

    /// Code-unit DTFT amplitudes and dimensionless ratios remain distinct in SI output mode.
    inline void recordTiltedSpectrumUnits(conduit_cpp::Node node, const input::RunInput& request) {
        node["amplitudes_units"]         = "dimensionless";
        node["face_amplitudes_units"]    = "normalized";
        node["raw_units"]                = "normalized";
        node["angular_frequency_units"]  = request.units_m.isSI() ? "rad/s" : "normalized";
        node["transform"]                = "dt sum w f exp(+i omega t)";
        node["raw_phase_origin"]         = "lifted_probe_coordinates";
        node["coefficient_phase_origin"] = "requested_face_anchor";
    }

    /// Finish one checked CSV stream; the caller handles communicator-wide error propagation.
    inline void finishTiltedFile(std::ofstream& file) {
        file.close();
        if (!file)
            throw std::runtime_error("Failed to finish tilted-plane output");
    }

    /// Publish signed roots beside the measured face fits so a reader need not infer a reflection.
    inline void writeTiltedBin(std::ostream& file, double omega, double frequencyFactor,
                               const TiltedBinResult& bin) {
        file << omega * frequencyFactor << ',' << omega << ',' << bin.incidentRoot_m.normalK_m
             << ',' << bin.reflectedRoot_m.normalK_m << ',' << bin.transmittedRoot_m.normalK_m
             << ',' << bin.incidentRoot_m.groupVelocity_m << ','
             << bin.reflectedRoot_m.groupVelocity_m << ',' << bin.transmittedRoot_m.groupVelocity_m;
        for (const auto value :
             {bin.incident_m, bin.reflected_m, bin.transmitted_m, bin.r_m, bin.t_m})
            file << ',' << value.real() << ',' << value.imag();
        file << ',' << bin.incidentResidual_m << ',' << bin.reflectedResidual_m << ','
             << bin.transmittedResidual_m << ',' << (bin.valid() ? 1 : 0) << ',' << bin.reason_m
             << '\n';
    }

    inline void writeTiltedSpectrum(const input::RunInput& request, const TiltedRow& r,
                                    const TiltedResult& result) {
        std::ofstream file(request.outputDirectory_m + "/spectra.csv");
        file << std::setprecision(17)
             << (request.units_m.isSI() ? "angular_frequency_rad_per_s" : "angular_frequency")
             << ",omega_code,incident_normal_k,reflected_normal_k,transmitted_normal_k,"
                "incident_normal_group_velocity,reflected_normal_group_velocity,transmitted_normal_"
                "group_velocity,"
                "incident_face_real,incident_face_imag,reflected_face_real,reflected_face_imag,"
                "transmitted_face_real,transmitted_face_imag,r_real,r_imag,t_real,t_imag,"
                "fit_residual_incident,fit_residual_r,fit_residual_t,valid,reason\n";
        for (size_t bin = 0; bin < result.bins_m.size(); ++bin)
            writeTiltedBin(file, r.plane_m.omega_m[bin], request.units_m.frequencyFactor(),
                           result.bins_m[bin]);
        finishTiltedFile(file);
    }

    /// Raw point spectra retain the wrapped coordinate and the explicit integer Cartesian lift.
    inline void writeTiltedPointSide(std::ostream& file, const TiltedRow& r,
                                     const TiltedResult& result, size_t bin, bool incident) {
        const auto& probe = incident ? result.incidentProbe_m : result.transmittedProbe_m;
        for (size_t point = 0; point < probe.coordinates_m.size(); ++point) {
            const auto coordinate = probe.coordinates_m[point];
            const auto uv         = tiltedCoordinates(r, coordinate);
            const auto incoming =
                incident ? result.analyticIncident_m[bin * probe.coordinates_m.size() + point]
                         : std::complex<double>{};
            const auto total = probe.value(point, 0, bin);
            file << bin << ',' << r.plane_m.omega_m[bin] << ','
                 << (incident ? "incident" : "transmitted") << ',' << point << ',' << coordinate[0]
                 << ',' << coordinate[1] << ',' << coordinate[2] << ',' << tiltedLift(r, 0) << ','
                 << tiltedLift(r, 2) << ',' << uv[0] << ',' << uv[1] << ',' << incoming.real()
                 << ',' << incoming.imag() << ',' << total.real() << ',' << total.imag() << '\n';
        }
    }

    inline void writeTiltedPoints(const input::RunInput& request, const TiltedRow& r,
                                  const TiltedResult& result) {
        std::ofstream file(request.outputDirectory_m + "/probe-spectra.csv");
        file << std::setprecision(17)
             << "bin,omega_code,side,point,x,y,z,lift_x,lift_z,u,v,incident_real,incident_imag,"
                "total_real,total_imag\n";
        for (size_t bin = 0; bin < result.bins_m.size(); ++bin) {
            writeTiltedPointSide(file, r, result, bin, true);
            writeTiltedPointSide(file, r, result, bin, false);
        }
        finishTiltedFile(file);
    }

    inline void writeTiltedModes(const input::RunInput& request, const TiltedRow& r) {
        std::ofstream file(request.outputDirectory_m + "/modes.csv");
        file << std::setprecision(17) << "kx,kz,omega_code,weight,normal_k\n";
        for (const auto& mode : r.modes_m)
            file << mode[0] << ',' << mode[1] << ',' << mode[2] << ',' << mode[3] << ',' << mode[4]
                 << '\n';
        finishTiltedFile(file);
    }

    /// Every rank learns a root-side file failure before the later collective run-record write.
    inline void writeTiltedResult(const input::RunInput& request, const TiltedRow& r,
                                  const TiltedResult& result) {
        if (request.outputDirectory_m.empty())
            return;
        int failed = 0;
        if (ippl::Comm->rank() == 0) {
            try {
                writeTiltedSpectrum(request, r, result);
                writeTiltedPoints(request, r, result);
                writeTiltedModes(request, r);
            } catch (const std::exception&) {
                failed = 1;
            }
        }
        MPI_Bcast(&failed, 1, MPI_INT, 0, ippl::Comm->getCommunicator());
        if (failed)
            throw std::runtime_error("Failed to write tilted-plane spectral output");
    }
}  // namespace chdr::apps::planarInterface
