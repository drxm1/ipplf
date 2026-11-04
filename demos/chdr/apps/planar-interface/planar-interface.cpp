/** @file planar-interface.cpp
 * @brief Existing periodic interface cases using common input, observations and run records.
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <vector>

#include "apps/planar-interface/InterfaceRun.h"
#include "apps/planar-interface/PacketCase.h"
#include "apps/planar-interface/PlaneInput.h"
#include "apps/planar-interface/PlanePulse.h"
#include "core/input/Application.h"
#include "core/input/RunInput.h"
#include "core/output/RunRecord.h"
#include "opt-in/output/FieldRunOutput.h"

namespace chdr::apps::planarInterface {
    /// Translate the existing command line into the same input used by YAML runs.
    input::RunInput legacyInput(const std::vector<std::string>& arguments) {
        if (arguments.empty())
            throw std::invalid_argument("Use CASE [--smoke] [--output DIRECTORY]");
        input::RunInput result;
        result.caseName_m = arguments.front();
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            if (arguments[i] == "--smoke")
                result.parameters_m["smoke"] = "true";
            else if (arguments[i] == "--output" && i + 1 < arguments.size())
                result.outputDirectory_m = arguments[++i];
            else
                throw std::invalid_argument("Use CASE [--smoke] [--output DIRECTORY]");
        }
        if (result.caseName_m == "plane" && arguments.size() != 1)
            throw std::invalid_argument("The legacy plane command accepts no additional options");
        return result;
    }

    /// Preserve renderer metadata and explicitly identify the existing packet's effective data.
    std::string packetMetadata(const PacketCase& selected) {
        const auto& box    = selected.box_m;
        const auto& packet = selected.start_m;
        std::ostringstream text;
        text << std::setprecision(17) << "case=" << selected.name_m
             << " plane=zx cells=" << box.cells_m[0] << ',' << box.cells_m[1] << ','
             << box.cells_m[2] << " steps=" << selected.steps_m << " h=" << box.spacing_m[0] << ','
             << box.spacing_m[1] << ',' << box.spacing_m[2]
             << " S=" << selected.solver_m.safetyFactor_m
             << " epsilon=2.13 slab=" << selected.slab_m.lower_m << ',' << selected.slab_m.upper_m
             << " centre=" << packet.centreX_m << ',' << packet.centreZ_m
             << " widths=" << packet.width_m << ',' << packet.transverseWidth_m
             << " k=" << packet.waveNumber_m << " direction=" << packet.directionX_m << ','
             << packet.directionZ_m << " TM=" << packet.tm_m;
        return text.str();
    }

    /// Save actual solver settings and resolved measurement schedules after successful output.
    void recordPacket(const input::RunInput& request, const PacketCase& selected, InterfaceRun& run,
                      const std::string& metadata, const output::FieldRunSettings& observations) {
        auto description = output::describeRun("planar-interface", selected.box_m, run.solver(),
                                               selected.steps_m, metadata);
        output::recordFieldRunOutput(description, request, observations);
        output::writeRunRecord(request.outputDirectory_m, request, description,
                               run.solver().A_n.getLayout().comm);
    }

    /// Keep output sampling in the manager callback and complete files after successful evolution.
    void runPacketOutput(const input::RunInput& request, const PacketCase& selected,
                         InterfaceRun& run) {
        if (request.outputDirectory_m.empty()) {
            run.advance(selected.steps_m, InterfaceRun::Observer{});
            return;
        }
        const auto observations = output::fieldRunSettings(request, selected.frameEvery_m);
        const auto metadata     = packetMetadata(selected);
        output::FieldRunOutput writer(request.outputDirectory_m, selected.box_m, run.solver(),
                                      run.material(), metadata, request.units_m);
        run.advance(selected.steps_m, [&](long long k, auto& state, const auto& owner) {
            writer.observe(k, state, owner, selected.steps_m, observations);
        });
        writer.finish();
        recordPacket(request, selected, run, metadata, observations);
    }

    /// Execute a packet preset without introducing a second initialization or stepping path.
    void runPacket(const input::RunInput& request) {
        auto selected =
            packetCase(request.caseName_m, input::booleanParameter(request, "smoke", false));
        InterfaceRun run(selected.box_m, material::TwoMaterialMap{selected.slab_m, 2.13, 1.0},
                         selected.solver_m);
        selected.steps_m = input::resolvedSteps(request, run.solver().getDt(), selected.steps_m);
        run.load(selected.start_m);
        runPacketOutput(request, selected, run);
    }

    /// Spectrum sampling remains at every completed state; diagnostic cadence is independent.
    int checkPlaneSchedules(const input::RunInput& request) {
        input::validateObservationNames(request, {"fourier", "diagnostics", "fields"});
        const auto fourier = input::observationSchedule(request, "fourier", {});
        const auto energy =
            input::observationSchedule(request, "diagnostics", {.includeFinal_m = true});
        if (!input::observationEnabled(request, "fourier") || fourier.every_m != 1
            || fourier.first_m != 0 || fourier.last_m != -1)
            throw std::invalid_argument("plane requires Fourier sampling at every completed step");
        if (!input::observationEnabled(request, "diagnostics") || energy.first_m != 0
            || energy.last_m != -1 || !energy.includeFinal_m)
            throw std::invalid_argument(
                "plane requires diagnostics from its initial through final state");
        return energy.every_m;
    }

    /// Plane frames are opt-in; their default cadence is independent of energy and Fourier data.
    output::FieldRunSettings planeOutputSettings(const input::RunInput& request, int cadence) {
        auto result            = output::fieldRunSettings(request, 100);
        result.diagnostics_m   = {.every_m = cadence, .includeFinal_m = true};
        result.fieldsEnabled_m = request.observations_m.contains("fields")
                                 && input::observationEnabled(request, "fields");
        if (result.fieldsEnabled_m && request.outputDirectory_m.empty())
            throw std::invalid_argument("Plane field frames require an output directory");
        return result;
    }

    /// A user duration may extend the evolution, but cannot truncate either resolved probe gate.
    int checkPlaneInput(const input::RunInput& request, PlaneRow& p) {
        p.steps_m = input::resolvedSteps(request, p.dt_m, p.steps_m);
        if (p.steps_m < std::max(p.gateIncident_m, p.gateOther_m))
            throw std::invalid_argument("plane duration must include both complete probe gates");
        return checkPlaneSchedules(request);
    }

    /// Print measured values; continuum and lattice references belong to the comparison utility.
    void printPlaneResult(const PlaneRow& p, const PlaneResult& result) {
        if (ippl::Comm->rank() != 0)
            return;
        std::cout << std::setprecision(17) << "r=" << result.r_m[0] << " t=" << result.t_m[0]
                  << " omega0=" << p.omega0_m << " angle=" << p.actualAngleDegrees_m
                  << " cells_per_wavelength=" << p.cellsPerWavelength_m
                  << " energy drift=" << result.drift_m << '\n';
        if (!std::cout)
            throw std::runtime_error("Failed to write plane measurement");
    }

    /// Publish raw values in code DTFT units; only the compatibility frequency column is scaled.
    void writePlaneBin(std::ostream& file, const input::RunInput& request,
                       const PlaneResult& result, size_t bin) {
        file << result.omega_m[bin] * request.units_m.frequencyFactor() << ','
             << result.r_m[bin].real() << ',' << result.r_m[bin].imag() << ','
             << result.t_m[bin].real() << ',' << result.t_m[bin].imag() << ','
             << result.omega_m[bin] << ',' << result.incident_m[bin].real() << ','
             << result.incident_m[bin].imag() << ',' << result.total_m[bin].real() << ','
             << result.total_m[bin].imag() << ',' << result.transmitted_m[bin].real() << ','
             << result.transmitted_m[bin].imag() << ',' << (result.valid(bin) ? 1 : 0) << ','
             << result.reason_m[bin] << '\n';
    }

    /// Retain the original scaled frequency and coefficient names while adding raw audit columns.
    void writePlaneResult(const input::RunInput& request, const PlaneResult& result) {
        if (request.outputDirectory_m.empty() || ippl::Comm->rank() != 0)
            return;
        std::ofstream file(request.outputDirectory_m + "/spectra.csv");
        file << std::setprecision(17)
             << (request.units_m.isSI() ? "angular_frequency_rad_per_s" : "angular_frequency")
             << ",r_real,r_imag,t_real,t_imag,omega_code,incident_real,incident_imag,"
                "total_real,total_imag,transmitted_real,transmitted_imag,valid,reason\n";
        for (size_t bin = 0; bin < result.omega_m.size(); ++bin)
            writePlaneBin(file, request, result, bin);
        file.close();
        if (!file)
            throw std::runtime_error("Failed to write plane spectra");
    }

    /// Record the requested and realized spectral parameters in normalized code units.
    void recordPlaneSpectrum(conduit_cpp::Node node, const PlaneRow& p) {
        node["polarization"]                   = p.polarization_m == Polarization::TE ? "TE" : "TM";
        node["incidence"]                      = p.incidentFromVacuum_m ? "vacuum" : "dielectric";
        node["epsilon"]                        = p.epsilon_m;
        node["requested_angle_degrees"]        = p.requestedAngleDegrees_m;
        node["actual_angle_degrees"]           = p.actualAngleDegrees_m;
        node["requested_cells_per_wavelength"] = p.requestedCellsPerWavelength_m;
        node["cells_per_wavelength"]           = p.cellsPerWavelength_m;
        node["omega0"]                         = p.omega0_m;
        node["kz"]                             = p.kz_m;
        node["kx0"]                            = p.kx0_m;
        node["sigma_k"]                        = p.sigmaK_m;
        node["tail_sigmas"]                    = p.tailSigmas_m;
        node["cutoff_fraction"]                = p.cutoffFraction_m;
        node["frequency_units"]                = "normalized";
        node["angle_definition"] = "continuum phase angle relative to interface normal";
    }

    /// Retain both the requested face and the transition produced by inclusive cell-centre
    /// sampling.
    void recordPlaneMaterial(conduit_cpp::Node& node, const PlaneRow& p) {
        node["face_position"]             = p.facePosition_m;
        node["represented_face_position"] = p.dielectricCells_m * p.h_m[0];
        node["face_position_error"]       = p.dielectricCells_m * p.h_m[0] - p.facePosition_m;
        node["material_sampling"]         = "cell_centres_inclusive_upper_face";
        node["cell_size"]                 = p.h_m[0];
        node["dielectric_cells"]          = p.dielectricCells_m;
    }

    /// Record component positions and the actual finite observation intervals.
    void recordPlaneGeometry(conduit_cpp::Node node, const PlaneRow& p) {
        node["geometry"]         = "aligned_plane";
        const bool tm            = p.polarization_m == Polarization::TM;
        node["coordinate_units"] = "normalized";
        node["incidence_plane"]  = "xz";
        node["component"]        = tm ? "H_y" : "E_y";
        recordPlaneMaterial(node, p);
        node["incident_probe_x"]  = p.incidentProbe_m;
        node["other_probe_x"]     = p.otherProbe_m;
        node["probe_z"]           = (tm ? 1 : 0.5) * p.h_m[2];
        node["gate_incident"]     = p.gateIncident_m;
        node["gate_other"]        = p.gateOther_m;
        node["timing_samples"]    = p.timingSamples_m;
        node["separation_factor"] = p.separationFactor_m;
    }

    /// Spectral schema and sample labels suffice to reproduce every raw-to-face conversion.
    void recordPlaneBand(conduit_cpp::Node node, const PlaneRow& p) {
        node["schema_version"] = 1;
        node["central_bin"]    = 0;
        node["angular_frequencies"].set(p.omega_m);
        node["frequency_grid"]       = "uniform_normal_wavenumber";
        node["band_sigmas"]          = p.bandSigmas_m;
        node["incident_floor"]       = p.incidentFloor_m;
        node["cutoff_k"]             = p.cutoffK_m;
        node["pulse_centre"]         = p.pulseCentre_m;
        node["envelope_half_width"]  = p.envelopeHalfWidth_m;
        node["mode_cutoff_sigmas"]   = 9;
        node["incident_probe_index"] = p.incidentProbeIndex_m;
        node["other_probe_index"]    = p.otherProbeIndex_m;
        node["probe_y"]              = (p.polarization_m == Polarization::TM ? 0.5 : 1) * p.h_m[1];
        node["gate_first"]           = 0;
        node["field_time_offset_steps"] = -0.5;
        node["window"]                  = "rectangular";
    }

    /// Raw and coefficient phase origins differ; record both explicitly in all output unit modes.
    void recordPlaneSpectrumUnits(conduit_cpp::Node node, const input::RunInput& request) {
        node["amplitudes_units"]         = "dimensionless";
        node["angular_frequency_units"]  = request.units_m.isSI() ? "rad/s" : "normalized";
        node["raw_units"]                = "normalized";
        node["transform"]                = "dt sum w f exp(+i omega t)";
        node["raw_phase_origin"]         = "probe_coordinates";
        node["coefficient_phase_origin"] = "requested_face";
    }

    /// VTK raw B has integer time; the magnetic Fourier channel instead averages adjacent times.
    void recordPlaneFields(conduit_cpp::Node node, const InterfaceSolver& solver,
                           const input::RunInput& request) {
        node["plane"]                      = "zx";
        node["electric_field"]             = "E_raw";
        node["magnetic_field"]             = "B_raw";
        node["electric_time_offset_steps"] = -0.5;
        node["magnetic_time_offset_steps"] = 0.0;
        node["coordinate_units"]           = request.units_m.isSI() ? "m" : "normalized";
        node["time_units"]                 = request.units_m.isSI() ? "s" : "normalized";
        node["component_offset_units"]     = "cell_fraction";
        node["component_order"]            = "E_x,E_y,E_z,B_x,B_y,B_z";
        const auto offsets                 = output::fieldOffsets(solver.levelContract());
        node["component_offsets"].set(std::vector<double>(offsets.begin(), offsets.end()));
    }

    /// Preserve the renderer's x-z view and select the y component for either polarization.
    std::string planeMetadata(const PlaneRow& p) {
        return std::string("case=plane plane=zx polarization=")
               + (p.polarization_m == Polarization::TM ? "TM" : "TE");
    }

    /// Describe the measured case while its manager still owns the solver and layout.
    output::RunRecordDescription describePlane(const PlaneRow& p, InterfaceSolver& solver,
                                               const input::RunInput& request,
                                               const output::FieldRunSettings& settings) {
        const simulation::RunBox box{{p.cellsX_m, p.cellsY_m, p.cellsZ_m}, p.h_m};
        auto result =
            output::describeRun("planar-interface", box, solver, p.steps_m, planeMetadata(p));
        recordPlaneSpectrum(result.effective_m["plane"], p);
        recordPlaneGeometry(result.effective_m["plane"], p);
        recordPlaneBand(result.effective_m["plane"], p);
        output::recordObservation(result, "fourier", {});
        output::recordFieldRunOutput(result, request, settings);
        recordPlaneSpectrumUnits(result.effective_m["output/spectra"], request);
        recordPlaneFields(result.effective_m["output/fields"], solver, request);
        return result;
    }

    /// Finish and destroy output before its borrowed run, also preserving that order on failure.
    PlaneResult measureRecordedPlane(const input::RunInput& request, const PlaneRow& p,
                                     const output::FieldRunSettings& settings,
                                     output::RunRecordDescription& description) {
        auto run = planeRun(p, settings.diagnostics_m.every_m, 0);
        std::optional<output::FieldRunOutput> writer;
        if (!request.outputDirectory_m.empty())
            writer.emplace(request.outputDirectory_m,
                           simulation::RunBox{{p.cellsX_m, p.cellsY_m, p.cellsZ_m}, p.h_m},
                           run.solver(), run.material(), planeMetadata(p), request.units_m);
        const auto result =
            measurePlaneRun(run, p, settings.diagnostics_m.every_m, {}, {},
                            [&](long long k, auto& state, const auto& owner) {
                                if (writer)
                                    writer->observe(k, state, owner, p.steps_m, settings);
                            });
        if (writer)
            writer->finish();
        description = describePlane(p, run.solver(), request, settings);
        return result;
    }

    /// Run the existing plane measurement and save its result through the common record writer.
    void runPlane(const input::RunInput& request) {
        auto p              = planeRow(request);
        const int cadence   = checkPlaneInput(request, p);
        const auto settings = planeOutputSettings(request, cadence);
        output::RunRecordDescription description;
        const auto result = measureRecordedPlane(request, p, settings, description);
        description.effective_m["plane/incident_peak"] = result.incidentPeak_m;
        description.effective_m["plane/mode_count"]    = static_cast<conduit_int64>(result.modes_m);
        printPlaneResult(p, result);
        writePlaneResult(request, result);
        output::writeRunRecord(request.outputDirectory_m, request, description,
                               ippl::Comm->getCommunicator());
    }

    /// Validate case-specific keys before constructing any distributed field state.
    void runConfigured(const input::RunInput& request) {
        if (request.caseName_m == "plane")
            return runPlane(request);
        input::validateParameters(request, {"smoke"});
        input::validateObservationNames(request, {"fields", "diagnostics"});
        if (request.outputDirectory_m.empty())
            for (const auto& [name, options] : request.observations_m)
                if (input::observationEnabled(request, name))
                    throw std::invalid_argument("Packet observation '" + name
                                                + "' requires an output directory");
        runPacket(request);
    }
}  // namespace chdr::apps::planarInterface

int main(int argc, char** argv) {
    return chdr::input::runApplication(argc, argv, chdr::apps::planarInterface::legacyInput,
                                       chdr::apps::planarInterface::runConfigured);
}
