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
        input::validateObservationNames(request, {"fourier", "diagnostics"});
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

    /// Save existing complex amplitudes; only their frequency coordinate has physical units.
    void writePlaneResult(const input::RunInput& request, const PlaneResult& result) {
        if (request.outputDirectory_m.empty() || ippl::Comm->rank() != 0)
            return;
        std::ofstream file(request.outputDirectory_m + "/spectra.csv");
        file << std::setprecision(17)
             << (request.units_m.isSI() ? "angular_frequency_rad_per_s" : "angular_frequency")
             << ",r_real,r_imag,t_real,t_imag\n";
        for (std::size_t f = 0; f < result.omega_m.size(); ++f)
            file << result.omega_m[f] * request.units_m.frequencyFactor() << ','
                 << result.r_m[f].real() << ',' << result.r_m[f].imag() << ','
                 << result.t_m[f].real() << ',' << result.t_m[f].imag() << '\n';
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

    /// Record component positions and the actual finite observation intervals.
    void recordPlaneGeometry(conduit_cpp::Node node, const PlaneRow& p) {
        node["geometry"]          = "aligned_plane";
        const bool tm             = p.polarization_m == Polarization::TM;
        node["coordinate_units"]  = "normalized";
        node["incidence_plane"]   = "xz";
        node["component"]         = tm ? "H_y" : "E_y";
        node["face_position"]     = p.dielectricCells_m * p.h_m[0];
        node["incident_probe_x"]  = p.incidentProbe_m;
        node["other_probe_x"]     = p.otherProbe_m;
        node["probe_z"]           = (tm ? 1 : 0.5) * p.h_m[2];
        node["gate_incident"]     = p.gateIncident_m;
        node["gate_other"]        = p.gateOther_m;
        node["timing_samples"]    = p.timingSamples_m;
        node["separation_factor"] = p.separationFactor_m;
    }

    /// Describe the measured case while its manager still owns the solver and layout.
    output::RunRecordDescription describePlane(const PlaneRow& p, InterfaceSolver& solver,
                                               const input::RunInput& request, int cadence) {
        const simulation::RunBox box{{p.cellsX_m, p.cellsY_m, p.cellsZ_m}, p.h_m};
        auto result =
            output::describeRun("planar-interface", box, solver, p.steps_m, "case=plane plane=zx");
        recordPlaneSpectrum(result.effective_m["plane"], p);
        recordPlaneGeometry(result.effective_m["plane"], p);
        output::recordObservation(result, "fourier", {});
        output::recordObservation(result, "diagnostics",
                                  {.every_m = cadence, .includeFinal_m = true});
        result.effective_m["output/spectra/amplitudes_units"] = "dimensionless";
        result.effective_m["output/spectra/angular_frequency_units"] =
            request.units_m.isSI() ? "rad/s" : "normalized";
        return result;
    }

    /// Run the existing plane measurement and save its result through the common record writer.
    void runPlane(const input::RunInput& request) {
        auto p            = planeRow(request);
        const int cadence = checkPlaneInput(request, p);
        output::prepareRunDirectory(request.outputDirectory_m, ippl::Comm->getCommunicator());
        std::optional<output::RunRecordDescription> description;
        const auto result = measurePlane(p, cadence, 0, {}, [&](auto& solver, const auto&) {
            description = describePlane(p, solver, request, cadence);
        });
        printPlaneResult(p, result);
        writePlaneResult(request, result);
        output::writeRunRecord(request.outputDirectory_m, request, *description,
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
