/** @file planar-interface.cpp
 * @brief Existing periodic interface cases using common input, observations and run records.
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

#include "apps/planar-interface/InterfaceRun.h"
#include "apps/planar-interface/PacketCase.h"
#include "apps/planar-interface/PlanePulse.h"
#include "core/input/Application.h"
#include "core/input/RunInput.h"
#include "core/output/RunRecord.h"
#include "opt-in/output/FieldRunOutput.h"

namespace chdr::apps::planarInterface {
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

    /// Preserve the fixed probe gates and reference case while accepting its full duration in SI.
    int checkPlaneInput(const input::RunInput& request, const NormalPlane& p) {
        input::validateObservationNames(request, {"fourier", "diagnostics"});
        if (input::booleanParameter(request, "smoke", false)
            || input::resolvedSteps(request, p.dt_m, p.steps_m) != p.steps_m)
            throw std::invalid_argument(
                "plane requires its existing mesh and complete probe gates");
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

    /// Preserve the existing dimensionless comparison printed by the plane case.
    void printPlaneResult(const NormalPlane& p, const PlaneResult& result) {
        if (ippl::Comm->rank() != 0)
            return;
        const double n = std::sqrt(p.epsilon_m), r = (1 - n) / (1 + n);
        const double kt = 2 * std::asin(n * std::sin(p.omega0_m * p.dt_m / 2) / p.dt_m);
        // Sine matching at the two adjacent cell centres gives the exact TE lattice
        // coefficients [yee1966MaxwellIsotropicMedia, Eqs. (5)-(6); deduction].
        std::cout << std::setprecision(17) << "r=" << result.r_m[0] << " t=" << result.t_m[0]
                  << " Fresnel=" << r << ',' << 1 + r
                  << " lattice=" << std::sin((p.kx0_m - kt) / 2) / std::sin((p.kx0_m + kt) / 2)
                  << ',' << std::sin(p.kx0_m) / std::sin((p.kx0_m + kt) / 2)
                  << " energy drift=" << result.drift_m << '\n';
        if (!std::cout)
            throw std::runtime_error("Failed to write plane comparison");
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

    /// Describe the measured case while its manager still owns the solver and its layout.
    output::RunRecordDescription describePlane(const NormalPlane& p, InterfaceSolver& solver,
                                               const input::RunInput& request, int cadence) {
        const simulation::RunBox box{{p.cellsX_m, p.cellsY_m, p.cellsZ_m}, p.h_m};
        auto result = output::describeRun("planar-interface", box, solver, p.steps_m,
                                          "case=plane polarization=TE epsilon=2.13");
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
        const NormalPlane p;
        const int cadence = checkPlaneInput(request, p);
        output::prepareRunDirectory(request.outputDirectory_m, ippl::Comm->getCommunicator());
        std::optional<output::RunRecordDescription> description;
        const auto result = measurePlane(p, false, cadence, 0, {}, [&](auto& solver, const auto&) {
            description = describePlane(p, solver, request, cadence);
        });
        printPlaneResult(p, result);
        writePlaneResult(request, result);
        output::writeRunRecord(request.outputDirectory_m, request, *description,
                               ippl::Comm->getCommunicator());
    }

    /// Validate case-specific keys before constructing any distributed field state.
    void runConfigured(const input::RunInput& request) {
        input::validateParameters(request, {"smoke"});
        if (request.caseName_m == "plane")
            return runPlane(request);
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
    return chdr::input::runApplication(argc, argv, chdr::apps::planarInterface::runConfigured);
}
