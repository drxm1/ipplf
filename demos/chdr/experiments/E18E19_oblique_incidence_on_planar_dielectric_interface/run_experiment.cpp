/** @file run_experiment.cpp
 * @brief Measure oblique reflection from the vacuum side of a planar dielectric interface.
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "core/input/Application.h"
#include "core/output/RunRecord.h"
#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19InitialSampler.h"
#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19Measurements.h"
#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19Output.h"
#include "experiments/E18E19_oblique_incidence_on_planar_dielectric_interface/E18E19Setup.h"

namespace chdr::experiments::e18e19 {
    /// Save each specified longitudinal profile at its reconstructed field time.
    void writeProfiles(Solver& solver, const Reference& reference, ProfileCsv<Reference>& profiles,
                       long long step) {
        const char* checkpoint = step == 0 ? "initial" : "final";
        const double target    = solver.labelTime(Quantity::VectorPotential, step);
        const auto write       = [&](auto& field, Quantity quantity, unsigned component) {
            profiles.write(field, quantity, component, checkpoint, target, step,
                                 solver.labelTime(quantity, step));
        };
        if (reference.tm_m) {
            write(solver.magneticField(), Quantity::Magnetic, 1);
            write(solver.electricField(), Quantity::Electric, 0);
            write(solver.electricField(), Quantity::Electric, 2);
        } else {
            write(solver.electricField(), Quantity::Electric, 1);
        }
    }

    /// Fourier magnitudes of the numerical field and, for the S cases, the continuum field.
    struct EndpointMagnitudes {
        std::array<double, 12> numerical_m{}, reference_m{};
    };

    /// Project both fields onto the incident or reflected direction at their native samples.
    EndpointMagnitudes sampleEndpoint(Solver& solver, const Reference& reference,
                                      ProfileCsv<Reference>& profiles, bool windowReference) {
        const auto step = solver.completedSolves();
        writeProfiles(solver, reference, profiles, step);
        auto& field         = reference.tm_m ? solver.magneticField() : solver.electricField();
        const auto quantity = reference.tm_m ? Quantity::Magnetic : Quantity::Electric;
        EndpointMagnitudes magnitudes;
        magnitudes.numerical_m = measureFourierMagnitude(
            field, solver.levelContract(), quantity, Reference::KX, reference.measuredWavenumbers(),
            Reference::INTERFACE_Z, step != 0);
        if (windowReference)
            magnitudes.reference_m =
                measureReferenceFourierMagnitude(field, solver.levelContract(), quantity, reference,
                                                 solver.labelTime(quantity, step), step != 0);
        return magnitudes;
    }

    /** @brief Sample both energies and the maximum at the periodic interface at each prescribed
     * step.
     * @details E18.meas03 divides the largest electric component seen on the two outer layers
     * up through this observation by the initial peak in the whole box. The final row thus
     * contains the maximum over the full measurement schedule.
     */
    void advance(Simulation& simulation, const Reference& reference, int cells, int steps,
                 double initialPeak, MeasurementCsv& history, FieldSnapshots& snapshots) {
        double secondFacePeak = 0;
        simulation.advance(steps, [&](long long step, Solver& solver, const auto& material) {
            if (step % cells == 0 || step == steps) {
                auto& electric    = solver.electricField();
                const auto energy = measureSideEnergy(electric, solver.magneticField(), material,
                                                      reference.tm_m, Reference::INTERFACE_Z);
                secondFacePeak    = std::max(secondFacePeak, measureElectricPeak(electric, true));
                history.row(step, solver.labelTime(Quantity::Electric, step), energy[0], energy[1],
                            secondFacePeak / initialPeak);
            }
            snapshots.observe(step, solver);
        });
    }

    /// Attach the spectrum, component conventions and normalization to the common solver record.
    void recordRun(Simulation& simulation, const input::ChdrYamlRunInput& request,
                   const simulation::ChdrCartesianGridSettings& grid, const Reference& reference,
                   const FieldSnapshots& snapshots, double initialPeak, int steps) {
        const int cells  = grid.cells_m[0];
        auto description = output::describeRun("chdr-experiment-e18-e19", grid, simulation.solver(),
                                               steps, request.caseName_m + " oblique incidence");
        auto& effective  = description.effective_m;
        reference.record(effective["analytic_reference"], request.caseName_m);
        effective["material/side2_cell_count"] = simulation.material().statistics().marked_m;
        effective["measurements/reflection"]   = reference.tm_m ? "B_y" : "E_y";
        const auto quantity = reference.tm_m ? Quantity::Magnetic : Quantity::Electric;
        effective["measurements/incident_time"]  = simulation.solver().labelTime(quantity, 0);
        effective["measurements/reflected_time"] = simulation.solver().labelTime(quantity, steps);
        effective["measurements/energy_every_steps"]    = cells;
        effective["measurements/initial_electric_peak"] = initialPeak;
        effective["profiles/line_i"]                    = 0;
        effective["profiles/line_j"]                    = 0;
        snapshots.record(description);
        output::writeRunRecord(request.outputDirectory_m, request, description,
                               simulation.solver().A_n.getLayout().comm);
    }

    /// Close all CSV and field writers before publishing the completed run record.
    void observeRun(Simulation& simulation, const input::ChdrYamlRunInput& request,
                    const simulation::ChdrCartesianGridSettings& grid, const Reference& reference,
                    int steps) {
        auto& solver             = simulation.solver();
        const auto comm          = solver.A_n.getLayout().comm;
        const auto& directory    = request.outputDirectory_m;
        const int cells          = grid.cells_m[0];
        const double initialPeak = measureElectricPeak(solver.electricField(), false);
        if (!(initialPeak > 0) || !std::isfinite(initialPeak))
            throw std::runtime_error("Initial electric field has no finite positive peak");
        output::prepareRunDirectory(directory, comm);
        FieldSnapshots snapshots(request, grid, solver, snapshotSteps(steps));
        MeasurementCsv history(directory, comm);
        ProfileCsv<Reference> profiles(directory, reference, solver.levelContract(), comm);
        writeModes(directory, reference, comm);
        const bool windowReference = usesSafetyFactor(request);
        const auto incident        = sampleEndpoint(solver, reference, profiles, windowReference);
        advance(simulation, reference, cells, steps, initialPeak, history, snapshots);
        const auto reflected = sampleEndpoint(solver, reference, profiles, windowReference);
        writeReflection(directory, reference, incident.numerical_m, reflected.numerical_m, comm,
                        windowReference ? &incident.reference_m : nullptr,
                        windowReference ? &reflected.reference_m : nullptr);
        profiles.finish();
        history.finish();
        snapshots.finish();
        recordRun(simulation, request, grid, reference, snapshots, initialPeak, steps);
    }

    /// Load both analytic potential histories at their native positions before observing step zero.
    void run(const input::ChdrYamlRunInput& request) {
        const auto grid = makeGrid(validateInput(request));
        const Reference reference(request.caseName_m == "E19" || request.caseName_m == "E19-S");
        Simulation simulation(grid, materialMap(), solverSettings(request, grid));
        const int steps = selectedSteps(request, simulation.solver(), grid.cells_m[0]);
        simulation.load(makeInitialSampler(reference, simulation.solver()));
        observeRun(simulation, request, grid, reference, steps);
    }
}  // namespace chdr::experiments::e18e19

int main(int argc, char** argv) {
    return chdr::input::runApplication(argc, argv, chdr::experiments::e18e19::run);
}
