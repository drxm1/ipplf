/** @file run_experiment.cpp
 * @brief E16 and E17 normal incidence at one planar dielectric face.
 */
#include <algorithm>
#include <array>
#include <iostream>
#include <string>

#include "core/input/Application.h"
#include "core/output/RunRecord.h"
#include "experiments/E16E17_normal_incidence_on_planar_dielectric_interface/E16E17Measurements.h"
#include "experiments/E16E17_normal_incidence_on_planar_dielectric_interface/E16E17Output.h"
#include "experiments/E16E17_normal_incidence_on_planar_dielectric_interface/E16E17Setup.h"

namespace chdr::experiments::e16e17 {
    /// Save the two energy sums defined by E16.meas03 from the reconstructed fields.
    void observeEnergy(Solver& solver, const material::MaterialOwner& material,
                       const Reference& reference, MeasurementCsv& csv, const char* checkpoint,
                       double target, long long step) {
        const auto energy = measureSideEnergy(solver.electricField(), solver.magneticField(),
                                              material.epsilon(), reference.zf_m);
        const double time = solver.labelTime(Quantity::Electric, step);
        csv.row("U_side1", checkpoint, target, step, time, energy[0]);
        csv.row("U_side2", checkpoint, target, step, time, energy[1]);
    }

    /// Preserve the unnormalized squared sums of E16.meas01 and E16.meas02.
    void observeErrors(Solver& solver, const Reference& reference, MeasurementCsv& csv,
                       const char* checkpoint, double target, long long step) {
        const std::array quantities{Quantity::Electric, Quantity::Magnetic};
        const std::array prefixes{"E", "B"};
        const std::array suffixes{"_error_side1", "_reference_side1", "_error_side2",
                                  "_reference_side2"};
        for (unsigned q = 0; q < quantities.size(); ++q) {
            auto& field       = q == 0 ? solver.electricField() : solver.magneticField();
            const double time = solver.labelTime(quantities[q], step);
            const auto sums   = measureSideFieldSums(field, reference, solver.levelContract(),
                                                     quantities[q], time, reference.zf_m);
            for (unsigned index = 0; index < sums.size(); ++index)
                csv.row((std::string(prefixes[q]) + suffixes[index]).c_str(), checkpoint, target,
                        step, time, sums[index]);
        }
    }

    /// Save the two native reconstructed profiles at their common R2 time.
    void observeFieldProfiles(Solver& solver, ProfileCsv<Reference>& profiles,
                              const char* checkpoint, double target, long long step) {
        profiles.write(solver.electricField(), Quantity::Electric, checkpoint, target, step,
                       solver.labelTime(Quantity::Electric, step));
        profiles.write(solver.magneticField(), Quantity::Magnetic, checkpoint, target, step,
                       solver.labelTime(Quantity::Magnetic, step));
    }

    /// Select potential and field samples independently at each prescribed target.
    void observeCheckpoint(Solver& solver, const material::MaterialOwner& material,
                           const Reference& reference, const Checkpoint& point, long long step,
                           MeasurementCsv& csv, ProfileCsv<Reference>& profiles) {
        if (step == point.potentialStep_m)
            profiles.write(solver.A_n, Quantity::VectorPotential, point.name_m, point.target_m,
                           step, solver.labelTime(Quantity::VectorPotential, step));
        if (step != point.fieldStep_m)
            return;
        observeErrors(solver, reference, csv, point.name_m, point.target_m, step);
        observeEnergy(solver, material, reference, csv, point.name_m, point.target_m, step);
        observeFieldProfiles(solver, profiles, point.name_m, point.target_m, step);
    }

    /** @brief Observe a state once per quantity and target, with energy at the fixed cadence.
     * @details The initial sample is explicitly step zero. At a target that also lies on the
     * energy cadence, the checkpoint rows supply the energy history without duplicate rows.
     */
    void observe(Solver& solver, const material::MaterialOwner& material,
                 const Reference& reference, const std::array<Checkpoint, 3>& points,
                 long long step, int stride, MeasurementCsv& csv, ProfileCsv<Reference>& profiles) {
        bool fieldObserved = false;
        for (const auto& point : points) {
            observeCheckpoint(solver, material, reference, point, step, csv, profiles);
            fieldObserved = fieldObserved || step == point.fieldStep_m;
        }
        if (!fieldObserved && step % stride == 0)
            observeEnergy(solver, material, reference, csv, "history",
                          solver.labelTime(Quantity::Electric, step), step);
    }

    /// Record the targets and actual sample times so analysis uses the solver's labels.
    void recordTargets(output::RunRecordDescription& description, const Solver& solver,
                       const std::array<Checkpoint, 3>& points, int stride) {
        auto schedule            = description.effective_m["observations"];
        schedule["energy_every"] = stride;
        schedule["protocol_steps"] =
            std::max(points.back().potentialStep_m, points.back().fieldStep_m);
        for (const auto& point : points) {
            auto target              = schedule[point.name_m];
            target["target_time"]    = point.target_m;
            target["potential_step"] = point.potentialStep_m;
            target["field_step"]     = point.fieldStep_m;
            target["potential_time"] =
                solver.labelTime(Quantity::VectorPotential, point.potentialStep_m);
            target["field_time"] = solver.labelTime(Quantity::Electric, point.fieldStep_m);
        }
    }

    /// Publish configuration after all CSV and field files have closed successfully.
    void recordRun(Simulation& simulation, const input::ChdrYamlRunInput& request,
                   const simulation::ChdrCartesianGridSettings& grid, const Reference& reference,
                   const std::array<Checkpoint, 3>& points, int steps, int stride,
                   const FieldSnapshots& snapshots) {
        auto description = output::describeRun("chdr-experiment-e16-e17", grid, simulation.solver(),
                                               steps, request.caseName_m + " normal incidence");
        reference.record(description.effective_m["analytic_reference"], request.caseName_m);
        description.effective_m["material/epsilon_side1"]   = reference.epsilon1_m;
        description.effective_m["material/epsilon_side2"]   = reference.epsilon2_m;
        description.effective_m["material/interface_z"]     = reference.zf_m;
        description.effective_m["measurements/file"]        = "measurements.csv";
        description.effective_m["measurements/identifiers"] = "E16.meas01,E16.meas02,E16.meas03";
        description.effective_m["profiles/file"]            = "profiles.csv";
        recordTargets(description, simulation.solver(), points, stride);
        snapshots.record(description);
        output::writeRunRecord(request.outputDirectory_m, request, description,
                               simulation.solver().A_n.getLayout().comm);
    }

    /// Execute the prescribed observations and write a completion record for the selected run.
    void evolve(Simulation& simulation, const input::ChdrYamlRunInput& request,
                const simulation::ChdrCartesianGridSettings& grid, const Reference& reference,
                int pulseCells) {
        auto& solver      = simulation.solver();
        const auto points = checkpoints(solver, reference);
        const int steps = selectedSteps(request, points), stride = pulseCells / 16;
        const auto comm = solver.A_n.getLayout().comm;
        output::prepareRunDirectory(request.outputDirectory_m, comm);
        MeasurementCsv csv(request.outputDirectory_m, comm);
        ProfileCsv<Reference> profiles(request.outputDirectory_m, reference, solver.levelContract(),
                                       comm);
        FieldSnapshots snapshots(request, grid, solver, snapshotSteps(points));
        simulation.advance(steps, [&](long long step, Solver& state, const auto& material) {
            observe(state, material, reference, points, step, stride, csv, profiles);
            snapshots.observe(step, state);
        });
        csv.finish();
        profiles.finish();
        snapshots.finish();
        recordRun(simulation, request, grid, reference, points, steps, stride, snapshots);
    }

    /// Load continuum potential histories using the solver's spatial and temporal contract.
    void run(const input::ChdrYamlRunInput& request) {
        validateInput(request);
        const auto reference = makeReference(request);
        const int pulseCells = cellsPerPulse(request);
        const auto grid      = makeGrid(pulseCells);
        Simulation simulation(grid, materialMap(reference), {.timeStep_m = grid.spacing_m[2] / 2});
        simulation.load(reference);
        evolve(simulation, request, grid, reference, pulseCells);
    }
}  // namespace chdr::experiments::e16e17

int main(int argc, char** argv) {
    return chdr::input::runApplication(argc, argv, chdr::experiments::e16e17::run);
}
