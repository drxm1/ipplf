/** @file E01FieldSnapshots.h
 * @brief Select E01 field frames and supply their arrays to the shared writer.
 */
#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <stdexcept>

#include "core/input/RunInput.h"
#include "core/output/FieldOutput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/SimulationConfig.h"
#include "core/solver/ChdrSolver.h"

namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric {
    using E01Solver = solver::ChdrSolver<solver::StandardSolverPeriodic>;

    /** @brief Schedule E, R2 B and raw B snapshots without changing the measured fields.
     * @details The shared writer derives positions and times from the solver's LevelContract.
     * Raw magnetic reconstruction uses independent storage through computeRawMagneticField;
     * the R2 magnetic output remains the array used by E01's error measurements.
     * @pre The run directory has been prepared collectively before construction.
     * All ranks construct, observe and finish this object in the same order.
     */
    class E01FieldSnapshots {
    public:
        /// Resolve the common YAML schedule over E01's refinement-dependent default interval.
        E01FieldSnapshots(const input::ChdrYamlRunInput& request,
                          const simulation::ChdrCartesianGridSettings& grid, E01Solver& state,
                          int defaultEvery)
            : schedule_m(input::observationSchedule(
                  request, "fields", {.every_m = defaultEvery, .includeFinal_m = true}))
            , enabled_m(input::observationEnabled(request, "fields",
                                                  !request.outputDirectory_m.empty())) {
            if (!enabled_m)
                return;
            if (request.outputDirectory_m.empty())
                throw std::invalid_argument("E01 field snapshots require an output directory");
            writer_m = std::make_unique<output::FieldOutput>(
                std::filesystem::path(request.outputDirectory_m) / "fields",
                output::OutputGeometry{grid.cells_m, {0, 0, 0}, grid.spacing_m},
                state.A_n.getLayout().comm, "fields", true,
                output::fieldOutputScales(request.units_m));
            rawMagnetic_m.initialize(state.A_n.get_mesh(), state.A_n.getLayout());
        }

        /// Reconstruct raw B only at selected states, then serialize all three native arrays.
        void observe(long long step, E01Solver& state, int finalStep) {
            if (!enabled_m || !schedule_m.selects(step, finalStep))
                return;
            state.computeRawMagneticField(rawMagnetic_m);
            const std::array fields{
                output::VectorFieldSnapshot{"E", state.electricField(), Quantity::Electric, 1.0},
                output::VectorFieldSnapshot{"B_R2", state.magneticField(), Quantity::Magnetic, 1.0},
                output::VectorFieldSnapshot{"B_raw", rawMagnetic_m, Quantity::RawMagnetic, 1.0}};
            const output::FrameRecord frame{static_cast<int>(step),
                                            state.labelTime(Quantity::VectorPotential, step),
                                            state.getDt()};
            writer_m->writeFrame(frame, state.levelContract(), fields);
        }

        /// Publish the complete accepted frame series before the run record is written.
        void finish() {
            if (writer_m)
                writer_m->writeSeries();
        }

        /// Record the effective schedule even when snapshots were explicitly disabled.
        void record(output::RunRecordDescription& description) const {
            output::recordObservation(description, "fields", schedule_m, enabled_m);
            if (!enabled_m)
                return;
            auto node                 = description.effective_m["field_snapshots"];
            node["series"]            = "fields/fields.pvti.series";
            node["arrays"]            = "E,B_R2,B_raw";
            node["sampling_metadata"] = "derived from solver LevelContract";
        }

    private:
        simulation::ObservationSchedule schedule_m;
        bool enabled_m;
        std::unique_ptr<output::FieldOutput> writer_m;
        EMFieldsType rawMagnetic_m;
    };
}  // namespace chdr::experiments::E01E02E03E04_wave_pulses_in_homogeneous_dielectric
