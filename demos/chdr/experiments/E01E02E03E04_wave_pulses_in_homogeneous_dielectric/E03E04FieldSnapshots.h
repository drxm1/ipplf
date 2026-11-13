/** @file E03E04FieldSnapshots.h
 * @brief Optional full-volume field frames of the spherical pulse runs, for visualization.
 */
#pragma once

#include <array>
#include <filesystem>
#include <memory>

#include "core/input/RunInput.h"
#include "core/output/FieldOutput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/SimulationConfig.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03E04Setup.h"

namespace chdr::experiments::e03 {
    /** @brief Write E and R2 B of selected states of a spherical pulse run as a ParaView series.
     * @details The frames are a visualization product. They do not enter any measurement,
     * and they are disabled unless the input enables the `fields` observation, because one
     * frame of the finest mesh holds \f$2\cdot3\cdot324^3\f$ values. The arrays are the
     * reconstructed fields that the measurements use, at their native staggered positions;
     * the shared writer records each array's time and component positions from the
     * solver's LevelContract. With `single_precision_fields: true` the values are stored as
     * 4-byte numbers, which halves the size of a frame.
     * @pre The run directory has been prepared collectively before construction.
     * All ranks construct, observe and finish this object in the same order.
     */
    class SphericalFieldSnapshots {
        using Solver = SphericalSimulation::Solver;

    public:
        /** @brief Resolve the `fields` schedule of the input.
         * @details A `fields` block in `observations` enables the frames unless it says
         * `enabled: false`. Its `every` defaults to one frame per step.
         */
        SphericalFieldSnapshots(const input::ChdrYamlRunInput& request,
                                const simulation::ChdrCartesianGridSettings& grid, Solver& state)
            : schedule_m(input::observationSchedule(request, "fields",
                                                    {.every_m = 1, .includeFinal_m = true}))
            , enabled_m(input::observationEnabled(request, "fields",
                                                  request.observations_m.count("fields") > 0))
            , singlePrecision_m(
                  input::booleanParameter(request, "single_precision_fields", false)) {
            if (!enabled_m)
                return;
            writer_m = std::make_unique<output::FieldOutput>(
                std::filesystem::path(request.outputDirectory_m) / "fields",
                output::OutputGeometry{grid.cells_m, {0, 0, 0}, grid.spacing_m},
                state.A_n.getLayout().comm, "fields", true,
                output::fieldOutputScales(request.units_m), singlePrecision_m);
        }

        /// Serialize the two reconstructed fields of a selected completed state.
        void observe(long long step, Solver& state, int finalStep) {
            if (!enabled_m || !schedule_m.selects(step, finalStep))
                return;
            const std::array fields{
                output::VectorFieldSnapshot{"E", state.electricField(), Quantity::Electric, 1.0},
                output::VectorFieldSnapshot{"B_R2", state.magneticField(), Quantity::Magnetic,
                                            1.0}};
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

        /// Record the effective schedule and storage, also when the frames are disabled.
        void record(output::RunRecordDescription& description) const {
            output::recordObservation(description, "fields", schedule_m, enabled_m);
            if (!enabled_m)
                return;
            auto node                 = description.effective_m["field_snapshots"];
            node["series"]            = "fields/fields.pvti.series";
            node["arrays"]            = "E,B_R2";
            node["value_type"]        = singlePrecision_m ? "Float32" : "Float64";
            node["sampling_metadata"] = "derived from solver LevelContract";
        }

    private:
        simulation::ObservationSchedule schedule_m;
        bool enabled_m, singlePrecision_m;
        std::unique_ptr<output::FieldOutput> writer_m;
    };
}  // namespace chdr::experiments::e03
