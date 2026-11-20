/** @file E03E04FieldSnapshots.h
 * @brief Optional full-volume frames and z=0 sections of the spherical pulse runs.
 */
#pragma once

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "core/input/RunInput.h"
#include "core/output/FieldOutput.h"
#include "core/output/RunRecord.h"
#include "core/simulation/SimulationConfig.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03E04Setup.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03Plane.h"

namespace chdr::experiments::e03 {
    /** @brief Write E and B of selected states of a spherical pulse run as a ParaView series.
     * @details The frames are a visualization product. They do not enter any measurement,
     * and they are disabled unless the input enables the `fields` observation, because one
     * frame of the finest mesh holds \f$2\cdot3\cdot324^3\f$ values. The arrays are the
     * reconstructed fields that the measurements use, at their native positions;
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

    /** @brief Write the z=0 sections of E and B at prescribed and scheduled states.
     * @details E03's `write_plane` prescribes the initial state and its final field target.
     * A `plane` block in `observations` adds a schedule for two-dimensional videos; its
     * `every` defaults to one section per step, and the final state is always included.
     * writeE03Plane appends each section to `plane-z0.csv`: \f$N^2\f$ rows instead of the
     * \f$N^3\f$ cells of a full frame. Like the full frames, the sections enter no measurement.
     * @pre The run directory has been prepared collectively. All ranks observe every step.
     */
    class SphericalPlaneFrames {
    public:
        SphericalPlaneFrames(const input::ChdrYamlRunInput& request,
                             std::vector<int> prescribedSteps = {})
            : directory_m(request.outputDirectory_m)
            , prescribedSteps_m(std::move(prescribedSteps))
            , schedule_m(input::observationSchedule(request, "plane",
                                                    {.every_m = 1, .includeFinal_m = true}))
            , enabled_m(input::observationEnabled(request, "plane",
                                                  request.observations_m.count("plane") > 0)) {}

        /// True if the run writes `plane-z0.csv`.
        bool writes() const { return enabled_m || !prescribedSteps_m.empty(); }

        /// Write the section of this state if it is prescribed or scheduled.
        template <class Solver, DipolePulseReference Reference>
        void observe(long long step, Solver& state, const Reference& ref, int finalStep) const {
            const bool prescribed =
                std::find(prescribedSteps_m.begin(), prescribedSteps_m.end(), step)
                != prescribedSteps_m.end();
            if (prescribed || (enabled_m && schedule_m.selects(step, finalStep)))
                writeE03Plane(state, ref, directory_m, step);
        }

    private:
        std::string directory_m;
        std::vector<int> prescribedSteps_m;
        simulation::ObservationSchedule schedule_m;
        bool enabled_m;
    };
}  // namespace chdr::experiments::e03
