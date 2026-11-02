/** @file RunRecord.h
 * @brief Record common configuration and effective execution settings after a successful run.
 */
#pragma once

#include <array>
#include <mpi.h>
#include <string>
#include <vector>

#include "core/input/RunInput.h"
#include "core/simulation/SimulationConfig.h"

namespace chdr::output {
    /** @brief Effective settings supplied by an application and its completed solver.
     * @details The structured common fields supplement the case's existing metadata.
     * Local bounds are gathered when every rank participates in writeRunRecord().
     */
    struct RunRecordDescription {
        std::string application_m, caseMetadata_m;
        conduit_cpp::Node effective_m;
        std::array<int, 6> localBounds_m{};
    };

    namespace run_record_detail {
        /// Store common box settings without assigning physical units to normalized values.
        void geometryRecord(conduit_cpp::Node node, const simulation::RunBox& box);
        /// Give the actual temporal reconstruction a stable, readable record label.
        const char* reconstructionName(EBReconstructionConvention convention);

        /// Copy solver-selected settings rather than constructor defaults into the record.
        template <class Solver>
        void solverRecord(conduit_cpp::Node node, const Solver& state) {
            node["time_step"]     = state.getDt();
            node["safety_factor"] = state.safetyFactor();
            node["stencil"] =
                Solver::stencil() == solver::Stencil::Standard ? "Standard" : "NonStandard";
            node["formulation"]    = state.formulationName();
            node["reconstruction"] = reconstructionName(state.ebReconstructionConvention());
            node["spatial"] = state.spatialDiscretization() == SpatialDiscretization::Staggered
                                  ? "staggered"
                                  : "collocated";
            if (state.hasMaterial())
                node["formulation_alpha"] = state.formulationParameters().alpha_m;
        }
    }  // namespace run_record_detail

    /** @brief Describe the actual completed run before its fields and solver are destroyed.
     * @details All numerical values in effective_m use normalized solver units. The
     * final record additionally names the selected output units and conversion factors.
     * This function performs no MPI operation and is called once on each rank.
     */
    template <class Solver>
    RunRecordDescription describeRun(const std::string& application, const simulation::RunBox& box,
                                     const Solver& solver, int steps,
                                     const std::string& caseMetadata = "") {
        RunRecordDescription result;
        result.application_m  = application;
        result.caseMetadata_m = caseMetadata;
        run_record_detail::geometryRecord(result.effective_m["box"], box);
        run_record_detail::solverRecord(result.effective_m["solver"], solver);
        result.effective_m["requested_steps"] = steps;
        result.effective_m["completed_solves"] =
            static_cast<conduit_int64>(solver.completedSolves());
        const auto local = solver.A_n.getLayout().getLocalNDIndex();
        for (unsigned d = 0; d < 3; ++d) {
            result.localBounds_m[2 * d]     = local[d].first();
            result.localBounds_m[2 * d + 1] = local[d].last();
        }
        return result;
    }

    /// Record an application's resolved schedule, including defaults and disabled observations.
    void recordObservation(RunRecordDescription& description, const std::string& name,
                           const simulation::ObservationSchedule& schedule, bool enabled = true);

    /** @brief Collectively prepare a new or empty output directory before a run starts.
     * @details Use this when no field writer performs the same preparation. Empty paths
     * disable output. Rank zero creates missing directories and rejects existing files or
     * nonempty directories; its errors are broadcast before any rank continues.
     */
    void prepareRunDirectory(const std::string& directory, MPI_Comm communicator);

    /** @brief Collectively write run.json after successful simulation and output finalization.
     * @details Empty directory disables this record for legacy stdout-only runs. A missing
     * directory is created; an existing run.json is never replaced. Construct field writers
     * first, because they require an empty output directory. All ranks call with identical
     * settings. Rank-zero filesystem failures are broadcast before any caller returns.
     * Source/build definitions are injected by CMake; absent definitions read "unrecorded".
     */
    void writeRunRecord(const std::string& directory, const input::RunInput& input,
                        const RunRecordDescription& description, MPI_Comm communicator);
}  // namespace chdr::output
