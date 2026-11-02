/** @file RunRecord.cpp
 * @brief Serialize the shared run record without coupling it to any particular diagnostic.
 */
#include <Kokkos_Core.hpp>

#include "core/output/RunRecord.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>

#if __has_include("ChdrSourceIdentity.h")
#include "ChdrSourceIdentity.h"
#endif

#ifndef CHDR_SOURCE_REVISION
#define CHDR_SOURCE_REVISION "unrecorded"
#endif
#ifndef CHDR_SOURCE_DIRTY
#define CHDR_SOURCE_DIRTY "unrecorded"
#endif
#ifndef CHDR_BUILD_TYPE
#define CHDR_BUILD_TYPE "unrecorded"
#endif
#ifndef CHDR_COMPILER
#define CHDR_COMPILER "unrecorded"
#endif

namespace chdr::output {
    namespace run_record_detail {
        using Node = conduit_cpp::Node;

        void geometryRecord(Node node, const simulation::RunBox& box) {
            node["cells"].set(std::vector<int>(box.cells_m.begin(), box.cells_m.end()));
            node["spacing"].set(std::vector<double>(box.spacing_m.begin(), box.spacing_m.end()));
            node["periodic"].set(std::vector<int>(box.periodic_m.begin(), box.periodic_m.end()));
            node["split_axis"]     = box.splitAxis_m;
            node["partition_seam"] = box.seam_m;
        }

        const char* reconstructionName(EBReconstructionConvention convention) {
            switch (convention) {
                case EBReconstructionConvention::R1:
                    return "R1";
                case EBReconstructionConvention::R2:
                    return "R2";
                case EBReconstructionConvention::R2Prime:
                    return "R2Prime";
            }
            throw std::invalid_argument("Unknown field reconstruction in run record");
        }

        /// Store the source identity captured when this binary was built.
        void buildRecord(Node node) {
            node["source_revision"]  = CHDR_SOURCE_REVISION;
            node["source_dirty"]     = CHDR_SOURCE_DIRTY;
            node["type"]             = CHDR_BUILD_TYPE;
            node["compiler"]         = CHDR_COMPILER;
            node["execution_space"]  = Kokkos::DefaultExecutionSpace::name();
            node["host_concurrency"] = Kokkos::DefaultHostExecutionSpace().concurrency();
        }

        /// Record standard unit factors; each output description identifies the factors it uses.
        void unitsRecord(Node node, const input::UnitSystem& units) {
            node["system"]                = units.system_m;
            node["length"]                = units.lengthFactor();
            node["time"]                  = units.timeFactor();
            node["electric_field"]        = units.electricFactor();
            node["magnetic_field"]        = units.magneticFactor();
            node["scalar_potential"]      = units.scalarPotentialFactor();
            node["vector_potential"]      = units.vectorPotentialFactor();
            node["physical_field_energy"] = units.energyFactor();
            node["velocity"]              = units.velocityFactor();
            node["frequency"]             = units.frequencyFactor();
        }

        /// Gather the actual owned index ranges, including automatic MPI decomposition.
        std::vector<int> rankBounds(const RunRecordDescription& description, MPI_Comm communicator,
                                    int rank, int ranks) {
            std::vector<int> result(rank == 0 ? 6 * static_cast<std::size_t>(ranks) : 0);
            MPI_Gather(description.localBounds_m.data(), 6, MPI_INT, result.data(), 6, MPI_INT, 0,
                       communicator);
            return result;
        }

        /// Record only selected execution variables; unrelated process environment is excluded.
        void executionRecord(Node node, int ranks, const std::vector<int>& bounds) {
            node["ranks"] = ranks;
            node["rank_bounds_order"] =
                "x_first,x_last,y_first,y_last,z_first,z_last; blocks ordered by MPI rank";
            node["rank_bounds"].set(bounds);
            char processor[MPI_MAX_PROCESSOR_NAME];
            int length = 0;
            MPI_Get_processor_name(processor, &length);
            node["rank_zero_host"] = std::string(processor, length);
            for (const char* name : {"OMP_NUM_THREADS", "OMP_PROC_BIND"})
                if (const char* value = std::getenv(name))
                    node[std::string("environment/") + name] = std::string(value);
        }

        /// Assemble a completed record; effective values remain explicitly in solver units.
        Node makeRecord(const input::RunInput& input, const RunRecordDescription& description,
                        int ranks, const std::vector<int>& bounds) {
            Node result;
            result["record_version"] = 1;
            result["status"]         = "completed";
            result["application"]    = description.application_m;
            result["input"].set(input::runInputRecord(input));
            result["effective"].set(description.effective_m);
            result["effective/value_system"]  = "normalized";
            result["effective/case_metadata"] = description.caseMetadata_m;
            buildRecord(result["build"]);
            unitsRecord(result["output_factors"], input.units_m);
            executionRecord(result["execution"], ranks, bounds);
            return result;
        }

        /// Read one finite floating value, respecting Conduit's offsets and strides.
        template <class Value>
        double floatingValue(const Node& node, conduit_index_t index) {
            Value value;
            std::memcpy(&value, node.element_ptr(index), sizeof(Value));
            if (!std::isfinite(value))
                throw std::runtime_error("Run record contains a non-finite value");
            return value;
        }

        /// Preserve floating leaves before JSON serialization; Conduit's formatter rounds them.
        nlohmann::json floatingJson(const Node& node) {
            const auto type = node.dtype();
            if (!type.endianness_matches_machine() || (!type.is_float32() && !type.is_float64()))
                throw std::runtime_error("Unsupported floating representation in run record");
            auto result = nlohmann::json::array();
            for (conduit_index_t i = 0; i < type.number_of_elements(); ++i)
                result.push_back(type.is_float64() ? floatingValue<conduit_float64>(node, i)
                                                   : floatingValue<conduit_float32>(node, i));
            return result.size() == 1 ? result.front() : result;
        }

        /// Convert the record's objects, lists and scalar/array leaves without rounding doubles.
        nlohmann::json recordJson(const Node& node) {
            const auto type = node.dtype();
            if (type.is_object() || type.is_list()) {
                auto result = type.is_object() ? nlohmann::json::object() : nlohmann::json::array();
                for (conduit_index_t i = 0; i < node.number_of_children(); ++i) {
                    const auto child = node.child(i);
                    if (type.is_object())
                        result[child.name()] = recordJson(child);
                    else
                        result.push_back(recordJson(child));
                }
                return result;
            }
            if (type.is_floating_point())
                return floatingJson(node);
            if (type.is_integer() || type.is_string() || type.is_empty())
                return nlohmann::json::parse(node.to_json());
            throw std::runtime_error("Unsupported value type in run record: " + type.name());
        }

        /// Write only a new record, using JSON's shortest decimal that round-trips each double.
        void saveRecord(const std::string& directory, const Node& record) {
            const auto contents = recordJson(record).dump(2);
            std::filesystem::create_directories(directory);
            const auto path = std::filesystem::path(directory) / "run.json";
            if (std::filesystem::exists(path))
                throw std::runtime_error("Refusing to replace an existing run record: "
                                         + path.string());
            std::ofstream output(path);
            output << contents << '\n';
            output.close();
            if (!output)
                throw std::runtime_error("Cannot write run record: " + path.string());
        }

        /// Propagate a rank-zero filesystem failure before returning to any caller.
        void broadcastFailure(std::string message, int rank, MPI_Comm communicator) {
            int length = static_cast<int>(message.size());
            MPI_Bcast(&length, 1, MPI_INT, 0, communicator);
            if (rank != 0)
                message.resize(length);
            if (length != 0) {
                MPI_Bcast(message.data(), length, MPI_CHAR, 0, communicator);
                throw std::runtime_error(message);
            }
        }
    }  // namespace run_record_detail

    void prepareRunDirectory(const std::string& directory, MPI_Comm communicator) {
        if (directory.empty())
            return;
        int rank = 0;
        MPI_Comm_rank(communicator, &rank);
        std::string failure;
        if (rank == 0) {
            try {
                if (std::filesystem::exists(directory) && !std::filesystem::is_empty(directory))
                    throw std::runtime_error("Output directory must be new or empty: " + directory);
                std::filesystem::create_directories(directory);
            } catch (const std::exception& error) {
                failure = error.what();
            }
        }
        run_record_detail::broadcastFailure(failure, rank, communicator);
    }

    void recordObservation(RunRecordDescription& description, const std::string& name,
                           const simulation::ObservationSchedule& schedule, bool enabled) {
        schedule.validate();
        auto node             = description.effective_m["observations/" + name];
        node["enabled"]       = enabled ? 1 : 0;
        node["first"]         = schedule.first_m;
        node["every"]         = schedule.every_m;
        node["last"]          = schedule.last_m;
        node["include_final"] = schedule.includeFinal_m ? 1 : 0;
    }

    void writeRunRecord(const std::string& directory, const input::RunInput& input,
                        const RunRecordDescription& description, MPI_Comm communicator) {
        if (directory.empty())
            return;
        int rank = 0, ranks = 0;
        MPI_Comm_rank(communicator, &rank);
        MPI_Comm_size(communicator, &ranks);
        const auto bounds = run_record_detail::rankBounds(description, communicator, rank, ranks);
        std::string failure;
        if (rank == 0) {
            try {
                run_record_detail::saveRecord(
                    directory, run_record_detail::makeRecord(input, description, ranks, bounds));
            } catch (const std::exception& error) {
                failure = error.what();
            }
        }
        run_record_detail::broadcastFailure(failure, rank, communicator);
    }
}  // namespace chdr::output
