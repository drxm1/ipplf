/** @file E01E02Study.h
 * @brief Run a named pulse refinement ladder and collect a portable analysis bundle.
 */
#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <mpi.h>
#include <nlohmann/json.hpp>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "core/input/RunInput.h"
#include "core/output/RunRecord.h"
#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Run.h"

namespace chdr::experiments::e01e02 {
    /** @brief Plain host-side command retained until the application has finalized IPPL.
     * @details Only rank zero receives a command. The caller launches the executable
     * with the script path as a separate argument, after simulation/MPI teardown.
     */
    struct PulseAnalysisCommand {
        std::string pythonExecutable_m;
        std::filesystem::path script_m;
    };

    namespace pulseStudyDetail {
        using Json = nlohmann::json;

        /// An ordinary pulse input file and its locations in the study bundle.
        struct StudyLevel {
            std::filesystem::path sourceInput_m, copiedInput_m;
            std::string directory_m;
            int cellsPerPulse_m, transverseCells_m;
        };

        /// Resolved study choices; level settings remain in their individual YAML files.
        struct StudyPlan {
            std::filesystem::path directory_m;
            std::vector<StudyLevel> levels_m;
            bool runAnalysis_m;
            std::string pythonExecutable_m;
            PulseProtocol protocol_m;
            std::filesystem::path e01Reference_m;
        };

        /// Inputs are paths resolved from the invocation directory, like ordinary project paths.
        inline Json inputFiles(const input::ChdrYamlRunInput& request,
                               const PulseProtocol& protocol) {
            if (!request.parameters_m.has_path("input_files"))
                throw std::invalid_argument("Pulse study requires parameters.input_files");
            const auto values = Json::parse(request.parameters_m["input_files"].to_json());
            if (!values.is_array() || values.size() < protocol.minimumStudyLevels_m
                || values.size() > protocol.resolutions_m.size())
                throw std::invalid_argument(
                    "Pulse study input count does not match the experiment refinement ladder");
            for (const auto& value : values)
                if (!value.is_string() || value.get<std::string>().empty())
                    throw std::invalid_argument(
                        "Pulse study input_files entries must be nonempty paths");
            return values;
        }

        /// Reserve the bundle's own files and require an immediate child output directory.
        inline void validateLevelDirectory(const std::string& name) {
            const std::filesystem::path path(name);
            const std::set<std::string> reserved{"inputs",
                                                 "figures",
                                                 "README.md",
                                                 "study.json",
                                                 "analyze_e01_e02.py",
                                                 "analysis.json",
                                                 "analysis.log",
                                                 "convergence.csv",
                                                 "study.json.tmp",
                                                 "e01-reference",
                                                 "residual_convergence.csv",
                                                 "e01_comparison.csv"};
            if (name.empty() || path.is_absolute() || path.has_parent_path() || name == "."
                || name == ".." || reserved.contains(name))
                throw std::invalid_argument(
                    "Pulse study level output.directory must be a distinct child-directory name");
        }

        /// Use the ordinary parser and validation without modifying the loaded request.
        inline input::ChdrYamlRunInput readLevelInput(const std::filesystem::path& path,
                                                      const PulseProtocol& protocol) {
            auto request = input::readRunInput(path.string());
            validatePulseInput(request, protocol);
            if (request.steps_m || request.duration_m)
                throw std::invalid_argument(
                    "Pulse study level inputs must specify full runs without run limits");
            validateLevelDirectory(request.outputDirectory_m);
            return request;
        }

        /// Describe one input using its actual settings, with the expected ladder position.
        inline StudyLevel inspectLevel(const std::string& inputFile, std::size_t index,
                                       const PulseProtocol& protocol) {
            const auto source   = std::filesystem::absolute(inputFile).lexically_normal();
            const auto selected = readLevelInput(source, protocol);
            const int cells     = pulseCellsPerPulse(selected, protocol);
            if (cells != protocol.resolutions_m[index])
                throw std::invalid_argument(
                    "Pulse study refinement levels must follow the agreed ladder");
            return {source, std::filesystem::path("inputs") / (std::to_string(index + 1) + ".yaml"),
                    selected.outputDirectory_m, cells, pulseTransverseCells(selected, protocol)};
        }

        /// Preflight every real input file before creating any study output.
        inline std::vector<StudyLevel> studyLevels(const input::ChdrYamlRunInput& request,
                                                   const PulseProtocol& protocol) {
            const auto values = inputFiles(request, protocol);
            std::set<std::string> directories;
            std::vector<StudyLevel> result;
            for (std::size_t i = 0; i < values.size(); ++i) {
                auto level = inspectLevel(values[i].get<std::string>(), i, protocol);
                if (!result.empty() && level.transverseCells_m != result.front().transverseCells_m)
                    throw std::invalid_argument("Pulse study must keep the transverse mesh fixed");
                if (!directories.insert(level.directory_m).second)
                    throw std::invalid_argument(
                        "Pulse study level output directories must be distinct");
                result.push_back(std::move(level));
            }
            return result;
        }

        /// A Python executable is one argv entry, not a command string interpreted by a shell.
        inline std::string pythonExecutable(const input::ChdrYamlRunInput& request) {
            if (!request.parameters_m.has_path("python_executable"))
                return "python3";
            const auto node = request.parameters_m["python_executable"];
            if (!node.dtype().is_string() || node.as_string().empty()
                || node.as_string().find('\0') != std::string::npos)
                throw std::invalid_argument(
                    "Pulse study python_executable must be a nonempty executable name or path");
            return node.as_string();
        }

        /// E02 comparisons use an existing completed E01 study, never a new control run.
        inline std::filesystem::path e01ReferenceDirectory(const input::ChdrYamlRunInput& request,
                                                           const PulseProtocol& protocol) {
            if (protocol.periodicZ_m)
                return {};
            if (!request.parameters_m.has_path("e01_reference_directory")
                || !request.parameters_m["e01_reference_directory"].dtype().is_string())
                throw std::invalid_argument(
                    "E02-study requires parameters.e01_reference_directory");
            const auto path = request.parameters_m["e01_reference_directory"].as_string();
            if (path.empty() || !std::filesystem::is_directory(path))
                throw std::invalid_argument("E02-study E01 reference directory does not exist");
            return std::filesystem::absolute(path).lexically_normal();
        }

        /// Reject incomplete studies and unsupported controls before creating output directories.
        inline StudyPlan studyPlan(const input::ChdrYamlRunInput& request,
                                   const PulseProtocol& protocol) {
            if (request.caseName_m != std::string(protocol.name_m) + "-study"
                || request.units_m.system_m != "normalized")
                throw std::invalid_argument(
                    "Pulse study requires its named study case and normalized units");
            if (request.outputDirectory_m.empty() || request.steps_m || request.duration_m)
                throw std::invalid_argument(
                    "Pulse study requires output.directory and full runs without run limits");
            if (protocol.periodicZ_m)
                input::validateParameters(request,
                                          {"input_files", "run_analysis", "python_executable"});
            else
                input::validateParameters(request,
                                          {"input_files", "run_analysis", "python_executable",
                                           "e01_reference_directory"});
            input::validateObservationNames(request, {});
            if (!request.parameters_m.has_path("run_analysis"))
                throw std::invalid_argument("Pulse study requires parameters.run_analysis");
            return {std::filesystem::absolute(request.outputDirectory_m).lexically_normal(),
                    studyLevels(request, protocol),
                    input::booleanParameter(request, "run_analysis", false),
                    pythonExecutable(request),
                    protocol,
                    e01ReferenceDirectory(request, protocol)};
        }

        /// Every rank receives a root I/O failure before entering another simulation collective.
        inline void broadcastFailure(std::string message, int rank) {
            int length = static_cast<int>(message.size());
            MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (rank != 0)
                message.resize(length);
            if (length != 0) {
                MPI_Bcast(message.data(), length, MPI_CHAR, 0, MPI_COMM_WORLD);
                throw std::runtime_error(message);
            }
        }

        /// Execute filesystem work once, then collectively propagate its success or failure.
        inline void onRoot(int rank, const std::function<void()>& action) {
            std::string failure;
            if (rank == 0) {
                try {
                    action();
                } catch (const std::exception& error) {
                    failure = error.what();
                }
            }
            broadcastFailure(failure, rank);
        }

        /// Close a complete temporary file before replacing the study's progress record.
        inline void writeText(const std::filesystem::path& path, const std::string& contents) {
            const auto temporary = path.string() + ".tmp";
            std::ofstream output;
            output.exceptions(std::ios::failbit | std::ios::badbit);
            output.open(temporary, std::ios::out | std::ios::trunc);
            output << contents << '\n';
            output.close();
            std::filesystem::rename(temporary, path);
        }

        /// Relative child paths keep the manifest usable after moving the completed study.
        inline Json studyManifest(const input::ChdrYamlRunInput& request, const StudyPlan& plan) {
            Json manifest{{"schema_version", 1},
                          {"case", request.caseName_m},
                          {"status", "running"},
                          {"run_analysis", plan.runAnalysis_m},
                          {"python_executable", plan.pythonExecutable_m},
                          {"analysis_script", "analyze_e01_e02.py"},
                          {"levels", Json::array()},
                          {"input", Json::parse(input::runInputRecord(request).to_json())}};
            if (!plan.e01Reference_m.empty())
                manifest["e01_reference"] = {{"source_directory", plan.e01Reference_m.string()},
                                             {"bundled_directory", "e01-reference"}};
            for (const auto& level : plan.levels_m)
                manifest["levels"].push_back(
                    {{"cells_per_pulse", level.cellsPerPulse_m},
                     {"transverse_cells_per_axis", level.transverseCells_m},
                     {"directory", level.directory_m},
                     {"input_file", level.copiedInput_m.generic_string()},
                     {"status", "pending"}});
            return manifest;
        }

        /// The copied analysis locates its data relative to its own file, independent of cwd.
        inline constexpr const char* StudyReadme = R"(# Homogeneous pulse refinement study

`study.json` records the requested levels and their completion status.
Each listed level directory contains its own `run.json`, `errors.csv`,
`profiles.csv` and any enabled field snapshots. The `inputs/` directory
contains exact copies of the individual YAML inputs used by those runs.
Each input was read with the ordinary YAML parser and executed unchanged,
with this study directory as its working directory.

To repeat the analysis from this directory, run:

```sh
python3 analyze_e01_e02.py
```

The script requires Python 3 and matplotlib. It reads only completed study
levels and writes `convergence.csv`, `analysis.json` and figures under `figures/`.
The script also accepts `--directory PATH` to analyze another study directory.
Use the Python executable recorded in `study.json` if it differs from python3.
The data and copied script can be moved together; their internal paths are relative.
)";

        /// Copy archived data unchanged; the analysis checks their mesh and native sample times.
        inline void copyReferenceLevel(const StudyPlan& plan, const Json& study, int resolution) {
            for (const auto& entry : study.at("levels")) {
                if (entry.at("cells_per_pulse") != resolution || entry.at("status") != "completed")
                    continue;
                const auto source = plan.e01Reference_m / entry.at("directory").get<std::string>();
                const auto destination =
                    plan.directory_m / "e01-reference" / ("n" + std::to_string(resolution));
                std::filesystem::create_directories(destination);
                for (const char* name : {"run.json", "profiles.csv"})
                    std::filesystem::copy_file(source / name, destination / name);
                return;
            }
            throw std::runtime_error("E01 reference is missing a completed E02 resolution");
        }

        /// Bundle only the existing E01 records and profiles needed for the matching resolutions.
        inline void copyE01Reference(const StudyPlan& plan) {
            if (plan.e01Reference_m.empty())
                return;
            std::ifstream file(plan.e01Reference_m / "study.json");
            const auto study = Json::parse(file);
            if (study.at("case") != "E01-study" || study.at("status") != "completed")
                throw std::runtime_error("E02 reference must be a completed E01 study");
            for (const auto& level : plan.levels_m)
                copyReferenceLevel(plan, study, level.cellsPerPulse_m);
        }

        /// Copy the reviewed analysis source before the first run; never generate Python code.
        inline void prepareBundle(const StudyPlan& plan, const Json& manifest, int rank) {
            output::prepareRunDirectory(plan.directory_m.string(), MPI_COMM_WORLD);
            onRoot(rank, [&] {
                std::filesystem::copy_file(CHDR_E01_E02_ANALYSIS_SOURCE,
                                           plan.directory_m / "analyze_e01_e02.py");
                std::filesystem::create_directory(plan.directory_m / "inputs");
                for (const auto& level : plan.levels_m)
                    std::filesystem::copy_file(level.sourceInput_m,
                                               plan.directory_m / level.copiedInput_m);
                copyE01Reference(plan);
                writeText(plan.directory_m / "README.md", StudyReadme);
                writeText(plan.directory_m / "study.json", manifest.dump(2));
            });
        }

        /// Reread the copied YAML on every rank; prevent a changed input from mislabelling output.
        inline input::ChdrYamlRunInput copiedLevelInput(const StudyPlan& plan,
                                                        const StudyLevel& level) {
            auto selected = readLevelInput(plan.directory_m / level.copiedInput_m, plan.protocol_m);
            if (pulseCellsPerPulse(selected, plan.protocol_m) != level.cellsPerPulse_m
                || pulseTransverseCells(selected, plan.protocol_m) != level.transverseCells_m
                || selected.outputDirectory_m != level.directory_m)
                throw std::runtime_error("Copied pulse input disagrees with its study manifest");
            return selected;
        }

        /// Record each transition only after the preceding collective operation has finished.
        inline void recordLevel(const StudyPlan& plan, Json& manifest, std::size_t index,
                                const char* status, int rank) {
            onRoot(rank, [&] {
                manifest["levels"][index]["status"] = status;
                writeText(plan.directory_m / "study.json", manifest.dump(2));
                std::cout << plan.protocol_m.name_m << " n" << plan.levels_m[index].cellsPerPulse_m
                          << ": " << status << std::endl;
            });
        }

        /// Successful level return means all simulation objects and output streams are closed.
        inline void runLevels(const StudyPlan& plan, Json& manifest, int rank,
                              void (*runLevel)(const input::ChdrYamlRunInput&)) {
            for (std::size_t i = 0; i < plan.levels_m.size(); ++i) {
                const auto selected = copiedLevelInput(plan, plan.levels_m[i]);
                recordLevel(plan, manifest, i, "running", rank);
                runLevel(selected);
                recordLevel(plan, manifest, i, "completed", rank);
            }
            onRoot(rank, [&] {
                manifest["status"] = "completed";
                writeText(plan.directory_m / "study.json", manifest.dump(2));
            });
        }

        /** @brief Retain ordinary relative output paths by running from the study directory.
         * @details Restore explicitly on success. During exception unwinding, attempt
         * restoration without replacing the original failure with another exception.
         */
        class WorkingDirectory {
        public:
            explicit WorkingDirectory(const std::filesystem::path& directory)
                : original_m(std::filesystem::current_path()) {
                std::filesystem::current_path(directory);
            }
            WorkingDirectory(const WorkingDirectory&)            = delete;
            WorkingDirectory& operator=(const WorkingDirectory&) = delete;
            ~WorkingDirectory() {
                if (active_m) {
                    std::error_code ignored;
                    std::filesystem::current_path(original_m, ignored);
                }
            }
            void restore() {
                std::filesystem::current_path(original_m);
                active_m = false;
            }

        private:
            std::filesystem::path original_m;
            bool active_m = true;
        };
    }  // namespace pulseStudyDetail

    /** @brief Run the selected experiment levels on every rank and return optional post-finalize
     * analysis.
     * @param runLevel Named single-level execution function; called identically on all
     * ranks.
     * @details IPPL is initialized once by runApplication. Each level owns and releases its
     * simulation before the next level begins. Every level is read from its own copied
     * YAML through readRunInput and passed unchanged to runLevel. Only the parent study
     * directory is prepared here; each child remains new/empty for ordinary pulse output.
     * No Python process is launched while this function or its MPI collectives are active.
     */
    inline std::optional<PulseAnalysisCommand> runPulseStudy(
        const input::ChdrYamlRunInput& request, const PulseProtocol& protocol,
        void (*runLevel)(const input::ChdrYamlRunInput&)) {
        if (!runLevel)
            throw std::invalid_argument("Pulse study requires a single-level execution function");
        const auto plan = pulseStudyDetail::studyPlan(request, protocol);
        int rank        = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        auto manifest = pulseStudyDetail::studyManifest(request, plan);
        pulseStudyDetail::prepareBundle(plan, manifest, rank);
        pulseStudyDetail::WorkingDirectory directory(plan.directory_m);
        pulseStudyDetail::runLevels(plan, manifest, rank, runLevel);
        directory.restore();
        if (rank != 0 || !plan.runAnalysis_m)
            return std::nullopt;
        return PulseAnalysisCommand{plan.pythonExecutable_m,
                                    plan.directory_m / "analyze_e01_e02.py"};
    }
}  // namespace chdr::experiments::e01e02
