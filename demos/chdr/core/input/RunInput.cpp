/** @file RunInput.cpp
 * @brief Decode common YAML options without changing experiment-specific semantics.
 */
#include "core/input/RunInput.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace chdr::input {
    namespace {
        using Node = conduit_cpp::Node;

        /// Include the configuration path in every validation error.
        [[noreturn]] void fail(const std::string& path, const std::string& reason) {
            throw std::invalid_argument("ChDR input '" + path + "': " + reason);
        }

        /// Require a mapping before using named children; null is not an empty mapping.
        void mapping(const Node& node, const std::string& path) {
            if (!node.dtype().is_object())
                fail(path, "expected a mapping");
        }

        /// Validate direct child names, without allowing Conduit path lookup to create keys.
        void keys(const Node& node, const std::string& path,
                  std::initializer_list<std::string> allowed) {
            mapping(node, path);
            for (conduit_index_t i = 0; i < node.number_of_children(); ++i) {
                const auto name = node.child(i).name();
                if (std::find(allowed.begin(), allowed.end(), name) == allowed.end())
                    fail(path + "/" + name, "unsupported key");
            }
        }

        /// Fetch a required value without silently constructing a missing node.
        Node required(const Node& node, const std::string& path) {
            if (!node.has_path(path))
                fail(path, "missing required value");
            return node[path];
        }

        /// Keep arrays and strings out of numeric scalar conversions.
        void scalar(const Node& node, const std::string& path) {
            if (!node.dtype().is_number() || node.number_of_elements() != 1)
                fail(path, "expected one numeric scalar");
        }

        /// Convert one numeric scalar to a finite double, preserving explicit validation.
        double number(const Node& node, const std::string& path) {
            scalar(node, path);
            const double value = node.to_float64();
            if (!std::isfinite(value))
                fail(path, "value must be finite");
            return value;
        }

        /// Avoid a floating conversion for integer YAML nodes and reject fractional indices.
        std::int64_t integer64(const Node& node, const std::string& path) {
            scalar(node, path);
            if (node.dtype().is_signed_integer())
                return node.to_int64();
            if (node.dtype().is_unsigned_integer()) {
                const auto value = node.to_uint64();
                if (value > std::uint64_t(std::numeric_limits<std::int64_t>::max()))
                    fail(path, "integer is outside the supported range");
                return static_cast<std::int64_t>(value);
            }
            const double value = number(node, path), limit = std::ldexp(1.0, 63);
            if (std::trunc(value) != value || value < -limit || value >= limit)
                fail(path, "expected an in-range integer");
            return static_cast<std::int64_t>(value);
        }

        /// Reject integer overflow before narrowing step counts or cadences.
        int integer(const Node& node, const std::string& path) {
            const auto value = integer64(node, path);
            if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
                fail(path, "integer is outside the supported range");
            return static_cast<int>(value);
        }

        /// Conduit's YAML versions represent Boolean values as strings or numeric zero/one.
        bool boolean(const Node& node, const std::string& path) {
            if (node.dtype().is_string()) {
                if (node.as_string() == "true")
                    return true;
                if (node.as_string() == "false")
                    return false;
            } else if (node.dtype().is_number()) {
                const double value = number(node, path);
                if (value == 0 || value == 1)
                    return value == 1;
            }
            fail(path, "expected true or false");
        }

        /// Preserve string typing; numerical case names or paths are not silently accepted.
        std::string string(const Node& node, const std::string& path) {
            if (!node.dtype().is_string())
                fail(path, "expected a string");
            return node.as_string();
        }

        /// Decode duration or step count; the application's numerical step is never adjusted.
        void readDuration(const Node& run, RunInput& input) {
            keys(run, "run", {"steps", "duration"});
            if (run.has_path("steps") && run.has_path("duration"))
                fail("run", "steps and duration are mutually exclusive");
            if (run.has_path("steps")) {
                input.steps_m = integer(run["steps"], "run/steps");
                if (*input.steps_m < 0)
                    fail("run/steps", "must be nonnegative");
            }
            if (run.has_path("duration")) {
                input.duration_m = number(run["duration"], "run/duration");
                if (*input.duration_m < 0)
                    fail("run/duration", "must be nonnegative");
            }
        }

        /// Decode only common run settings; case-specific limits remain with the application.
        void readRun(const Node& root, RunInput& input) {
            if (root.has_path("run"))
                readDuration(root["run"], input);
            if (root.has_path("output")) {
                const auto output = root["output"];
                keys(output, "output", {"directory"});
                input.outputDirectory_m = string(required(output, "directory"), "output/directory");
                if (input.outputDirectory_m.empty())
                    fail("output/directory", "must not be empty; omit output for stdout-only runs");
            }
        }

        /// Require SI scales only when values are explicitly requested in SI units.
        UnitSystem readUnits(const Node& node) {
            UnitSystem result;
            result.system_m = string(required(node, "system"), "units/system");
            if (result.isSI()) {
                keys(node, "units", {"system", "length_m", "charge_density_C_m3"});
                result.length_m = number(required(node, "length_m"), "units/length_m");
                result.chargeDensityScale_m =
                    number(required(node, "charge_density_C_m3"), "units/charge_density_C_m3");
            } else
                keys(node, "units", {"system"});
            result.validate();
            return result;
        }

        /// Reject intrinsically invalid overrides even if their observation is disabled.
        void validateObservation(const ObservationOptions& options) {
            simulation::ObservationSchedule schedule;
            schedule.first_m        = options.first_m.value_or(0);
            schedule.every_m        = options.every_m.value_or(1);
            schedule.last_m         = options.last_m.value_or(-1);
            schedule.includeFinal_m = options.includeFinal_m.value_or(false);
            schedule.validate();
        }

        /// Parse optional schedule fields independently, keeping unspecified defaults intact.
        ObservationOptions readObservation(const Node& node, const std::string& path) {
            keys(node, path, {"first", "every", "last", "include_final", "enabled"});
            ObservationOptions result;
            if (node.has_path("first"))
                result.first_m = integer64(node["first"], path + "/first");
            if (node.has_path("every"))
                result.every_m = integer(node["every"], path + "/every");
            if (node.has_path("last"))
                result.last_m = integer64(node["last"], path + "/last");
            if (node.has_path("include_final"))
                result.includeFinal_m = boolean(node["include_final"], path + "/include_final");
            if (node.has_path("enabled"))
                result.enabled_m = boolean(node["enabled"], path + "/enabled");
            validateObservation(result);
            return result;
        }

        /// Preserve observation names as direct identifiers rather than Conduit paths.
        void readObservations(const Node& node, RunInput& input) {
            mapping(node, "observations");
            for (conduit_index_t i = 0; i < node.number_of_children(); ++i) {
                const auto child = node.child(i);
                const auto name  = child.name();
                if (name.empty() || name.find('/') != std::string::npos)
                    fail("observations", "observation names must be nonempty direct keys");
                input.observations_m.emplace(name, readObservation(child, "observations/" + name));
            }
        }

        /// Serialize exactly the explicit schedule overrides, including a false Boolean.
        void observationRecord(Node node, const ObservationOptions& options) {
            if (options.first_m)
                node["first"] = *options.first_m;
            if (options.every_m)
                node["every"] = *options.every_m;
            if (options.last_m)
                node["last"] = *options.last_m;
            if (options.includeFinal_m)
                node["include_final"] = *options.includeFinal_m ? 1 : 0;
            if (options.enabled_m)
                node["enabled"] = *options.enabled_m ? 1 : 0;
        }

        /// Decode the identity and unit convention shared by every application.
        void readIdentity(const Node& root, RunInput& input) {
            keys(
                root, "root",
                {"schema_version", "case", "parameters", "run", "output", "observations", "units"});
            if (integer(required(root, "schema_version"), "schema_version") != 1)
                fail("schema_version", "only version 1 is supported");
            input.caseName_m = string(required(root, "case"), "case");
            if (input.caseName_m.empty())
                fail("case", "must not be empty");
            input.units_m = readUnits(required(root, "units"));
        }

        /// Convert a physical duration to an integer count without changing the numerical step.
        int durationSteps(double duration, double timeFactor, double codeTimeStep) {
            if (!std::isfinite(codeTimeStep) || codeTimeStep <= 0)
                fail("run/duration", "time step must be finite and positive");
            const double ratio   = duration / timeFactor / codeTimeStep;
            const double rounded = std::round(ratio);
            const double tolerance =
                32 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(ratio));
            if (!std::isfinite(ratio) || ratio < 0 || rounded > std::numeric_limits<int>::max()
                || std::abs(ratio - rounded) > tolerance)
                fail("run/duration",
                     "duration must be an integer number of time steps within roundoff");
            return static_cast<int>(rounded);
        }

        /// Preserve explicit SI scales without assigning them to normalized-only input.
        void unitInputRecord(Node node, const UnitSystem& units) {
            node["system"] = units.system_m;
            if (units.isSI()) {
                node["length_m"]            = units.length_m;
                node["charge_density_C_m3"] = units.chargeDensityScale_m;
            }
        }
    }  // namespace

    RunInput readRunInput(const std::string& path) {
        Node root;
        conduit_node_load(conduit_cpp::c_node(&root), path.c_str(), "yaml");
        RunInput input;
        readIdentity(root, input);
        if (root.has_path("parameters")) {
            mapping(root["parameters"], "parameters");
            input.parameters_m.set(root["parameters"]);
        }
        readRun(root, input);
        if (root.has_path("observations"))
            readObservations(root["observations"], input);
        return input;
    }

    void validateParameters(const RunInput& input, std::initializer_list<std::string> allowed) {
        if (!input.parameters_m.dtype().is_empty())
            keys(input.parameters_m, "parameters", allowed);
    }

    void validateObservationNames(const RunInput& input,
                                  std::initializer_list<std::string> allowed) {
        for (const auto& [name, options] : input.observations_m)
            if (std::find(allowed.begin(), allowed.end(), name) == allowed.end())
                fail("observations/" + name, "unsupported observation for this application");
    }

    bool booleanParameter(const RunInput& input, const std::string& name, bool fallback) {
        return input.parameters_m.has_path(name)
                   ? boolean(input.parameters_m[name], "parameters/" + name)
                   : fallback;
    }

    bool observationEnabled(const RunInput& input, const std::string& name, bool fallback) {
        const auto found = input.observations_m.find(name);
        return found == input.observations_m.end() ? fallback
                                                   : found->second.enabled_m.value_or(fallback);
    }

    simulation::ObservationSchedule observationSchedule(const RunInput& input,
                                                        const std::string& name,
                                                        simulation::ObservationSchedule fallback) {
        if (const auto found = input.observations_m.find(name);
            found != input.observations_m.end()) {
            const auto& options     = found->second;
            fallback.first_m        = options.first_m.value_or(fallback.first_m);
            fallback.every_m        = options.every_m.value_or(fallback.every_m);
            fallback.last_m         = options.last_m.value_or(fallback.last_m);
            fallback.includeFinal_m = options.includeFinal_m.value_or(fallback.includeFinal_m);
        }
        fallback.validate();
        return fallback;
    }

    int resolvedSteps(const RunInput& input, double codeTimeStep, int defaultSteps) {
        if (input.steps_m && input.duration_m)
            fail("run", "steps and duration are mutually exclusive");
        if (!input.duration_m) {
            const int result = input.steps_m.value_or(defaultSteps);
            if (result < 0)
                fail("run/steps", "must be nonnegative");
            return result;
        }
        return durationSteps(*input.duration_m, input.units_m.timeFactor(), codeTimeStep);
    }

    conduit_cpp::Node runInputRecord(const RunInput& input) {
        Node result;
        result["schema_version"] = 1;
        result["case"]           = input.caseName_m;
        unitInputRecord(result["units"], input.units_m);
        if (!input.parameters_m.dtype().is_empty())
            result["parameters"].set(input.parameters_m);
        if (input.steps_m)
            result["run/steps"] = *input.steps_m;
        if (input.duration_m)
            result["run/duration"] = *input.duration_m;
        if (!input.outputDirectory_m.empty())
            result["output/directory"] = input.outputDirectory_m;
        for (const auto& [name, options] : input.observations_m)
            observationRecord(result["observations/" + name], options);
        return result;
    }
}  // namespace chdr::input
