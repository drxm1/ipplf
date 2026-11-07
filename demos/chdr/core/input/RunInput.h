/** @file RunInput.h
 * @brief Common run options and strict YAML input shared by field applications.
 */
#pragma once

#include <catalyst_conduit.hpp>
#include <cstdint>
#include <initializer_list>
#include <map>
#include <optional>
#include <string>

#include "core/input/UnitSystem.h"
#include "core/simulation/ObservationSchedule.h"

namespace chdr::input {
    /** @brief Optional overrides of an application's documented observation defaults.
     * @details An omitted key preserves its default independently of other keys.
     * Disabling an observation therefore does not change its sampling schedule.
     */
    struct ObservationOptions {
        std::optional<int> every_m;
        std::optional<std::int64_t> first_m, last_m;
        std::optional<bool> includeFinal_m, enabled_m;
    };

    /** @brief Common run configuration read from the shared YAML schema.
     * @details Each application validates its case name, supported parameter keys and
     * observation names before allocating fields. Empty outputDirectory_m retains
     * stdout-only runs. parameters_m owns its data independently of the input document.
     */
    struct ChdrYamlRunInput {
        std::string caseName_m, outputDirectory_m;
        std::optional<int> steps_m;
        std::optional<double> duration_m;
        std::map<std::string, ObservationOptions> observations_m;
        UnitSystem units_m;
        conduit_cpp::Node parameters_m;
    };

    /** @brief Read schema version 1 without accepting unknown common configuration keys.
     * @details Required keys are schema_version, case and units. Optional parameters
     * and observations are mappings. Application-specific names are checked separately.
     * All ranks read the same file before simulation setup; this function is host-only.
     * @throws std::invalid_argument For malformed or unsupported settings.
     */
    ChdrYamlRunInput readRunInput(const std::string& path);

    /// Reject parameter names not implemented by this application, including misspellings.
    void validateParameters(const ChdrYamlRunInput& input,
                            std::initializer_list<std::string> allowed);
    /// Reject observation names not implemented by this application, even when disabled.
    void validateObservationNames(const ChdrYamlRunInput& input,
                                  std::initializer_list<std::string> allowed);

    /// Read a Boolean parameter in Conduit's true/false or zero/one representation.
    bool booleanParameter(const ChdrYamlRunInput& input, const std::string& name, bool fallback);
    /// Return the selected enabled state, retaining the caller's default when omitted.
    bool observationEnabled(const ChdrYamlRunInput& input, const std::string& name,
                            bool fallback = true);
    /** @brief Merge explicitly supplied schedule fields into the application's defaults.
     * @details The combined schedule is validated, including interactions between a
     * supplied endpoint and an endpoint inherited from the application.
     */
    simulation::ObservationSchedule observationSchedule(const ChdrYamlRunInput& input,
                                                        const std::string& name,
                                                        simulation::ObservationSchedule fallback);

    /** @brief Resolve steps or duration without changing the selected numerical time step.
     * @details Duration is expressed in the selected time units. Its ratio to the actual
     * time step must be an integer within floating-point roundoff; no truncation is used.
     */
    int resolvedSteps(const ChdrYamlRunInput& input, double codeTimeStep, int defaultSteps);

    /** @brief Serialize supplied options for the common run record.
     * @details This preserves explicit overrides. The application's effective run
     * description records defaults resolved from its case and actual solver state.
     */
    conduit_cpp::Node runInputRecord(const ChdrYamlRunInput& input);
}  // namespace chdr::input
