/** @file Application.h
 * @brief Common input selection, IPPL lifetime and failure handling for field experiments.
 */
#pragma once

#include <string>
#include <vector>

#include "core/input/RunInput.h"

namespace chdr::input {
    /// Configure and run one experiment from the supplied input; throw on failure.
    using ExperimentFunction = void (*)(const RunInput&);

    /** @brief Remove IPPL options already consumed during initialization from application input.
     * @details Retain the experiment arguments for the common YAML selector.
     * A missing IPPL option value is still rejected explicitly.
     */
    std::vector<std::string> applicationArguments(int argc, char** argv);

    /** @brief Read the YAML file selected by exactly --input FILE.
     * @throws std::invalid_argument For missing, additional or malformed arguments.
     */
    RunInput applicationInput(const std::vector<std::string>& args);

    /** @brief Run an experiment with YAML input selected by --input FILE.
     * @param execute Function taking const RunInput& and returning void. It validates the
     * experiment settings, constructs the simulation and performs its measurements/output.
     * @details Initialize IPPL before reading input and finalize it after execute returns.
     * Simulation and output objects must be destroyed before execute returns.
     * Input and execution exceptions use the shared MPI-aware fatal-error handler.
     * All ranks must enter collective operations in the same order.
     * @pre execute is non-null.
     */
    int runApplication(int argc, char** argv, ExperimentFunction execute);

}  // namespace chdr::input
