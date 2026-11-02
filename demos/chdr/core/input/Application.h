/** @file Application.h
 * @brief Common input selection, IPPL lifetime and failure handling for field experiments.
 */
#pragma once

#include <algorithm>
#include <iostream>
#include <locale>
#include <string>
#include <vector>

#include "Utility/IpplException.h"

#include "core/base/Util.h"
#include "core/input/RunInput.h"

namespace chdr::input {
    /** @brief Remove IPPL options already consumed during initialization from application input.
     * @details Retain all experiment arguments for either the common YAML selector or the
     * compatibility parser. A missing IPPL option value is still rejected explicitly.
     */
    inline std::vector<std::string> applicationArguments(int argc, char** argv) {
        std::vector<std::string> result;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--info" || arg == "-i" || arg == "--overallocate" || arg == "-b"
                || arg == "--timer-fences") {
                if (++i >= argc)
                    throw std::invalid_argument("Missing value after " + arg);
            } else if (arg != "--debug" && arg != "-g") {
                result.push_back(arg);
            }
        }
        return result;
    }

    /** @brief Select a strict YAML file or translate the application's historical arguments.
     * @details YAML selection is exactly --input FILE; mixing it with legacy overrides is
     * rejected so a saved input cannot silently differ from the simulation that was run.
     */
    template <class LegacyParser>
    RunInput applicationInput(const std::vector<std::string>& args, LegacyParser parseLegacy) {
        if (std::find(args.begin(), args.end(), "--input") == args.end())
            return parseLegacy(args);
        if (args.size() != 2 || args.front() != "--input")
            throw std::invalid_argument("Use --input FILE without legacy case overrides");
        return readRunInput(args[1]);
    }

    /// Execute one input while retaining IPPL for rank-aware fatal-error reporting.
    template <class LegacyParser, class Execute>
    int executeApplication(int argc, char** argv, LegacyParser parseLegacy, Execute execute) {
        try {
            std::cout.imbue(std::locale::classic());
            execute(applicationInput(applicationArguments(argc, argv), parseLegacy));
            return 0;
        } catch (const IpplException& error) {
            return util::abortFailure(error.where(), error.what());
        } catch (const std::exception& error) {
            return util::abortFailure("CHDR application", error.what());
        }
    }

    /** @brief Run an experiment with the common input and resource-lifetime rules.
     * @param parseLegacy Translate existing arguments into RunInput; it does not run a case.
     * @param execute Configure and run the case from RunInput, including its observations.
     * @details Simulation and output objects created inside execute are destroyed before
     * IPPL finalization. Collective callbacks must be entered in the same order on all ranks.
     */
    template <class LegacyParser, class Execute>
    int runApplication(int argc, char** argv, LegacyParser parseLegacy, Execute execute) {
        ippl::initialize(argc, argv);
        const int result = executeApplication(argc, argv, parseLegacy, execute);
        ippl::finalize();
        return result;
    }
}  // namespace chdr::input
