/** @file Application.cpp
 * @brief YAML selection and shared IPPL lifetime for field experiments.
 */
#include "core/input/Application.h"

#include <iostream>
#include <locale>
#include <stdexcept>

#include "Utility/IpplException.h"

#include "core/base/Util.h"

namespace chdr::input {
    namespace {
        /** @brief Read input and execute while IPPL is available for fatal-error reporting.
         * @details Input and simulation objects are released
         * before returning to the function that finalizes IPPL.
         */
        int executeApplication(int argc, char** argv, ExperimentFunction execute) {
            try {
                std::cout.imbue(std::locale::classic());
                const auto input = applicationInput(applicationArguments(argc, argv));
                execute(input);
                return 0;
            } catch (const IpplException& error) {
                return util::abortFailure(error.where(), error.what());
            } catch (const std::exception& error) {
                return util::abortFailure("CHDR application", error.what());
            }
        }

    }  // namespace

    std::vector<std::string> applicationArguments(int argc, char** argv) {
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

    RunInput applicationInput(const std::vector<std::string>& args) {
        if (args.size() != 2 || args.front() != "--input")
            throw std::invalid_argument("Use --input FILE without additional arguments");
        return readRunInput(args[1]);
    }

    int runApplication(int argc, char** argv, ExperimentFunction execute) {
        ippl::initialize(argc, argv);
        const int result = executeApplication(argc, argv, execute);
        ippl::finalize();
        return result;
    }
}  // namespace chdr::input
