/** @file E01E02AnalysisProcess.h
 * @brief Invoke the bundled analysis once, after the application has finalized IPPL.
 */
#pragma once

#include <array>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <iostream>
#include <spawn.h>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>

#include "experiments/E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01E02Study.h"

extern char** environ;

namespace chdr::experiments::e01e02 {
    namespace pulseAnalysisProcessDetail {
        /// Convert POSIX spawn setup errors, whose error number is the return value.
        inline void requireSuccess(int error) {
            if (error != 0)
                throw std::runtime_error(std::string("Cannot launch Pulse analysis: ")
                                         + std::strerror(error));
        }

        /// Own the child-only redirection; the parent keeps its original console streams.
        class AnalysisLog {
        public:
            explicit AnalysisLog(const std::filesystem::path& path) {
                requireSuccess(posix_spawn_file_actions_init(&actions_m));
                try {
                    requireSuccess(
                        posix_spawn_file_actions_addopen(&actions_m, STDOUT_FILENO, path.c_str(),
                                                         O_WRONLY | O_CREAT | O_EXCL, 0644));
                    requireSuccess(
                        posix_spawn_file_actions_adddup2(&actions_m, STDOUT_FILENO, STDERR_FILENO));
                } catch (...) {
                    posix_spawn_file_actions_destroy(&actions_m);
                    throw;
                }
            }
            ~AnalysisLog() { posix_spawn_file_actions_destroy(&actions_m); }
            AnalysisLog(const AnalysisLog&)            = delete;
            AnalysisLog& operator=(const AnalysisLog&) = delete;
            const posix_spawn_file_actions_t* actions() const { return &actions_m; }

        private:
            posix_spawn_file_actions_t actions_m;
        };

        /// Pass executable and script as separate arguments, without command-string evaluation.
        inline pid_t start(const PulseAnalysisCommand& command) {
            std::string python = command.pythonExecutable_m;
            std::string script = command.script_m.string();
            std::array<char*, 3> arguments{python.data(), script.data(), nullptr};
            AnalysisLog log(command.script_m.parent_path() / "analysis.log");
            pid_t child{};
            requireSuccess(posix_spawnp(&child, python.c_str(), log.actions(), nullptr,
                                        arguments.data(), environ));
            return child;
        }

        /// Wait through interrupting signals and preserve the analysis program's exit status.
        inline int wait(pid_t child) {
            int status{};
            while (waitpid(child, &status, 0) < 0) {
                if (errno != EINTR)
                    throw std::runtime_error(std::string("Cannot wait for Pulse analysis: ")
                                             + std::strerror(errno));
            }
            if (WIFEXITED(status))
                return WEXITSTATUS(status);
            return WIFSIGNALED(status) ? 128 + WTERMSIG(status) : 1;
        }
    }  // namespace pulseAnalysisProcessDetail

    /** @brief Run Python on the former rank-zero process after IPPL/MPI finalization.
     * @details A failed analysis does not remove completed simulation results. Its output
     * is retained in analysis.log, and a nonzero exit status reaches the calling process.
     * @pre All simulation objects are destroyed and input::runApplication has returned.
     */
    inline int runPulseAnalysis(const PulseAnalysisCommand& command) {
        try {
            const int status =
                pulseAnalysisProcessDetail::wait(pulseAnalysisProcessDetail::start(command));
            std::cout << "Pulse analysis " << (status == 0 ? "completed" : "failed") << "; see "
                      << command.script_m.parent_path() / "analysis.log" << '\n';
            return status;
        } catch (const std::exception& error) {
            std::cerr << error.what() << "; completed study results are retained in "
                      << command.script_m.parent_path() << '\n';
            return 1;
        }
    }
}  // namespace chdr::experiments::e01e02
