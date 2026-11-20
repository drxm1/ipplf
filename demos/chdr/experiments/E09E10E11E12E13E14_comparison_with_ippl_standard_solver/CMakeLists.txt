add_executable(experiments-comparison-with-ippl-standard-solver run_experiment.cpp)
target_link_libraries(experiments-comparison-with-ippl-standard-solver
                      PRIVATE chdr-application nlohmann_json::nlohmann_json)
