add_executable(experiment-e15 run_experiment.cpp)
target_link_libraries(experiment-e15 PRIVATE chdr-application nlohmann_json::nlohmann_json)
