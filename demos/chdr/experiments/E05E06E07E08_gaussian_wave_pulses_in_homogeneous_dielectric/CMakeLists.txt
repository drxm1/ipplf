add_executable(experiments-gaussian-wave-pulses-in-homogeneous-dielectric run_experiment.cpp)
target_link_libraries(experiments-gaussian-wave-pulses-in-homogeneous-dielectric
                      PRIVATE chdr-application nlohmann_json::nlohmann_json)
