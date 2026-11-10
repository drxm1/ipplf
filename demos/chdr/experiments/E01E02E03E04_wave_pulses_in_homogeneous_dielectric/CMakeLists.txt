add_executable(experiments-wave-pulses-in-homogeneous-dielectric run_experiment.cpp E03.cpp E04.cpp)
target_link_libraries(experiments-wave-pulses-in-homogeneous-dielectric
                      PRIVATE chdr-application nlohmann_json::nlohmann_json)
