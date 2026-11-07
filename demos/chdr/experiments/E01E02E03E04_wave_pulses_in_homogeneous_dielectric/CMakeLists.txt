add_executable(experiments-wave-pulses-in-homogeneous-dielectric run_experiment.cpp)
target_link_libraries(experiments-wave-pulses-in-homogeneous-dielectric
                      PRIVATE chdr-application nlohmann_json::nlohmann_json)

# Bundle the reviewed Python source with this build, then copy it into each study.
configure_file(analyze_e01.py analyze_e01.py COPYONLY)
target_compile_definitions(
  experiments-wave-pulses-in-homogeneous-dielectric
  PRIVATE CHDR_E01_ANALYSIS_SOURCE="${CMAKE_CURRENT_BINARY_DIR}/analyze_e01.py")
