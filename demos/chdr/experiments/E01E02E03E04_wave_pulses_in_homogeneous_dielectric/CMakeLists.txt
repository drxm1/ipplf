add_executable(experiments-wave-pulses-in-homogeneous-dielectric run_experiment.cpp E03.cpp)
target_link_libraries(experiments-wave-pulses-in-homogeneous-dielectric
                      PRIVATE chdr-application nlohmann_json::nlohmann_json)

# Bundle the reviewed Python source with this build, then copy it into each study.
configure_file(analyze_e01_e02.py analyze_e01_e02.py COPYONLY)
target_compile_definitions(
  experiments-wave-pulses-in-homogeneous-dielectric
  PRIVATE CHDR_E01_E02_ANALYSIS_SOURCE="${CMAKE_CURRENT_BINARY_DIR}/analyze_e01_e02.py")

if(BUILD_TESTING AND IPPL_ENABLE_UNIT_TESTS)
  include(AddIpplTest)
  add_ippl_test(test-chdr-e03-reference COMPILE_ONLY SOURCES
                tests/TestE03AnalyticContinuumReference.cpp)
  target_link_libraries(test-chdr-e03-reference PRIVATE chdr-testing)
  add_test(NAME chdr.experiments.e03.reference COMMAND test-chdr-e03-reference)
  set_tests_properties(chdr.experiments.e03.reference PROPERTIES TIMEOUT 30 LABELS "chdr;unit"
                                                                 ENVIRONMENT "OMP_NUM_THREADS=1")
endif()
