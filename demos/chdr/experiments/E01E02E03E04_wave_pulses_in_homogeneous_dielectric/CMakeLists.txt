add_executable(experiments-wave-pulses-in-homogeneous-dielectric run_experiment.cpp E03.cpp E04.cpp)
target_link_libraries(experiments-wave-pulses-in-homogeneous-dielectric
                      PRIVATE chdr-application nlohmann_json::nlohmann_json)

if(BUILD_TESTING AND IPPL_ENABLE_UNIT_TESTS)
  include(AddIpplTest)
  add_ippl_test(test-chdr-e03-reference COMPILE_ONLY SOURCES
                tests/TestE03AnalyticContinuumReference.cpp)
  target_link_libraries(test-chdr-e03-reference PRIVATE chdr-testing)
  add_test(NAME chdr.experiments.e03.reference COMMAND test-chdr-e03-reference)
  set_tests_properties(chdr.experiments.e03.reference PROPERTIES TIMEOUT 30 LABELS "chdr;unit"
                                                                 ENVIRONMENT "OMP_NUM_THREADS=1")
  add_ippl_test(test-chdr-e04-monitor COMPILE_ONLY SOURCES tests/TestE04Monitor.cpp)
  target_link_libraries(test-chdr-e04-monitor PRIVATE chdr-application chdr-testing)
  add_test(NAME chdr.experiments.e04.monitor.serial COMMAND test-chdr-e04-monitor)
  add_test(NAME chdr.experiments.e04.monitor.mpi2
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 --bind-to none
                   $<TARGET_FILE:test-chdr-e04-monitor>)
  set_tests_properties(
    chdr.experiments.e04.monitor.serial chdr.experiments.e04.monitor.mpi2
    PROPERTIES TIMEOUT 30 LABELS "chdr;unit" ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  set_tests_properties(chdr.experiments.e04.monitor.mpi2 PROPERTIES PROCESSORS 2)
endif()
