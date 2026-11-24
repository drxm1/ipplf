add_executable(chdr-experiment-e18-e19 run_experiment.cpp)
target_link_libraries(chdr-experiment-e18-e19 PRIVATE chdr-application)

if(BUILD_TESTING)
  add_executable(test-e18-e19-reference test_analytic_reference.cpp)
  target_link_libraries(test-e18-e19-reference PRIVATE chdr-application)
  add_test(NAME chdr.experiments.e18-e19.reference COMMAND test-e18-e19-reference)
  set_tests_properties(chdr.experiments.e18-e19.reference PROPERTIES LABELS
                                                                     "chdr;experiment;record")

  add_executable(test-e18-e19-measurements test_measurements.cpp)
  target_link_libraries(test-e18-e19-measurements PRIVATE chdr-application)
  add_test(NAME chdr.experiments.e18-e19.measurements COMMAND test-e18-e19-measurements)
  set_tests_properties(chdr.experiments.e18-e19.measurements PROPERTIES LABELS
                                                                        "chdr;experiment;unit")
  set_tests_properties(chdr.experiments.e18-e19.reference chdr.experiments.e18-e19.measurements
                       PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  if(MPIEXEC_EXECUTABLE)
    add_test(NAME chdr.experiments.e18-e19.measurements.mpi2
             COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                     $<TARGET_FILE:test-e18-e19-measurements> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
      chdr.experiments.e18-e19.measurements.mpi2
      PROPERTIES LABELS "chdr;experiment;unit;mpi" PROCESSORS 2 ENVIRONMENT
                 "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  endif()
endif()
