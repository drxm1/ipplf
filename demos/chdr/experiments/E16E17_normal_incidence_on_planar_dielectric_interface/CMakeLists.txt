add_executable(chdr-experiment-e16-e17 run_experiment.cpp)
target_link_libraries(chdr-experiment-e16-e17 PRIVATE chdr-application)

if(BUILD_TESTING)
  add_executable(test-e16-e17-reference test_analytic_reference.cpp)
  target_link_libraries(test-e16-e17-reference PRIVATE chdr-application)
  add_test(NAME chdr.experiments.e16-e17.reference COMMAND test-e16-e17-reference)
  set_tests_properties(chdr.experiments.e16-e17.reference PROPERTIES LABELS record)
endif()
