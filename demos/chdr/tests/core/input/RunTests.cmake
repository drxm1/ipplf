find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(
    NAME chdr.input.runs.si
    COMMAND
      "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_runs.py"
      $<TARGET_FILE:chdr-vacuum-wave> $<TARGET_FILE:interface-movies>
      $<TARGET_FILE:planar-interface>)
  set_tests_properties(
    chdr.input.runs.si
    PROPERTIES LABELS
               "chdr;input;integration"
               TIMEOUT
               240
               PROCESSORS
               1
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;PYTHONDONTWRITEBYTECODE=1")
else()
  add_test(NAME chdr.input.runs.si COMMAND ${CMAKE_COMMAND} -E echo
                                           "SKIP: SI run comparison requires Python3")
  set_tests_properties(chdr.input.runs.si PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:" LABELS
                                                     "chdr;input;integration;unavailable")
endif()
