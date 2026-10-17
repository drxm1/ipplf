add_test(
  NAME chdr.interface.unit.smoke
  COMMAND
    "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:interface-demo>"
    "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/smoke-frames" -P
    "${CMAKE_CURRENT_SOURCE_DIR}/smoke.cmake")
set_tests_properties(
  chdr.interface.unit.smoke PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false" TIMEOUT
                                       60 LABELS "chdr;interface;unit")
