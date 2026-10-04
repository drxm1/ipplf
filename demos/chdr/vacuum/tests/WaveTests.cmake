# Propagation and CLI checks; field-output and real-reader acceptance are deferred.
if(NOT Python3_Interpreter_FOUND OR NOT MPIEXEC_EXECUTABLE)
  if(NOT Python3_Interpreter_FOUND)
    set(_wave_unavailable "Python3 interpreter unavailable")
  else()
    set(_wave_unavailable "MPI launcher unavailable; the full three-case benchmark cannot run")
  endif()
  foreach(_test IN ITEMS propagation cli)
    chdr_vacuum_skip_test("chdr.vacuum.${_test}" "${_wave_unavailable}" integration)
  endforeach()
  return()
endif()
set(_wave_results "${CHDR_VACUUM_TEST_OUTPUT_DIR}/wave")
set(_wave_mpi_arguments "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}")
foreach(_flag IN LISTS MPIEXEC_PREFLAGS)
  list(APPEND _wave_mpi_arguments "--mpi-before=${_flag}")
endforeach()
foreach(_flag IN LISTS MPIEXEC_POSTFLAGS)
  list(APPEND _wave_mpi_arguments "--mpi-after=${_flag}")
endforeach()

add_test(
  NAME chdr.vacuum.propagation
  COMMAND
    "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_wave_checks.py" --binary
    $<TARGET_FILE:chdr-vacuum-wave> --results "${_wave_results}" --mpiexec "${MPIEXEC_EXECUTABLE}"
    ${_wave_mpi_arguments})
set_tests_properties(chdr.vacuum.propagation PROPERTIES LABELS "chdr;vacuum;integration" TIMEOUT
                                                        360 PROCESSORS 2)

add_test(
  NAME chdr.vacuum.cli
  COMMAND
    "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_cli_checks.py" --binary
    $<TARGET_FILE:chdr-vacuum-wave> --results "${CHDR_VACUUM_TEST_OUTPUT_DIR}/cli" --mpiexec
    "${MPIEXEC_EXECUTABLE}" ${_wave_mpi_arguments})
set_tests_properties(chdr.vacuum.cli PROPERTIES LABELS "chdr;vacuum;integration" TIMEOUT 180
                                                PROCESSORS 2)

set_tests_properties(chdr.vacuum.propagation chdr.vacuum.cli PROPERTIES ENVIRONMENT
                                                                        "PYTHONDONTWRITEBYTECODE=1")
