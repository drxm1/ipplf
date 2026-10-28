# Full-period and partial-run output are checked by separate actual-reader tests.
if(NOT Python3_Interpreter_FOUND OR NOT MPIEXEC_EXECUTABLE)
  if(NOT Python3_Interpreter_FOUND)
    set(_wave_unavailable "Python3 interpreter unavailable")
  else()
    set(_wave_unavailable "MPI launcher unavailable; the full three-case benchmark cannot run")
  endif()
  foreach(_test IN ITEMS propagation cli)
    chdr_vacuum_skip_test("chdr.vacuum.${_test}" "${_wave_unavailable}" integration)
  endforeach()
  foreach(_test IN ITEMS readback readback.partial)
    chdr_vacuum_skip_test("chdr.vacuum.${_test}" "${_wave_unavailable}" reader)
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

# Fixtures force fresh successful runs before either readback consumes its manifest.
set_tests_properties(chdr.vacuum.propagation PROPERTIES FIXTURES_SETUP ChdrVacuumWaveOutput)
set_tests_properties(chdr.vacuum.cli PROPERTIES FIXTURES_SETUP ChdrVacuumCliOutput)

add_test(
  NAME chdr.vacuum.readback
  COMMAND
    "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_optional_reader.py"
    "${CHDR_VACUUM_PVPYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/verify_wave_output.py" --manifest
    "${_wave_results}/wave-runs.json" --summary "${_wave_results}/wave-readback.json")
add_test(
  NAME chdr.vacuum.readback.partial
  COMMAND
    "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_optional_reader.py"
    "${CHDR_VACUUM_PVPYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/verify_partial_output.py" --manifest
    "${CHDR_VACUUM_TEST_OUTPUT_DIR}/cli/cli-checks.json" --summary
    "${CHDR_VACUUM_TEST_OUTPUT_DIR}/cli/partial-readback.json")
set_tests_properties(
  chdr.vacuum.readback chdr.vacuum.readback.partial
  PROPERTIES LABELS
             "chdr;vacuum;reader"
             TIMEOUT
             120
             PROCESSORS
             1
             SKIP_RETURN_CODE
             77
             ENVIRONMENT
             "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;PYTHONDONTWRITEBYTECODE=1")
if(CHDR_VACUUM_PVPYTHON AND EXISTS "${CHDR_VACUUM_PVPYTHON}")
  set_tests_properties(chdr.vacuum.readback PROPERTIES FIXTURES_REQUIRED ChdrVacuumWaveOutput)
  set_tests_properties(chdr.vacuum.readback.partial PROPERTIES FIXTURES_REQUIRED
                                                               ChdrVacuumCliOutput)
endif()
