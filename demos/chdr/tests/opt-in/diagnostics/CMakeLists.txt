include(AddIpplTest)
add_ippl_test(
  test-chdr-diagnostics
  COMPILE_ONLY
  SOURCES
  test_interface_energy.cpp
  test_monitors.cpp
  test_monitor_records.cpp
  test_probe.cpp
  test_mur_layer.cpp)
# SerialReference.h (borrowed from tests/core/solver) brings in RankRecords.h, which uses JSON.
target_link_libraries(test-chdr-diagnostics PRIVATE chdr-core chdr-diagnostics chdr-testing
                                                    nlohmann_json::nlohmann_json)
target_compile_options(test-chdr-diagnostics PRIVATE -UNDEBUG)
find_package(Python3 QUIET COMPONENTS Interpreter)
set(CHDR_DIAGNOSTICS_RECORD_ROOT "${CMAKE_CURRENT_BINARY_DIR}/runs"
    CACHE PATH "Directory for the records of the diagnostics tests")

function(chdr_diagnostics_skip name reason)
  add_test(NAME chdr.diagnostics.${name} COMMAND ${CMAKE_COMMAND} -E echo "SKIP: ${reason}")
  set_tests_properties(chdr.diagnostics.${name} PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:" LABELS
                                                           "chdr;diagnostics;unit;unavailable")
endfunction()

# Register both variants by default; SERIAL_ONLY or MPI2_ONLY keeps the specified rank count.
function(chdr_diagnostics_test suite filter)
  cmake_parse_arguments(DIAGNOSTICS "SERIAL_ONLY;MPI2_ONLY" "" "" ${ARGN})
  if(NOT DIAGNOSTICS_MPI2_ONLY)
    add_test(NAME chdr.diagnostics.${suite}.serial COMMAND test-chdr-diagnostics
                                                           --gtest_filter=${filter})
    set_tests_properties(
      chdr.diagnostics.${suite}.serial
      PROPERTIES TIMEOUT 60 LABELS "chdr;diagnostics;unit" ENVIRONMENT
                 "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  endif()
  if(MPIEXEC_EXECUTABLE AND NOT DIAGNOSTICS_SERIAL_ONLY)
    add_test(
      NAME chdr.diagnostics.${suite}.mpi2
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
              $<TARGET_FILE:test-chdr-diagnostics> ${MPIEXEC_POSTFLAGS} --gtest_filter=${filter})
    set_tests_properties(
      chdr.diagnostics.${suite}.mpi2
      PROPERTIES TIMEOUT
                 60
                 PROCESSORS
                 2
                 LABELS
                 "chdr;diagnostics;unit;mpi"
                 ENVIRONMENT
                 "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  endif()
endfunction()

# NormalSplit and TangentialSplit use the normal and tangential decomposition axes. Both compare
# with a one-rank evolution on MPI_COMM_SELF, so the two-rank entry also checks serial evolution.
chdr_diagnostics_test(energy "InterfaceEnergy.*" MPI2_ONLY)
chdr_diagnostics_test(monitors "ChdrMonitors.*")
chdr_diagnostics_test(monitor-records "ChdrMonitorRecords.*")
chdr_diagnostics_test(probe "FourierProbe.*")
chdr_diagnostics_test(mur-layer-refusals "MurLayer.ProbeAndMonitorRefusals" SERIAL_ONLY)
include(ProbeTests.cmake)
