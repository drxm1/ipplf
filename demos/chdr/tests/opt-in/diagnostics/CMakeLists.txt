include(AddIpplTest)
# test_probe.cpp and test_probe_full.cpp build their periodic fields with the vacuum-wave setup, as
# test_solver_contract.cpp of the solver tests does. VacuumSetup.cpp stays an app file until that
# setup moves onto BoxSetup.
add_ippl_test(
  test-chdr-diagnostics
  COMPILE_ONLY
  SOURCES
  ../../../apps/vacuum-wave/VacuumSetup.cpp
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

# One serial and one two-rank entry per suite; an optional third argument narrows the serial filter.
function(chdr_diagnostics_test suite filter)
  set(serial_filter "${filter}")
  if(ARGC GREATER 2)
    set(serial_filter "${ARGV2}")
  endif()
  add_test(NAME chdr.diagnostics.${suite}.serial COMMAND test-chdr-diagnostics
                                                         --gtest_filter=${serial_filter})
  set_tests_properties(
    chdr.diagnostics.${suite}.serial PROPERTIES TIMEOUT 60 LABELS "chdr;diagnostics;unit"
                                                ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  if(MPIEXEC_EXECUTABLE)
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

# InterfaceEnergy.NormalSplit and .TangentialSplit differ only in the MPI decomposition axis: along
# the slab normal x (the rank seam lies on the slab face) or along y. One rank has no decomposition,
# so its entry runs NormalSplit only; two ranks run both.
chdr_diagnostics_test(energy "InterfaceEnergy.*" "InterfaceEnergy.NormalSplit")
chdr_diagnostics_test(monitors "ChdrMonitors.*")
chdr_diagnostics_test(monitor-records "ChdrMonitorRecords.*")
chdr_diagnostics_test(probe "FourierProbe.*")
chdr_diagnostics_test(mur-layer-refusals "MurLayer.ProbeAndMonitorRefusals")
include(ProbeTests.cmake)
