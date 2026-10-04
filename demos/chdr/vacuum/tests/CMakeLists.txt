# One standard GoogleTest executable; run the same suites on one and two ranks.
include(AddIpplTest)

# Only CSV/CLI integration checks need Python; no field reader is needed in this unit.
find_package(Python3 QUIET COMPONENTS Interpreter)
set(CHDR_VACUUM_TEST_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/runs"
    CACHE PATH "Directory for bounded CHDR vacuum test logs and CSV")

function(chdr_vacuum_skip_test test_name reason)
  add_test(NAME "${test_name}" COMMAND "${CMAKE_COMMAND}" -E echo "SKIP: ${reason}")
  set_tests_properties("${test_name}" PROPERTIES LABELS "chdr;vacuum;unavailable;${ARGN}" TIMEOUT
                                                 10 SKIP_REGULAR_EXPRESSION "^SKIP:")
endfunction()

add_ippl_test(
  test-chdr-vacuum
  COMPILE_ONLY
  SOURCES
  TestMain.cpp
  test_first_step.cpp
  test_setup.cpp
  test_diagnostics.cpp
  test_policy.cpp
  ../VacuumSetup.cpp
  INCLUDE_DIRS
  "${CMAKE_CURRENT_SOURCE_DIR}/.."
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../fel")

# Exercise the actual debug guards in this test target, including VacuumSetup.cpp. This does not
# change the production executable's build configuration.
target_compile_options(
  test-chdr-vacuum PRIVATE $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-UNDEBUG>
                           $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/UNDEBUG>)

add_test(NAME chdr.vacuum.unit.serial COMMAND $<TARGET_FILE:test-chdr-vacuum>)
set_tests_properties(
  chdr.vacuum.unit.serial
  PROPERTIES LABELS
             "unit;chdr;vacuum"
             TIMEOUT
             60
             PROCESSORS
             1
             ENVIRONMENT
             "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")

# Smoke checks exercise the real driver; the GoogleTests establish numerical accuracy.
add_test(NAME chdr.vacuum.driver.serial COMMAND $<TARGET_FILE:chdr-vacuum-wave> --steps 1)
set_tests_properties(
  chdr.vacuum.driver.serial
  PROPERTIES LABELS
             "smoke;chdr;vacuum"
             TIMEOUT
             60
             PROCESSORS
             1
             ENVIRONMENT
             "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")

if(MPIEXEC_EXECUTABLE)
  add_test(NAME chdr.vacuum.unit.mpi2
           COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:test-chdr-vacuum> ${MPIEXEC_POSTFLAGS})
  add_test(
    NAME chdr.vacuum.unit.mpi2.shuffle
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test-chdr-vacuum> ${MPIEXEC_POSTFLAGS} --gtest_shuffle --gtest_random_seed=0)
  set_tests_properties(
    chdr.vacuum.unit.mpi2 chdr.vacuum.unit.mpi2.shuffle
    PROPERTIES LABELS
               "unit;chdr;vacuum;mpi"
               TIMEOUT
               60
               PROCESSORS
               2
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  add_test(NAME chdr.vacuum.driver.mpi2
           COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:chdr-vacuum-wave> ${MPIEXEC_POSTFLAGS} --steps 1)
  set_tests_properties(
    chdr.vacuum.driver.mpi2
    PROPERTIES LABELS
               "smoke;chdr;vacuum;mpi"
               TIMEOUT
               60
               PROCESSORS
               2
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/WaveTests.cmake")
