include(AddIpplTest)
add_ippl_test(
  test-chdr-solver
  COMPILE_ONLY
  SOURCES
  ../../vacuum/tests/TestMain.cpp
  INCLUDE_DIRS
  "${CMAKE_CURRENT_SOURCE_DIR}/.."
  "${CMAKE_CURRENT_SOURCE_DIR}/../../vacuum"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../fel")
target_link_libraries(test-chdr-solver PRIVATE chdr-solver)
target_compile_options(test-chdr-solver PRIVATE -UNDEBUG)
target_link_libraries(test-chdr-solver PRIVATE nlohmann_json::nlohmann_json)
find_package(Python3 QUIET COMPONENTS Interpreter)
set(CHDR_SOLVER_RECORD_ROOT "${CMAKE_CURRENT_BINARY_DIR}/runs"
    CACHE PATH "Directory containing solver test records and rank fixtures")

function(chdr_solver_skip name reason)
  add_test(NAME chdr.solver.${name} COMMAND ${CMAKE_COMMAND} -E echo "SKIP: ${reason}")
  set_tests_properties(chdr.solver.${name} PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:" LABELS
                                                      "chdr;solver;unit;unavailable")
endfunction()

function(chdr_solver_test suite filter)
  add_test(NAME chdr.solver.${suite}.serial COMMAND test-chdr-solver --gtest_filter=${filter})
  set_tests_properties(
    chdr.solver.${suite}.serial PROPERTIES TIMEOUT 60 LABELS "chdr;solver;unit" ENVIRONMENT
                                           "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  if(MPIEXEC_EXECUTABLE)
    add_test(NAME chdr.solver.${suite}.mpi2
             COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                     $<TARGET_FILE:test-chdr-solver> ${MPIEXEC_POSTFLAGS} --gtest_filter=${filter})
    set_tests_properties(
      chdr.solver.${suite}.mpi2
      PROPERTIES TIMEOUT
                 60
                 PROCESSORS
                 2
                 LABELS
                 "chdr;solver;unit;mpi"
                 ENVIRONMENT
                 "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  endif()
endfunction()
target_sources(test-chdr-solver PRIVATE test_setup.cpp)
chdr_solver_test(setup "PeriodicSetup.*")
