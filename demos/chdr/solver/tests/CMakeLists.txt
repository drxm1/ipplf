include(AddIpplTest)
add_ippl_test(
  test-chdr-solver
  COMPILE_ONLY
  SOURCES
  ../../vacuum/tests/TestMain.cpp
  test_shell.cpp
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
target_sources(test-chdr-solver PRIVATE test_monitors.cpp)
target_sources(test-chdr-solver PRIVATE test_monitor_records.cpp)
target_sources(test-chdr-solver PRIVATE test_compact.cpp)
target_sources(test-chdr-solver PRIVATE test_staggered.cpp)

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
chdr_solver_test(shell "SolverShell.*")
chdr_solver_test(monitors "ChdrMonitors.*")
chdr_solver_test(monitor-records "ChdrMonitorRecords.*")
chdr_solver_test(compact "CompactOperators.*")
chdr_solver_test(staggered "ChdrStaggered.*")
target_sources(test-chdr-solver PRIVATE test_setup.cpp)
chdr_solver_test(setup "BoxSetup.*")
target_sources(test-chdr-solver PRIVATE test_operators.cpp)
chdr_solver_test(operators "Operators.*")
target_sources(test-chdr-solver PRIVATE test_timestep.cpp)
chdr_solver_test(timestep "TimeStep.*")
foreach(case StandardBounded AnisotropicBounded NonStandardBounded)
  chdr_solver_test(long.timestep.${case} "TimeStepLong.${case}")
  set_tests_properties(chdr.solver.long.timestep.${case}.serial PROPERTIES TIMEOUT 1800 LABELS
                                                                           "chdr;solver;long")
  if(MPIEXEC_EXECUTABLE)
    set_tests_properties(chdr.solver.long.timestep.${case}.mpi2 PROPERTIES TIMEOUT 1800 LABELS
                                                                           "chdr;solver;long;mpi")
  endif()
endforeach()

function(chdr_open_record name filter)
  add_test(NAME chdr.solver.open.record.${name} COMMAND test-chdr-solver --gtest_filter=${filter})
  set_tests_properties(
    chdr.solver.open.record.${name} PROPERTIES LABELS "chdr;solver;record" TIMEOUT 1800 ENVIRONMENT
                                               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endfunction()
chdr_open_record(zero-ghost OpenOuterLayerRecords.ZeroGhostMagnitude)
target_sources(test-chdr-solver PRIVATE test_gather.cpp)
chdr_solver_test(gather "RankGather.*")
target_sources(test-chdr-solver PRIVATE test_manufactured.cpp)
chdr_solver_test(manufactured "Manufactured.*")
target_sources(test-chdr-solver PRIVATE test_ranks.cpp)
chdr_solver_test(rank_cases "RankCases.*")
include(RankTests.cmake)
include(ConvergenceTests.cmake)
chdr_solver_test(wave_norms "WaveNorms.*")
include(LongTests.cmake)
include(DispersionTests.cmake)
target_sources(test-chdr-solver PRIVATE test_endpoint_records.cpp)
chdr_solver_test(record.endpoint "TimeStepRecords.*")
target_sources(test-chdr-solver PRIVATE test_serial_reference.cpp)
chdr_solver_test(serial_reference "SerialReference.*")
target_sources(test-chdr-solver PRIVATE test_statistics.cpp)
chdr_solver_test(statistics "StableStatistics.*")
target_sources(test-chdr-solver PRIVATE test_levels.cpp)
chdr_solver_test(levels "Levels.*")
target_sources(test-chdr-solver PRIVATE test_equivalence.cpp)
chdr_solver_test(equivalence "LevelEquivalence.*")
foreach(case ObliqueStandard ObliqueNonStandard AxialNonStandardStable
             AxialNonStandardBoundaryRecord)
  chdr_solver_test(gauge.${case} "GaugeWave.${case}")
endforeach()

target_sources(test-chdr-solver PRIVATE test_probe.cpp)
chdr_solver_test(probe "FourierProbe.*")
include(ProbeTests.cmake)
target_sources(test-chdr-solver PRIVATE test_comoving.cpp)
chdr_solver_test(comoving "ComovingT7.*")
foreach(variant serial mpi2)
  if(TEST chdr.solver.comoving.${variant})
    set_tests_properties(
      chdr.solver.comoving.${variant}
      PROPERTIES
        ENVIRONMENT
        "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;CHDR_SOLVER_RECORD_DIR=${CHDR_SOLVER_RECORD_ROOT}/comoving-${variant}"
    )
  endif()
endforeach()
include(StaggeredTests.cmake)
include(MediumTests.cmake)
target_sources(test-chdr-solver PRIVATE test_open_reflection.cpp test_open_field_reflection.cpp
                                        test_open_refinement.cpp)
chdr_solver_test(open.reflection "OpenReflection.ExactTwoWaveRows")
chdr_solver_test(open.field-reflection "OpenFieldReflection.SignedElectricAndMagnetic")
chdr_solver_test(open.refinement "OpenReflectionRefinement.CentrePlaneOrder")
foreach(variant serial mpi2)
  if(TEST chdr.solver.open.records.${variant})
    set_tests_properties(chdr.solver.open.records.${variant} PROPERTIES LABELS "chdr;solver;record")
  endif()
endforeach()
target_sources(test-chdr-solver PRIVATE test_open_shell.cpp test_open_medium.cpp
                                        test_open_probe.cpp)
set(randomPotentialHistoryFixtureDirectory
    "${CMAKE_CURRENT_BINARY_DIR}/random-potential-history-fixtures")
add_custom_command(
  OUTPUT "${randomPotentialHistoryFixtureDirectory}/medium.bin"
  COMMAND
    "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/GenerateRandomPotentialHistoryFixtures.py"
    "${randomPotentialHistoryFixtureDirectory}"
  DEPENDS GenerateRandomPotentialHistoryFixtures.py
  VERBATIM)
add_custom_target(chdr-random-potential-history-fixtures
                  DEPENDS "${randomPotentialHistoryFixtureDirectory}/medium.bin")
add_dependencies(test-chdr-solver chdr-random-potential-history-fixtures)
target_compile_definitions(
  test-chdr-solver
  PRIVATE CHDR_RANDOM_POTENTIAL_HISTORY_FIXTURE_DIR="${randomPotentialHistoryFixtureDirectory}")
chdr_solver_test(open.shell "OpenShell.*")
target_include_directories(
  test-chdr-solver PRIVATE "${PROJECT_SOURCE_DIR}/demos/chdr/boundaries/unit_tests/MaxwellSolvers")
target_sources(test-chdr-solver PRIVATE test_outer_layer.cpp)
chdr_solver_test(open.outer-layer "OpenOuterLayer.*")
