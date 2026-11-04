include(AddIpplTest)
# Targets live with their sources; registration stays here to preserve test working directories.
add_subdirectory(support)
add_subdirectory(operators)
add_subdirectory(accuracy)
add_subdirectory(mur)
find_package(Python3 QUIET COMPONENTS Interpreter)
set(CHDR_SOLVER_RECORD_ROOT "${CMAKE_CURRENT_BINARY_DIR}/runs"
    CACHE PATH "Directory containing solver test records and rank fixtures")

function(chdr_solver_skip name reason)
  add_test(NAME chdr.solver.${name} COMMAND ${CMAKE_COMMAND} -E echo "SKIP: ${reason}")
  set_tests_properties(chdr.solver.${name} PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:" LABELS
                                                      "chdr;solver;unit;unavailable")
endfunction()

# Register both rank counts unless a suite requests only one of them.
function(chdr_solver_test testTarget suite filter)
  cmake_parse_arguments(CASE "SERIAL_ONLY;MPI2_ONLY" "" "" ${ARGN})
  if(NOT CASE_MPI2_ONLY)
    add_test(NAME chdr.solver.${suite}.serial COMMAND ${testTarget} --gtest_filter=${filter})
    set_tests_properties(
      chdr.solver.${suite}.serial PROPERTIES TIMEOUT 60 LABELS "chdr;solver;unit" ENVIRONMENT
                                             "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  endif()
  if(MPIEXEC_EXECUTABLE AND NOT CASE_SERIAL_ONLY)
    add_test(NAME chdr.solver.${suite}.mpi2
             COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                     $<TARGET_FILE:${testTarget}> ${MPIEXEC_POSTFLAGS} --gtest_filter=${filter})
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
chdr_solver_test(test-chdr-solver-operators contract "ChdrSolverContract.*")
chdr_solver_test(test-chdr-solver-operators compact "CompactOperators.*" SERIAL_ONLY)
chdr_solver_test(test-chdr-solver-operators staggered "ChdrStaggered.*")
chdr_solver_test(test-chdr-solver-operators setup "BoxSetup.*")
chdr_solver_test(test-chdr-solver-operators operators "Operators.*")
chdr_solver_test(test-chdr-solver-operators timestep "TimeStep.*")
foreach(case StandardBounded AnisotropicBounded NonStandardBounded)
  chdr_solver_test(test-chdr-solver-operators long.timestep.${case} "TimeStepLong.${case}")
  set_tests_properties(chdr.solver.long.timestep.${case}.serial PROPERTIES TIMEOUT 1800 LABELS
                                                                           "chdr;solver;long")
  if(MPIEXEC_EXECUTABLE)
    set_tests_properties(chdr.solver.long.timestep.${case}.mpi2 PROPERTIES TIMEOUT 1800 LABELS
                                                                           "chdr;solver;long;mpi")
  endif()
endforeach()

function(chdr_open_record name filter)
  add_test(NAME chdr.solver.open.record.${name} COMMAND test-chdr-solver-mur
                                                        --gtest_filter=${filter})
  set_tests_properties(
    chdr.solver.open.record.${name} PROPERTIES LABELS "chdr;solver;record" TIMEOUT 1800 ENVIRONMENT
                                               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endfunction()
chdr_open_record(gauge OpenGaugeRecords.VacuumRepresentatives)
chdr_open_record(zero-ghost MurOuterLayerRecords.ZeroGhostMagnitude)
chdr_open_record(all-cell-drift MurStabilityRecords.AllCellDrift)
chdr_solver_test(test-chdr-solver-support gather "RankGather.*")
chdr_solver_test(test-chdr-solver-accuracy manufactured "Manufactured.*")
chdr_solver_test(test-chdr-solver-accuracy rank_cases "RankCases.*")
include(support/RankTests.cmake)
include(accuracy/ConvergenceTests.cmake)
chdr_solver_test(test-chdr-solver-accuracy wave_norms "WaveNorms.*")
include(accuracy/LongTests.cmake)
include(accuracy/DispersionTests.cmake)
chdr_solver_test(test-chdr-solver-operators record.endpoint "TimeStepRecords.*")
chdr_solver_test(test-chdr-solver-support serial_reference "SerialReference.*")
chdr_solver_test(test-chdr-solver-support statistics "StableStatistics.*")
chdr_solver_test(test-chdr-solver-accuracy levels "Levels.*")
chdr_solver_test(test-chdr-solver-accuracy equivalence "LevelEquivalence.*")
foreach(case ObliqueStandard ObliqueNonStandard AxialNonStandardStable
             AxialNonStandardBoundaryRecord)
  chdr_solver_test(test-chdr-solver-accuracy gauge.${case} "GaugeCancellation.${case}")
endforeach()

chdr_solver_test(test-chdr-solver-accuracy comoving "Comoving.*")
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
include(accuracy/StaggeredTests.cmake)
include(accuracy/MediumTests.cmake)
chdr_solver_test(test-chdr-solver-mur mur.reflection "MurReflection.ExactTwoWaveRows")
foreach(variant serial mpi2)
  if(TEST chdr.solver.mur.reflection.${variant})
    set_tests_properties(chdr.solver.mur.reflection.${variant} PROPERTIES TIMEOUT 180)
  endif()
endforeach()
chdr_solver_test(test-chdr-solver-mur mur.field-reflection
                 "MurFieldReflection.SignedElectricAndMagnetic")
chdr_solver_test(test-chdr-solver-mur mur.refinement "MurReflectionRefinement.CentrePlaneOrder"
                 SERIAL_ONLY)
chdr_solver_test(test-chdr-solver-mur mur.records "MurReflectionRefinement.CodeConfirmedRows"
                 SERIAL_ONLY)
set_tests_properties(chdr.solver.mur.refinement.serial chdr.solver.mur.records.serial
                     PROPERTIES LABELS "chdr;solver;record")
chdr_solver_test(test-chdr-solver-mur mur.anisotropic "MurReflection.AnisotropicMedium" MPI2_ONLY)
set(randomPotentialHistoryFixtureDirectory
    "${CMAKE_CURRENT_BINARY_DIR}/random-potential-history-fixtures")
add_custom_command(
  OUTPUT "${randomPotentialHistoryFixtureDirectory}/medium.bin"
         "${randomPotentialHistoryFixtureDirectory}/stability-open.bin"
         "${randomPotentialHistoryFixtureDirectory}/stability-mixed.bin"
  COMMAND
    "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/support/GenerateRandomPotentialHistoryFixtures.py"
    "${randomPotentialHistoryFixtureDirectory}"
  DEPENDS support/GenerateRandomPotentialHistoryFixtures.py
  VERBATIM)
add_custom_target(
  chdr-random-potential-history-fixtures
  DEPENDS "${randomPotentialHistoryFixtureDirectory}/medium.bin"
          "${randomPotentialHistoryFixtureDirectory}/stability-open.bin"
          "${randomPotentialHistoryFixtureDirectory}/stability-mixed.bin")
add_dependencies(test-chdr-solver-mur chdr-random-potential-history-fixtures)
target_compile_definitions(
  test-chdr-solver-mur
  PRIVATE CHDR_RANDOM_POTENTIAL_HISTORY_FIXTURE_DIR="${randomPotentialHistoryFixtureDirectory}")
chdr_solver_test(test-chdr-solver-mur mur.faces "MurFaces.*")
chdr_solver_test(test-chdr-solver-mur mur.medium "MurMedium.*")
chdr_solver_test(test-chdr-solver-mur mur.outer-layer "MurOuterLayer.*")
chdr_solver_test(test-chdr-solver-mur mur.families "MurEdgeCornerFamilies.*")
chdr_solver_test(test-chdr-solver-mur mur.stability "MurStability.InnerDataEnvelope")
chdr_solver_test(test-chdr-solver-mur mur.stability-control "MurStability.MismatchedSpeedControl"
                 SERIAL_ONLY)
foreach(suite stability stability-control)
  foreach(variant serial mpi2)
    if(TEST chdr.solver.mur.${suite}.${variant})
      set_tests_properties(chdr.solver.mur.${suite}.${variant} PROPERTIES TIMEOUT 1800 LABELS
                                                                          "chdr;solver;long")
      if(variant STREQUAL "mpi2")
        set_property(TEST chdr.solver.mur.${suite}.${variant} APPEND PROPERTY LABELS mpi)
      endif()
    endif()
  endforeach()
endforeach()
chdr_solver_test(test-chdr-solver-mur open.beam-ranks "OpenBeam.ReducedRankWitness" MPI2_ONLY)
chdr_solver_test(test-chdr-solver-mur open.beam-full "OpenBeam.FullBoxAgainstPeriodicWitness"
                 MPI2_ONLY)
foreach(suite beam-ranks beam-full)
  foreach(variant serial mpi2)
    if(TEST chdr.solver.open.${suite}.${variant})
      set_tests_properties(chdr.solver.open.${suite}.${variant} PROPERTIES TIMEOUT 1800 LABELS
                                                                           "chdr;solver;long")
      if(variant STREQUAL "mpi2")
        set_property(TEST chdr.solver.open.${suite}.${variant} APPEND PROPERTY LABELS mpi)
      endif()
    endif()
  endforeach()
endforeach()
