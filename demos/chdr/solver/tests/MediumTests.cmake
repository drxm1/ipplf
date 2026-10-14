# Development cases use explicit filters; the independently held cases are never registered.
target_sources(test-chdr-solver PRIVATE test_medium_guard.cpp)
function(chdr_solver_medium_test suite filter)
  chdr_solver_test(medium.${suite} "${filter}")
  foreach(variant serial mpi2)
    if(TEST chdr.solver.medium.${suite}.${variant})
      set_tests_properties(
        chdr.solver.medium.${suite}.${variant}
        PROPERTIES
          ENVIRONMENT
          "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;CHDR_SOLVER_RECORD_DIR=${CHDR_SOLVER_RECORD_ROOT}/medium-${suite}-${variant}"
      )
      if(suite MATCHES "^(long|record)[.]")
        set_tests_properties(chdr.solver.medium.${suite}.${variant} PROPERTIES TIMEOUT 1800 LABELS
                                                                               "chdr;solver;long")
        if(variant STREQUAL "mpi2")
          set_property(TEST chdr.solver.medium.${suite}.${variant} APPEND PROPERTY LABELS mpi)
        endif()
      endif()
    endif()
  endforeach()
endfunction()
foreach(
  suite
  Validity
  Guard
  Equivalence
  Manufactured
  TimeStep
  Dispersion
  Source)
  chdr_solver_medium_test(${suite} "Medium${suite}.*")
endforeach()
foreach(reading Collocated Staggered)
  chdr_solver_medium_test(Convergence.${reading} "MediumConvergence.${reading}")
endforeach()
chdr_solver_medium_test(long.Dispersion "MediumDispersionLong.*")
chdr_solver_medium_test(record.Dispersion "MediumDispersionRecord.*")
chdr_solver_medium_test(long.TimeStep "MediumTimeStepLong.*")
foreach(reading Collocated Staggered)
  chdr_solver_medium_test(long.Convergence.${reading} "MediumConvergenceLong.${reading}")
  chdr_solver_medium_test(long.Source.${reading} "MediumSourceLong.${reading}")
endforeach()

# The existing reader compares manufactured histories and measured mode phases. Continuum cases
# perform their communicator-self comparison inside the executable. The generation fixture clears
# old outputs before either rank count runs. Both producer suites must finish before comparison,
# which checks the complete expected file set as well as values.
if(Python3_Interpreter_FOUND AND MPIEXEC_EXECUTABLE)
  set(mediumRoot "${CHDR_SOLVER_RECORD_ROOT}/medium-ranks")
  add_test(NAME chdr.solver.medium.ranks.prepare
           COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/prepare_rank_records.py"
                   "${mediumRoot}")
  set_tests_properties(
    chdr.solver.medium.ranks.prepare PROPERTIES FIXTURES_SETUP chdr_medium_generation TIMEOUT 10
                                                LABELS "chdr;solver;unit;ranks")
  foreach(suite Manufactured Dispersion)
    foreach(variant serial mpi2)
      set_tests_properties(
        chdr.solver.medium.${suite}.${variant}
        PROPERTIES
          FIXTURES_SETUP chdr_medium_ranks FIXTURES_REQUIRED chdr_medium_generation ENVIRONMENT
          "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;CHDR_SOLVER_RECORD_DIR=${mediumRoot}/${variant}")
    endforeach()
  endforeach()
  add_test(NAME chdr.solver.medium.ranks
           COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/compare_ranks.py"
                   "${mediumRoot}/serial" "${mediumRoot}/mpi2" --medium)
  set_tests_properties(
    chdr.solver.medium.ranks PROPERTIES FIXTURES_REQUIRED chdr_medium_ranks TIMEOUT 60 LABELS
                                        "chdr;solver;unit;ranks;mpi")
endif()
