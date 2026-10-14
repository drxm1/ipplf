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
foreach(suite Validity Guard Equivalence Manufactured Dispersion TimeStep)
  chdr_solver_medium_test(${suite} "Medium${suite}.*")
endforeach()
chdr_solver_medium_test(long.Dispersion "MediumDispersionLong.*")
chdr_solver_medium_test(record.Dispersion "MediumDispersionRecord.*")
chdr_solver_medium_test(long.TimeStep "MediumTimeStepLong.*")
