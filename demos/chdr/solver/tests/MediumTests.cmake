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
    endif()
  endforeach()
endfunction()
foreach(suite Validity Equivalence)
  chdr_solver_medium_test(${suite} "Medium${suite}.*")
endforeach()
