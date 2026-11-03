# Exercise the face/edge reconstruction with the same analytic waves and error checks as the
# collocated reading; keep its long ladder out of the quick suite.
chdr_solver_test(test-chdr-solver-accuracy convergence.staggered "StaggeredConvergence.*")
chdr_solver_rank_producer(convergence.staggered)
chdr_solver_test(test-chdr-solver-accuracy long.staggered "LongStaggeredConvergence.*")
foreach(variant IN ITEMS serial mpi2)
  if(TEST chdr.solver.long.staggered.${variant})
    set_tests_properties(
      chdr.solver.long.staggered.${variant}
      PROPERTIES
        TIMEOUT
        1800
        LABELS
        "chdr;solver;long"
        ENVIRONMENT
        "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;CHDR_SOLVER_RECORD_DIR=${CHDR_SOLVER_RECORD_ROOT}/long-${variant}"
    )
    if(variant STREQUAL "mpi2")
      set_property(TEST chdr.solver.long.staggered.${variant} APPEND PROPERTY LABELS mpi)
    endif()
  endif()
endforeach()
