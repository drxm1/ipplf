# Single-root phase checks and the full direction/safety-factor records.
target_sources(test-chdr-solver PRIVATE test_dispersion.cpp)
chdr_solver_test(dispersion "Dispersion.*")
chdr_solver_rank_producer(dispersion)
chdr_solver_test(record.dispersion "RecordDispersion.*")
foreach(variant IN ITEMS serial mpi2)
  if(TEST chdr.solver.record.dispersion.${variant})
    set_tests_properties(
      chdr.solver.record.dispersion.${variant}
      PROPERTIES
        TIMEOUT
        600
        LABELS
        "chdr;solver;record"
        ENVIRONMENT
        "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;CHDR_SOLVER_RECORD_DIR=${CHDR_SOLVER_RECORD_ROOT}/record-${variant}"
    )
    if(variant STREQUAL "mpi2")
      set_property(TEST chdr.solver.record.dispersion.${variant} APPEND PROPERTY LABELS mpi)
    endif()
  endif()
endforeach()
