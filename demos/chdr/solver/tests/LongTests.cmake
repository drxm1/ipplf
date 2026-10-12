# Full ladders and ten-period gauge histories; the independent anisotropic case remains disabled.
target_sources(test-chdr-solver PRIVATE test_gauge_growth.cpp test_convergence_long.cpp)
chdr_solver_test(long.GaugeGrowth.StaggeredRecord "LongGaugeGrowth.StaggeredRecord")
foreach(stencil IN ITEMS Standard NonStandard StaggeredRecord)
  foreach(family IN ITEMS Convergence GaugeGrowth)
    if(stencil STREQUAL "StaggeredRecord" AND family STREQUAL "Convergence")
      continue()
    endif()
    if(NOT stencil STREQUAL "StaggeredRecord")
      chdr_solver_test(long.${family}.${stencil} "Long${family}.${stencil}")
    endif()
    foreach(variant IN ITEMS serial mpi2)
      if(TEST chdr.solver.long.${family}.${stencil}.${variant})
        set_tests_properties(
          chdr.solver.long.${family}.${stencil}.${variant}
          PROPERTIES
            TIMEOUT
            1800
            LABELS
            "chdr;solver;long"
            ENVIRONMENT
            "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;CHDR_SOLVER_RECORD_DIR=${CHDR_SOLVER_RECORD_ROOT}/long-${variant}"
        )
        if(variant STREQUAL "mpi2")
          set_property(TEST chdr.solver.long.${family}.${stencil}.${variant} APPEND PROPERTY LABELS
                                                                                             mpi)
        endif()
      endif()
    endforeach()
  endforeach()
endforeach()
