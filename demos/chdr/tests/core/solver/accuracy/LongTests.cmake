# Full ladders, anisotropic cases and ten-period gauge histories.
chdr_solver_test(test-chdr-solver-accuracy long.GaugeGrowth.StaggeredRecord
                 "LongGaugeGrowth.StaggeredRecord")
foreach(stencil IN ITEMS Standard NonStandard StaggeredRecord)
  foreach(family IN ITEMS Convergence GaugeGrowth)
    set(costLabel long)
    if(stencil STREQUAL "StaggeredRecord")
      set(costLabel record)
    endif()
    if(stencil STREQUAL "StaggeredRecord" AND family STREQUAL "Convergence")
      continue()
    endif()
    if(NOT stencil STREQUAL "StaggeredRecord")
      chdr_solver_test(test-chdr-solver-accuracy long.${family}.${stencil}
                       "Long${family}.${stencil}")
    endif()
    foreach(variant IN ITEMS serial mpi2)
      if(TEST chdr.solver.long.${family}.${stencil}.${variant})
        set_tests_properties(
          chdr.solver.long.${family}.${stencil}.${variant}
          PROPERTIES
            TIMEOUT
            1800
            LABELS
            "chdr;solver;${costLabel}"
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
foreach(reading Collocated StaggeredRecord)
  chdr_solver_test(test-chdr-solver-accuracy long.anisotropic.${reading}
                   "AnisotropicConvergence.${reading}" MPI2_ONLY)
  if(TEST chdr.solver.long.anisotropic.${reading}.mpi2)
    set_tests_properties(
      chdr.solver.long.anisotropic.${reading}.mpi2
      PROPERTIES
        TIMEOUT
        1800
        LABELS
        "chdr;solver;long;mpi"
        ENVIRONMENT
        "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;CHDR_SOLVER_RECORD_DIR=${CHDR_SOLVER_RECORD_ROOT}/anisotropic-${reading}-mpi2"
    )
  endif()
endforeach()
