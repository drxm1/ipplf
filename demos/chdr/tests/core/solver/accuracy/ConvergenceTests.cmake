# Bounded continuum ladders with records for decomposition comparisons.
foreach(stencil IN ITEMS Standard NonStandard)
  foreach(family Mixed Sourced)
    foreach(reconstruction R1 R2 R2Prime)
      set(suite convergence.${stencil}.${family}.${reconstruction})
      chdr_solver_test(test-chdr-solver-accuracy ${suite}
                       "Convergence.${stencil}${family}${reconstruction}")
      chdr_solver_rank_producer(${suite})
    endforeach()
  endforeach()
endforeach()
