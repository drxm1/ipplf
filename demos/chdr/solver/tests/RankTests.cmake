# Comparison fixture: every producer completes before the record reader.
if(Python3_Interpreter_FOUND)
  add_test(NAME chdr.solver.ranks.prepare
           COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/prepare_rank_records.py"
                   "${CHDR_SOLVER_RECORD_ROOT}")
  set_tests_properties(
    chdr.solver.ranks.prepare PROPERTIES FIXTURES_SETUP chdr_solver_rank_generation TIMEOUT 10
                                         LABELS "chdr;solver;unit;ranks")
else()
  chdr_solver_skip(ranks.prepare "Rank record generation requires Python3")
endif()

function(chdr_solver_rank_producer suite)
  if(NOT Python3_Interpreter_FOUND)
    return()
  endif()
  foreach(variant IN ITEMS serial mpi2)
    if(TEST chdr.solver.${suite}.${variant})
      set_tests_properties(
        chdr.solver.${suite}.${variant}
        PROPERTIES
          ENVIRONMENT
          "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;CHDR_SOLVER_RECORD_DIR=${CHDR_SOLVER_RECORD_ROOT}/${variant}"
          FIXTURES_SETUP
          chdr_solver_ranks
          FIXTURES_REQUIRED
          chdr_solver_rank_generation)
    endif()
  endforeach()
endfunction()

chdr_solver_rank_producer(manufactured)
