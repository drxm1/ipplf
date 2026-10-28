# Tests of the IPPL library changes made for ChDR, one directory per library component. Each
# executable is built by add_ippl_test(... COMPILE_ONLY) and registered by chdr_src_test under the
# ChDR naming scheme, so that ctest -R '^chdr\.' and ctest -L chdr select them with the rest.

# Register chdr.src.<name>.serial (the executable alone) and chdr.src.<name>.mpi2 (two ranks). ARGS
# go to both entries, SERIAL_ARGS and MPI_ARGS to one of them; TIMEOUT defaults to 60 s.
function(chdr_src_test name)
  cmake_parse_arguments(SRC "" "TIMEOUT" "ARGS;SERIAL_ARGS;MPI_ARGS" ${ARGN})
  if(NOT SRC_TIMEOUT)
    set(SRC_TIMEOUT 60)
  endif()
  add_test(NAME chdr.src.${name}.serial COMMAND $<TARGET_FILE:${name}> ${SRC_ARGS}
                                                ${SRC_SERIAL_ARGS})
  set_tests_properties(
    chdr.src.${name}.serial
    PROPERTIES TIMEOUT
               ${SRC_TIMEOUT}
               PROCESSORS
               1
               LABELS
               "chdr;src;unit"
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  if(MPIEXEC_EXECUTABLE)
    add_test(NAME chdr.src.${name}.mpi2
             COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                     $<TARGET_FILE:${name}> ${MPIEXEC_POSTFLAGS} ${SRC_ARGS} ${SRC_MPI_ARGS})
    set_tests_properties(
      chdr.src.${name}.mpi2
      PROPERTIES TIMEOUT
                 ${SRC_TIMEOUT}
                 PROCESSORS
                 2
                 LABELS
                 "chdr;src;unit;mpi"
                 ENVIRONMENT
                 "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  endif()
endfunction()

add_subdirectory(FieldLayout)
add_subdirectory(MaxwellSolvers)
