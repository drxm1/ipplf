add_ippl_test(OpenFaceClosure NUM_PROCS 1 PROPERTIES PROCESSORS 1)
add_test(NAME OpenFaceClosureMPI COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:OpenFaceClosure>)
set_tests_properties(
  OpenFaceClosure OpenFaceClosureMPI PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false"
                                                LABELS "unit;maxwell")
set_tests_properties(OpenFaceClosureMPI PROPERTIES PROCESSORS 2 TIMEOUT 60)
