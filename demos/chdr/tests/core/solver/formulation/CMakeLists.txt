include(AddIpplTest)
add_ippl_test(test-chdr-formulation COMPILE_ONLY SOURCES test_o1vs_regression.cpp)
target_link_libraries(test-chdr-formulation PRIVATE chdr-core chdr-diagnostics chdr-testing)
target_compile_options(test-chdr-formulation PRIVATE -UNDEBUG)
add_test(NAME chdr.solver.formulation.o1vs.serial COMMAND test-chdr-formulation)
set_tests_properties(
  chdr.solver.formulation.o1vs.serial
  PROPERTIES TIMEOUT 60 LABELS "chdr;solver;formulation;unit" ENVIRONMENT
             "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
if(MPIEXEC_EXECUTABLE)
  add_test(NAME chdr.solver.formulation.o1vs.mpi2
           COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:test-chdr-formulation> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    chdr.solver.formulation.o1vs.mpi2
    PROPERTIES TIMEOUT
               60
               PROCESSORS
               2
               LABELS
               "chdr;solver;formulation;unit;mpi"
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endif()
