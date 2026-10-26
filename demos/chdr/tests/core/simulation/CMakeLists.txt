include(AddIpplTest)
add_ippl_test(test-chdr-simulation COMPILE_ONLY SOURCES test_simulation.cpp)
target_link_libraries(test-chdr-simulation PRIVATE chdr-core chdr-testing)
target_compile_options(test-chdr-simulation PRIVATE -UNDEBUG)
add_test(NAME chdr.simulation.unit.serial COMMAND test-chdr-simulation)
set_tests_properties(
  chdr.simulation.unit.serial PROPERTIES TIMEOUT 60 LABELS "chdr;simulation;unit" ENVIRONMENT
                                         "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
if(MPIEXEC_EXECUTABLE)
  add_test(NAME chdr.simulation.unit.mpi2
           COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:test-chdr-simulation> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    chdr.simulation.unit.mpi2
    PROPERTIES TIMEOUT
               60
               PROCESSORS
               2
               LABELS
               "chdr;simulation;unit;mpi"
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endif()
