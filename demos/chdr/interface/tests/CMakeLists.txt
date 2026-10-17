include(AddIpplTest)
add_ippl_test(test-plane-interface COMPILE_ONLY SOURCES ../../vacuum/tests/TestMain.cpp
              test_plane_interface.cpp)
target_link_libraries(test-plane-interface PRIVATE chdr-solver)
add_test(NAME chdr.interface.normal COMMAND test-plane-interface
                                            --gtest_filter=PlaneInterface.Normal)
add_test(NAME chdr.interface.long.brewster COMMAND test-plane-interface
                                                   --gtest_filter=PlaneInterface.Brewster)
set_tests_properties(
  chdr.interface.normal chdr.interface.long.brewster
  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false" TIMEOUT 1800 LABELS
             "chdr;interface")
set_tests_properties(chdr.interface.long.brewster PROPERTIES LABELS "chdr;interface;long")
add_test(
  NAME chdr.interface.unit.smoke
  COMMAND
    "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:interface-demo>"
    "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/smoke-frames" -P
    "${CMAKE_CURRENT_SOURCE_DIR}/smoke.cmake")
set_tests_properties(
  chdr.interface.unit.smoke PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false" TIMEOUT
                                       60 LABELS "chdr;interface;unit")
if(MPIEXEC_EXECUTABLE)
  set(interface_mpi "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                    $<TARGET_FILE:test-plane-interface> ${MPIEXEC_POSTFLAGS})
  add_test(NAME chdr.interface.normal.mpi2 COMMAND ${interface_mpi}
                                                   --gtest_filter=PlaneInterface.NormalSeam)
  set_tests_properties(
    chdr.interface.normal.mpi2
    PROPERTIES PROCESSORS
               2
               TIMEOUT
               60
               LABELS
               "chdr;interface;mpi"
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endif()
