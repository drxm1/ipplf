include(AddIpplTest)
add_ippl_test(test-plane-interface COMPILE_ONLY SOURCES test_plane_interface.cpp)
target_link_libraries(test-plane-interface PRIVATE chdr-core chdr-diagnostics chdr-testing)
add_test(NAME chdr.interface.normal COMMAND test-plane-interface
                                            --gtest_filter=PlaneInterface.Normal)
add_test(NAME chdr.interface.long.brewster COMMAND test-plane-interface
                                                   --gtest_filter=PlaneInterface.Brewster)
set_tests_properties(
  chdr.interface.normal chdr.interface.long.brewster
  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false" TIMEOUT 1800 LABELS
             "chdr;interface")
set_tests_properties(chdr.interface.long.brewster PROPERTIES LABELS "chdr;interface;long")
# planar-interface runs its half-space case for four steps on a small box.
add_test(
  NAME chdr.interface.smoke
  COMMAND
    "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:planar-interface>"
    "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/smoke-frames" -P
    "${CMAKE_CURRENT_SOURCE_DIR}/smoke.cmake")
set_tests_properties(
  chdr.interface.smoke PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false" TIMEOUT 60
                                  LABELS "chdr;interface;smoke")
if(MPIEXEC_EXECUTABLE)
  set(interface_mpi "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                    $<TARGET_FILE:test-plane-interface> ${MPIEXEC_POSTFLAGS})
  add_test(NAME chdr.interface.normal.mpi2 COMMAND ${interface_mpi}
                                                   --gtest_filter=PlaneInterface.NormalSeam)
  add_test(NAME chdr.interface.reserved COMMAND ${interface_mpi} --gtest_also_run_disabled_tests
                                                --gtest_filter=PlaneInterface.DISABLED_Reserved)
  set_tests_properties(chdr.interface.normal.mpi2 chdr.interface.reserved
                       PROPERTIES PROCESSORS 2 ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  set_tests_properties(chdr.interface.normal.mpi2 PROPERTIES TIMEOUT 60 LABELS "chdr;interface;mpi")
  set_tests_properties(chdr.interface.reserved PROPERTIES DISABLED TRUE TIMEOUT 1800 LABELS
                                                          "chdr;interface;reserved;long")
endif()
