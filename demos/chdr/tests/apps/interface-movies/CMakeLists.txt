include(AddIpplTest)
add_ippl_test(test-chdr-interface-movies COMPILE_ONLY SOURCES test_movie_start.cpp)
target_link_libraries(test-chdr-interface-movies PRIVATE chdr-core chdr-diagnostics chdr-testing)
target_compile_options(test-chdr-interface-movies PRIVATE -UNDEBUG)
add_test(NAME chdr.interface.movies.unit.serial COMMAND test-chdr-interface-movies)
set_tests_properties(
  chdr.interface.movies.unit.serial PROPERTIES TIMEOUT 60 LABELS "chdr;interface;unit" ENVIRONMENT
                                               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
# interface-movies runs all eight TE/TM cases for four steps on a small box.
add_test(
  NAME chdr.interface.movies.smoke
  COMMAND
    "${CMAKE_COMMAND}" "-DPROGRAM=$<TARGET_FILE:interface-movies>"
    "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/movie-smoke-frames" -P
    "${CMAKE_CURRENT_SOURCE_DIR}/movies-smoke.cmake")
set_tests_properties(
  chdr.interface.movies.smoke PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false"
                                         TIMEOUT 60 LABELS "chdr;interface;smoke")
if(MPIEXEC_EXECUTABLE)
  add_test(NAME chdr.interface.movies.unit.mpi2
           COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:test-chdr-interface-movies> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    chdr.interface.movies.unit.mpi2
    PROPERTIES TIMEOUT
               60
               PROCESSORS
               2
               LABELS
               "chdr;interface;unit;mpi"
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endif()
