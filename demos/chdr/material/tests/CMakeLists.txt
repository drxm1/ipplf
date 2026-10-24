find_package(GTest REQUIRED)
include(AddIpplTest)
add_ippl_test(
  test-chdr-material
  COMPILE_ONLY
  SOURCES
  test_owner.cpp
  test_halo.cpp
  test_derived.cpp
  INCLUDE_DIRS
  "${CMAKE_CURRENT_SOURCE_DIR}/.."
  "${CMAKE_CURRENT_SOURCE_DIR}/../../solver"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../fel")
target_link_libraries(test-chdr-material PRIVATE chdr-material chdr-testing)
target_compile_options(test-chdr-material PRIVATE -UNDEBUG)
foreach(suite Owner Halo Derived)
  add_test(NAME chdr.material.${suite}.serial COMMAND test-chdr-material
                                                      --gtest_filter=Material${suite}.*)
  set_tests_properties(
    chdr.material.${suite}.serial PROPERTIES TIMEOUT 60 LABELS "chdr;material;unit" ENVIRONMENT
                                             "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  if(MPIEXEC_EXECUTABLE)
    add_test(
      NAME chdr.material.${suite}.mpi2
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test-chdr-material> ${MPIEXEC_POSTFLAGS} --gtest_filter=Material${suite}.*)
    set_tests_properties(
      chdr.material.${suite}.mpi2
      PROPERTIES TIMEOUT
                 60
                 PROCESSORS
                 2
                 LABELS
                 "chdr;material;unit;mpi"
                 ENVIRONMENT
                 "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
  endif()
endforeach()
