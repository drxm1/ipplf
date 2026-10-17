include(AddIpplTest)
add_ippl_test(
  test-chdr-interface-unit
  COMPILE_ONLY
  SOURCES
  ../../vacuum/tests/TestMain.cpp
  test_o1vs_regression.cpp
  INCLUDE_DIRS
  "${CMAKE_CURRENT_SOURCE_DIR}/../../solver/tests")
target_link_libraries(test-chdr-interface-unit PRIVATE chdr-solver nlohmann_json::nlohmann_json)
target_compile_options(test-chdr-interface-unit PRIVATE -UNDEBUG)
foreach(suite O1VS)
  add_test(NAME chdr.interface.unit.${suite}.serial COMMAND test-chdr-interface-unit
                                                            --gtest_filter=${suite}.*)
  set_tests_properties(
    chdr.interface.unit.${suite}.serial
    PROPERTIES TIMEOUT 60 LABELS "chdr;interface;unit" ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endforeach()
