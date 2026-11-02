include(AddIpplTest)
add_ippl_test(test-chdr-input COMPILE_ONLY SOURCES test_input.cpp)
target_link_libraries(test-chdr-input PRIVATE chdr-application chdr-testing)
target_compile_options(test-chdr-input PRIVATE -UNDEBUG)
add_test(NAME chdr.input.unit.serial COMMAND test-chdr-input)
set_tests_properties(
  chdr.input.unit.serial PROPERTIES TIMEOUT 30 LABELS "chdr;input;unit" ENVIRONMENT
                                    "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
