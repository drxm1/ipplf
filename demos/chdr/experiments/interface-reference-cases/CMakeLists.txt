# This reference calculator uses only standard C++ and the independent host formulas.
add_executable(interface-reference reference.cpp)
target_compile_features(interface-reference PRIVATE cxx_std_20)
target_include_directories(
  interface-reference PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../.."
                              "${CMAKE_CURRENT_BINARY_DIR}/../../core/generated")
add_dependencies(interface-reference chdr-source-identity)

if(BUILD_TESTING AND IPPL_ENABLE_UNIT_TESTS)
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    add_test(NAME chdr.interface.reference-tools
             COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_tools.py"
                     --reference "$<TARGET_FILE:interface-reference>")
    set_tests_properties(chdr.interface.reference-tools PROPERTIES TIMEOUT 60 LABELS
                                                                   "chdr;interface;analysis")
  endif()
endif()
