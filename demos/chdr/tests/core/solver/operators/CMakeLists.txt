add_ippl_test(
  test-chdr-solver-operators
  COMPILE_ONLY
  SOURCES
  test_solver_contract.cpp
  test_compact.cpp
  test_staggered.cpp
  test_setup.cpp
  test_operators.cpp
  test_timestep.cpp
  test_endpoint_records.cpp)
target_link_libraries(test-chdr-solver-operators PRIVATE chdr-core chdr-diagnostics chdr-testing
                                                         nlohmann_json::nlohmann_json)
target_compile_options(test-chdr-solver-operators PRIVATE -UNDEBUG)
