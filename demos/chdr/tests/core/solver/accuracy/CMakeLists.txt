add_ippl_test(
  test-chdr-solver-accuracy
  COMPILE_ONLY
  SOURCES
  test_manufactured.cpp
  test_convergence.cpp
  test_convergence_long.cpp
  test_staggered_convergence.cpp
  test_dispersion.cpp
  test_levels.cpp
  test_equivalence.cpp
  test_gauge_growth.cpp
  test_ranks.cpp
  test_comoving.cpp
  test_medium_guard.cpp)
target_link_libraries(test-chdr-solver-accuracy PRIVATE chdr-core chdr-diagnostics chdr-testing
                                                        nlohmann_json::nlohmann_json)
target_compile_options(test-chdr-solver-accuracy PRIVATE -UNDEBUG)
