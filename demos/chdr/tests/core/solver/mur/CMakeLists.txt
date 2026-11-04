add_ippl_test(
  test-chdr-solver-mur
  COMPILE_ONLY
  SOURCES
  test_mur_faces.cpp
  test_mur_medium.cpp
  test_mur_outer_layer.cpp
  test_mur_stability.cpp
  test_mur_families.cpp
  test_mur_reflection.cpp
  test_mur_field_reflection.cpp
  test_mur_refinement.cpp
  test_open_beam.cpp
  test_open_gauge_records.cpp)
target_link_libraries(test-chdr-solver-mur PRIVATE chdr-core chdr-diagnostics chdr-testing
                                                   nlohmann_json::nlohmann_json)
target_compile_options(test-chdr-solver-mur PRIVATE -UNDEBUG)
