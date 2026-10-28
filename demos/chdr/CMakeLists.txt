# ChDR: the potential solver core, the opt-in diagnostics, the programs under apps/, the mesh
# geometry for meep comparisons and the tests. README.md maps the directory.
message(STATUS "Configuring demos/chdr/")

# Shared test support: every ChDR test executable links chdr-testing, which compiles the common
# GoogleTest main into it and provides the utilities under testing/ as "testing/<file>".
add_library(chdr-testing INTERFACE)
target_sources(chdr-testing INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/testing/TestMain.cpp")
target_include_directories(chdr-testing INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

# The field container and field types come from the FEL demo; one target names that dependency.
add_library(chdr-fel-headers INTERFACE)
target_include_directories(chdr-fel-headers INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/../fel")

add_subdirectory(core)
add_subdirectory(mesh)
add_subdirectory(opt-in/diagnostics)
add_subdirectory(opt-in/output)
add_subdirectory(apps/vacuum-wave)
add_subdirectory(apps/planar-interface)
add_subdirectory(apps/interface-movies)

if(BUILD_TESTING AND IPPL_ENABLE_UNIT_TESTS)
  add_subdirectory(tests)
endif()
