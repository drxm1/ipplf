# ChDR geometry and vacuum-wave demos.
message(STATUS "Configuring demos/chdr/")

# Shared test support: every ChDR test executable links chdr-testing, which compiles the common
# GoogleTest main into it and provides the utilities under testing/ as "testing/<file>".
add_library(chdr-testing INTERFACE)
target_sources(chdr-testing INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/testing/TestMain.cpp")
target_include_directories(chdr-testing INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

# The field container and field types come from the FEL demo; one target names that dependency.
add_library(chdr-fel-headers INTERFACE)
target_include_directories(chdr-fel-headers INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/../fel")

add_subdirectory(material)
add_subdirectory(mesh)
add_subdirectory(vacuum)
add_subdirectory(solver)
add_subdirectory(opt-in/diagnostics)
add_subdirectory(opt-in/output)
add_subdirectory(interface)

if(BUILD_TESTING AND IPPL_ENABLE_UNIT_TESTS)
  add_subdirectory(tests)
endif()
