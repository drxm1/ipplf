# ChDR geometry and vacuum-wave demos.
message(STATUS "Configuring demos/chdr/")

# One GoogleTest main for every ChDR test executable: IPPL lifetime and rank-global exit status.
set(CHDR_TEST_MAIN "${CMAKE_CURRENT_SOURCE_DIR}/testing/TestMain.cpp")

add_subdirectory(material)
add_subdirectory(mesh)
add_subdirectory(vacuum)
add_subdirectory(solver)
add_subdirectory(interface)
