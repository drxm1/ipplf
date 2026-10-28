# The chdr-vacuum-wave executable. Its tests live in demos/chdr/tests/apps/vacuum-wave.
add_executable(chdr-vacuum-wave vacuum-wave.cpp VacuumSetup.cpp)
target_link_libraries(chdr-vacuum-wave PRIVATE chdr-core IPPL::ippl)
