add_executable(chdr-vacuum-wave vacuum-wave.cpp VacuumSetup.cpp)
target_link_libraries(chdr-vacuum-wave PRIVATE chdr-solver IPPL::ippl)

if(BUILD_TESTING AND IPPL_ENABLE_UNIT_TESTS)
  add_subdirectory(tests)
endif()
