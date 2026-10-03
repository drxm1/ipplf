add_executable(chdr-vacuum-wave vacuum-wave.cpp VacuumSetup.cpp)
target_link_libraries(chdr-vacuum-wave PRIVATE IPPL::ippl)
target_include_directories(chdr-vacuum-wave PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../fel")

if(BUILD_TESTING AND IPPL_ENABLE_UNIT_TESTS)
  add_subdirectory(tests)
endif()
