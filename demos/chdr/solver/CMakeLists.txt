# Periodic potential solver and audited reconstruction utilities.
add_library(chdr-solver STATIC PeriodicSetup.cpp ../vacuum/VacuumSetup.cpp)
target_include_directories(
  chdr-solver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/../vacuum"
                     "${CMAKE_CURRENT_SOURCE_DIR}/../../fel")
target_link_libraries(chdr-solver PUBLIC ippl)
if(BUILD_TESTING AND IPPL_ENABLE_UNIT_TESTS)
  add_subdirectory(tests)
endif()
