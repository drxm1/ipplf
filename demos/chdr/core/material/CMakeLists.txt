add_library(chdr-material INTERFACE)
target_include_directories(chdr-material INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(chdr-material INTERFACE IPPL::ippl)
