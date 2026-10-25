add_library(chdr-material INTERFACE)
# The demos/chdr root: material headers are included as "core/material/<File>.h".
target_include_directories(chdr-material INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/../..")
target_link_libraries(chdr-material INTERFACE IPPL::ippl)
