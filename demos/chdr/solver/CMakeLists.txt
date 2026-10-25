# Periodic potential solver and audited reconstruction utilities.
add_library(chdr-solver STATIC BoxSetup.cpp FieldContainer.cpp)
# The demos/chdr root is exported so that consumers include shared headers by their path.
target_include_directories(chdr-solver PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
                                              "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(chdr-solver PUBLIC ippl chdr-material chdr-fel-headers)
