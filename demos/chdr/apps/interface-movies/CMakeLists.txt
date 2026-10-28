# The open-box dielectric interface movies. Tests live in demos/chdr/tests/apps/interface-movies.
add_executable(interface-movies interface-movies.cpp)
target_link_libraries(interface-movies PRIVATE chdr-core chdr-diagnostics chdr-run-output)
