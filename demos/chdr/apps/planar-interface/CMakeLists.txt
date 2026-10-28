# One planar dielectric in a periodic box: the plane-pulse reflection measurement and the packet
# cases. Tests live in demos/chdr/tests/apps/planar-interface.
add_executable(planar-interface planar-interface.cpp)
target_link_libraries(planar-interface PRIVATE chdr-core chdr-diagnostics chdr-run-output)
