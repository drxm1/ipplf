add_ippl_test(test-chdr-solver-support COMPILE_ONLY SOURCES test_gather.cpp test_statistics.cpp
              test_serial_reference.cpp)
target_link_libraries(test-chdr-solver-support PRIVATE chdr-core chdr-diagnostics chdr-testing
                                                       nlohmann_json::nlohmann_json)
target_compile_options(test-chdr-solver-support PRIVATE -UNDEBUG)
