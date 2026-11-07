# The ChDR core: the base types shared by everything, the potential solver with its operators and
# formulations, the material owner, the simulation manager and the generic writer. Nothing
# experiment-, test- or opt-in-related lives here. Consumers include core headers by their path from
# demos/chdr, e.g. "core/solver/ChdrSolver.h".
add_subdirectory(material)

add_library(chdr-core STATIC base/FieldContainer.cpp solver/BoxSetup.cpp)
# Only the demos/chdr root is exported: every core header is included by its path, and the compiler
# rejects a bare include.
target_include_directories(chdr-core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(chdr-core PUBLIC ippl chdr-material chdr-fel-headers)

# Common input and run records remain optional for low-level solver/unit-test consumers.
add_library(chdr-application STATIC input/Application.cpp input/RunInput.cpp output/RunRecord.cpp)
target_link_libraries(chdr-application PUBLIC chdr-core catalyst::catalyst
                      PRIVATE nlohmann_json::nlohmann_json)
# Refresh provenance for each build; unchanged content keeps RunRecord.cpp up to date.
set(chdr_source_header "${CMAKE_CURRENT_BINARY_DIR}/generated/ChdrSourceIdentity.h")
add_custom_target(
  chdr-source-identity
  COMMAND
    "${CMAKE_COMMAND}" "-DCHDR_SOURCE_DIRECTORY=${CMAKE_CURRENT_SOURCE_DIR}"
    "-DCHDR_OUTPUT_HEADER=${chdr_source_header}" -P
    "${CMAKE_CURRENT_SOURCE_DIR}/output/WriteSourceIdentity.cmake"
  BYPRODUCTS "${chdr_source_header}"
  VERBATIM)
add_dependencies(chdr-application chdr-source-identity)
target_include_directories(chdr-application PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/generated")
target_compile_definitions(
  chdr-application PRIVATE CHDR_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
                           CHDR_COMPILER="${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")
