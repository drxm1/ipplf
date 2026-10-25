# The ChDR core: field types, material owner, operators, formulations, the potential solver, the
# simulation manager and the generic writer. Nothing experiment-, test- or opt-in-related lives
# here. Consumers include core headers by their path from demos/chdr, e.g.
# "core/solver/ChdrSolver.h".
add_subdirectory(material)

add_library(chdr-core STATIC base/FieldContainer.cpp solver/BoxSetup.cpp)
# Only the demos/chdr root is exported: every core header is included by its path, and the compiler
# rejects a bare include.
target_include_directories(chdr-core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(chdr-core PUBLIC ippl chdr-material chdr-fel-headers)
