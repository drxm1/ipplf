# The ChDR core: field types, material owner, operators, formulations, the potential solver, the
# simulation manager and the generic writer. Nothing experiment-, test- or opt-in-related lives
# here. Consumers include core headers by their path from demos/chdr, e.g.
# "core/solver/ChdrSolver.h".
add_subdirectory(material)

add_library(chdr-core STATIC base/FieldContainer.cpp solver/BoxSetup.cpp)
# The demos/chdr root carries the qualified includes. The unit directories stay exported until the
# include conversion commit; they let the bare includes of the old layout resolve.
target_include_directories(
  chdr-core
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/.."
         "${CMAKE_CURRENT_SOURCE_DIR}/base"
         "${CMAKE_CURRENT_SOURCE_DIR}/operators"
         "${CMAKE_CURRENT_SOURCE_DIR}/formulation"
         "${CMAKE_CURRENT_SOURCE_DIR}/solver"
         "${CMAKE_CURRENT_SOURCE_DIR}/simulation"
         "${CMAKE_CURRENT_SOURCE_DIR}/output")
target_link_libraries(chdr-core PUBLIC ippl chdr-material chdr-fel-headers)
