# Opt-in measurements of a running field problem; header-only, built on the solver core.
add_library(chdr-diagnostics INTERFACE)
target_link_libraries(chdr-diagnostics INTERFACE chdr-solver)
