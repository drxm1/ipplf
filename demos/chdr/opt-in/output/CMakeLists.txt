# Frames plus energy CSV for one kind of run (O1-VS, staggered R2); header-only, built on the
# generic writer of the core and the opt-in diagnostics.
add_library(chdr-run-output INTERFACE)
target_link_libraries(chdr-run-output INTERFACE chdr-diagnostics)
