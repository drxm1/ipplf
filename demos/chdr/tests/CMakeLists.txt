# Tests of the ChDR code, one directory per unit of the code they test.
add_subdirectory(src) #< src/ holds the tests of the IPPL library changes made for ChDR (per-axis
                      #< periodic layouts, Mur faces).
add_subdirectory(core/material)
add_subdirectory(core/solver)
add_subdirectory(core/formulation)
add_subdirectory(core/simulation)
add_subdirectory(opt-in/diagnostics)
add_subdirectory(apps/vacuum-wave)
add_subdirectory(apps/planar-interface)
add_subdirectory(apps/interface-movies)
add_subdirectory(1d)
