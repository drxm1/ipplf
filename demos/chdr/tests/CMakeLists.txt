# Tests of the ChDR code, one directory per unit of the code they test.
#
# src/ holds the tests of the IPPL library changes made for ChDR (per-axis periodic layouts, Mur
# faces). They are written with IPPL's own test support (add_ippl_test and TestUtils.h from
# unit_tests/), which IPPL provides only when it is the top-level project, so they are entered only
# then, as before their move.
if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(src)
endif()
add_subdirectory(core/material)
add_subdirectory(core/solver)
add_subdirectory(core/solver/formulation)
add_subdirectory(core/simulation)
add_subdirectory(opt-in/diagnostics)
add_subdirectory(apps/vacuum-wave)
add_subdirectory(apps/planar-interface)
add_subdirectory(apps/interface-movies)
add_subdirectory(1d)
