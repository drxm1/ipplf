# -----------------------------------------------------------------------------
# demos/CMakeLists.txt
#
# Top-level entry point for all demo applications bundled with IPPL.
# -----------------------------------------------------------------------------
message(STATUS "Configuring demos/")

# ------------------------------------------------------------------------------
# Include optional source directories for demos
# ------------------------------------------------------------------------------

if(IPPL_ENABLE_ALPINE)
  add_subdirectory(alpine)
endif()

if(IPPL_ENABLE_ESPIF AND IPPL_ENABLE_FFT)
  add_subdirectory(electrostaticPIF)
endif()

if(IPPL_ENABLE_COSMOLOGY)
  add_subdirectory(cosmology)
endif()

if(IPPL_ENABLE_FEL)
  add_subdirectory(fel)
endif()

if(IPPL_ENABLE_COLLISIONS)
  add_subdirectory(collisions)
endif()

# Boundary regressions also build on the standalone library branches.
if(PROJECT_IS_TOP_LEVEL AND BUILD_TESTING AND IPPL_ENABLE_UNIT_TESTS)
  add_subdirectory(chdr/boundaries/unit_tests)
endif()
