# Register channel-wiring checks and the serial/two-rank spectrum comparison. The basic point and
# configuration tests are registered by the CMakeLists.txt that includes this file.
target_sources(test-chdr-diagnostics PRIVATE test_probe_full.cpp)
chdr_diagnostics_test(
  probe_wiring
  "FourierProbeSpectra.IdenticalSamplesGiveTheLabelShiftFactor:FourierProbeSpectra.SolverLabelsAndChannelCoordinates:FourierProbeSpectra.StaggeredPlaneClosedForm:FourierProbeSpectra.LiveSolverChannels"
)
chdr_diagnostics_test(probe_gates "FourierProbeSpectra.FrozenGatesAndNaNGhosts")

if(Python3_Interpreter_FOUND AND MPIEXEC_EXECUTABLE)
  set(probe_mpi_flags)
  foreach(flag IN LISTS MPIEXEC_PREFLAGS)
    list(APPEND probe_mpi_flags "--mpi-preflag=${flag}")
  endforeach()
  foreach(flag IN LISTS MPIEXEC_POSTFLAGS)
    list(APPEND probe_mpi_flags "--mpi-postflag=${flag}")
  endforeach()
  add_test(
    NAME chdr.diagnostics.probe.ranks
    COMMAND
      "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/compare_probe_ranks.py" --executable
      $<TARGET_FILE:test-chdr-diagnostics> --mpiexec "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}" ${probe_mpi_flags} --output-dir
      "${CHDR_DIAGNOSTICS_RECORD_ROOT}/probe-ranks")
  # The script runs FourierProbeSpectra.FrozenGatesAndNaNGhosts on one rank and on two ranks, each
  # run bounded by 60 s inside the script, and compares the two record files.
  set_tests_properties(
    chdr.diagnostics.probe.ranks
    PROPERTIES LABELS
               "chdr;diagnostics;unit;probe;mpi"
               TIMEOUT
               150
               PROCESSORS
               2
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
else()
  chdr_diagnostics_skip(probe.ranks "Probe rank comparison requires Python3 and an MPI launcher")
endif()
