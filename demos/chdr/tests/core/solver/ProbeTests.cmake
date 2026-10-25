# Register channel-wiring checks and the serial/two-rank spectrum comparison. The basic point and
# configuration tests are registered by the parent CMakeLists.txt.
target_sources(test-chdr-solver PRIVATE test_probe_full.cpp)
chdr_solver_test(
  probe_wiring
  "FourierProbeT8.IdenticalSamplesGiveTheLabelShiftFactor:FourierProbeT8.ShellLabelsAndChannelCoordinates:FourierProbeT8.StaggeredPlaneClosedForm:FourierProbeT8.LiveSolverChannels"
)
chdr_solver_test(probe_gates "FourierProbeT8.FrozenGatesAndNaNGhosts")

if(Python3_Interpreter_FOUND AND MPIEXEC_EXECUTABLE)
  set(probe_mpi_flags)
  foreach(flag IN LISTS MPIEXEC_PREFLAGS)
    list(APPEND probe_mpi_flags "--mpi-preflag=${flag}")
  endforeach()
  foreach(flag IN LISTS MPIEXEC_POSTFLAGS)
    list(APPEND probe_mpi_flags "--mpi-postflag=${flag}")
  endforeach()
  add_test(
    NAME chdr.solver.probe.ranks
    COMMAND
      "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/compare_probe_ranks.py" --executable
      $<TARGET_FILE:test-chdr-solver> --mpiexec "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}" ${probe_mpi_flags} --output-dir
      "${CHDR_SOLVER_RECORD_ROOT}/probe-ranks")
  set_tests_properties(
    chdr.solver.probe.ranks
    PROPERTIES LABELS
               "chdr;solver;unit;probe;mpi"
               TIMEOUT
               60
               PROCESSORS
               2
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
else()
  chdr_solver_skip(probe.ranks "Probe rank comparison requires Python3 and an MPI launcher")
endif()
