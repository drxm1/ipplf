# Register channel-wiring checks and the serial/two-rank spectrum comparison. The basic point and
# configuration tests are registered by the parent CMakeLists.txt.
target_sources(test-chdr-solver PRIVATE test_probe_full.cpp)
chdr_solver_test(
  probe_wiring
  "FourierProbeT8.IdenticalSamplesGiveTheLabelShiftFactor:FourierProbeT8.ShellLabelsAndChannelCoordinates:FourierProbeT8.StaggeredPlaneClosedForm:FourierProbeT8.LiveSolverChannels"
)
chdr_solver_test(probe_gates "FourierProbeT8.FrozenGatesAndNaNGhosts")
