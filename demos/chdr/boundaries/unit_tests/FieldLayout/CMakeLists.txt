add_ippl_test(PeriodicAxes PROPERTIES PROCESSORS 2)
add_test(NAME PeriodicAxesSerial COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:PeriodicAxes>)
set_tests_properties(PeriodicAxesSerial PROPERTIES PROCESSORS 1 LABELS unit TIMEOUT 60)
