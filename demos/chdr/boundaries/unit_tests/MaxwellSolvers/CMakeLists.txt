add_ippl_test(HaloOrder NUM_PROCS 2 PROPERTIES PROCESSORS 2)
add_test(NAME HaloOrderSerial COMMAND HaloOrder --gtest_filter=HaloOrder.*Drift)
set_tests_properties(HaloOrderSerial PROPERTIES PROCESSORS 1 TIMEOUT 60 LABELS unit)
set_tests_properties(HaloOrder HaloOrderSerial PROPERTIES ENVIRONMENT
                                                          "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
