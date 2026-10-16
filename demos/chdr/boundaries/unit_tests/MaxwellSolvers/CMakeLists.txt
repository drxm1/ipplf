add_ippl_test(HaloOrder NUM_PROCS 2 PROPERTIES PROCESSORS 2)
target_sources(HaloOrder PRIVATE HaloOrderPulse.cpp)
add_test(NAME HaloOrderSerial COMMAND HaloOrder --gtest_filter=HaloOrder.*Drift)
set_tests_properties(HaloOrderSerial PROPERTIES PROCESSORS 1 TIMEOUT 60 LABELS unit)
set_tests_properties(HaloOrder HaloOrderSerial PROPERTIES ENVIRONMENT
                                                          "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")

add_ippl_test(
  OpenBoundaryRegression
  NUM_PROCS
  2
  TIMEOUT
  180
  ARGS
  --gtest_filter=OpenBoundaryRegression.OriginalOrder
  PROPERTIES
  PROCESSORS
  2)
add_test(NAME OpenBoundaryRegressionSerial
         COMMAND OpenBoundaryRegression --gtest_filter=OpenBoundaryRegression.OriginalOrder)
set_tests_properties(OpenBoundaryRegressionSerial PROPERTIES PROCESSORS 1 TIMEOUT 180 LABELS unit)
set_tests_properties(OpenBoundaryRegression OpenBoundaryRegressionSerial
                     PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
