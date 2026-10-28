add_ippl_test(HaloOrder COMPILE_ONLY SOURCES HaloOrder.cpp HaloOrderPulse.cpp HaloOrderMixed.cpp)
# One rank runs the drift checks only; two ranks run every HaloOrder case.
chdr_src_test(HaloOrder SERIAL_ARGS --gtest_filter=HaloOrder.*Drift)

add_ippl_test(OpenBoundaryRegression COMPILE_ONLY)
chdr_src_test(OpenBoundaryRegression TIMEOUT 180 ARGS
              --gtest_filter=OpenBoundaryRegression.OriginalOrder)

add_ippl_test(OpenFaceClosure COMPILE_ONLY)
chdr_src_test(OpenFaceClosure)

foreach(name FDTDBoundary MurSpeed InteriorBounds MurDispatch)
  add_ippl_test(${name} COMPILE_ONLY)
  chdr_src_test(${name})
endforeach()
