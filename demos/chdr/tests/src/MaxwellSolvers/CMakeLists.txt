add_ippl_test(HaloOrder COMPILE_ONLY SOURCES HaloOrder.cpp HaloOrderPulse.cpp HaloOrderMixed.cpp)
# The two-rank suite also runs the drift checks on MPI_COMM_SELF.
chdr_src_test(HaloOrder MPI2_ONLY)

add_ippl_test(OpenBoundaryRegression COMPILE_ONLY)
chdr_src_test(OpenBoundaryRegression TIMEOUT 180)

add_ippl_test(OpenFaceClosure COMPILE_ONLY)
chdr_src_test(OpenFaceClosure)

add_ippl_test(MurSpeed COMPILE_ONLY)
chdr_src_test(MurSpeed SERIAL_ONLY)

foreach(name FDTDBoundary InteriorBounds MurDispatch)
  add_ippl_test(${name} COMPILE_ONLY)
  chdr_src_test(${name})
endforeach()
