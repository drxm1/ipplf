include(AddIpplTest)
add_ippl_test(test-chdr-diagnostics COMPILE_ONLY SOURCES test_interface_energy.cpp)
# SerialReference.h (borrowed from tests/core/solver) brings in RankRecords.h, which uses JSON.
target_link_libraries(test-chdr-diagnostics PRIVATE chdr-core chdr-diagnostics chdr-testing
                                                    nlohmann_json::nlohmann_json)
target_compile_options(test-chdr-diagnostics PRIVATE -UNDEBUG)
# As before the split: one rank runs the normal split only, two ranks run both splits.
add_test(NAME chdr.diagnostics.energy.serial COMMAND test-chdr-diagnostics
                                                     --gtest_filter=InterfaceEnergy.NormalSplit)
set_tests_properties(
  chdr.diagnostics.energy.serial PROPERTIES TIMEOUT 60 LABELS "chdr;diagnostics;unit" ENVIRONMENT
                                            "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
if(MPIEXEC_EXECUTABLE)
  add_test(NAME chdr.diagnostics.energy.mpi2
           COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:test-chdr-diagnostics> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(
    chdr.diagnostics.energy.mpi2
    PROPERTIES TIMEOUT
               60
               PROCESSORS
               2
               LABELS
               "chdr;diagnostics;unit;mpi"
               ENVIRONMENT
               "OMP_NUM_THREADS=1;OMP_PROC_BIND=false")
endif()
