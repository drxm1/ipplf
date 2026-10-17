add_test(NAME chdr.interface.unit.smoke COMMAND interface-demo half-space --smoke)
set_tests_properties(
  chdr.interface.unit.smoke PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1;OMP_PROC_BIND=false" TIMEOUT
                                       60 LABELS "chdr;interface;unit")
