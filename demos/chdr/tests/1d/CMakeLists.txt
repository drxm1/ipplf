# The one-dimensional spectral prototype (demos/chdr/1d, Python) carries 13 unittest checks of its
# physical identities and its command line. They need NumPy, SciPy and Matplotlib; CHDR_1D_PYTHON
# names an interpreter that has them (default: the Python3 found by the build). Without them the
# entry is registered as skipped.
find_package(Python3 QUIET COMPONENTS Interpreter)
set(CHDR_1D_PYTHON "${Python3_EXECUTABLE}"
    CACHE FILEPATH
          "Python interpreter with NumPy, SciPy and Matplotlib for the tests of demos/chdr/1d")
set(chdr_1d_source_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../1d")
set(chdr_1d_ready FALSE)
if(CHDR_1D_PYTHON)
  execute_process(COMMAND "${CHDR_1D_PYTHON}" -c "import numpy, scipy, matplotlib"
                  RESULT_VARIABLE chdr_1d_import_status OUTPUT_QUIET ERROR_QUIET)
  if(chdr_1d_import_status EQUAL 0)
    set(chdr_1d_ready TRUE)
  endif()
endif()
if(chdr_1d_ready)
  add_test(NAME chdr.1d.unittest COMMAND "${CHDR_1D_PYTHON}" -m unittest -v test_chdr_1d
           WORKING_DIRECTORY "${chdr_1d_source_dir}")
  set_tests_properties(chdr.1d.unittest PROPERTIES TIMEOUT 120 LABELS "chdr;1d;unit" ENVIRONMENT
                                                   "PYTHONDONTWRITEBYTECODE=1;MPLBACKEND=Agg")
else()
  add_test(NAME chdr.1d.unittest
           COMMAND "${CMAKE_COMMAND}" -E echo
                   "SKIP: NumPy, SciPy or Matplotlib missing in CHDR_1D_PYTHON=${CHDR_1D_PYTHON}")
  set_tests_properties(chdr.1d.unittest PROPERTIES SKIP_REGULAR_EXPRESSION "^SKIP:" TIMEOUT 10
                                                   LABELS "chdr;1d;unit;unavailable")
endif()
