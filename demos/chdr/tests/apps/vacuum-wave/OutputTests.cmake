# Included from this directory's CMakeLists.txt while BUILD_TESTING is enabled.
add_executable(test-chdr-vacuum-output output_fixture.cpp)
target_compile_features(test-chdr-vacuum-output PRIVATE cxx_std_20)
# The fixture includes the core writer by its path from the demos/chdr root.
target_include_directories(test-chdr-vacuum-output PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../..")
target_link_libraries(test-chdr-vacuum-output PRIVATE chdr-core)
target_compile_options(test-chdr-vacuum-output
                       PRIVATE $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Werror=sign-compare>)

if(NOT Python3_Interpreter_FOUND)
  chdr_vacuum_skip_test("chdr.vacuum.output.single" "Python3 interpreter unavailable" output reader)
  foreach(_case IN ITEMS serial mpi2 failure.piece failure.wrapper failure.series failure.nonfinite)
    chdr_vacuum_skip_test("chdr.vacuum.output.${_case}" "Python3 interpreter unavailable" output
                          reader)
  endforeach()
  return()
endif()

# The arguments after output_case are the fixture command. The optional keyword ENCODING selects the
# stored value format; without it the launcher chooses by the rank count.
function(chdr_add_output_test test_name ranks output_case)
  cmake_parse_arguments(PARSE_ARGV 3 _output "" "ENCODING" "")
  set(_encoding_arguments)
  if(_output_ENCODING)
    set(_encoding_arguments --encoding "${_output_ENCODING}")
  endif()
  add_test(
    NAME "${test_name}"
    COMMAND
      "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_output_fixture.py" --reader
      "${CHDR_VACUUM_PVPYTHON}" --output-parent "${CHDR_VACUUM_TEST_OUTPUT_DIR}/writer" --ranks
      "${ranks}" --case "${output_case}" ${_encoding_arguments} -- ${_output_UNPARSED_ARGUMENTS})
  set_tests_properties(
    "${test_name}"
    PROPERTIES LABELS
               "chdr;vacuum;output;reader"
               TIMEOUT
               180
               PROCESSORS
               "${ranks}"
               SKIP_RETURN_CODE
               77)
  if(ranks GREATER 1)
    set_property(TEST "${test_name}" APPEND PROPERTY LABELS mpi)
  endif()
endfunction()

chdr_add_output_test(chdr.vacuum.output.serial 1 valid $<TARGET_FILE:test-chdr-vacuum-output>)
# The same sentinel frames with 4-byte values: the reader must return them rounded to binary32.
chdr_add_output_test(chdr.vacuum.output.single 1 valid ENCODING binary32
                     $<TARGET_FILE:test-chdr-vacuum-output>)
if(MPIEXEC_EXECUTABLE)
  set(_chdr_output_mpi "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
                       $<TARGET_FILE:test-chdr-vacuum-output> ${MPIEXEC_POSTFLAGS})
  chdr_add_output_test(chdr.vacuum.output.mpi2 2 valid ${_chdr_output_mpi})
  foreach(_case IN ITEMS piece wrapper series nonfinite)
    chdr_add_output_test("chdr.vacuum.output.failure.${_case}" 2 "${_case}" ${_chdr_output_mpi})
  endforeach()
else()
  foreach(_case IN ITEMS mpi2 failure.piece failure.wrapper failure.series failure.nonfinite)
    chdr_vacuum_skip_test("chdr.vacuum.output.${_case}" "MPI launcher unavailable" output reader)
  endforeach()
endif()
