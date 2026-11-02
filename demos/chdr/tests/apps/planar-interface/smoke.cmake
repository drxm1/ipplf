# Replace only this build's dedicated smoke-output directory before invoking the writer.
file(REMOVE_RECURSE "${OUTPUT}")
file(
  WRITE "${OUTPUT}.yaml"
  "schema_version: 1\ncase: half-space\nunits: {system: normalized}\nparameters: {smoke: true}\noutput: {directory: '${OUTPUT}'}\nobservations:\n  fields: {every: 2}\n  diagnostics: {every: 3}\n"
)
execute_process(COMMAND "${PROGRAM}" --input "${OUTPUT}.yaml" COMMAND_ERROR_IS_FATAL ANY)
file(READ "${OUTPUT}/wave.pvti.series" series)
string(JSON count LENGTH "${series}" files)
if(NOT count EQUAL 3)
  message(FATAL_ERROR "Expected field frames at states 0, 2 and 4")
endif()
file(STRINGS "${OUTPUT}/energy.csv" rows REGEX "^[0-9]+,")
list(LENGTH rows count)
if(NOT count EQUAL 3)
  message(FATAL_ERROR "Expected energy records at states 0, 3 and 4")
endif()
list(GET rows 1 middle)
if(NOT middle MATCHES "^3,")
  message(FATAL_ERROR "Energy cadence must be independent of field output")
endif()
file(READ "${OUTPUT}/run.json" record)
string(JSON count GET "${record}" effective completed_solves)
if(NOT count EQUAL 4)
  message(FATAL_ERROR "Completed run record is missing the four solver steps")
endif()

# The historical command runs the same four-step case with its original output cadence.
set(legacy "${OUTPUT}-legacy")
file(REMOVE_RECURSE "${legacy}")
execute_process(COMMAND "${PROGRAM}" half-space --smoke --output "${legacy}" COMMAND_ERROR_IS_FATAL
                        ANY)
file(READ "${legacy}/run.json" legacy_record)
foreach(key solver box completed_solves)
  string(JSON settings GET "${record}" effective "${key}")
  string(JSON legacy_settings GET "${legacy_record}" effective "${key}")
  if(NOT "${settings}" STREQUAL "${legacy_settings}")
    message(FATAL_ERROR "Legacy and YAML packet settings differ: ${key}")
  endif()
endforeach()
foreach(step 0 4)
  execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${OUTPUT}/wave-${step}-rank-0.vti"
                          "${legacy}/wave-${step}-rank-0.vti" COMMAND_ERROR_IS_FATAL ANY)
endforeach()
