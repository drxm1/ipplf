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
foreach(index RANGE 0 2)
  math(EXPR step "2 * ${index}")
  string(
    JSON
    name
    GET
    "${series}"
    files
    ${index}
    name)
  if(NOT name STREQUAL "wave-${step}.pvti")
    message(FATAL_ERROR "Unexpected field frame at series index ${index}: ${name}")
  endif()
  if(NOT EXISTS "${OUTPUT}/${name}" OR NOT EXISTS "${OUTPUT}/wave-${step}-rank-0.vti")
    message(FATAL_ERROR "Field frame at state ${step} is incomplete")
  endif()
endforeach()
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
