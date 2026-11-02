# A build-time snapshot of Git identity. Failure must never look like a clean checkout.
set(chdr_source_revision "unrecorded")
set(chdr_source_dirty "unrecorded")
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY "${CHDR_SOURCE_DIRECTORY}"
  RESULT_VARIABLE chdr_revision_result
  OUTPUT_VARIABLE chdr_revision
  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
execute_process(
  COMMAND git status --porcelain
  WORKING_DIRECTORY "${CHDR_SOURCE_DIRECTORY}"
  RESULT_VARIABLE chdr_status_result
  OUTPUT_VARIABLE chdr_changes
  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(chdr_revision_result STREQUAL "0" AND chdr_status_result STREQUAL "0" AND chdr_revision MATCHES
                                                                             "^[0-9a-fA-F]+$")
  set(chdr_source_revision "${chdr_revision}")
  if(chdr_changes)
    set(chdr_source_dirty "true")
  else()
    set(chdr_source_dirty "false")
  endif()
endif()

set(chdr_contents
    "#pragma once\n#define CHDR_SOURCE_REVISION \"${chdr_source_revision}\"\n#define CHDR_SOURCE_DIRTY \"${chdr_source_dirty}\"\n"
)
set(chdr_previous "")
if(EXISTS "${CHDR_OUTPUT_HEADER}")
  file(READ "${CHDR_OUTPUT_HEADER}" chdr_previous)
endif()
if(NOT chdr_contents STREQUAL chdr_previous)
  get_filename_component(chdr_header_directory "${CHDR_OUTPUT_HEADER}" DIRECTORY)
  file(MAKE_DIRECTORY "${chdr_header_directory}")
  file(WRITE "${CHDR_OUTPUT_HEADER}" "${chdr_contents}")
endif()
