# Replace only this build's dedicated smoke-output directory before invoking the writer.
file(REMOVE_RECURSE "${OUTPUT}")
execute_process(COMMAND "${PROGRAM}" half-space --smoke --output "${OUTPUT}" COMMAND_ERROR_IS_FATAL
                        ANY)
