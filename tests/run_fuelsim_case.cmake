if(NOT DEFINED FUELSIM_EXECUTABLE OR FUELSIM_EXECUTABLE STREQUAL "")
    message(FATAL_ERROR "FUELSIM_EXECUTABLE is required")
endif()
if(NOT DEFINED INPUT_FILE OR INPUT_FILE STREQUAL "")
    message(FATAL_ERROR "INPUT_FILE is required")
endif()
if(NOT DEFINED CHECK_COMMAND OR CHECK_COMMAND STREQUAL "")
    message(FATAL_ERROR "CHECK_COMMAND is required")
endif()

set(output_files)
if(DEFINED OUTPUT_FILE AND NOT OUTPUT_FILE STREQUAL "")
    list(APPEND output_files "${OUTPUT_FILE}")
endif()
if(DEFINED OUTPUT_FILES AND NOT OUTPUT_FILES STREQUAL "")
    list(APPEND output_files ${OUTPUT_FILES})
endif()
foreach(output_file IN LISTS output_files)
    file(REMOVE "${output_file}")
endforeach()

execute_process(
    COMMAND "${FUELSIM_EXECUTABLE}" -i "${INPUT_FILE}"
    RESULT_VARIABLE fuelsim_status
    OUTPUT_VARIABLE fuelsim_stdout
    ERROR_VARIABLE fuelsim_stderr
)
if(NOT fuelsim_status EQUAL 0)
    message(FATAL_ERROR
        "fuelsim failed with status ${fuelsim_status}\n"
        "stdout:\n${fuelsim_stdout}\n"
        "stderr:\n${fuelsim_stderr}"
    )
endif()
foreach(output_file IN LISTS output_files)
    if(NOT EXISTS "${output_file}")
        message(FATAL_ERROR "fuelsim completed without producing '${output_file}'")
    endif()
endforeach()

set(check_arguments)
if(DEFINED CHECK_ARGUMENTS AND NOT CHECK_ARGUMENTS STREQUAL "")
    list(APPEND check_arguments ${CHECK_ARGUMENTS})
endif()
execute_process(
    COMMAND "${CHECK_COMMAND}" ${check_arguments}
    RESULT_VARIABLE check_status
    OUTPUT_VARIABLE check_stdout
    ERROR_VARIABLE check_stderr
)
if(NOT check_status EQUAL 0)
    message(FATAL_ERROR
        "Result verification failed with status ${check_status}\n"
        "fuelsim stdout:\n${fuelsim_stdout}\n"
        "fuelsim stderr:\n${fuelsim_stderr}\n"
        "verification stdout:\n${check_stdout}\n"
        "verification stderr:\n${check_stderr}"
    )
endif()

message(STATUS "fuelsim stdout:\n${fuelsim_stdout}")
message(STATUS "verification stdout:\n${check_stdout}")
