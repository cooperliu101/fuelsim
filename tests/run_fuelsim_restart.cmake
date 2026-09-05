# Both input cards are complete, checked-in files. This script only schedules runs.
execute_process(
    COMMAND "${CMAKE_COMMAND}"
        "-DFUELSIM_EXECUTABLE=${FUELSIM_EXECUTABLE}"
        "-DINPUT_FILE=${CHECKPOINT_INPUT}"
        "-DOUTPUT_FILE=${CHECKPOINT_FILE}"
        "-DCHECK_COMMAND=${CMAKE_COMMAND}"
        "-DCHECK_ARGUMENTS=-E;true"
        -P "${CMAKE_CURRENT_LIST_DIR}/run_fuelsim_case.cmake"
    RESULT_VARIABLE checkpoint_status
    OUTPUT_VARIABLE checkpoint_stdout
    ERROR_VARIABLE checkpoint_stderr
)
if(NOT checkpoint_status EQUAL 0)
    message(FATAL_ERROR "Checkpoint run failed:\n${checkpoint_stdout}\n${checkpoint_stderr}")
endif()
include("${CMAKE_CURRENT_LIST_DIR}/run_fuelsim_case.cmake")
