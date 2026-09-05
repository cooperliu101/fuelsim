# Keep one coarse mesh and one distorted mesh in the unified regression. The
# regular refined mesh remains a manually runnable sensitivity case.
foreach(branch IN ITEMS coarse distorted)
    set(input "${CASE_DIRECTORY}/transient_b59_${branch}.fsi")
    set(result "${CASE_DIRECTORY}/transient_b59_${branch}_results.e")
    set(summary "${CASE_DIRECTORY}/transient_b59_${branch}_summary.csv")
    file(REMOVE "${result}" "${summary}")
    execute_process(COMMAND "${FUELSIM_EXECUTABLE}" -i "${input}"
        RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(NOT status EQUAL 0 OR NOT EXISTS "${result}" OR NOT EXISTS "${summary}")
        message(FATAL_ERROR "B5.9 ${branch} production run failed: ${status}\n${output}\n${error}")
    endif()
endforeach()
execute_process(COMMAND "${CHECK_COMMAND}" hex8-multimaterial
    "${CASE_DIRECTORY}/transient_b59_coarse_results.e"
    "${CASE_DIRECTORY}/transient_b59_distorted_results.e"
    "${REFERENCE_DIRECTORY}/b59_hex8_c3d8t_multimaterial_nodal.csv"
    "${REFERENCE_DIRECTORY}/b59_hex8_c3d8t_multimaterial_integration.csv"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "B5.9 production result comparison failed: ${status}\n${output}\n${error}")
endif()
message(STATUS "${output}")
# Only shared-node ownership and the original cooled steady assembly audit remain.
execute_process(COMMAND "${INTERNAL_COMMAND}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "B5.9 internal contracts failed: ${status}\n${output}\n${error}")
endif()
message(STATUS "${output}")
