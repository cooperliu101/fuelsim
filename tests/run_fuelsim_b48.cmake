# Complete and split/restarted paths use explicit, unmodified production cards.
foreach(stem IN ITEMS transient_hex8_multi_contact_path_abaqus transient_b48_checkpoint transient_b48_restart)
    file(REMOVE "${CASE_DIRECTORY}/${stem}_results.e" "${CASE_DIRECTORY}/${stem}_summary.csv")
endforeach()
file(REMOVE "${CASE_DIRECTORY}/b48_full.checkpoint" "${CASE_DIRECTORY}/b48_split.checkpoint"
    "${CASE_DIRECTORY}/b48_restart.checkpoint" "${CASE_DIRECTORY}/transient_b48_restart_results.part1.e")
foreach(stem IN ITEMS transient_hex8_multi_contact_path_abaqus transient_b48_checkpoint transient_b48_restart)
    execute_process(COMMAND "${FUELSIM_EXECUTABLE}" -i "${CASE_DIRECTORY}/${stem}.fsi"
        RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR "B4.8 ${stem} production solve failed: ${status}\n${output}\n${error}")
    endif()
endforeach()
execute_process(COMMAND "${CHECK_COMMAND}" hex8-multi-contact-path
    "${CASE_DIRECTORY}/transient_hex8_multi_contact_path_abaqus_results.e"
    "${CASE_DIRECTORY}/transient_hex8_multi_contact_path_abaqus_summary.csv" "${REFERENCE_DIRECTORY}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "B4.8 Abaqus comparison failed: ${status}\n${output}\n${error}")
endif()
message(STATUS "${output}")
execute_process(COMMAND "${CHECK_COMMAND}" b48-restart
    "${CASE_DIRECTORY}/transient_b48_restart_results.part1.e" "${CASE_DIRECTORY}/transient_b48_restart_summary.csv"
    "${CASE_DIRECTORY}/transient_hex8_multi_contact_path_abaqus_results.e"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "B4.8 restart output comparison failed: ${status}\n${output}\n${error}")
endif()
message(STATUS "${output}")
execute_process(COMMAND "${INTERNAL_COMMAND}" "${CASE_DIRECTORY}/transient_hex8_multi_contact_path_abaqus.fsi"
    "${CASE_DIRECTORY}/b48_full.checkpoint" "${CASE_DIRECTORY}/b48_restart.checkpoint"
    "${CASE_DIRECTORY}/b48_split.checkpoint"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "B4.8 retained internal contracts failed: ${status}\n${output}\n${error}")
endif()
message(STATUS "${output}")
