if(NOT DEFINED ROOT)
    message(FATAL_ERROR "ROOT is required")
endif()

set(stage_b_baseline 16216)
set(maximum_source_lines 12972)
file(GLOB_RECURSE source_files
    "${ROOT}/src/*.cpp"
    "${ROOT}/src/*.hpp"
    "${ROOT}/include/*.hpp"
)

set(source_lines 0)
foreach(source_file IN LISTS source_files)
    file(READ "${source_file}" contents)
    string(REGEX REPLACE "[^\n]" "" newlines "${contents}")
    string(LENGTH "${newlines}" file_lines)
    math(EXPR source_lines "${source_lines} + ${file_lines}")
endforeach()

math(EXPR removed_lines "${stage_b_baseline} - ${source_lines}")
math(EXPR reduction_per_mille "1000 * ${removed_lines} / ${stage_b_baseline}")
if(source_lines GREATER maximum_source_lines)
    message(FATAL_ERROR
        "src/include physical lines grew to ${source_lines}; the post-stage-B limit is ${maximum_source_lines} "
        "from the ${stage_b_baseline}-line baseline"
    )
endif()

message(STATUS
    "src/include physical lines: ${source_lines}; removed ${removed_lines} from ${stage_b_baseline}; "
    "reduction ${reduction_per_mille} per mille"
)
