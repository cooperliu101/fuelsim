if(NOT DEFINED ROOT)
    message(FATAL_ERROR "ROOT is required")
endif()

set(nonempty_line_baseline 11007)
set(maximum_nonempty_lines 9906)

file(READ "${ROOT}/.clang-format" format_configuration)
string(FIND "${format_configuration}" "ColumnLimit: 120" column_limit)
if(column_limit EQUAL -1)
    message(FATAL_ERROR "The repository clang-format column limit must remain 120")
endif()
foreach(required_style "AlignAfterOpenBracket: DontAlign" "RemoveBracesLLVM: true" "MaxEmptyLinesToKeep: 0")
    string(FIND "${format_configuration}" "${required_style}" style_location)
    if(style_location EQUAL -1)
        message(FATAL_ERROR "The compact 120-column source style is missing: ${required_style}")
    endif()
endforeach()
file(GLOB_RECURSE source_files
    "${ROOT}/src/*.cpp"
    "${ROOT}/src/*.hpp"
    "${ROOT}/include/*.hpp"
)

execute_process(
    COMMAND awk "NF { count++ } END { print count }" ${source_files}
    RESULT_VARIABLE count_status
    OUTPUT_VARIABLE nonempty_lines
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT count_status EQUAL 0 OR NOT nonempty_lines MATCHES "^[0-9]+$")
    message(FATAL_ERROR "Could not count src/include nonempty lines")
endif()

math(EXPR removed_lines "${nonempty_line_baseline} - ${nonempty_lines}")
math(EXPR reduction_per_mille "1000 * ${removed_lines} / ${nonempty_line_baseline}")
if(nonempty_lines GREATER maximum_nonempty_lines)
    message(FATAL_ERROR
        "src/include nonempty lines grew to ${nonempty_lines}; the 10 percent refactor limit is "
        "${maximum_nonempty_lines} from the ${nonempty_line_baseline}-line baseline"
    )
endif()

message(STATUS
    "120-column src/include nonempty lines: ${nonempty_lines}; removed ${removed_lines} from "
    "${nonempty_line_baseline}; "
    "reduction ${reduction_per_mille} per mille"
)
