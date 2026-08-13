if(NOT DEFINED ROOT)
    message(FATAL_ERROR "ROOT is required")
endif()

set(coexistence_baseline 12774)
set(maximum_source_lines 11496)

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

set(source_lines 0)
foreach(source_file IN LISTS source_files)
    file(READ "${source_file}" contents)
    string(REGEX REPLACE "[^\n]" "" newlines "${contents}")
    string(LENGTH "${newlines}" file_lines)
    math(EXPR source_lines "${source_lines} + ${file_lines}")
endforeach()

math(EXPR removed_lines "${coexistence_baseline} - ${source_lines}")
math(EXPR reduction_per_mille "1000 * ${removed_lines} / ${coexistence_baseline}")
if(source_lines GREATER maximum_source_lines)
    message(FATAL_ERROR
        "src/include physical lines grew to ${source_lines}; the 10 percent coexistence-refactor limit is "
        "${maximum_source_lines} from the ${coexistence_baseline}-line baseline"
    )
endif()

message(STATUS
    "120-column src/include physical lines: ${source_lines}; removed ${removed_lines} from ${coexistence_baseline}; "
    "reduction ${reduction_per_mille} per mille"
)
