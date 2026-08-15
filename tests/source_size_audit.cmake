if(NOT DEFINED ROOT)
    message(FATAL_ERROR "ROOT is required")
endif()

set(nonempty_line_baseline 11007)
set(maximum_nonempty_lines 9906)
set(hex8_inelastic_capability_lines 115)
set(hex8_finite_strain_capability_lines 276)

file(READ "${ROOT}/.clang-format" format_configuration)
string(FIND "${format_configuration}" "ColumnLimit: 120" column_limit)
if(column_limit EQUAL -1)
    message(FATAL_ERROR "The repository clang-format column limit must remain 120")
endif()
foreach(required_style
    "AlignAfterOpenBracket: DontAlign"
    "RemoveBracesLLVM: true"
    "MaxEmptyLinesToKeep: 1"
    "SeparateDefinitionBlocks: Always"
    "EmptyLineBeforeAccessModifier: LogicalBlock"
)
    string(FIND "${format_configuration}" "${required_style}" style_location)
    if(style_location EQUAL -1)
        message(FATAL_ERROR "The repository 120-column source style is missing: ${required_style}")
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
math(EXPR capability_adjusted_limit
    "${maximum_nonempty_lines} + ${hex8_inelastic_capability_lines} + ${hex8_finite_strain_capability_lines}"
)
if(nonempty_lines GREATER capability_adjusted_limit)
    message(FATAL_ERROR
        "src/include nonempty lines grew to ${nonempty_lines}; the original 10 percent refactor limit is "
        "${maximum_nonempty_lines}, with ${hex8_inelastic_capability_lines} additional lines allowed for "
        "three-dimensional inelasticity and ${hex8_finite_strain_capability_lines} additional lines allowed for "
        "three-dimensional finite strain"
    )
endif()

message(STATUS
    "120-column src/include nonempty lines: ${nonempty_lines}; removed ${removed_lines} from "
    "${nonempty_line_baseline}; "
    "reduction ${reduction_per_mille} per mille; original refactor limit ${maximum_nonempty_lines} plus "
    "${hex8_inelastic_capability_lines} three-dimensional inelastic capability lines plus "
    "${hex8_finite_strain_capability_lines} three-dimensional finite-strain capability lines"
)
