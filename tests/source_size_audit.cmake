if(NOT DEFINED ROOT)
    message(FATAL_ERROR "ROOT is required")
endif()

set(nonempty_line_baseline 11007)
set(maximum_nonempty_lines 9906)
set(hex8_inelastic_capability_lines 115)
set(hex8_finite_strain_capability_lines 276)
set(hex8_contact_capability_lines 1070)
set(dynamic_contact_assembly_performance_lines 171)
set(large_sliding_contact_search_lines 534)
set(cartesian_current_pressure_input_lines 6)
set(m58_two_process_efficiency_lines 98)
set(m58_four_process_efficiency_lines 234)
set(memory_statistics_lines 73)
set(shared_node_capability_lines 120)
set(hex8_narrow_constitutive_ad_lines 177)
set(rz_narrow_constitutive_ad_lines 179)
set(source_layout_refactor_lines 24)
set(optional_input_sections_lines 5)
set(bound_material_function_lines 109)
set(thermal_time_term_option_lines 16)
set(pressure_configuration_selection_lines 31)
set(hex20_u2_t1_capability_lines 1521)
set(hex20_face_quadrature_lines 38)
set(hex20_contact_capability_lines 1150)
set(hex20_surface_mechanical_contact_lines 136)
set(hex20_small_sliding_surface_contact_lines 238)
set(hex20_abaqus_averaged_contact_lines 407)
set(hex20_curved_abaqus_contact_lines 27)
set(hex20_abaqus_averaged_friction_lines 143)
set(hex20_abaqus_contact_observables_lines 6)
set(hex20_objective_friction_history_lines 69)
set(hex20_finite_sliding_contact_lines 191)
set(hex8_abaqus_small_sliding_surface_contact_lines 317)
set(hex8_finite_sliding_contact_lines 350)
set(hex8_finite_averaged_sts_lines 315)

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
    "${maximum_nonempty_lines} + ${hex8_inelastic_capability_lines} + ${hex8_finite_strain_capability_lines} + ${hex8_contact_capability_lines} + ${dynamic_contact_assembly_performance_lines} + ${large_sliding_contact_search_lines} + ${cartesian_current_pressure_input_lines} + ${m58_two_process_efficiency_lines} + ${m58_four_process_efficiency_lines} + ${memory_statistics_lines} + ${shared_node_capability_lines} + ${hex8_narrow_constitutive_ad_lines} + ${rz_narrow_constitutive_ad_lines} + ${source_layout_refactor_lines} + ${optional_input_sections_lines} + ${bound_material_function_lines} + ${thermal_time_term_option_lines} + ${pressure_configuration_selection_lines} + ${hex20_u2_t1_capability_lines} + ${hex20_face_quadrature_lines} + ${hex20_contact_capability_lines} + ${hex20_surface_mechanical_contact_lines} + ${hex20_small_sliding_surface_contact_lines} + ${hex20_abaqus_averaged_contact_lines} + ${hex20_curved_abaqus_contact_lines} + ${hex20_abaqus_averaged_friction_lines} + ${hex20_abaqus_contact_observables_lines} + ${hex20_objective_friction_history_lines} + ${hex20_finite_sliding_contact_lines} + ${hex8_abaqus_small_sliding_surface_contact_lines} + ${hex8_finite_sliding_contact_lines} + ${hex8_finite_averaged_sts_lines}"
)
if(nonempty_lines GREATER capability_adjusted_limit)
    message(FATAL_ERROR
        "src/include nonempty lines grew to ${nonempty_lines}; the original 10 percent refactor limit is "
        "${maximum_nonempty_lines}, with ${hex8_inelastic_capability_lines} additional lines allowed for "
        "three-dimensional inelasticity and ${hex8_finite_strain_capability_lines} additional lines allowed for "
        "three-dimensional finite strain, plus ${hex8_contact_capability_lines} additional lines allowed for "
        "three-dimensional contact and ${dynamic_contact_assembly_performance_lines} additional lines allowed for "
        "dynamic active-contact assembly performance, plus ${large_sliding_contact_search_lines} additional lines "
        "allowed for exact large-sliding contact search, plus ${cartesian_current_pressure_input_lines} additional "
        "lines allowed for Cartesian current-configuration pressure input, plus ${m58_two_process_efficiency_lines} "
        "lines allowed for M5.8 two-process parallel efficiency, plus ${m58_four_process_efficiency_lines} "
        "additional lines allowed for M5.8 four-process parallel efficiency, plus ${memory_statistics_lines} "
        "additional lines allowed for solver memory statistics, plus ${shared_node_capability_lines} additional lines "
        "allowed for native shared-node material interfaces, plus ${hex8_narrow_constitutive_ad_lines} additional "
        "lines allowed for narrow two-level constitutive automatic differentiation, plus "
        "${rz_narrow_constitutive_ad_lines} additional lines allowed for narrow two-level RZ constitutive "
        "automatic differentiation, plus ${source_layout_refactor_lines} lines for explicit module boundaries "
        "introduced by the source-layout refactor, plus ${optional_input_sections_lines} lines for optional input "
        "sections, plus ${bound_material_function_lines} lines for two-stage bound material evaluators, plus "
        "${thermal_time_term_option_lines} lines for the transient thermal time-term option, plus "
        "${pressure_configuration_selection_lines} lines for selectable pressure-boundary configurations, plus "
        "${hex20_u2_t1_capability_lines} lines for the non-contact mixed-order HEX20-U2/T1 capability, plus "
        "${hex20_face_quadrature_lines} lines for independent HEX20 face quadrature rules, plus "
        "${hex20_contact_capability_lines} lines for HEX20 thermal and mechanical contact, plus "
        "${hex20_surface_mechanical_contact_lines} lines for HEX20 surface-to-surface mechanical contact, plus "
        "${hex20_small_sliding_surface_contact_lines} lines for reference-segmented HEX20 small-sliding contact, plus "
        "${hex20_abaqus_averaged_contact_lines} lines for Abaqus-style averaged HEX20 contact, plus "
        "${hex20_curved_abaqus_contact_lines} lines for curved Abaqus-style HEX20 contact, plus "
        "${hex20_abaqus_averaged_friction_lines} lines for Abaqus-style averaged HEX20 Coulomb friction, plus "
        "${hex20_abaqus_contact_observables_lines} lines for signed force-vector and total-slip observables, plus "
        "${hex20_objective_friction_history_lines} lines for objective HEX20 friction-history transport, plus "
        "${hex20_finite_sliding_contact_lines} lines for HEX20 finite-sliding current-configuration search and "
        "history, plus ${hex8_abaqus_small_sliding_surface_contact_lines} lines for Abaqus-style HEX8 "
        "small-sliding surface contact, plus ${hex8_finite_sliding_contact_lines} lines for HEX8 finite-sliding "
        "current-configuration search and objective friction history, plus ${hex8_finite_averaged_sts_lines} lines "
        "for the identified HEX8 finite-sliding averaged surface operator"
    )
endif()

message(STATUS
    "120-column src/include nonempty lines: ${nonempty_lines}; removed ${removed_lines} from "
    "${nonempty_line_baseline}; "
    "reduction ${reduction_per_mille} per mille; original refactor limit ${maximum_nonempty_lines} plus "
    "${hex8_inelastic_capability_lines} three-dimensional inelastic capability lines plus "
    "${hex8_finite_strain_capability_lines} three-dimensional finite-strain capability lines plus "
    "${hex8_contact_capability_lines} three-dimensional contact capability lines plus "
    "${dynamic_contact_assembly_performance_lines} dynamic active-contact assembly performance lines plus "
    "${large_sliding_contact_search_lines} exact large-sliding contact-search capability lines plus "
    "${cartesian_current_pressure_input_lines} Cartesian current-pressure input lines plus "
    "${m58_two_process_efficiency_lines} M5.8 two-process parallel-efficiency lines plus "
    "${m58_four_process_efficiency_lines} M5.8 four-process parallel-efficiency lines plus "
    "${memory_statistics_lines} solver-memory-statistics lines plus ${shared_node_capability_lines} native "
    "shared-node material-interface lines plus ${hex8_narrow_constitutive_ad_lines} narrow two-level constitutive-"
    "differentiation lines plus ${rz_narrow_constitutive_ad_lines} narrow two-level RZ constitutive-differentiation "
    "lines plus ${source_layout_refactor_lines} source-layout refactor boundary lines plus "
    "${optional_input_sections_lines} optional-input-section lines plus ${bound_material_function_lines} "
    "two-stage bound-material-evaluator lines plus ${thermal_time_term_option_lines} transient-thermal-time-term "
    "option lines plus ${pressure_configuration_selection_lines} selectable-pressure-configuration lines plus "
    "${hex20_u2_t1_capability_lines} non-contact mixed-order HEX20-U2/T1 capability lines plus "
    "${hex20_face_quadrature_lines} independent HEX20 face-quadrature lines plus "
    "${hex20_contact_capability_lines} HEX20 thermal and mechanical contact lines plus "
    "${hex20_surface_mechanical_contact_lines} HEX20 surface-to-surface mechanical-contact lines plus "
    "${hex20_small_sliding_surface_contact_lines} reference-segmented HEX20 small-sliding contact lines plus "
    "${hex20_abaqus_averaged_contact_lines} Abaqus-style averaged HEX20 contact lines plus "
    "${hex20_curved_abaqus_contact_lines} curved Abaqus-style HEX20 contact lines plus "
    "${hex20_abaqus_averaged_friction_lines} Abaqus-style averaged HEX20 Coulomb-friction lines plus "
    "${hex20_abaqus_contact_observables_lines} signed force-vector and total-slip observable lines plus "
    "${hex20_objective_friction_history_lines} objective HEX20 friction-history transport lines plus "
    "${hex20_finite_sliding_contact_lines} HEX20 finite-sliding current-configuration search and history lines plus "
    "${hex8_abaqus_small_sliding_surface_contact_lines} Abaqus-style HEX8 small-sliding surface-contact lines plus "
    "${hex8_finite_sliding_contact_lines} HEX8 finite-sliding current-configuration search and history lines"
)
