cmake_minimum_required(VERSION 3.20)

project(
    fuelsim
    VERSION 0.1.0
    DESCRIPTION "Minimal axisymmetric nuclear-fuel thermo-mechanics solver"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(FUELSIM_ENABLE_SANITIZERS
       "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(FUELSIM_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(FUELSIM_ENABLE_RELEASE_IPO "Enable interprocedural optimization for Release builds" ON)

if(FUELSIM_ENABLE_RELEASE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT fuelsim_ipo_supported OUTPUT fuelsim_ipo_error LANGUAGES CXX)
    if(NOT fuelsim_ipo_supported)
        message(FATAL_ERROR "Release interprocedural optimization is unavailable: ${fuelsim_ipo_error}")
    endif()
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

find_package(adlite 0.2.3 CONFIG REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(SEACASExodus CONFIG REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc petsc)

function(fuelsim_enable_warnings target)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options("${target}" PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion
            -Wsign-conversion
        )
        if(FUELSIM_WARNINGS_AS_ERRORS)
            target_compile_options("${target}" PRIVATE -Werror)
        endif()
        if(FUELSIM_ENABLE_SANITIZERS)
            target_compile_options("${target}" PRIVATE
                -fsanitize=address,undefined
                -fno-omit-frame-pointer
            )
            target_link_options("${target}" PRIVATE
                -fsanitize=address,undefined
            )
        endif()
    elseif(MSVC)
        target_compile_options("${target}" PRIVATE /W4 /permissive-)
        if(FUELSIM_WARNINGS_AS_ERRORS)
            target_compile_options("${target}" PRIVATE /WX)
        endif()
    endif()
endfunction()

function(fuelsim_link_petsc target)
    target_link_libraries("${target}" PRIVATE PkgConfig::PETSC)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND PETSC_STATIC_LIBRARY_DIRS)
        list(JOIN PETSC_STATIC_LIBRARY_DIRS ":" petsc_rpath_link_dirs)
        target_link_options("${target}" PUBLIC
            "LINKER:-rpath-link,${petsc_rpath_link_dirs}"
        )
    endif()
endfunction()

include(CTest)
add_subdirectory(elements)

add_library(fuelsim_core
    src/core/radial_assembly.cpp
    src/core/rz_assembly.cpp
    src/core/cartesian3d_assembly.cpp
    src/core/contact_search.cpp
    src/core/rz8_assembly.cpp
    src/core/rz8_contact.cpp
    src/core/mesh.cpp
    src/core/mesh_bar2.cpp
    src/core/mesh_quad8.cpp
    src/core/nonlinear_problem.cpp
    src/core/spatial_problem.cpp
)
add_library(fuelsim::core ALIAS fuelsim_core)

target_include_directories(fuelsim_core
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_link_libraries(fuelsim_core PUBLIC fuelsim::elements)
fuelsim_enable_warnings(fuelsim_core)

add_library(fuelsim_io
    src/io/case_input.cpp
    src/io/checkpoint.cpp
    src/io/problem_signature.cpp
    src/io/results_io.cpp
    src/io/hex8_result_fields.cpp
)
add_library(fuelsim::io ALIAS fuelsim_io)

target_include_directories(fuelsim_io
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_link_libraries(fuelsim_io
    PUBLIC fuelsim::core
    PRIVATE SEACASExodus::exodus
)
fuelsim_enable_warnings(fuelsim_io)

add_library(fuelsim_solver
    src/solver/solver.cpp
    src/solver/solve_workflows.cpp
)
add_library(fuelsim::solver ALIAS fuelsim_solver)

target_include_directories(fuelsim_solver
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_link_libraries(fuelsim_solver
    PUBLIC fuelsim::core
)
fuelsim_link_petsc(fuelsim_solver)
fuelsim_enable_warnings(fuelsim_solver)

add_executable(fuelsim
    src/app/main.cpp
)
target_include_directories(fuelsim PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(fuelsim
    PRIVATE fuelsim::io fuelsim::solver
)
fuelsim_enable_warnings(fuelsim)

include(CTest)

if(BUILD_TESTING)
    add_test(NAME fuelsim_cax2t_gps_tests COMMAND fuelsim_cax2t_gps_tests)
    add_test(NAME fuelsim_ring_gps_tests COMMAND fuelsim_ring_gps_tests)
    add_executable(fuelsim_bar2_mesh_tests tests/bar2_mesh_tests.cpp)
    target_link_libraries(fuelsim_bar2_mesh_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_bar2_mesh_tests)
    add_test(NAME fuelsim_bar2_mesh_tests COMMAND fuelsim_bar2_mesh_tests)
    add_executable(fuelsim_bar2_io_tests tests/bar2_io_tests.cpp)
    target_link_libraries(fuelsim_bar2_io_tests PRIVATE fuelsim::io SEACASExodus::exodus)
    fuelsim_enable_warnings(fuelsim_bar2_io_tests)
    add_test(NAME fuelsim_bar2_io_tests COMMAND fuelsim_bar2_io_tests "${CMAKE_CURRENT_BINARY_DIR}/bar2_io")
    add_test(NAME fuelsim_cax4t_tests COMMAND fuelsim_cax4t_tests)
    add_test(NAME fuelsim_cax4rt_tests COMMAND fuelsim_cax4rt_tests)
    add_test(NAME fuelsim_cax8t_tests COMMAND fuelsim_cax8t_tests)
    add_test(NAME fuelsim_cax8rt_tests COMMAND fuelsim_cax8rt_tests)
    add_test(NAME fuelsim_c3d8t_tests COMMAND fuelsim_c3d8t_tests)
    add_test(NAME fuelsim_c3d8rt_tests COMMAND fuelsim_c3d8rt_tests)
    add_test(NAME fuelsim_c3d20t_tests COMMAND fuelsim_c3d20t_tests)
    add_test(NAME fuelsim_c3d20rt_tests COMMAND fuelsim_c3d20rt_tests)
    add_test(NAME fuelsim_line2_rz_tests COMMAND fuelsim_line2_rz_tests)
    add_test(NAME fuelsim_line3_rz_tests COMMAND fuelsim_line3_rz_tests)
    add_test(NAME fuelsim_quad4_face_tests COMMAND fuelsim_quad4_face_tests)
    add_test(NAME fuelsim_quad8_face_tests COMMAND fuelsim_quad8_face_tests)
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/elements/tests")

    add_library(fuelsim_test_mesh_fixture
        tests/support/mesh_fixture.cpp
    )
    target_include_directories(fuelsim_test_mesh_fixture
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    )
    target_link_libraries(fuelsim_test_mesh_fixture PUBLIC fuelsim::core)
    fuelsim_enable_warnings(fuelsim_test_mesh_fixture)

    add_library(fuelsim_exodus_result_test_support
        tests/support/exodus_result_reader.cpp
    )
    target_include_directories(fuelsim_exodus_result_test_support
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    )
    target_link_libraries(fuelsim_exodus_result_test_support
        PRIVATE SEACASExodus::exodus
    )
    fuelsim_enable_warnings(fuelsim_exodus_result_test_support)

    add_library(fuelsim_field_error_metrics tests/support/field_error_metrics.cpp)
    target_include_directories(fuelsim_field_error_metrics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_field_error_metrics)

    add_executable(fuelsim_radial_assembly_tests tests/radial_assembly_tests.cpp)
    target_link_libraries(fuelsim_radial_assembly_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_radial_assembly_tests)
    add_test(NAME fuelsim_radial_assembly_tests COMMAND fuelsim_radial_assembly_tests)
    add_executable(fuelsim_radial_transaction_tests tests/radial_transaction_tests.cpp)
    target_link_libraries(fuelsim_radial_transaction_tests PRIVATE fuelsim::io)
    fuelsim_enable_warnings(fuelsim_radial_transaction_tests)
    add_test(NAME fuelsim_radial_transaction_tests COMMAND fuelsim_radial_transaction_tests
        "${CMAKE_CURRENT_BINARY_DIR}/radial_transaction.checkpoint")
    add_executable(fuelsim_production_radial tests/production_radial.cpp)
    target_link_libraries(fuelsim_production_radial PRIVATE
        fuelsim_exodus_result_test_support fuelsim_field_error_metrics SEACASExodus::exodus)
    fuelsim_enable_warnings(fuelsim_production_radial)
    add_executable(fuelsim_production_radial_equivalence tests/production_radial_equivalence.cpp)
    target_link_libraries(fuelsim_production_radial_equivalence PRIVATE
        fuelsim_exodus_result_test_support SEACASExodus::exodus)
    fuelsim_enable_warnings(fuelsim_production_radial_equivalence)
    add_executable(fuelsim_production_radial_chain tests/production_radial_chain.cpp)
    target_link_libraries(fuelsim_production_radial_chain PRIVATE
        fuelsim_exodus_result_test_support fuelsim_field_error_metrics)
    fuelsim_enable_warnings(fuelsim_production_radial_chain)
    add_executable(fuelsim_production_radial_inelastic tests/production_radial_inelastic.cpp)
    target_link_libraries(fuelsim_production_radial_inelastic PRIVATE
        fuelsim_exodus_result_test_support fuelsim_field_error_metrics)
    fuelsim_enable_warnings(fuelsim_production_radial_inelastic)

    add_executable(fuelsim_production_pcmi_abaqus tests/production_pcmi_abaqus.cpp)
    target_link_libraries(fuelsim_production_pcmi_abaqus PRIVATE
        fuelsim_exodus_result_test_support fuelsim_field_error_metrics)
    fuelsim_enable_warnings(fuelsim_production_pcmi_abaqus)
    add_executable(fuelsim_production_nts_abaqus tests/production_nts_abaqus.cpp)
    target_link_libraries(fuelsim_production_nts_abaqus PRIVATE
        fuelsim_exodus_result_test_support fuelsim_field_error_metrics)
    fuelsim_enable_warnings(fuelsim_production_nts_abaqus)

    add_executable(fuelsim_production_result_tests
        tests/production_result_tests.cpp
        tests/production_rz_abaqus.cpp
        tests/production_cax4t_material_temperature.cpp
        tests/production_hex20_contact.cpp
        tests/production_hex8_contact.cpp
        tests/production_hex8_plate.cpp
        tests/production_hex8_sticking.cpp
        tests/production_hex20_curved.cpp
        tests/production_hex20_curved_friction.cpp
        tests/production_hex8_sliding.cpp
        tests/production_hex8_biaxial_friction.cpp
        tests/production_hex20_partial_contact.cpp
        tests/production_hex20_friction_path_33.cpp
        tests/production_hex20_friction_path_36.cpp
        tests/production_hex20_integrated.cpp
        tests/production_hex8_b55.cpp
        tests/production_hex8_b58.cpp
        tests/production_hex8_b531.cpp
        tests/production_hex8_inelastic_abaqus.cpp
        tests/production_hex8_norton_abaqus.cpp
        tests/production_hex8_multimaterial.cpp
        tests/production_hex8_finite_contact.cpp
        tests/production_hex8_multi_contact.cpp
        tests/support/production_hex8_full_field.cpp
    )
    target_link_libraries(fuelsim_production_result_tests
        PRIVATE fuelsim_exodus_result_test_support fuelsim_field_error_metrics
    )
    fuelsim_enable_warnings(fuelsim_production_result_tests)

    function(fuelsim_stage_blackbox_case output_variable case_name input_name)
        set(case_root "${CMAKE_CURRENT_BINARY_DIR}/blackbox/${case_name}/verification")
        set(case_directory "${case_root}/fuelsim")
        file(MAKE_DIRECTORY "${case_directory}")
        foreach(reference_directory IN ITEMS moose abaqus meshes)
            if(NOT EXISTS "${case_root}/${reference_directory}")
                file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/verification/${reference_directory}"
                    "${case_root}/${reference_directory}" SYMBOLIC COPY_ON_ERROR)
            endif()
        endforeach()
        configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/${input_name}"
            "${case_directory}/${input_name}"
            COPYONLY
        )
        set(${output_variable} "${case_directory}/${input_name}" PARENT_SCOPE)
    endfunction()

    function(fuelsim_add_hex20_blackbox_test test_name case_name input_name
        temperature_reference displacement_reference tolerance field_selection test_labels)
        fuelsim_stage_blackbox_case(staged_input "${case_name}" "${input_name}")
        get_filename_component(case_directory "${staged_input}" DIRECTORY)
        get_filename_component(case_stem "${input_name}" NAME_WE)
        set(results "${case_directory}/${case_stem}_results.e")
        add_test(NAME "${test_name}"
            COMMAND "${CMAKE_COMMAND}"
                "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
                "-DINPUT_FILE=${staged_input}"
                "-DOUTPUT_FILE=${results}"
                "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
                "-DCHECK_ARGUMENTS=hex20-fields;${results};${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/${temperature_reference};${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/${displacement_reference};${tolerance};${field_selection}"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        set_tests_properties("${test_name}" PROPERTIES LABELS "${test_labels}")
    endfunction()

    foreach(material_case IN ITEMS elastic_moduli elastic_moduli_finite
        plastic_yield plastic_yield_finite plastic_hardening plastic_hardening_finite
        creep_coefficient_finite creep_exponent_finite)
        set(test_name "fuelsim_cax4t_material_${material_case}_abaqus")
        set(input_stem "transient_cax4t_${material_case}")
        fuelsim_stage_blackbox_case(staged_input "${test_name}" "${input_stem}.fsi")
        get_filename_component(case_directory "${staged_input}" DIRECTORY)
        set(results "${case_directory}/${input_stem}_results.e")
        add_test(NAME "${test_name}"
            COMMAND "${CMAKE_COMMAND}"
                "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
                "-DINPUT_FILE=${staged_input}"
                "-DOUTPUT_FILE=${results}"
                "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
                "-DCHECK_ARGUMENTS=cax4t-material-temperature;${results};${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/cax4t_material_temperature_probe;${material_case}"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        set_tests_properties("${test_name}" PROPERTIES LABELS "abaqus;production;axisymmetric")
    endforeach()

    function(fuelsim_add_radial_case test_name input_stem mode)
        fuelsim_stage_blackbox_case(staged_input "${test_name}" "${input_stem}.fsi")
        get_filename_component(case_directory "${staged_input}" DIRECTORY)
        set(results "${case_directory}/${input_stem}_results.e")
        set(summary "${case_directory}/${input_stem}_summary.csv")
        add_test(NAME "${test_name}" COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DINPUT_FILE=${staged_input}"
            "-DOUTPUT_FILES=${results};${summary};${case_directory}/${input_stem}.checkpoint"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_radial>"
            "-DCHECK_ARGUMENTS=${mode};${results};${summary};${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/gps_1d"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
    endfunction()
    fuelsim_add_radial_case(fuelsim_gps_uniform_small_abaqus transient_gps_uniform_small uniform-small)
    fuelsim_add_radial_case(fuelsim_gps_uniform_finite_abaqus transient_gps_uniform_finite uniform-finite)
    fuelsim_add_radial_case(fuelsim_gps_nonuniform_finite_abaqus transient_gps_nonuniform_finite nonuniform-finite)
    fuelsim_add_radial_case(fuelsim_gps_material_temperature_finite_abaqus
        transient_gps_material_temperature_finite material-temperature-finite)
    fuelsim_add_radial_case(fuelsim_gps_contact_abaqus transient_gps_contact_small contact-small)
    fuelsim_stage_blackbox_case(gps_steady_input gps_steady steady_gps_uniform_finite.fsi)
    get_filename_component(gps_steady_directory "${gps_steady_input}" DIRECTORY)
    set(gps_steady_result "${gps_steady_directory}/steady_gps_uniform_finite_results.e")
    set(gps_steady_summary "${gps_steady_directory}/steady_gps_uniform_finite_summary.csv")
    add_test(NAME fuelsim_gps_steady_finite COMMAND "${CMAKE_COMMAND}"
        "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
        "-DINPUT_FILE=${gps_steady_input}"
        "-DOUTPUT_FILES=${gps_steady_result};${gps_steady_summary}"
        "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_radial>"
        "-DCHECK_ARGUMENTS=steady-finite;${gps_steady_result};${gps_steady_summary};${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/gps_1d"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
    set_tests_properties(fuelsim_gps_contact_abaqus PROPERTIES FIXTURES_SETUP gps_contact_reference)
    set(gps_contact_reference "${CMAKE_CURRENT_BINARY_DIR}/blackbox/fuelsim_gps_contact_abaqus/verification/fuelsim/transient_gps_contact_small_results.e")
    foreach(formulation IN ITEMS small finite)
        set(case_name "gps_chain_${formulation}")
        set(input_stem "transient_gps_chain_${formulation}")
        fuelsim_stage_blackbox_case(staged_input "${case_name}" "${input_stem}.fsi")
        get_filename_component(case_directory "${staged_input}" DIRECTORY)
        set(results "${case_directory}/${input_stem}_results.e")
        set(summary "${case_directory}/${input_stem}_summary.csv")
        add_test(NAME "fuelsim_${case_name}" COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DINPUT_FILE=${staged_input}"
            "-DOUTPUT_FILES=${results};${summary};${case_directory}/${input_stem}.checkpoint"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_radial_chain>"
            "-DCHECK_ARGUMENTS=${formulation};${results};${summary}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        if(formulation STREQUAL "small")
            set_tests_properties(fuelsim_gps_chain_small PROPERTIES FIXTURES_SETUP gps_small_chain_reference)
            add_test(NAME fuelsim_gps_chain_abaqus COMMAND "$<TARGET_FILE:fuelsim_production_radial>"
                chain-small "${results}" "${summary}" "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/gps_1d")
            set_tests_properties(fuelsim_gps_chain_abaqus PROPERTIES FIXTURES_REQUIRED gps_small_chain_reference)
        endif()
    endforeach()
    set_tests_properties(fuelsim_gps_chain_finite PROPERTIES FIXTURES_SETUP gps_finite_reference)
    fuelsim_stage_blackbox_case(gps_inelastic_input gps_inelastic transient_gps_inelastic_finite.fsi)
    get_filename_component(gps_inelastic_directory "${gps_inelastic_input}" DIRECTORY)
    set(gps_inelastic_result "${gps_inelastic_directory}/transient_gps_inelastic_finite_results.e")
    set(gps_inelastic_summary "${gps_inelastic_directory}/transient_gps_inelastic_finite_summary.csv")
    add_test(NAME fuelsim_gps_inelastic_finite COMMAND "${CMAKE_COMMAND}"
        "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
        "-DINPUT_FILE=${gps_inelastic_input}"
        "-DOUTPUT_FILES=${gps_inelastic_result};${gps_inelastic_summary};${gps_inelastic_directory}/transient_gps_inelastic_finite.checkpoint"
        "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_radial_inelastic>"
        "-DCHECK_ARGUMENTS=${gps_inelastic_result};${gps_inelastic_summary}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
    fuelsim_stage_blackbox_case(gps_checkpoint_input gps_restart transient_gps_contact_checkpoint.fsi)
    fuelsim_stage_blackbox_case(gps_restart_input gps_restart transient_gps_contact_restart.fsi)
    get_filename_component(gps_restart_directory "${gps_restart_input}" DIRECTORY)
    set(gps_restart_result "${gps_restart_directory}/transient_gps_contact_restart_results.part1.e")
    add_test(NAME fuelsim_gps_contact_restart COMMAND "${CMAKE_COMMAND}"
        "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
        "-DCHECKPOINT_INPUT=${gps_checkpoint_input}"
        "-DCHECKPOINT_FILE=${gps_restart_directory}/transient_gps_contact_checkpoint.checkpoint"
        "-DINPUT_FILE=${gps_restart_input}"
        "-DOUTPUT_FILES=${gps_restart_result};${gps_restart_directory}/transient_gps_contact_restart.checkpoint"
        "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_radial_equivalence>"
        "-DCHECK_ARGUMENTS=${gps_contact_reference};${gps_restart_result}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_restart.cmake")
    set_tests_properties(fuelsim_gps_contact_restart PROPERTIES FIXTURES_REQUIRED gps_contact_reference)

    # Test scheduling only: every physical input card is copied byte-for-byte.
    function(fuelsim_add_production_test test_name case_name input_stem check_mode)
        fuelsim_stage_blackbox_case(staged_input "${case_name}" "${input_stem}.fsi")
        get_filename_component(case_directory "${staged_input}" DIRECTORY)
        set(results "${case_directory}/${input_stem}_results.e")
        set(summary "${case_directory}/${input_stem}_summary.csv")
        set(output_files "${results};${summary}")
        if(check_mode STREQUAL "hex8-inelastic" OR check_mode STREQUAL "b34-transient" OR check_mode STREQUAL "b40"
            OR check_mode STREQUAL "b549" OR check_mode STREQUAL "b550" OR check_mode STREQUAL "b555"
            OR check_mode STREQUAL "b531" OR check_mode STREQUAL "b534" OR check_mode STREQUAL "rz-sliding-abaqus"
            OR input_stem STREQUAL "transient_b510_small_j2" OR input_stem STREQUAL "transient_b515_finite_j2")
            list(APPEND output_files "${case_directory}/${input_stem}.checkpoint")
        endif()
        add_test(NAME "${test_name}"
            COMMAND "${CMAKE_COMMAND}"
                "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
                "-DINPUT_FILE=${staged_input}"
                "-DOUTPUT_FILES=${output_files}"
                "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
                "-DCHECK_ARGUMENTS=${check_mode};${results};${summary};${ARGN}"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        set_tests_properties("${test_name}" PROPERTIES LABELS "acceptance;production" TIMEOUT 120)
    endfunction()

    add_executable(fuelsim_hex8_result_fields_tests tests/hex8_result_fields_tests.cpp)
    target_link_libraries(fuelsim_hex8_result_fields_tests PRIVATE fuelsim::io)
    target_include_directories(fuelsim_hex8_result_fields_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    fuelsim_enable_warnings(fuelsim_hex8_result_fields_tests)
    add_test(NAME fuelsim_hex8_result_fields_tests COMMAND fuelsim_hex8_result_fields_tests)
    set_tests_properties(fuelsim_hex8_result_fields_tests PROPERTIES LABELS "unit;io;hex8;analytic")
    add_executable(fuelsim_test_mesh_export EXCLUDE_FROM_ALL benchmarks/test_mesh_export.cpp)
    target_link_libraries(fuelsim_test_mesh_export PRIVATE fuelsim::io)
    fuelsim_enable_warnings(fuelsim_test_mesh_export)

    add_library(fuelsim_abaqus_hex8_test_support
        tests/support/abaqus_hex8_full_field.cpp
    )
    target_include_directories(fuelsim_abaqus_hex8_test_support
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    )
    target_link_libraries(fuelsim_abaqus_hex8_test_support
        PUBLIC fuelsim_field_error_metrics fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_abaqus_hex8_test_support)

    add_executable(fuelsim_input_tests tests/input_tests.cpp)
    target_link_libraries(fuelsim_input_tests PRIVATE fuelsim::io)
    fuelsim_enable_warnings(fuelsim_input_tests)
    add_test(NAME fuelsim_input_tests
        COMMAND fuelsim_input_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_c3d8rt_thermoelastic.fsi"
            "${CMAKE_CURRENT_BINARY_DIR}/malformed_input.fsi"
    )
    set_tests_properties(fuelsim_input_tests PROPERTIES LABELS "unit")

    add_executable(fuelsim_exodus_mesh_io_tests
        tests/exodus_mesh_io_tests.cpp
    )
    target_link_libraries(fuelsim_exodus_mesh_io_tests
        PRIVATE fuelsim::io
    )
    fuelsim_enable_warnings(fuelsim_exodus_mesh_io_tests)
    add_test(NAME fuelsim_exodus_mesh_io_tests
        COMMAND fuelsim_exodus_mesh_io_tests
            "${CMAKE_CURRENT_BINARY_DIR}/fuelsim_two_quad_smoke.exo"
    )
    set_tests_properties(fuelsim_exodus_mesh_io_tests PROPERTIES
        LABELS "dependency"
    )

    add_executable(fuelsim_checkpoint_results_tests
        tests/checkpoint_results_tests.cpp
    )
    target_link_libraries(fuelsim_checkpoint_results_tests
        PRIVATE
            fuelsim::io
            fuelsim::solver
            SEACASExodus::exodus
    )
    fuelsim_enable_warnings(fuelsim_checkpoint_results_tests)
    add_test(NAME fuelsim_checkpoint_results_tests
        COMMAND fuelsim_checkpoint_results_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
            "${CMAKE_CURRENT_BINARY_DIR}/m3_checkpoint.bin"
            "${CMAKE_CURRENT_BINARY_DIR}/m3_results.e"
    )
    set_tests_properties(fuelsim_checkpoint_results_tests PROPERTIES
        LABELS "integration"
    )

    add_executable(fuelsim_m3_load_boundary_tests
        tests/m3_load_boundary_tests.cpp
    )
    target_link_libraries(fuelsim_m3_load_boundary_tests
        PRIVATE
            fuelsim_field_error_metrics
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m3_load_boundary_tests)
    add_test(NAME fuelsim_m3_load_boundary_tests
        COMMAND fuelsim_m3_load_boundary_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
    )
    set_tests_properties(fuelsim_m3_load_boundary_tests PROPERTIES
        LABELS "unit;time-control;transaction"
    )

    fuelsim_add_production_test(fuelsim_pressure_cylinder_analytic_tests pressure_cylinder steady_pressure_boundary pressure-cylinder)

    add_executable(fuelsim_rz_contact_projection_contract_tests tests/rz_contact_projection_contract_tests.cpp)
    target_link_libraries(fuelsim_rz_contact_projection_contract_tests PRIVATE fuelsim::io)
    target_include_directories(fuelsim_rz_contact_projection_contract_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_rz_contact_projection_contract_tests)
    add_test(NAME fuelsim_rz_contact_projection_contract COMMAND fuelsim_rz_contact_projection_contract_tests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi")
    set_tests_properties(fuelsim_rz_contact_projection_contract PROPERTIES LABELS "unit;contact;projection")

    add_executable(fuelsim_m58_integrated_hex8_benchmark
        benchmarks/m58_integrated_hex8.cpp
    )
    target_link_libraries(fuelsim_m58_integrated_hex8_benchmark
        PRIVATE
            fuelsim_field_error_metrics
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m58_integrated_hex8_benchmark)
    add_executable(fuelsim_b60_fuel_plate_c3d8rt_benchmark
        benchmarks/b60_fuel_plate_c3d8rt.cpp
    )
    target_link_libraries(fuelsim_b60_fuel_plate_c3d8rt_benchmark
        PRIVATE fuelsim_field_error_metrics fuelsim::io fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_b60_fuel_plate_c3d8rt_benchmark)
    add_executable(fuelsim_b61_c3d20t_material_replay
        benchmarks/b61_c3d20t_material_replay.cpp
    )
    target_link_libraries(fuelsim_b61_c3d20t_material_replay
        PRIVATE fuelsim_field_error_metrics fuelsim::io fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_b61_c3d20t_material_replay)
    fuelsim_add_production_test(fuelsim_b60_c3d8rt_finite_steady_bending_tests b60_steady
        steady_b60_fuel_plate_c3d8rt_finite_bending b60 steady
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b60_fuel_plate_c3d8rt_finite_steady_bending_nodal.csv")
    set_tests_properties(fuelsim_b60_c3d8rt_finite_steady_bending_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;finite-strain;steady;thermal;mechanical;bending"
        TIMEOUT 120
    )
    fuelsim_add_production_test(fuelsim_b60_c3d8rt_finite_transient_bending_tests b60_transient
        transient_b60_fuel_plate_c3d8rt_finite_ramped_bending b60 transient
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b60_fuel_plate_c3d8rt_finite_ramped_bending_nodal.csv")
    set_tests_properties(fuelsim_b60_c3d8rt_finite_transient_bending_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;finite-strain;transient;thermal;mechanical;bending"
        TIMEOUT 120
    )
    fuelsim_add_production_test(fuelsim_b61_fuel_plate_c3d8rt_finite_inelastic_abaqus_tests b61
        transient_b61_fuel_plate_c3d8rt_finite_inelastic_bending hex8-bulk-abaqus b61
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b61_fuel_plate_c3d8rt_finite_inelastic_bending")
    set_tests_properties(fuelsim_b61_fuel_plate_c3d8rt_finite_inelastic_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;thermal;mechanical;finite-strain;bending;plasticity;creep;abaqus"
        TIMEOUT 120
    )
    add_executable(fuelsim_m58_hex8_to_hex20
        benchmarks/m58_hex8_to_hex20.cpp
    )
    target_link_libraries(fuelsim_m58_hex8_to_hex20
        PRIVATE fuelsim::io
    )
    fuelsim_enable_warnings(fuelsim_m58_hex8_to_hex20)
    add_executable(fuelsim_m58_integrated_hex20_results
        benchmarks/m58_integrated_hex20_results.cpp
    )
    target_link_libraries(fuelsim_m58_integrated_hex20_results
        PRIVATE fuelsim_field_error_metrics fuelsim::io fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m58_integrated_hex20_results)
    add_executable(fuelsim_b549_small_c3d20t_mesh
        benchmarks/b549_small_c3d20t_mesh.cpp
    )
    target_link_libraries(fuelsim_b549_small_c3d20t_mesh
        PRIVATE fuelsim::io
    )
    fuelsim_enable_warnings(fuelsim_b549_small_c3d20t_mesh)
    fuelsim_add_production_test(fuelsim_m52_transient_reference m52_transient
        transient_large_sliding_contact m52-transient)
    set_tests_properties(fuelsim_m52_transient_reference PROPERTIES FIXTURES_SETUP m52_transient_production)
    fuelsim_stage_blackbox_case(m52_checkpoint_input m52_restart transient_large_sliding_contact_checkpoint.fsi)
    fuelsim_stage_blackbox_case(m52_restart_input m52_restart transient_large_sliding_contact_restart.fsi)
    get_filename_component(m52_restart_directory "${m52_restart_input}" DIRECTORY)
    set(m52_full_directory "${CMAKE_CURRENT_BINARY_DIR}/blackbox/m52_transient/verification/fuelsim")
    set(m52_restart_result "${m52_restart_directory}/transient_large_sliding_contact_restart_results.part1.e")
    set(m52_restart_summary "${m52_restart_directory}/transient_large_sliding_contact_restart_summary.csv")
    add_test(NAME fuelsim_m52_large_sliding_restart_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DCHECKPOINT_INPUT=${m52_checkpoint_input}"
            "-DCHECKPOINT_FILE=${m52_restart_directory}/transient_large_sliding_contact_split.checkpoint"
            "-DINPUT_FILE=${m52_restart_input}"
            "-DOUTPUT_FILES=${m52_restart_result};${m52_restart_summary}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=m52-restart;${m52_restart_result};${m52_restart_summary};${m52_full_directory}/transient_large_sliding_contact_results.e"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_restart.cmake")
    set_tests_properties(fuelsim_m52_large_sliding_restart_tests PROPERTIES
        FIXTURES_REQUIRED m52_transient_production FIXTURES_SETUP m52_restart_production
        LABELS "acceptance;production;restart;contact" TIMEOUT 120)
    add_executable(fuelsim_rz_restart_step_contract_tests tests/rz_restart_step_contract_tests.cpp)
    target_link_libraries(fuelsim_rz_restart_step_contract_tests PRIVATE fuelsim::io)
    fuelsim_enable_warnings(fuelsim_rz_restart_step_contract_tests)
    add_test(NAME fuelsim_m52_restart_step_contract
        COMMAND fuelsim_rz_restart_step_contract_tests "${m52_restart_input}"
            "${m52_restart_directory}/transient_large_sliding_contact_split.checkpoint")
    set_tests_properties(fuelsim_m52_restart_step_contract PROPERTIES
        FIXTURES_REQUIRED m52_restart_production LABELS "unit;checkpoint;transaction")

    add_executable(fuelsim_m54_transient_augmented_contact_tests
        tests/m54_transient_augmented_contact_tests.cpp
    )
    target_link_libraries(fuelsim_m54_transient_augmented_contact_tests
        PRIVATE
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m54_transient_augmented_contact_tests)
    add_test(NAME fuelsim_m54_transient_augmented_contact_tests
        COMMAND fuelsim_m54_transient_augmented_contact_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
    )
    set_tests_properties(fuelsim_m54_transient_augmented_contact_tests
        PROPERTIES LABELS "integration"
    )

    add_executable(fuelsim_core_tests tests/core_tests.cpp)
    target_link_libraries(fuelsim_core_tests
        PRIVATE fuelsim_test_mesh_fixture fuelsim::core
    )
    fuelsim_enable_warnings(fuelsim_core_tests)
    add_test(NAME fuelsim_core_tests COMMAND fuelsim_core_tests)
    set_tests_properties(fuelsim_core_tests PROPERTIES LABELS "unit")


    add_executable(fuelsim_c3d20rt_probe_results tests/production_c3d20rt_probe.cpp)
    add_executable(fuelsim_c3d20rt_material_results tests/production_c3d20rt_material.cpp)
    target_link_libraries(fuelsim_c3d20rt_material_results PRIVATE
        fuelsim_exodus_result_test_support fuelsim_field_error_metrics)
    fuelsim_enable_warnings(fuelsim_c3d20rt_material_results)
    add_executable(fuelsim_c3d20rt_contact_results tests/production_c3d20rt_contact.cpp)
    target_link_libraries(fuelsim_c3d20rt_contact_results PRIVATE
        fuelsim_exodus_result_test_support fuelsim_field_error_metrics)
    fuelsim_enable_warnings(fuelsim_c3d20rt_contact_results)
    # Manual prescribed-state diagnostic; it does not solve an end-to-end case.
    add_executable(fuelsim_c3d20rt_friction_history_replay EXCLUDE_FROM_ALL
        tests/c3d20rt_friction_history_replay.cpp)
    target_include_directories(fuelsim_c3d20rt_friction_history_replay PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_link_libraries(fuelsim_c3d20rt_friction_history_replay PRIVATE fuelsim::io)
    fuelsim_enable_warnings(fuelsim_c3d20rt_friction_history_replay)
    fuelsim_stage_blackbox_case(c3d20rt_contact_input c3d20rt_contact_friction
        transient_c3d20rt_contact_friction.fsi)
    get_filename_component(c3d20rt_contact_directory "${c3d20rt_contact_input}" DIRECTORY)
    set(c3d20rt_contact_output "${c3d20rt_contact_directory}/transient_c3d20rt_contact_friction_results.e")
    add_test(NAME fuelsim_c3d20rt_contact_friction_abaqus_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DINPUT_FILE=${c3d20rt_contact_input}"
            "-DOUTPUT_FILE=${c3d20rt_contact_output}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_c3d20rt_contact_results>"
            "-DCHECK_ARGUMENTS=${c3d20rt_contact_output};${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/c3d20rt_contact_friction_nodes.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/c3d20rt_contact_friction_contact.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/c3d20rt_contact_friction_points.csv;70"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
    set_tests_properties(fuelsim_c3d20rt_contact_friction_abaqus_tests PROPERTIES
        LABELS "three-dimensional;quadratic;abaqus;contact;thermal;friction;transient" TIMEOUT 180)
    # A short zero-friction branch complements the full manual friction-reversal case.
    fuelsim_stage_blackbox_case(c3d20rt_frictionless_input c3d20rt_finite_frictionless
        transient_c3d20rt_finite_contact_frictionless.fsi)
    get_filename_component(c3d20rt_frictionless_directory "${c3d20rt_frictionless_input}" DIRECTORY)
    set(c3d20rt_frictionless_output "${c3d20rt_frictionless_directory}/transient_c3d20rt_finite_contact_frictionless_results.e")
    add_test(NAME fuelsim_c3d20rt_finite_frictionless_abaqus_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DINPUT_FILE=${c3d20rt_frictionless_input}"
            "-DOUTPUT_FILE=${c3d20rt_frictionless_output}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_c3d20rt_contact_results>"
            "-DCHECK_ARGUMENTS=${c3d20rt_frictionless_output};${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/c3d20rt_finite_contact_frictionless_short_nodes.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/c3d20rt_finite_contact_frictionless_short_contact.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/c3d20rt_finite_contact_frictionless_short_points.csv;80;8;--qualified-finite-frictionless"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
    set_tests_properties(fuelsim_c3d20rt_finite_frictionless_abaqus_tests PROPERTIES
        LABELS "three-dimensional;quadratic;abaqus;contact;thermal;transient" TIMEOUT 180)
    foreach(c3d20rt_material_case IN ITEMS small_coupled finite_coupled small_plastic finite_plastic finite_creep small_creep_hold)
        string(REGEX REPLACE "^(small|finite)_" "" c3d20rt_material_mode "${c3d20rt_material_case}")
        set(c3d20rt_material_name "c3d20rt_${c3d20rt_material_case}")
        fuelsim_stage_blackbox_case(c3d20rt_material_input "${c3d20rt_material_name}"
            "transient_${c3d20rt_material_name}.fsi")
        get_filename_component(c3d20rt_material_directory "${c3d20rt_material_input}" DIRECTORY)
        set(c3d20rt_material_output "${c3d20rt_material_directory}/transient_${c3d20rt_material_name}_results.e")
        add_test(NAME "fuelsim_${c3d20rt_material_name}_abaqus_tests"
            COMMAND "${CMAKE_COMMAND}"
                "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
                "-DINPUT_FILE=${c3d20rt_material_input}"
                "-DOUTPUT_FILE=${c3d20rt_material_output}"
                "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_c3d20rt_material_results>"
                "-DCHECK_ARGUMENTS=${c3d20rt_material_output};${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${c3d20rt_material_name}_nodes.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${c3d20rt_material_name}_points.csv;${c3d20rt_material_mode}"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        set_tests_properties("fuelsim_${c3d20rt_material_name}_abaqus_tests" PROPERTIES
            LABELS "three-dimensional;quadratic;abaqus;inelastic;transient")
    endforeach()
    target_link_libraries(fuelsim_c3d20rt_probe_results PRIVATE
        fuelsim_exodus_result_test_support fuelsim_field_error_metrics)
    fuelsim_enable_warnings(fuelsim_c3d20rt_probe_results)
    foreach(c3d20rt_probe_name IN ITEMS c3d20rt_thermal_probe c3d20rt_nonaffine_probe c3d20rt_small_probe)
        fuelsim_stage_blackbox_case(c3d20rt_probe_input ${c3d20rt_probe_name} transient_${c3d20rt_probe_name}.fsi)
        get_filename_component(c3d20rt_probe_directory "${c3d20rt_probe_input}" DIRECTORY)
        set(c3d20rt_probe_results "${c3d20rt_probe_directory}/transient_${c3d20rt_probe_name}_results.e")
        add_test(NAME fuelsim_${c3d20rt_probe_name}_abaqus_tests
            COMMAND "${CMAKE_COMMAND}"
                "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
                "-DINPUT_FILE=${c3d20rt_probe_input}"
                "-DOUTPUT_FILE=${c3d20rt_probe_results}"
                "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_c3d20rt_probe_results>"
                "-DCHECK_ARGUMENTS=${c3d20rt_probe_results};${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${c3d20rt_probe_name}_nodes.csv"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        set_tests_properties(fuelsim_${c3d20rt_probe_name}_abaqus_tests PROPERTIES
            LABELS "three-dimensional;quadratic;abaqus;thermal;transient")
    endforeach()


    add_executable(fuelsim_hex8_solver_tests tests/hex8_solver_tests.cpp)
    target_link_libraries(fuelsim_hex8_solver_tests PRIVATE fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_hex8_solver_tests)
    add_test(NAME fuelsim_hex8_solver_tests
        COMMAND fuelsim_hex8_solver_tests
            "${CMAKE_CURRENT_BINARY_DIR}/hex8_steady_results.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex8_transient_results.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex8_checkpoint.bin")
    set_tests_properties(fuelsim_hex8_solver_tests PROPERTIES LABELS "integration;three-dimensional")

    add_executable(fuelsim_b38_hex8_sts_identification_tests
        tests/b38_hex8_sts_identification_tests.cpp)
    target_link_libraries(fuelsim_b38_hex8_sts_identification_tests
        PRIVATE fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b38_hex8_sts_identification_tests)
    add_test(NAME fuelsim_b38_hex8_sts_identification_tests
        COMMAND fuelsim_b38_hex8_sts_identification_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b38_hex8_sts_operator.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b38_hex8_sts_contact_summary.txt")
    set_tests_properties(fuelsim_b38_hex8_sts_identification_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;abaqus;identification")

    add_executable(fuelsim_b49_hex8_c3d8t_operator_abaqus_tests
        tests/b49_hex8_c3d8t_operator_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b49_hex8_c3d8t_operator_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b49_hex8_c3d8t_operator_abaqus_tests)
    add_test(NAME fuelsim_b49_hex8_c3d8t_operator_abaqus_tests
        COMMAND fuelsim_b49_hex8_c3d8t_operator_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b49_hex8_c3d8t_operator_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b49_hex8_c3d8t_operator_integration_points.csv")
    set_tests_properties(fuelsim_b49_hex8_c3d8t_operator_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;thermo-mechanical;abaqus;identification")

    add_executable(fuelsim_b50_hex8_c3d8t_capacity_abaqus_tests
        tests/b50_hex8_c3d8t_capacity_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b50_hex8_c3d8t_capacity_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b50_hex8_c3d8t_capacity_abaqus_tests)
    add_test(NAME fuelsim_b50_hex8_c3d8t_capacity_abaqus_tests
        COMMAND fuelsim_b50_hex8_c3d8t_capacity_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b50_hex8_c3d8t_capacity.csv")
    set_tests_properties(fuelsim_b50_hex8_c3d8t_capacity_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;abaqus;identification")

    add_executable(fuelsim_b51_hex8_c3d8t_finite_heat_abaqus_tests
        tests/b51_hex8_c3d8t_finite_heat_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b51_hex8_c3d8t_finite_heat_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b51_hex8_c3d8t_finite_heat_abaqus_tests)
    add_test(NAME fuelsim_b51_hex8_c3d8t_finite_heat_abaqus_tests
        COMMAND fuelsim_b51_hex8_c3d8t_finite_heat_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b51_hex8_c3d8t_finite_heat_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b51_hex8_c3d8t_finite_heat_integration_points.csv")
    set_tests_properties(fuelsim_b51_hex8_c3d8t_finite_heat_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;finite-strain;thermal;abaqus;identification")

    add_executable(fuelsim_b52_hex8_c3d8t_thermal_contact_abaqus_tests
        tests/b52_hex8_c3d8t_thermal_contact_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b52_hex8_c3d8t_thermal_contact_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b52_hex8_c3d8t_thermal_contact_abaqus_tests)
    add_test(NAME fuelsim_b52_hex8_c3d8t_thermal_contact_abaqus_tests
        COMMAND fuelsim_b52_hex8_c3d8t_thermal_contact_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b52_hex8_c3d8t_thermal_contact_nodal.csv")
    set_tests_properties(fuelsim_b52_hex8_c3d8t_thermal_contact_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;thermal;contact;abaqus")

    add_executable(fuelsim_b520_hex8_c3d8t_gap_conductance_abaqus_tests
        tests/b520_hex8_c3d8t_gap_conductance_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b520_hex8_c3d8t_gap_conductance_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b520_hex8_c3d8t_gap_conductance_abaqus_tests)
    add_test(NAME fuelsim_b520_hex8_c3d8t_gap_conductance_abaqus_tests
        COMMAND fuelsim_b520_hex8_c3d8t_gap_conductance_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b520_hex8_c3d8t_gap_conductance_clearance_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b520_hex8_c3d8t_gap_conductance_pressure_nodal.csv")
    set_tests_properties(fuelsim_b520_hex8_c3d8t_gap_conductance_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;thermo-mechanical;contact;abaqus;identification")

    add_executable(fuelsim_b521_hex8_c3d8t_thermal_contact_path_abaqus_tests
        tests/b521_hex8_c3d8t_thermal_contact_path_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b521_hex8_c3d8t_thermal_contact_path_abaqus_tests
        PRIVATE fuelsim::io fuelsim_field_error_metrics)
    fuelsim_enable_warnings(fuelsim_b521_hex8_c3d8t_thermal_contact_path_abaqus_tests)
    add_test(NAME fuelsim_b521_hex8_c3d8t_thermal_contact_path_abaqus_tests
        COMMAND fuelsim_b521_hex8_c3d8t_thermal_contact_path_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b521_hex8_c3d8t_thermal_contact_path_nodal.csv"
            "${CMAKE_CURRENT_BINARY_DIR}/b521_hex8_c3d8t_thermal_contact_checkpoint.bin")
    set_tests_properties(fuelsim_b521_hex8_c3d8t_thermal_contact_path_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;thermo-mechanical;contact;abaqus")

    add_executable(fuelsim_b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests
        tests/b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests
        PRIVATE fuelsim::io fuelsim_field_error_metrics)
    fuelsim_enable_warnings(fuelsim_b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests)
    add_test(NAME fuelsim_b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests
        COMMAND fuelsim_b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_f6_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_f12_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_swapped_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_tight_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_f6_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_f12_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_swapped_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b522_hex8_c3d8t_faceted_thermal_contact_tight_contact.csv")
    set_tests_properties(fuelsim_b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;thermo-mechanical;contact;curved;abaqus")

    add_executable(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests EXCLUDE_FROM_ALL
        benchmarks/abaqus_hex8_contact_scans.cpp)
    target_include_directories(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    target_link_libraries(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests
        PRIVATE fuelsim_abaqus_hex8_test_support fuelsim_field_error_metrics fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests)
    fuelsim_add_production_test(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests b523
        transient_b523_integrated_contact hex8-contact-path b523
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path")
    set_tests_properties(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;finite-strain;plasticity;creep;friction;contact;abaqus")

    foreach(b526_case IN ITEMS b526_contact_cycle b526_friction_reversal)
        string(REPLACE "b526_" "" b526_suffix "${b526_case}")
        set(b526_test "fuelsim_b526_hex8_c3d8t_${b526_suffix}_abaqus_tests")
        fuelsim_add_production_test(${b526_test} ${b526_case}
            transient_${b526_case} hex8-contact-path ${b526_case}
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${b526_case}")
        set_tests_properties(${b526_test} PROPERTIES
            LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;finite-strain;friction;contact;abaqus;transition")
    endforeach()

    foreach(b53_contact_case IN ITEMS b538_contact_cycle b539_friction_reversal)
        string(REPLACE "b538_" "" b53_contact_suffix "${b53_contact_case}")
        string(REPLACE "b539_" "" b53_contact_suffix "${b53_contact_suffix}")
        set(b53_contact_test "fuelsim_${b53_contact_case}_c3d8rt_abaqus_tests")
        fuelsim_add_production_test(${b53_contact_test} ${b53_contact_case}
            transient_${b53_contact_case} hex8-contact-path ${b53_contact_case}
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${b53_contact_case}")
        set_tests_properties(${b53_contact_test} PROPERTIES
            LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;thermal;mechanical;finite-strain;friction;contact;abaqus;transition")
    endforeach()

    fuelsim_add_production_test(fuelsim_b540_c3d8rt_nonmatching_contact_abaqus_tests b540
        transient_b540_nonmatching_contact_cycle hex8-contact-path b540_nonmatching_contact_cycle
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b540_nonmatching_contact_cycle")
    set_tests_properties(fuelsim_b540_c3d8rt_nonmatching_contact_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;thermal;mechanical;finite-strain;friction;contact;nonmatching;abaqus;transition")

    fuelsim_add_production_test(fuelsim_b544_hex8_c3d8rt_distorted_bending_abaqus_tests b544
        transient_b544_distorted_bending hex8-bulk-abaqus b544
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b544_hex8_c3d8rt_distorted_bending")
    set_tests_properties(fuelsim_b544_hex8_c3d8rt_distorted_bending_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;thermal;mechanical;finite-strain;bending;distorted-mesh;temperature-dependent;abaqus")

    add_executable(fuelsim_b527_hex8_engineering_fuel_clad_tests EXCLUDE_FROM_ALL
        benchmarks/abaqus_hex8_engineering_scans.cpp)
    target_include_directories(fuelsim_b527_hex8_engineering_fuel_clad_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    target_link_libraries(fuelsim_b527_hex8_engineering_fuel_clad_tests
        PRIVATE fuelsim_abaqus_hex8_test_support fuelsim_field_error_metrics fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b527_hex8_engineering_fuel_clad_tests)
    fuelsim_add_production_test(fuelsim_b527_medium_abaqus_tests b527
        transient_b527_engineering_contact hex8-contact-path b527_medium
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b527_medium")
    set_tests_properties(fuelsim_b527_medium_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;finite-strain;friction;contact;engineering;abaqus")

    add_executable(fuelsim_b552_hex8_c3d8t_sts_cross_face_abaqus_tests
        tests/b552_hex8_sts_cross_face_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b552_hex8_c3d8t_sts_cross_face_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b552_hex8_c3d8t_sts_cross_face_abaqus_tests)
    add_test(NAME fuelsim_b552_hex8_c3d8t_sts_cross_face_abaqus_tests
        COMMAND fuelsim_b552_hex8_c3d8t_sts_cross_face_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus")
    set_tests_properties(fuelsim_b552_hex8_c3d8t_sts_cross_face_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8t;steady;mechanical;finite-strain;contact;surface-to-surface;finite-sliding;operator-identification;abaqus")

    add_executable(fuelsim_b53_hex8_c3d8t_thermal_load_abaqus_tests
        tests/b53_hex8_c3d8t_thermal_load_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b53_hex8_c3d8t_thermal_load_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b53_hex8_c3d8t_thermal_load_abaqus_tests)
    add_test(NAME fuelsim_b53_hex8_c3d8t_thermal_load_abaqus_tests
        COMMAND fuelsim_b53_hex8_c3d8t_thermal_load_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b53_hex8_c3d8t_thermal_load.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b53_hex8_c3d8t_thermal_load_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b53_hex8_c3d8t_surface_geometry.csv")
    set_tests_properties(fuelsim_b53_hex8_c3d8t_thermal_load_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;thermal;abaqus;identification")

    add_executable(fuelsim_b54_hex8_c3d8t_finite_thermal_load_abaqus_tests
        tests/b54_hex8_c3d8t_finite_thermal_load_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b54_hex8_c3d8t_finite_thermal_load_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b54_hex8_c3d8t_finite_thermal_load_abaqus_tests)
    add_test(NAME fuelsim_b54_hex8_c3d8t_finite_thermal_load_abaqus_tests
        COMMAND fuelsim_b54_hex8_c3d8t_finite_thermal_load_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b54_hex8_c3d8t_finite_capacity.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b54_hex8_c3d8t_finite_thermal_load.csv")
    set_tests_properties(fuelsim_b54_hex8_c3d8t_finite_thermal_load_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;finite-strain;transient;thermal;abaqus;identification")

    fuelsim_add_production_test(fuelsim_b55_hex8_c3d8t_transient_abaqus_tests b55
            transient_b55_hex8_c3d8t_thermal b55
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b55_hex8_c3d8t_transient_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b55_hex8_c3d8t_transient_integration.csv")
    set_tests_properties(fuelsim_b55_hex8_c3d8t_transient_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;abaqus")

    add_test(NAME fuelsim_b56_hex8_c3d8t_temperature_operator_abaqus_tests
        COMMAND fuelsim_b49_hex8_c3d8t_operator_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b56_hex8_c3d8t_temperature_operator_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b56_hex8_c3d8t_temperature_operator_integration.csv"
            temperature_dependent)
    set_tests_properties(fuelsim_b56_hex8_c3d8t_temperature_operator_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;thermal;mechanical;temperature-dependent;abaqus;identification")

    add_executable(fuelsim_b57_hex8_c3d8t_temperature_capacity_abaqus_tests
        tests/b57_hex8_c3d8t_temperature_capacity_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b57_hex8_c3d8t_temperature_capacity_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b57_hex8_c3d8t_temperature_capacity_abaqus_tests)
    add_test(NAME fuelsim_b57_hex8_c3d8t_temperature_capacity_abaqus_tests
        COMMAND fuelsim_b57_hex8_c3d8t_temperature_capacity_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b57_hex8_c3d8t_temperature_capacity.csv")
    set_tests_properties(fuelsim_b57_hex8_c3d8t_temperature_capacity_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;temperature-dependent;abaqus;identification")

    fuelsim_add_production_test(fuelsim_b58_hex8_c3d8t_multistep_abaqus_tests b58
            transient_b58_hex8_c3d8t_multistep b58
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b58_hex8_c3d8t_multistep_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b58_hex8_c3d8t_multistep_integration.csv")
    set_tests_properties(fuelsim_b58_hex8_c3d8t_multistep_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;temperature-dependent;abaqus")

    add_executable(fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests
        tests/b59_hex8_c3d8t_multimaterial_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests)
    fuelsim_stage_blackbox_case(b59_coarse_input b59 transient_b59_coarse.fsi)
    fuelsim_stage_blackbox_case(b59_distorted_input b59 transient_b59_distorted.fsi)
    get_filename_component(b59_case_directory "${b59_coarse_input}" DIRECTORY)
    add_test(NAME fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DCASE_DIRECTORY=${b59_case_directory}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DINTERNAL_COMMAND=$<TARGET_FILE:fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests>"
            "-DREFERENCE_DIRECTORY=${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_b59.cmake")
    set_tests_properties(fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;multimaterial;distorted;abaqus")

    fuelsim_add_production_test(fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests b510
        transient_b510_small_j2 hex8-inelastic-abaqus plastic
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b510_hex8_c3d8t_small_j2_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b510_hex8_c3d8t_small_j2_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b510_hex8_c3d8t_small_j2_energy.csv")
    set_tests_properties(fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests PROPERTIES
        FIXTURES_SETUP b510_production
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;abaqus")

    fuelsim_add_production_test(fuelsim_b512_hex8_c3d8t_small_coupled_abaqus_tests b512
        transient_b512_small_coupled hex8-inelastic-abaqus coupled
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b512_hex8_c3d8t_small_coupled_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b512_hex8_c3d8t_small_coupled_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b512_hex8_c3d8t_small_coupled_energy.csv")
    set_tests_properties(fuelsim_b512_hex8_c3d8t_small_coupled_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;creep;abaqus")

    fuelsim_add_production_test(fuelsim_b513_hex8_c3d8t_small_noncoaxial_abaqus_tests b513
        transient_b513_small_noncoaxial hex8-inelastic-abaqus noncoaxial
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b513_hex8_c3d8t_small_noncoaxial_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b513_hex8_c3d8t_small_noncoaxial_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b513_hex8_c3d8t_small_noncoaxial_energy.csv")
    set_tests_properties(fuelsim_b513_hex8_c3d8t_small_noncoaxial_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;creep;noncoaxial;reversal;abaqus")

    add_executable(fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests
        tests/b514_hex8_c3d8t_finite_elastic_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests)
    add_test(NAME fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests
        COMMAND fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b514_hex8_c3d8t_finite_elastic_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b514_hex8_c3d8t_finite_elastic_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b514_hex8_c3d8t_finite_elastic_energy.csv")
    set_tests_properties(fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;elasticity;finite-strain;rotation;abaqus")

    fuelsim_add_production_test(fuelsim_b515_hex8_c3d8t_finite_j2_abaqus_tests b515
        transient_b515_finite_j2 hex8-inelastic-abaqus finite_plastic
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b515_hex8_c3d8t_finite_j2_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b515_hex8_c3d8t_finite_j2_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b515_hex8_c3d8t_finite_j2_energy.csv")
    set_tests_properties(fuelsim_b515_hex8_c3d8t_finite_j2_abaqus_tests PROPERTIES
        FIXTURES_SETUP b515_production
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;finite-strain;reversal;abaqus")

    fuelsim_add_production_test(fuelsim_b517_hex8_c3d8t_finite_coupled_abaqus_tests b517
        transient_b517_finite_coupled hex8-inelastic-abaqus finite_coupled
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b517_hex8_c3d8t_finite_coupled_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b517_hex8_c3d8t_finite_coupled_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b517_hex8_c3d8t_finite_coupled_energy.csv")
    set_tests_properties(fuelsim_b517_hex8_c3d8t_finite_coupled_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;creep;finite-strain;abaqus")

    fuelsim_add_production_test(fuelsim_b535_hex8_c3d8rt_finite_j2_abaqus_tests b535
        transient_b535_finite_c3d8rt_j2 hex8-inelastic-abaqus finite_reduced_plastic
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b535_hex8_c3d8rt_finite_j2_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b535_hex8_c3d8rt_finite_j2_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b535_hex8_c3d8rt_finite_j2_energy.csv")
    set_tests_properties(fuelsim_b535_hex8_c3d8rt_finite_j2_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;mechanical;plasticity;finite-strain;reversal;abaqus")

    fuelsim_add_production_test(fuelsim_b537_hex8_c3d8rt_finite_coupled_abaqus_tests b537
        transient_b537_finite_c3d8rt_coupled hex8-inelastic-abaqus finite_reduced_coupled
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b537_hex8_c3d8rt_finite_coupled_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b537_hex8_c3d8rt_finite_coupled_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b537_hex8_c3d8rt_finite_coupled_energy.csv")
    set_tests_properties(fuelsim_b537_hex8_c3d8rt_finite_coupled_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;mechanical;plasticity;creep;finite-strain;abaqus")

    fuelsim_add_production_test(fuelsim_b541_hex8_c3d8rt_small_j2_abaqus_tests b541
        transient_b541_small_c3d8rt_j2 hex8-inelastic-abaqus reduced_plastic
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b541_hex8_c3d8rt_small_j2_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b541_hex8_c3d8rt_small_j2_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b541_hex8_c3d8rt_small_j2_energy.csv")
    set_tests_properties(fuelsim_b541_hex8_c3d8rt_small_j2_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;mechanical;plasticity;small-strain;reversal;abaqus")

    fuelsim_add_production_test(fuelsim_b543_hex8_c3d8rt_small_coupled_abaqus_tests b543
        transient_b543_small_c3d8rt_coupled hex8-inelastic-abaqus reduced_coupled
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b543_hex8_c3d8rt_small_coupled_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b543_hex8_c3d8rt_small_coupled_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b543_hex8_c3d8rt_small_coupled_energy.csv")
    set_tests_properties(fuelsim_b543_hex8_c3d8rt_small_coupled_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;mechanical;plasticity;creep;small-strain;abaqus")

    fuelsim_add_production_test(fuelsim_b545_hex8_c3d8rt_finite_noncoaxial_abaqus_tests b545
        transient_b545_finite_c3d8rt_noncoaxial hex8-inelastic-abaqus finite_reduced_noncoaxial
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b545_hex8_c3d8rt_finite_noncoaxial_100step_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b545_hex8_c3d8rt_finite_noncoaxial_100step_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b545_hex8_c3d8rt_finite_noncoaxial_100step_energy.csv")
    set_tests_properties(fuelsim_b545_hex8_c3d8rt_finite_noncoaxial_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;mechanical;plasticity;creep;finite-strain;noncoaxial;large-rotation;abaqus")

    fuelsim_add_production_test(fuelsim_b518_hex8_c3d8t_finite_noncoaxial_abaqus_tests b518
        transient_b518_finite_noncoaxial hex8-inelastic-abaqus finite_noncoaxial
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b518_hex8_c3d8t_finite_noncoaxial_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b518_hex8_c3d8t_finite_noncoaxial_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b518_hex8_c3d8t_finite_noncoaxial_energy.csv")
    set_tests_properties(fuelsim_b518_hex8_c3d8t_finite_noncoaxial_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;creep;finite-strain;noncoaxial;reversal;abaqus")

    add_executable(fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests
        tests/b519_hex8_c3d8t_finite_selective_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests)
    add_test(NAME fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests
        COMMAND fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b519_hex8_c3d8t_finite_selective_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b519_hex8_c3d8t_finite_selective_integration.csv")
    set_tests_properties(fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;mechanical;elasticity;finite-strain;selective-integration;abaqus")

    fuelsim_add_production_test(fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests b511
        transient_b511_small_norton hex8-norton-abaqus small
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b511_hex8_c3d8t_small_norton_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b511_hex8_c3d8t_small_norton_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b511_hex8_c3d8t_small_norton_energy.csv")
    set_tests_properties(fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;creep;abaqus")

    fuelsim_add_production_test(fuelsim_b516_hex8_c3d8t_finite_norton_abaqus_tests b516
        transient_b516_finite_norton hex8-norton-abaqus finite
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b516_hex8_c3d8t_finite_norton_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b516_hex8_c3d8t_finite_norton_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b516_hex8_c3d8t_finite_norton_energy.csv")
    set_tests_properties(fuelsim_b516_hex8_c3d8t_finite_norton_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;creep;finite-strain;abaqus")

    fuelsim_add_production_test(fuelsim_b536_hex8_c3d8rt_finite_norton_abaqus_tests b536
        transient_b536_finite_c3d8rt_norton hex8-norton-abaqus finite_reduced
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b536_hex8_c3d8rt_finite_norton_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b536_hex8_c3d8rt_finite_norton_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b536_hex8_c3d8rt_finite_norton_energy.csv")
    set_tests_properties(fuelsim_b536_hex8_c3d8rt_finite_norton_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;mechanical;creep;finite-strain;abaqus")

    fuelsim_add_production_test(fuelsim_b542_hex8_c3d8rt_small_norton_abaqus_tests b542
        transient_b542_small_c3d8rt_norton hex8-norton-abaqus reduced
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b542_hex8_c3d8rt_small_norton_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b542_hex8_c3d8rt_small_norton_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b542_hex8_c3d8rt_small_norton_energy.csv")
    set_tests_properties(fuelsim_b542_hex8_c3d8rt_small_norton_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;c3d8rt;reduced-integration;transient;mechanical;creep;small-strain;abaqus")

    add_executable(fuelsim_hex8_sts_contract_tests tests/hex8_sts_contract_tests.cpp)
    target_link_libraries(fuelsim_hex8_sts_contract_tests PRIVATE fuelsim::io fuelsim_exodus_result_test_support)
    target_include_directories(fuelsim_hex8_sts_contract_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_hex8_sts_contract_tests)
    foreach(b39_case IN ITEMS matching nonmatching tilted_gap)
        set(b39_test "fuelsim_b39_hex8_sts_${b39_case}_abaqus_tests")
        set(b39_stem "steady_hex8_sts_${b39_case}_abaqus")
        fuelsim_add_production_test(${b39_test} b39_${b39_case} ${b39_stem} planar-sts-hex8
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_${b39_case}_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_${b39_case}_force.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_${b39_case}_reaction.csv")
        set_tests_properties(${b39_test} PROPERTIES FIXTURES_SETUP "b39_${b39_case}_production"
            LABELS "acceptance;production;three-dimensional;hex8;contact;abaqus")
        set(b39_directory "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b39_${b39_case}/verification/fuelsim")
        add_test(NAME "fuelsim_b39_${b39_case}_contact_contract"
            COMMAND fuelsim_hex8_sts_contract_tests "${b39_directory}/${b39_stem}.fsi"
                "${b39_directory}/${b39_stem}_results.e")
        set_tests_properties("fuelsim_b39_${b39_case}_contact_contract" PROPERTIES
            FIXTURES_REQUIRED "b39_${b39_case}_production" LABELS "unit;contact;jacobian;conservation")
    endforeach()

    fuelsim_add_production_test(fuelsim_b40_hex8_sts_friction_abaqus_tests b40
        transient_hex8_sts_friction_abaqus b40
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b40_hex8_sts_friction_contact.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b40_hex8_sts_friction_reaction.csv")
    set_tests_properties(fuelsim_b40_hex8_sts_friction_abaqus_tests
        PROPERTIES FIXTURES_SETUP b40_production LABELS "acceptance;production;hex8;contact;friction;abaqus")
    fuelsim_stage_blackbox_case(b40_checkpoint_input b40_restart transient_hex8_sts_friction_checkpoint.fsi)
    fuelsim_stage_blackbox_case(b40_restart_input b40_restart transient_hex8_sts_friction_restart.fsi)
    get_filename_component(b40_restart_directory "${b40_restart_input}" DIRECTORY)
    set(b40_full_directory "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b40/verification/fuelsim")
    set(b40_restart_result "${b40_restart_directory}/transient_hex8_sts_friction_restart_results.part1.e")
    set(b40_restart_summary "${b40_restart_directory}/transient_hex8_sts_friction_restart_summary.csv")
    add_test(NAME fuelsim_b40_production_restart
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DCHECKPOINT_INPUT=${b40_checkpoint_input}"
            "-DCHECKPOINT_FILE=${b40_restart_directory}/transient_hex8_sts_friction_split.checkpoint"
            "-DINPUT_FILE=${b40_restart_input}"
            "-DOUTPUT_FILES=${b40_restart_result};${b40_restart_summary};${b40_restart_directory}/transient_hex8_sts_friction_restart.checkpoint"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=b40-restart;${b40_restart_result};${b40_restart_summary};${b40_full_directory}/transient_hex8_sts_friction_abaqus_results.e"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_restart.cmake")
    set_tests_properties(fuelsim_b40_production_restart PROPERTIES FIXTURES_REQUIRED b40_production
        FIXTURES_SETUP b40_restart_production LABELS "acceptance;production;restart;friction" TIMEOUT 120)
    add_executable(fuelsim_hex8_biaxial_friction_contract_tests tests/hex8_biaxial_friction_contract_tests.cpp)
    target_link_libraries(fuelsim_hex8_biaxial_friction_contract_tests PRIVATE fuelsim::io)
    fuelsim_enable_warnings(fuelsim_hex8_biaxial_friction_contract_tests)
    add_test(NAME fuelsim_b40_friction_history_contract
        COMMAND fuelsim_hex8_biaxial_friction_contract_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_sts_friction_abaqus.fsi"
            "${b40_full_directory}/transient_hex8_sts_friction_abaqus.checkpoint"
            "${b40_restart_directory}/transient_hex8_sts_friction_restart.checkpoint"
            "${b40_restart_directory}/transient_hex8_sts_friction_split.checkpoint")
    set_tests_properties(fuelsim_b40_friction_history_contract PROPERTIES
        FIXTURES_REQUIRED "b40_production;b40_restart_production" LABELS "unit;hex8;friction;jacobian;checkpoint")

    add_executable(fuelsim_b41_hex8_finite_sliding_abaqus_tests
        tests/b41_hex8_finite_sliding_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b41_hex8_finite_sliding_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b41_hex8_finite_sliding_abaqus_tests)
    add_test(NAME fuelsim_b41_hex8_finite_sliding_abaqus_tests
        COMMAND fuelsim_b41_hex8_finite_sliding_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b41_hex8_finite_sliding_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b41_hex8_finite_sliding_small.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b41_hex8_finite_sliding_finite.csv")
    set_tests_properties(fuelsim_b41_hex8_finite_sliding_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;friction;finite-sliding;abaqus")
    add_executable(fuelsim_b42_hex8_friction_objectivity_abaqus_tests
        tests/b42_hex8_friction_objectivity_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b42_hex8_friction_objectivity_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b42_hex8_friction_objectivity_abaqus_tests)
    add_test(NAME fuelsim_b42_hex8_friction_objectivity_abaqus_tests
        COMMAND fuelsim_b42_hex8_friction_objectivity_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b40_hex8_sts_friction_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b42_hex8_sts_friction_objectivity.csv")
    set_tests_properties(fuelsim_b42_hex8_friction_objectivity_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;friction;finite-sliding;objectivity;abaqus")

    fuelsim_add_production_test(fuelsim_b43_hex8_finite_strain_contact_abaqus_tests b43
        transient_hex8_finite_strain_contact_abaqus hex8-finite-contact
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b43_hex8_finite_strain_contact_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b43_hex8_finite_strain_contact.csv")
    set_tests_properties(fuelsim_b43_hex8_finite_strain_contact_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;finite-sliding;finite-strain;abaqus")

    add_executable(fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests
        tests/b44_hex8_nonmatching_partial_contact_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests)
    add_test(NAME fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests
        COMMAND fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b44_hex8_nonmatching_partial_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b44_hex8_nonmatching_partial_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b44_hex8_nonmatching_partial_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b44_hex8_nonmatching_partial_swapped_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b44_hex8_nonmatching_partial_swapped_contact.csv")
    set_tests_properties(fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;finite-sliding;finite-strain;nonmatching;abaqus")

    add_executable(fuelsim_b45_hex8_nonplanar_sts_abaqus_tests
        tests/b45_hex8_nonplanar_sts_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b45_hex8_nonplanar_sts_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::io)
    fuelsim_enable_warnings(fuelsim_b45_hex8_nonplanar_sts_abaqus_tests)
    add_test(NAME fuelsim_b45_hex8_nonplanar_sts_abaqus_tests
        COMMAND fuelsim_b45_hex8_nonplanar_sts_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_faceted_cylinder_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_faceted_cylinder_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_faceted_cylinder_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_warped_bilinear_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_warped_bilinear_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_warped_bilinear_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_warped_bilinear_linear_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_warped_bilinear_linear_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_warped_bilinear_quadratic_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b45_hex8_warped_bilinear_quadratic_contact.csv")
    set_tests_properties(fuelsim_b45_hex8_nonplanar_sts_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;finite-sliding;curved;warped-geometry;abaqus")

    add_executable(fuelsim_b46_hex8_release_recontact_abaqus_tests
        tests/b46_hex8_release_recontact_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b46_hex8_release_recontact_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::io)
    fuelsim_enable_warnings(fuelsim_b46_hex8_release_recontact_abaqus_tests)
    add_test(NAME fuelsim_b46_hex8_release_recontact_abaqus_tests
        COMMAND fuelsim_b46_hex8_release_recontact_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b46_hex8_release_recontact_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b46_hex8_release_recontact_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b46_hex8_release_recontact_contact.csv")
    set_tests_properties(fuelsim_b46_hex8_release_recontact_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;finite-sliding;release;recontact;abaqus")

    add_executable(fuelsim_b47_hex8_scale_parameter_abaqus_tests
        tests/b47_hex8_scale_parameter_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b47_hex8_scale_parameter_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::io)
    fuelsim_enable_warnings(fuelsim_b47_hex8_scale_parameter_abaqus_tests)
    add_test(NAME fuelsim_b47_hex8_scale_parameter_abaqus_tests
        COMMAND fuelsim_b47_hex8_scale_parameter_abaqus_tests
            unit
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_unit_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_unit_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_unit_contact.csv"
            aspect16
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_aspect16_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_aspect16_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_aspect16_contact.csv"
            area4
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_area4_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_area4_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_area4_contact.csv"
            area_quarter
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_area_quarter_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_area_quarter_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b47_hex8_scale_area_quarter_contact.csv")
    set_tests_properties(fuelsim_b47_hex8_scale_parameter_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;friction;small-sliding;scale;abaqus")

    add_executable(fuelsim_b48_multi_contact_contract tests/b48_multi_contact_contract.cpp)
    target_link_libraries(fuelsim_b48_multi_contact_contract PRIVATE fuelsim::io)
    target_include_directories(fuelsim_b48_multi_contact_contract PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    fuelsim_enable_warnings(fuelsim_b48_multi_contact_contract)
    fuelsim_stage_blackbox_case(b48_input b48 transient_hex8_multi_contact_path_abaqus.fsi)
    fuelsim_stage_blackbox_case(b48_split b48 transient_b48_checkpoint.fsi)
    fuelsim_stage_blackbox_case(b48_restart b48 transient_b48_restart.fsi)
    get_filename_component(b48_directory "${b48_input}" DIRECTORY)
    add_test(NAME fuelsim_b48_hex8_multi_contact_path_abaqus_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DINTERNAL_COMMAND=$<TARGET_FILE:fuelsim_b48_multi_contact_contract>"
            "-DCASE_DIRECTORY=${b48_directory}"
            "-DREFERENCE_DIRECTORY=${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_b48.cmake")
    set_tests_properties(fuelsim_b48_hex8_multi_contact_path_abaqus_tests
        PROPERTIES
            FIXTURES_SETUP b48_hex8_multi_contact_mpi_reference
            LABELS "acceptance;three-dimensional;hex8;contact;friction;finite-sliding;multi-contact;restart;abaqus")

    add_executable(fuelsim_hex20_solver_tests tests/hex20_solver_tests.cpp)
    target_link_libraries(fuelsim_hex20_solver_tests PRIVATE fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_hex20_solver_tests)
    add_test(NAME fuelsim_hex20_solver_tests
        COMMAND fuelsim_hex20_solver_tests
            "${CMAKE_CURRENT_BINARY_DIR}/hex20_mesh.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex20_steady_results.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex20_checkpoint.bin")
    set_tests_properties(fuelsim_hex20_solver_tests
        PROPERTIES LABELS "integration;three-dimensional;quadratic;checkpoint;contact;multiblock")

    fuelsim_stage_blackbox_case(
        FUELSIM_B6_BLACKBOX_INPUT b6 steady_hex20_u2_t1_moose.fsi)
    get_filename_component(FUELSIM_B6_BLACKBOX_DIRECTORY
        "${FUELSIM_B6_BLACKBOX_INPUT}" DIRECTORY)
    set(FUELSIM_B6_BLACKBOX_RESULTS
        "${FUELSIM_B6_BLACKBOX_DIRECTORY}/steady_hex20_u2_t1_moose_results.e")
    add_test(NAME fuelsim_b6_hex20_moose_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DINPUT_FILE=${FUELSIM_B6_BLACKBOX_INPUT}"
            "-DOUTPUT_FILE=${FUELSIM_B6_BLACKBOX_RESULTS}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=b6;${FUELSIM_B6_BLACKBOX_RESULTS};${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b6_hex20_u2_t1_temperature_nodes_0001.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b6_hex20_u2_t1_displacement_nodes_0001.csv"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
    set_tests_properties(fuelsim_b6_hex20_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic")

    fuelsim_add_hex20_blackbox_test(fuelsim_h20_02_convection_moose_tests h20_02
        steady_hex20_convection_moose.fsi
        h20_02_convection_temperature_nodes_0001.csv
        h20_02_convection_displacement_nodes_0001.csv
        1e-3 temperature-only "acceptance;three-dimensional;quadratic")
    fuelsim_add_hex20_blackbox_test(fuelsim_h20_03_pressure_moose_tests h20_03
        steady_hex20_pressure_moose.fsi
        h20_03_pressure_temperature_nodes_0001.csv
        h20_03_pressure_displacement_nodes_0001.csv
        1e-5 all "acceptance;three-dimensional;quadratic")
    fuelsim_add_hex20_blackbox_test(fuelsim_h20_04_traction_moose_tests h20_04
        steady_hex20_traction_moose.fsi
        h20_04_traction_temperature_nodes_0001.csv
        h20_04_traction_displacement_nodes_0001.csv
        1e-5 all "acceptance;three-dimensional;quadratic")
    fuelsim_add_hex20_blackbox_test(fuelsim_h20_05_thermo_mechanical_moose_tests h20_05
        steady_hex20_thermo_mechanical_moose.fsi
        h20_05_thermo_mechanical_temperature_nodes_0001.csv
        h20_05_thermo_mechanical_displacement_nodes_0001.csv
        1e-5 all "acceptance;three-dimensional;quadratic")

    fuelsim_stage_blackbox_case(
        FUELSIM_H20_06_BLACKBOX_INPUT h20_06 transient_hex20_moose.fsi)
    get_filename_component(FUELSIM_H20_06_BLACKBOX_DIRECTORY
        "${FUELSIM_H20_06_BLACKBOX_INPUT}" DIRECTORY)
    set(FUELSIM_H20_06_BLACKBOX_RESULTS
        "${FUELSIM_H20_06_BLACKBOX_DIRECTORY}/transient_hex20_moose_results.e")
    set(FUELSIM_H20_06_BLACKBOX_SUMMARY
        "${FUELSIM_H20_06_BLACKBOX_DIRECTORY}/transient_hex20_moose_summary.csv")
    add_test(NAME fuelsim_h20_06_transient_moose_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DINPUT_FILE=${FUELSIM_H20_06_BLACKBOX_INPUT}"
            "-DOUTPUT_FILES=${FUELSIM_H20_06_BLACKBOX_RESULTS};${FUELSIM_H20_06_BLACKBOX_SUMMARY}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=hex20-transient;${FUELSIM_H20_06_BLACKBOX_RESULTS};${FUELSIM_H20_06_BLACKBOX_SUMMARY};${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_06_transient_temperature_nodes_0001.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_06_transient_displacement_nodes_0001.csv;1e-3"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
    set_tests_properties(fuelsim_h20_06_transient_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic;transient")

    add_executable(fuelsim_hex20_history_contract_tests tests/hex20_history_contract_tests.cpp)
    target_link_libraries(fuelsim_hex20_history_contract_tests PRIVATE fuelsim::io)
    target_include_directories(fuelsim_hex20_history_contract_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_hex20_history_contract_tests)
    add_test(NAME fuelsim_hex20_history_contract_tests COMMAND fuelsim_hex20_history_contract_tests
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_inelastic_plastic_moose.fsi"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_inelastic_creep_moose.fsi"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_inelastic_coupled_moose.fsi")
    set_tests_properties(fuelsim_hex20_history_contract_tests PROPERTIES LABELS "unit;quadratic;transaction")
    foreach(branch IN ITEMS plastic creep coupled)
        fuelsim_add_production_test(fuelsim_h20_inelastic_${branch}_moose_tests h20_inelastic_${branch}
            transient_hex20_inelastic_${branch}_moose hex20-inelastic
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_${branch}_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_${branch}_displacement_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_${branch}_element_state_0001.csv")
    endforeach()
    fuelsim_add_hex20_blackbox_test(fuelsim_h20_10_finite_strain_moose_tests h20_10
        steady_hex20_finite_moose.fsi
        h20_10_finite_temperature_nodes_0001.csv
        h20_10_finite_displacement_nodes_0001.csv
        1e-3 all "acceptance;three-dimensional;quadratic;finite-strain")
    add_executable(fuelsim_h20_24_hex20_nonmatching_contact_tests tests/h20_nonmatching_contact_tests.cpp)
    target_link_libraries(fuelsim_h20_24_hex20_nonmatching_contact_tests
        PRIVATE fuelsim::io fuelsim_exodus_result_test_support)
    target_include_directories(fuelsim_h20_24_hex20_nonmatching_contact_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_h20_24_hex20_nonmatching_contact_tests)
    add_test(NAME fuelsim_h20_24_hex20_nonmatching_contact_tests
        COMMAND fuelsim_h20_24_hex20_nonmatching_contact_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_nonmatching_contact_abaqus.fsi"
            "${CMAKE_CURRENT_BINARY_DIR}/blackbox/h20_24_nonmatching/verification/fuelsim/steady_hex20_nonmatching_contact_abaqus_results.e")
    set_tests_properties(fuelsim_h20_24_hex20_nonmatching_contact_tests
        PROPERTIES FIXTURES_REQUIRED h20_24_production LABELS "unit;three-dimensional;quadratic;contact;nonmatching;jacobian;conservation")

    fuelsim_stage_blackbox_case(
        FUELSIM_B3_BLACKBOX_INPUT b3 steady_hex8_thermoelastic.fsi)
    get_filename_component(FUELSIM_B3_BLACKBOX_DIRECTORY
        "${FUELSIM_B3_BLACKBOX_INPUT}" DIRECTORY)
    set(FUELSIM_B3_BLACKBOX_RESULTS
        "${FUELSIM_B3_BLACKBOX_DIRECTORY}/steady_hex8_thermoelastic_results.e")
    add_test(NAME fuelsim_b3_hex8_moose_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DINPUT_FILE=${FUELSIM_B3_BLACKBOX_INPUT}"
            "-DOUTPUT_FILE=${FUELSIM_B3_BLACKBOX_RESULTS}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=b3;${FUELSIM_B3_BLACKBOX_RESULTS};${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_thermoelastic_all_nodes_final.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_thermoelastic_element_stress_final.csv"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
    set_tests_properties(fuelsim_b3_hex8_moose_tests PROPERTIES LABELS "acceptance;three-dimensional")

    fuelsim_stage_blackbox_case(FUELSIM_B35_BLACKBOX_INPUT b35
        steady_hex8_shared_meat_clad_moose.fsi)
    get_filename_component(FUELSIM_B35_BLACKBOX_DIRECTORY
        "${FUELSIM_B35_BLACKBOX_INPUT}" DIRECTORY)
    set(FUELSIM_B35_BLACKBOX_RESULTS
        "${FUELSIM_B35_BLACKBOX_DIRECTORY}/steady_hex8_shared_meat_clad_moose_results.e")
    add_test(NAME fuelsim_b35_hex8_shared_nodes_moose_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DINPUT_FILE=${FUELSIM_B35_BLACKBOX_INPUT}"
            "-DOUTPUT_FILE=${FUELSIM_B35_BLACKBOX_RESULTS}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=b35;${FUELSIM_B35_BLACKBOX_RESULTS};${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b35_hex8_shared_meat_clad_all_nodes_0001.csv"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
    set_tests_properties(fuelsim_b35_hex8_shared_nodes_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional")

    fuelsim_add_production_test(fuelsim_b36_plate_meat_clad_moose_tests b36
            transient_plate_meat_clad_moose b36
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b36_plate_meat_clad_all_nodes_0040.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b36_plate_meat_clad_element_state_0040.csv")
    set_tests_properties(fuelsim_b36_plate_meat_clad_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;inelastic")

    add_executable(fuelsim_hex8_shared_nodes_contract_tests tests/hex8_shared_nodes_contract_tests.cpp)
    target_link_libraries(fuelsim_hex8_shared_nodes_contract_tests
        PRIVATE fuelsim::io fuelsim_exodus_result_test_support)
    target_include_directories(fuelsim_hex8_shared_nodes_contract_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_hex8_shared_nodes_contract_tests)
    set_tests_properties(fuelsim_b35_hex8_shared_nodes_moose_tests PROPERTIES FIXTURES_SETUP b35_production)
    set_tests_properties(fuelsim_b36_plate_meat_clad_moose_tests PROPERTIES FIXTURES_SETUP b36_production)
    foreach(plate_case IN ITEMS b35 b36)
        if(plate_case STREQUAL "b35")
            set(plate_stem steady_hex8_shared_meat_clad_moose)
        else()
            set(plate_stem transient_plate_meat_clad_moose)
        endif()
        add_test(NAME "fuelsim_${plate_case}_shared_nodes_contract"
            COMMAND fuelsim_hex8_shared_nodes_contract_tests
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/${plate_stem}.fsi"
                "${CMAKE_CURRENT_BINARY_DIR}/blackbox/${plate_case}/verification/fuelsim/${plate_stem}_results.e"
                "${plate_case}")
        set_tests_properties("fuelsim_${plate_case}_shared_nodes_contract" PROPERTIES
            FIXTURES_REQUIRED "${plate_case}_production" LABELS "unit;shared-nodes;conservation")
    endforeach()

    fuelsim_add_production_test(fuelsim_b33_frictionless_baseline b33_frictionless
        steady_hex8_contact_frictionless_baseline frictionless interface)
    set_tests_properties(fuelsim_b33_frictionless_baseline PROPERTIES FIXTURES_SETUP b33_frictionless)
    fuelsim_add_production_test(fuelsim_b33_hex8_contact_moose_tests b33
            steady_hex8_contact_moose b33
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b33_hex8_thermal_contact_all_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b33_hex8_contact_all_nodes_0010.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b33_hex8_contact_secondary_contact_0010.csv"
            "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b33_frictionless/verification/fuelsim/steady_hex8_contact_frictionless_baseline_results.e")
    set_tests_properties(fuelsim_b33_hex8_contact_moose_tests PROPERTIES
        FIXTURES_REQUIRED b33_frictionless FIXTURES_SETUP b33_production
        LABELS "acceptance;production;three-dimensional;contact;friction")

    fuelsim_add_production_test(fuelsim_b34_hex8_sliding_contact_abaqus_tests b34
            steady_hex8_sliding_contact_abaqus b34
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b34_hex8_sliding_contact_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b34_hex8_sliding_contact_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b34_hex8_sliding_contact_reaction.csv")
    set_tests_properties(fuelsim_b34_hex8_sliding_contact_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;friction;restart;abaqus")
    fuelsim_add_production_test(fuelsim_b34_transient_sliding_checkpoint b34_transient
        transient_hex8_sliding_contact_checkpoint b34-transient
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b34_hex8_sliding_contact_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b34_hex8_sliding_contact_contact.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b34_hex8_sliding_contact_reaction.csv")
    set_tests_properties(fuelsim_b34_transient_sliding_checkpoint PROPERTIES FIXTURES_SETUP b34_transient_production)
    add_executable(fuelsim_hex8_friction_checkpoint_contract_tests tests/hex8_friction_checkpoint_contract_tests.cpp)
    target_link_libraries(fuelsim_hex8_friction_checkpoint_contract_tests
        PRIVATE fuelsim::io fuelsim_exodus_result_test_support)
    target_include_directories(fuelsim_hex8_friction_checkpoint_contract_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_hex8_friction_checkpoint_contract_tests)
    add_test(NAME fuelsim_b34_friction_checkpoint_contract
        COMMAND fuelsim_hex8_friction_checkpoint_contract_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_sliding_contact_checkpoint.fsi"
            "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b34_transient/verification/fuelsim/transient_hex8_sliding_contact_checkpoint.checkpoint"
            "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b34_transient/verification/fuelsim/transient_hex8_sliding_contact_checkpoint_results.e"
            "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b34_transient/verification/fuelsim/roundtrip.checkpoint")
    set_tests_properties(fuelsim_b34_friction_checkpoint_contract PROPERTIES
        FIXTURES_REQUIRED b34_transient_production LABELS "unit;hex8;friction;checkpoint;transaction")

    fuelsim_add_production_test(fuelsim_h20_hex20_thermal_contact_reference_tests h20_16
            steady_hex20_thermal_contact_moose hex20-thermal-contact
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_16_hex20_thermal_contact_temperature_nodes_0001.csv")
    set_tests_properties(fuelsim_h20_hex20_thermal_contact_reference_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;thermal")

    # Archived HEX20 node-to-surface diagnostic sources and references remain
    # available for study. Do not build obsolete example executables for a
    # mechanical-contact path rejected by the production program.

    # These tools analyze archived Abaqus probes. They do not exercise the
    # current Fuelsim implementation, so keep them available for focused
    # studies without registering or building them in the unified regression.
    add_executable(fuelsim_h20_26_abaqus_sts_identification_tests EXCLUDE_FROM_ALL
        tests/h20_abaqus_sts_identification_tests.cpp)
    fuelsim_enable_warnings(fuelsim_h20_26_abaqus_sts_identification_tests)

    add_executable(fuelsim_h20_27_abaqus_sts_nonmatching_identification_tests EXCLUDE_FROM_ALL
        tests/h20_abaqus_sts_nonmatching_identification_tests.cpp)
    fuelsim_enable_warnings(fuelsim_h20_27_abaqus_sts_nonmatching_identification_tests)

    add_executable(fuelsim_h20_28_abaqus_sts_smoothing_identification_tests EXCLUDE_FROM_ALL
        tests/h20_abaqus_sts_smoothing_identification_tests.cpp)
    fuelsim_enable_warnings(fuelsim_h20_28_abaqus_sts_smoothing_identification_tests)

    foreach(h20_29_case IN ITEMS biaxial tilted graded initial_gap)
        set(h20_29_test "fuelsim_h20_29_hex20_sts_${h20_29_case}_abaqus_tests")
        fuelsim_add_production_test(${h20_29_test} h20_29_${h20_29_case}
                steady_hex20_sts_${h20_29_case}_abaqus planar-sts
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_29_hex20_sts_${h20_29_case}_displacement.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_29_hex20_sts_${h20_29_case}_force.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_29_hex20_sts_${h20_29_case}_reaction.csv")
        set_tests_properties(${h20_29_test}
            PROPERTIES LABELS "acceptance;production;three-dimensional;hex20;contact;nonmatching;abaqus;split-case")
    endforeach()

    add_executable(fuelsim_hex20_curved_contact_contract_tests tests/hex20_curved_contact_contract_tests.cpp)
    target_link_libraries(fuelsim_hex20_curved_contact_contract_tests
        PRIVATE fuelsim::io fuelsim_exodus_result_test_support)
    target_include_directories(fuelsim_hex20_curved_contact_contract_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_hex20_curved_contact_contract_tests)
    foreach(curved_case IN ITEMS faceted_cylinder quadratic_cylinder)
        set(curved_stem "steady_hex20_sts_${curved_case}_abaqus")
        set(curved_test "fuelsim_h20_30_${curved_case}_abaqus_tests")
        fuelsim_add_production_test("${curved_test}" "h20_30_${curved_case}"
            "${curved_stem}" hex20-curved "${curved_case}"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_${curved_case}_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_${curved_case}_force.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_${curved_case}_resultant.csv")
        set_tests_properties("${curved_test}" PROPERTIES FIXTURES_SETUP "h20_30_${curved_case}_production"
            LABELS "acceptance;production;three-dimensional;hex20;contact;curved;abaqus")
        add_test(NAME "fuelsim_h20_30_${curved_case}_contract"
            COMMAND fuelsim_hex20_curved_contact_contract_tests
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/${curved_stem}.fsi"
                "${CMAKE_CURRENT_BINARY_DIR}/blackbox/h20_30_${curved_case}/verification/fuelsim/${curved_stem}_results.e"
                "${curved_case}")
        set_tests_properties("fuelsim_h20_30_${curved_case}_contract" PROPERTIES
            FIXTURES_REQUIRED "h20_30_${curved_case}_production" LABELS "unit;hex20;contact;jacobian;conservation")
    endforeach()
    fuelsim_add_production_test(fuelsim_h20_friction_path_abaqus_tests h20_33
            transient_hex20_nonmatching_friction_abaqus hex20-friction-path-33
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_33_hex20_nonmatching_friction_path_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_33_hex20_nonmatching_friction_path_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_33_hex20_nonmatching_friction_path_reaction.csv")
    set_tests_properties(fuelsim_h20_friction_path_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction;nonmatching;abaqus;restart")

    fuelsim_add_production_test(fuelsim_h20_35_curved_friction_abaqus_tests h20_35
            steady_hex20_sts_quadratic_cylinder_friction_abaqus hex20-curved-friction
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_35_hex20_sts_quadratic_cylinder_friction_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_35_hex20_sts_quadratic_cylinder_friction_contact.csv")
    set_tests_properties(fuelsim_h20_35_curved_friction_abaqus_tests
        PROPERTIES FIXTURES_SETUP h20_35_production LABELS "acceptance;production;three-dimensional;hex20;contact;friction;curved;abaqus")
    add_test(NAME fuelsim_h20_35_curved_friction_contract
        COMMAND fuelsim_hex20_curved_contact_contract_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_sts_quadratic_cylinder_friction_abaqus.fsi"
            "${CMAKE_CURRENT_BINARY_DIR}/blackbox/h20_35/verification/fuelsim/steady_hex20_sts_quadratic_cylinder_friction_abaqus_results.e"
            quadratic_friction)
    set_tests_properties(fuelsim_h20_35_curved_friction_contract PROPERTIES
        FIXTURES_REQUIRED h20_35_production LABELS "unit;hex20;contact;friction;jacobian;conservation")

    fuelsim_add_production_test(fuelsim_h20_36_curved_friction_path_abaqus_tests h20_36
            transient_hex20_sts_quadratic_cylinder_friction_path_abaqus hex20-friction-path-36
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_36_hex20_sts_quadratic_cylinder_friction_path_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_36_hex20_sts_quadratic_cylinder_friction_path_contact.csv")
    set_tests_properties(fuelsim_h20_36_curved_friction_path_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction;curved;abaqus;restart")

    foreach(friction_path IN ITEMS 33 36)
        add_executable(fuelsim_hex20_friction_transaction_${friction_path}_tests
            tests/hex20_friction_transaction_${friction_path}_tests.cpp)
        target_link_libraries(fuelsim_hex20_friction_transaction_${friction_path}_tests
            PRIVATE fuelsim::io fuelsim::solver)
        fuelsim_enable_warnings(fuelsim_hex20_friction_transaction_${friction_path}_tests)
    endforeach()
    add_test(NAME fuelsim_h20_33_friction_transaction_contract
        COMMAND fuelsim_hex20_friction_transaction_33_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_nonmatching_friction_abaqus.fsi"
            "${CMAKE_CURRENT_BINARY_DIR}/h20_33_transaction.checkpoint")
    add_test(NAME fuelsim_h20_36_friction_transaction_contract
        COMMAND fuelsim_hex20_friction_transaction_36_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_sts_quadratic_cylinder_friction_path_abaqus.fsi"
            "${CMAKE_CURRENT_BINARY_DIR}/h20_36_transaction.checkpoint")
    set_tests_properties(fuelsim_h20_33_friction_transaction_contract fuelsim_h20_36_friction_transaction_contract
        PROPERTIES LABELS "unit;hex20;friction;transaction;checkpoint" TIMEOUT 120)

    add_executable(fuelsim_h20_38_friction_objectivity_abaqus_tests
        tests/h20_friction_objectivity_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_38_friction_objectivity_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::io)
    fuelsim_enable_warnings(fuelsim_h20_38_friction_objectivity_abaqus_tests)
    add_test(NAME fuelsim_h20_38_friction_objectivity_abaqus_tests
        COMMAND fuelsim_h20_38_friction_objectivity_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_38_hex20_sts_friction_objectivity_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_38_hex20_sts_friction_objectivity.inp"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_38_hex20_sts_friction_objectivity.csv")
    set_tests_properties(fuelsim_h20_38_friction_objectivity_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction;sliding;objectivity;abaqus")

    add_executable(fuelsim_h20_39_finite_sliding_abaqus_tests
        tests/h20_finite_sliding_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_39_finite_sliding_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_39_finite_sliding_abaqus_tests)
    add_test(NAME fuelsim_h20_39_finite_sliding_abaqus_tests
        COMMAND fuelsim_h20_39_finite_sliding_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_39_hex20_finite_sliding_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_39_hex20_finite_sliding_small.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_39_hex20_finite_sliding_finite.csv")
    set_tests_properties(fuelsim_h20_39_finite_sliding_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction;finite-sliding;abaqus")

    add_executable(fuelsim_h20_40_curved_finite_sliding_abaqus_tests
        tests/h20_curved_finite_sliding_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_40_curved_finite_sliding_abaqus_tests
        PRIVATE fuelsim_field_error_metrics fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_40_curved_finite_sliding_abaqus_tests)
    add_test(NAME fuelsim_h20_40_curved_finite_sliding_finite_abaqus_tests
        COMMAND fuelsim_h20_40_curved_finite_sliding_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_40_hex20_curved_finite_sliding_mesh.e"
            finite
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_40_hex20_curved_finite_sliding_finite.csv")
    set_tests_properties(fuelsim_h20_40_curved_finite_sliding_finite_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction;curved;finite-sliding;abaqus")

    fuelsim_add_production_test(fuelsim_h20_41_hex20_partial_contact_abaqus_tests h20_41
        steady_hex20_finite_sliding_partial_contact_abaqus hex20-partial-contact
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_41_hex20_finite_sliding_partial_contact_contact.csv")
    set_tests_properties(fuelsim_h20_41_hex20_partial_contact_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;partial-contact;finite-sliding;abaqus")

    fuelsim_add_production_test(fuelsim_b549_c3d20t_abaqus_tests b549 transient_b549_small_c3d20t b549
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b549_small_c3d20t_temperature.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b549_small_c3d20t_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b549_small_c3d20t_material.csv")
    set_tests_properties(fuelsim_b549_c3d20t_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex20;c3d20t;transient;thermal;mechanical;finite-strain;plasticity;creep;abaqus")
    fuelsim_add_production_test(fuelsim_b550_c3d20t_contact_abaqus_tests b550 transient_b550_small_c3d20t_contact b550
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b550_small_c3d20t_contact_temperature.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b550_small_c3d20t_contact_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b550_small_c3d20t_contact_material.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b550_small_c3d20t_contact.csv")
    set_tests_properties(fuelsim_b550_c3d20t_contact_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex20;c3d20t;transient;thermal;mechanical;contact;finite-sliding;finite-strain;plasticity;creep;abaqus")
    fuelsim_add_production_test(fuelsim_b555_c3d20t_friction_abaqus_tests b555 transient_b555_small_c3d20t_friction b555
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b555_small_c3d20t_friction_temperature.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b555_small_c3d20t_friction_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b555_small_c3d20t_friction_material.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b555_small_c3d20t_friction.csv")
    set_tests_properties(fuelsim_b555_c3d20t_friction_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex20;c3d20t;transient;thermal;mechanical;contact;friction;finite-sliding;finite-strain;plasticity;creep;abaqus")
    add_executable(fuelsim_hex20_integrated_output_contract_tests tests/hex20_integrated_output_contract_tests.cpp)
    target_link_libraries(fuelsim_hex20_integrated_output_contract_tests
        PRIVATE fuelsim::io fuelsim_exodus_result_test_support)
    fuelsim_enable_warnings(fuelsim_hex20_integrated_output_contract_tests)
    set_tests_properties(fuelsim_b549_c3d20t_abaqus_tests PROPERTIES FIXTURES_SETUP b549_production)
    set_tests_properties(fuelsim_b550_c3d20t_contact_abaqus_tests PROPERTIES FIXTURES_SETUP b550_production)
    set_tests_properties(fuelsim_b555_c3d20t_friction_abaqus_tests PROPERTIES FIXTURES_SETUP b555_production)
    foreach(integrated_case IN ITEMS b549 b550 b555)
        if(integrated_case STREQUAL "b549")
            set(integrated_stem transient_b549_small_c3d20t)
        elseif(integrated_case STREQUAL "b550")
            set(integrated_stem transient_b550_small_c3d20t_contact)
        else()
            set(integrated_stem transient_b555_small_c3d20t_friction)
        endif()
        add_test(NAME fuelsim_${integrated_case}_output_contract
            COMMAND fuelsim_hex20_integrated_output_contract_tests
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/${integrated_stem}.fsi"
                "${CMAKE_CURRENT_BINARY_DIR}/blackbox/${integrated_case}/verification/fuelsim/${integrated_stem}.checkpoint"
                "${CMAKE_CURRENT_BINARY_DIR}/blackbox/${integrated_case}/verification/fuelsim/${integrated_stem}_results.e")
        set_tests_properties(fuelsim_${integrated_case}_output_contract PROPERTIES
            FIXTURES_REQUIRED ${integrated_case}_production LABELS "unit;hex20;io;jacobian;contact")
    endforeach()

    add_executable(fuelsim_h20_31_hex20_curved_sts_identification_tests
        tests/h20_curved_sts_identification_tests.cpp)
    target_link_libraries(fuelsim_h20_31_hex20_curved_sts_identification_tests
        PRIVATE fuelsim_io)
    fuelsim_enable_warnings(fuelsim_h20_31_hex20_curved_sts_identification_tests)
    add_test(NAME fuelsim_h20_31_hex20_curved_sts_identification_tests
        COMMAND fuelsim_h20_31_hex20_curved_sts_identification_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_31_hex20_curved_sts_operator.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_quadratic_cylinder_mesh.e")
    set_tests_properties(fuelsim_h20_31_hex20_curved_sts_identification_tests
        PROPERTIES LABELS "identification;three-dimensional;hex20;contact;curved;abaqus")

    fuelsim_add_production_test(fuelsim_h20_24_hex20_nonmatching_abaqus_tests h20_24_nonmatching
            steady_hex20_nonmatching_contact_abaqus hex20-nonmatching
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_24_hex20_nonmatching_surface_contact_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_24_hex20_nonmatching_surface_contact_reaction.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_24_hex20_nonmatching_surface_contact_pressure.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_24_hex20_nonmatching_mortar_all_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_24_hex20_nonmatching_mortar.csv")
    set_tests_properties(fuelsim_h20_24_hex20_nonmatching_abaqus_tests
        PROPERTIES FIXTURES_SETUP h20_24_production LABELS "acceptance;three-dimensional;hex20;contact;nonmatching;abaqus;moose")

    fuelsim_add_production_test(fuelsim_b37_hex8_multi_contact_moose_tests b37_multi_contact
            steady_hex8_multi_contact_friction_moose hex8-multi-contact
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b37_hex8_multi_contact_all_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b37_hex8_multi_contact_pair_a_contact_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b37_hex8_multi_contact_pair_b_contact_0001.csv")
    set_tests_properties(fuelsim_b37_hex8_multi_contact_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;contact;friction")

    add_executable(fuelsim_hex8_history_contract_tests tests/hex8_history_contract_tests.cpp)
    target_link_libraries(fuelsim_hex8_history_contract_tests PRIVATE fuelsim::io fuelsim_exodus_result_test_support)
    target_include_directories(fuelsim_hex8_history_contract_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_hex8_history_contract_tests)
    foreach(strain IN ITEMS small finite)
        if(strain STREQUAL "small")
            set(hex8_prefix b3)
            set(hex8_strain "")
        else()
            set(hex8_prefix b32)
            set(hex8_strain "finite_")
        endif()
        foreach(branch IN ITEMS creep coupled)
            set(hex8_stem "transient_hex8_${hex8_strain}${branch}_moose")
            set(hex8_case "${hex8_prefix}_${branch}")
            set(hex8_test "fuelsim_${hex8_prefix}_hex8_${branch}_moose_tests")
            fuelsim_add_production_test(${hex8_test} ${hex8_case} ${hex8_stem} hex8-inelastic
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/${hex8_prefix}_hex8_${hex8_strain}${branch}_all_nodes_final.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/${hex8_prefix}_hex8_${hex8_strain}${branch}_element_state_final.csv")
            set_tests_properties(${hex8_test} PROPERTIES FIXTURES_SETUP "${hex8_case}_production")
            set(hex8_directory "${CMAKE_CURRENT_BINARY_DIR}/blackbox/${hex8_case}/verification/fuelsim")
            add_test(NAME "fuelsim_${hex8_case}_history_contract"
                COMMAND fuelsim_hex8_history_contract_tests
                    "${hex8_directory}/${hex8_stem}.fsi"
                    "${hex8_directory}/${hex8_stem}.checkpoint"
                    "${hex8_directory}/history_roundtrip.checkpoint"
                    "${hex8_directory}/${hex8_stem}_results.e")
            set_tests_properties("fuelsim_${hex8_case}_history_contract" PROPERTIES
                FIXTURES_REQUIRED "${hex8_case}_production" LABELS "unit;transaction;checkpoint")
        endforeach()
    endforeach()
    foreach(plastic_case IN ITEMS b510 b515)
        if(plastic_case STREQUAL "b510")
            set(plastic_stem transient_b510_small_j2)
        else()
            set(plastic_stem transient_b515_finite_j2)
        endif()
        set(plastic_directory "${CMAKE_CURRENT_BINARY_DIR}/blackbox/${plastic_case}/verification/fuelsim")
        add_test(NAME "fuelsim_${plastic_case}_history_contract"
            COMMAND fuelsim_hex8_history_contract_tests
                "${plastic_directory}/${plastic_stem}.fsi"
                "${plastic_directory}/${plastic_stem}.checkpoint"
                "${plastic_directory}/history_roundtrip.checkpoint"
                "${plastic_directory}/${plastic_stem}_results.e")
        set_tests_properties("fuelsim_${plastic_case}_history_contract" PROPERTIES
            FIXTURES_REQUIRED "${plastic_case}_production" LABELS "unit;transaction;checkpoint")
    endforeach()

    add_executable(fuelsim_multi_region_tests tests/multi_region_tests.cpp)
    target_link_libraries(fuelsim_multi_region_tests
        PRIVATE fuelsim::core fuelsim::solver
    )

    add_executable(fuelsim_b528_hex8_c3d8rt_operator_abaqus_tests
        tests/b528_hex8_c3d8rt_operator_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b528_hex8_c3d8rt_operator_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b528_hex8_c3d8rt_operator_abaqus_tests)
    add_test(NAME fuelsim_b528_hex8_c3d8rt_operator_abaqus_tests
        COMMAND fuelsim_b528_hex8_c3d8rt_operator_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b528_hex8_c3d8rt_operator_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b528_hex8_c3d8rt_warped_operator_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b528_hex8_c3d8rt_holdout_operator_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b532_hex8_c3d8rt_finite_regular_operator_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b532_hex8_c3d8rt_finite_warped_operator_nodal.csv")
    set_tests_properties(fuelsim_b528_hex8_c3d8rt_operator_abaqus_tests PROPERTIES
        LABELS "acceptance;abaqus;three-dimensional;reduced-integration")

    add_executable(fuelsim_b529_hex8_c3d8rt_capacity_abaqus_tests
        tests/b529_hex8_c3d8rt_capacity_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b529_hex8_c3d8rt_capacity_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b529_hex8_c3d8rt_capacity_abaqus_tests)
    add_test(NAME fuelsim_b529_hex8_c3d8rt_capacity_abaqus_tests
        COMMAND fuelsim_b529_hex8_c3d8rt_capacity_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b529_hex8_c3d8rt_regular_capacity.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b529_hex8_c3d8rt_warped_capacity.csv")
    set_tests_properties(fuelsim_b529_hex8_c3d8rt_capacity_abaqus_tests PROPERTIES
        LABELS "acceptance;abaqus;three-dimensional;reduced-integration;transient;thermal")

    add_executable(fuelsim_b530_hex8_c3d8rt_thermal_load_abaqus_tests
        tests/b530_hex8_c3d8rt_thermal_load_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b530_hex8_c3d8rt_thermal_load_abaqus_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b530_hex8_c3d8rt_thermal_load_abaqus_tests)
    add_test(NAME fuelsim_b530_hex8_c3d8rt_thermal_load_abaqus_tests
        COMMAND fuelsim_b530_hex8_c3d8rt_thermal_load_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b530_hex8_c3d8rt_thermal_load.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b530_hex8_c3d8rt_thermal_load_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b530_hex8_c3d8rt_surface_geometry.csv")
    set_tests_properties(fuelsim_b530_hex8_c3d8rt_thermal_load_abaqus_tests PROPERTIES
        LABELS "acceptance;abaqus;three-dimensional;reduced-integration;thermal")

    fuelsim_add_production_test(fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests b531
            transient_b531_hex8_c3d8rt b531
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b531_hex8_c3d8rt_transient_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b531_hex8_c3d8rt_transient_integration.csv")
    set_tests_properties(fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests PROPERTIES
        LABELS "acceptance;abaqus;three-dimensional;reduced-integration;transient;thermo-mechanical")
    fuelsim_add_production_test(fuelsim_b534_hex8_c3d8rt_finite_transient_abaqus_tests b534
            transient_b534_hex8_c3d8rt_finite b534
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b534_hex8_c3d8rt_finite_transient_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b534_hex8_c3d8rt_finite_transient_integration.csv")
    set_tests_properties(fuelsim_b534_hex8_c3d8rt_finite_transient_abaqus_tests PROPERTIES
        LABELS "acceptance;abaqus;three-dimensional;reduced-integration;finite-strain;transient;thermo-mechanical")
    add_executable(fuelsim_hex8_reduced_history_contract_tests tests/hex8_reduced_history_contract_tests.cpp)
    target_link_libraries(fuelsim_hex8_reduced_history_contract_tests PRIVATE fuelsim::io fuelsim_exodus_result_test_support)
    fuelsim_enable_warnings(fuelsim_hex8_reduced_history_contract_tests)
    set_tests_properties(fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests PROPERTIES FIXTURES_SETUP b531_production)
    set_tests_properties(fuelsim_b534_hex8_c3d8rt_finite_transient_abaqus_tests PROPERTIES FIXTURES_SETUP b534_production)
    foreach(reduced_case IN ITEMS b531 b534)
        if(reduced_case STREQUAL "b531")
            set(reduced_stem transient_b531_hex8_c3d8rt)
        else()
            set(reduced_stem transient_b534_hex8_c3d8rt_finite)
        endif()
        set(reduced_directory "${CMAKE_CURRENT_BINARY_DIR}/blackbox/${reduced_case}/verification/fuelsim")
        add_test(NAME fuelsim_${reduced_case}_history_contract
            COMMAND fuelsim_hex8_reduced_history_contract_tests
                "${reduced_directory}/${reduced_stem}.fsi"
                "${reduced_directory}/${reduced_stem}.checkpoint"
                "${reduced_directory}/${reduced_stem}_results.e")
        set_tests_properties(fuelsim_${reduced_case}_history_contract PROPERTIES
            FIXTURES_REQUIRED ${reduced_case}_production LABELS "unit;transaction;reduced-integration")
    endforeach()
    fuelsim_enable_warnings(fuelsim_multi_region_tests)
    add_test(NAME fuelsim_multi_region_tests COMMAND fuelsim_multi_region_tests)
    set_tests_properties(fuelsim_multi_region_tests PROPERTIES LABELS "acceptance")

    add_executable(fuelsim_m2_core_tests tests/m2_core_tests.cpp)
    target_link_libraries(fuelsim_m2_core_tests
        PRIVATE fuelsim_test_mesh_fixture fuelsim::core
    )
    fuelsim_enable_warnings(fuelsim_m2_core_tests)
    add_test(NAME fuelsim_m2_core_tests COMMAND fuelsim_m2_core_tests)
    set_tests_properties(fuelsim_m2_core_tests PROPERTIES LABELS "unit")

    add_executable(fuelsim_solver_tests tests/solver_tests.cpp)
    target_link_libraries(fuelsim_solver_tests
        PRIVATE fuelsim_test_mesh_fixture fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_solver_tests)
    add_test(NAME fuelsim_solver_tests COMMAND fuelsim_solver_tests)
    set_tests_properties(fuelsim_solver_tests PROPERTIES LABELS "acceptance")

    find_program(FUELSIM_MPIEXEC_EXECUTABLE
        NAMES mpiexec
        HINTS "${PETSC_PREFIX}/bin"
        NO_DEFAULT_PATH
    )
    set(h20_production_reference "${CMAKE_CURRENT_BINARY_DIR}/blackbox/h20_inelastic_coupled/verification/fuelsim/transient_hex20_inelastic_coupled_moose_results.e")
    set_tests_properties(fuelsim_h20_inelastic_coupled_moose_tests PROPERTIES
        FIXTURES_SETUP h20_production_reference)
    fuelsim_stage_blackbox_case(h20_checkpoint_input h20_restart transient_hex20_coupled_checkpoint.fsi)
    fuelsim_stage_blackbox_case(h20_restart_input h20_restart transient_hex20_coupled_restart.fsi)
    get_filename_component(h20_restart_directory "${h20_restart_input}" DIRECTORY)
    set(h20_restart_result "${h20_restart_directory}/transient_hex20_coupled_restart_results.part1.e")
    set(h20_restart_summary "${h20_restart_directory}/transient_hex20_coupled_restart_summary.csv")
    add_test(NAME fuelsim_production_hex20_restart
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DCHECKPOINT_INPUT=${h20_checkpoint_input}"
            "-DCHECKPOINT_FILE=${h20_restart_directory}/transient_hex20_coupled.checkpoint"
            "-DINPUT_FILE=${h20_restart_input}"
            "-DOUTPUT_FILES=${h20_restart_result};${h20_restart_summary}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=restart;${h20_restart_result};${h20_restart_summary};${h20_production_reference}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_restart.cmake")
    set_tests_properties(fuelsim_production_hex20_restart PROPERTIES
        FIXTURES_REQUIRED h20_production_reference LABELS "acceptance;production;restart" TIMEOUT 120)
    if(FUELSIM_MPIEXEC_EXECUTABLE)
        fuelsim_stage_blackbox_case(gps_mpi_input gps_mpi transient_gps_chain_finite.fsi)
        get_filename_component(gps_mpi_directory "${gps_mpi_input}" DIRECTORY)
        set(gps_mpi_result "${gps_mpi_directory}/transient_gps_chain_finite_results.e")
        add_test(NAME fuelsim_gps_finite_mpi_equivalence COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DFUELSIM_LAUNCHER=${FUELSIM_MPIEXEC_EXECUTABLE};-n;2"
            "-DINPUT_FILE=${gps_mpi_input}"
            "-DOUTPUT_FILES=${gps_mpi_result};${gps_mpi_directory}/transient_gps_chain_finite.checkpoint"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_radial_equivalence>"
            "-DCHECK_ARGUMENTS=${CMAKE_CURRENT_BINARY_DIR}/blackbox/gps_chain_finite/verification/fuelsim/transient_gps_chain_finite_results.e;${gps_mpi_result}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        set_tests_properties(fuelsim_gps_finite_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED gps_finite_reference PROCESSORS 2 TIMEOUT 120)
        fuelsim_stage_blackbox_case(h20_mpi_input h20_mpi transient_hex20_inelastic_coupled_moose.fsi)
        get_filename_component(h20_mpi_directory "${h20_mpi_input}" DIRECTORY)
        set(h20_mpi_result "${h20_mpi_directory}/transient_hex20_inelastic_coupled_moose_results.e")
        set(h20_mpi_summary "${h20_mpi_directory}/transient_hex20_inelastic_coupled_moose_summary.csv")
        add_test(NAME fuelsim_production_hex20_mpi_equivalence
            COMMAND "${CMAKE_COMMAND}"
                "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
                "-DFUELSIM_LAUNCHER=${FUELSIM_MPIEXEC_EXECUTABLE};-n;2"
                "-DINPUT_FILE=${h20_mpi_input}"
                "-DOUTPUT_FILES=${h20_mpi_result};${h20_mpi_summary}"
                "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
                "-DCHECK_ARGUMENTS=mpi-equivalence;${h20_mpi_result};${h20_mpi_summary};${h20_production_reference}"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        set_tests_properties(fuelsim_production_hex20_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED h20_production_reference LABELS "acceptance;production;mpi"
            PROCESSORS 2 TIMEOUT 120)
        set(c3d8rt_mpi_reference "${CMAKE_CURRENT_BINARY_DIR}/c3d8rt_contact_one_rank.txt")
        add_test(NAME fuelsim_c3d8rt_contact_one_rank_reference
            COMMAND fuelsim_mpi_equivalence_benchmark
                write
                "${c3d8rt_mpi_reference}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_c3d8rt_contact_mpi.fsi")
        set_tests_properties(fuelsim_c3d8rt_contact_one_rank_reference PROPERTIES
            FIXTURES_SETUP c3d8rt_contact_reference
            LABELS "integration;mpi-reference;three-dimensional;hex8;c3d8rt;finite-strain;contact")
        add_test(NAME fuelsim_c3d8rt_contact_two_rank_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${c3d8rt_mpi_reference}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_c3d8rt_contact_mpi.fsi")
        set_tests_properties(fuelsim_c3d8rt_contact_two_rank_equivalence PROPERTIES
            FIXTURES_REQUIRED c3d8rt_contact_reference
            PROCESSORS 2
            LABELS "integration;mpi;three-dimensional;hex8;c3d8rt;finite-strain;contact;equivalence")

        add_test(NAME fuelsim_mpi_hex8_solver_tests
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                    "$<TARGET_FILE:fuelsim_hex8_solver_tests>"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex8_steady_results.e"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex8_transient_results.e"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex8_checkpoint.bin"
                    --mpi-only
                    -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps)
        set_tests_properties(fuelsim_mpi_hex8_solver_tests PROPERTIES
            PROCESSORS 2
            LABELS "integration;mpi;three-dimensional")

        add_test(NAME fuelsim_mpi_hex20_solver_tests
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                    "$<TARGET_FILE:fuelsim_hex20_solver_tests>"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex20_mesh.e"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex20_steady_results.e"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex20_checkpoint.bin"
                    --mpi-only
                    -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps)
        set_tests_properties(fuelsim_mpi_hex20_solver_tests PROPERTIES
            PROCESSORS 2
            LABELS "integration;mpi;three-dimensional;quadratic")

        add_test(NAME fuelsim_mpi_solver_tests
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_solver_tests>
                --case runtime-layout
                -ksp_type preonly
                -pc_type lu
                -pc_factor_mat_solver_type mumps
        )
        set_tests_properties(fuelsim_mpi_solver_tests PROPERTIES
            LABELS "acceptance;parallel"
            PROCESSORS 2
        )

        add_test(NAME fuelsim_mpi_runtime_layout_field_split
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_solver_tests>
                --case runtime-layout
        )
        set_tests_properties(fuelsim_mpi_runtime_layout_field_split PROPERTIES
            LABELS "acceptance;parallel;architecture"
            PROCESSORS 2
            TIMEOUT 60
        )

        add_executable(fuelsim_mpi_equivalence_benchmark
            benchmarks/mpi_equivalence.cpp
        )
        target_link_libraries(fuelsim_mpi_equivalence_benchmark
            PRIVATE fuelsim::io fuelsim::solver fuelsim_exodus_result_test_support
        )
        target_include_directories(fuelsim_mpi_equivalence_benchmark PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
        fuelsim_enable_warnings(fuelsim_mpi_equivalence_benchmark)
        add_test(NAME fuelsim_rz_contact_mpi_reference
            COMMAND fuelsim_mpi_equivalence_benchmark write
                "${CMAKE_CURRENT_BINARY_DIR}/rz_contact_one_rank.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi")
        set_tests_properties(fuelsim_rz_contact_mpi_reference PROPERTIES
            FIXTURES_SETUP rz_contact_one_rank LABELS "unit;mpi-reference")
        add_test(NAME fuelsim_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${CMAKE_CURRENT_BINARY_DIR}/rz_contact_one_rank.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
        )
        set_tests_properties(fuelsim_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED rz_contact_one_rank
            LABELS "acceptance;parallel"
            PROCESSORS 2
        )

        add_test(NAME fuelsim_b33_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b33/verification/fuelsim/steady_hex8_contact_moose_results.e"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_b33_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED b33_production
            LABELS "acceptance;parallel;three-dimensional;contact"
            PROCESSORS 2
        )
        add_test(NAME fuelsim_b521_thermal_contact_mpi_reference
            COMMAND fuelsim_mpi_equivalence_benchmark
                write
                "${CMAKE_CURRENT_BINARY_DIR}/b521_thermal_contact_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_pressure_gap_conductance_abaqus.fsi"
        )
        set_tests_properties(fuelsim_b521_thermal_contact_mpi_reference PROPERTIES
            FIXTURES_SETUP b521_thermal_contact_mpi_reference
            LABELS "integration;parallel;three-dimensional;hex8;thermo-mechanical;contact;abaqus"
        )
        add_test(NAME fuelsim_b521_thermal_contact_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${CMAKE_CURRENT_BINARY_DIR}/b521_thermal_contact_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_pressure_gap_conductance_abaqus.fsi"
        )
        set_tests_properties(fuelsim_b521_thermal_contact_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED b521_thermal_contact_mpi_reference
            LABELS "acceptance;parallel;three-dimensional;hex8;thermo-mechanical;contact;abaqus"
            PROCESSORS 2
        )
        add_test(NAME fuelsim_b39_hex8_sts_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b39_nonmatching/verification/fuelsim/steady_hex8_sts_nonmatching_abaqus_results.e"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_sts_nonmatching_abaqus.fsi"
        )
        set_tests_properties(fuelsim_b39_hex8_sts_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED b39_nonmatching_production
            LABELS "acceptance;parallel;three-dimensional;hex8;contact;surface-to-surface"
            PROCESSORS 2
        )
        add_test(NAME fuelsim_b41_hex8_finite_sliding_mpi_reference
            COMMAND fuelsim_mpi_equivalence_benchmark
                write
                "${CMAKE_CURRENT_BINARY_DIR}/b41_hex8_finite_sliding_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_finite_sliding_abaqus.fsi"
        )
        set_tests_properties(fuelsim_b41_hex8_finite_sliding_mpi_reference PROPERTIES
            FIXTURES_SETUP b41_hex8_finite_sliding_mpi_reference
            LABELS "integration;parallel;three-dimensional;hex8;contact;surface-to-surface;finite-sliding"
        )
        add_test(NAME fuelsim_b41_hex8_finite_sliding_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${CMAKE_CURRENT_BINARY_DIR}/b41_hex8_finite_sliding_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_finite_sliding_abaqus.fsi"
        )
        set_tests_properties(fuelsim_b41_hex8_finite_sliding_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED b41_hex8_finite_sliding_mpi_reference
            LABELS "acceptance;parallel;three-dimensional;hex8;contact;surface-to-surface;finite-sliding"
            PROCESSORS 2
        )
        fuelsim_stage_blackbox_case(b48_mpi_input b48_mpi transient_hex8_multi_contact_path_abaqus.fsi)
        get_filename_component(b48_mpi_directory "${b48_mpi_input}" DIRECTORY)
        add_test(NAME fuelsim_b48_hex8_multi_contact_path_mpi_equivalence
            COMMAND "${CMAKE_COMMAND}"
                "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
                "-DFUELSIM_LAUNCHER=${FUELSIM_MPIEXEC_EXECUTABLE};-n;2"
                "-DINPUT_FILE=${b48_mpi_input}"
                "-DOUTPUT_FILES=${b48_mpi_directory}/transient_hex8_multi_contact_path_abaqus_results.e;${b48_mpi_directory}/transient_hex8_multi_contact_path_abaqus_summary.csv"
                "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
                "-DCHECK_ARGUMENTS=b48-mpi;${b48_mpi_directory}/transient_hex8_multi_contact_path_abaqus_results.e;${b48_mpi_directory}/transient_hex8_multi_contact_path_abaqus_summary.csv;${b48_directory}/transient_hex8_multi_contact_path_abaqus_results.e"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        set_tests_properties(fuelsim_b48_hex8_multi_contact_path_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED b48_hex8_multi_contact_mpi_reference
            LABELS "acceptance;parallel;three-dimensional;hex8;contact;friction;surface-to-surface;finite-sliding;multi-contact"
            PROCESSORS 2
        )
        add_test(NAME fuelsim_mpi_field_split
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare_field_split
                "${CMAKE_CURRENT_BINARY_DIR}/rz_contact_one_rank.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
        )
        set_tests_properties(fuelsim_mpi_field_split PROPERTIES
            FIXTURES_REQUIRED rz_contact_one_rank
            LABELS "acceptance;parallel"
            PROCESSORS 2
            TIMEOUT 60
        )
        add_test(NAME fuelsim_mpi_block_jacobi
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare_block_jacobi
                "${CMAKE_CURRENT_BINARY_DIR}/rz_contact_one_rank.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
        )
        set_tests_properties(fuelsim_mpi_block_jacobi PROPERTIES
            FIXTURES_REQUIRED rz_contact_one_rank
            LABELS "integration;parallel"
            PROCESSORS 2
            TIMEOUT 60
        )
        add_test(NAME fuelsim_mpi_hypre
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare_hypre
                "${CMAKE_CURRENT_BINARY_DIR}/rz_contact_one_rank.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
        )
        set_tests_properties(fuelsim_mpi_hypre PROPERTIES
            FIXTURES_REQUIRED rz_contact_one_rank
            LABELS "integration;parallel"
            PROCESSORS 2
            TIMEOUT 120
        )
        add_test(NAME fuelsim_transient_mpi_reference
            COMMAND fuelsim_mpi_equivalence_benchmark
                write_transient
                "${CMAKE_CURRENT_BINARY_DIR}/m35_transient_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
        )
        set_tests_properties(fuelsim_transient_mpi_reference PROPERTIES
            FIXTURES_SETUP m35_transient_mpi_reference
            LABELS "integration;parallel"
        )
        add_test(NAME fuelsim_transient_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare_transient
                "${CMAKE_CURRENT_BINARY_DIR}/m35_transient_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
        )
        set_tests_properties(fuelsim_transient_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED m35_transient_mpi_reference
            LABELS "acceptance;parallel"
            PROCESSORS 2
            TIMEOUT 60
        )
        add_test(NAME fuelsim_mpi_root_io_failure
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                test_io_failure
                unused-reference
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
        )
        set_tests_properties(fuelsim_mpi_root_io_failure PROPERTIES
            LABELS "integration;parallel"
            PROCESSORS 2
            TIMEOUT 30
        )
    endif()

    add_executable(fuelsim_m21_manufactured_heat_tests
        tests/m21_manufactured_heat_tests.cpp
    )
    target_link_libraries(fuelsim_m21_manufactured_heat_tests
        PRIVATE fuelsim_test_mesh_fixture fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m21_manufactured_heat_tests)
    add_test(NAME fuelsim_m21_manufactured_heat_tests
        COMMAND fuelsim_m21_manufactured_heat_tests
    )
    set_tests_properties(fuelsim_m21_manufactured_heat_tests
        PROPERTIES LABELS "acceptance"
    )

    # Primary axisymmetric Abaqus qualification. Only scheduling metadata is generated;
    # every physical input is a complete, tracked card copied without substitutions.
    foreach(rz_case IN ITEMS b70_rz_small_coupled b71_rz_finite_coupled
        b72_rz_small_creep b73_rz_finite_creep b74_rz_small_plastic b75_rz_finite_plastic
        b76_rz_small_thermal b77_rz_finite_thermal b78_rz_small_contact b79_rz_finite_contact)
        string(REGEX REPLACE ".*_" "" rz_mechanism "${rz_case}")
        if(rz_case STREQUAL "b72_rz_small_creep")
            set(rz_mechanism "creep-integration")
        endif()
        if(rz_mechanism STREQUAL "contact")
            set(rz_input "steady_${rz_case}")
        else()
            set(rz_input "transient_${rz_case}")
        endif()
        fuelsim_add_production_test("fuelsim_${rz_case}_abaqus_tests" "${rz_case}" "${rz_input}" rz-abaqus
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${rz_case}_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${rz_case}_points.csv" "${rz_mechanism}")
    endforeach()

    fuelsim_add_production_test(fuelsim_b80_rz_small_friction_abaqus_tests b80_rz_small_friction
        transient_b80_rz_small_friction rz-friction-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b80_rz_small_friction_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b80_rz_small_friction_points.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b80_rz_small_friction_contact.csv")
    fuelsim_add_production_test(fuelsim_b81_rz_finite_sliding_abaqus_tests b81_rz_finite_sliding
        transient_b81_rz_finite_sliding rz-sliding-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b81_rz_finite_sliding_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b81_rz_finite_sliding_points.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b81_rz_finite_sliding_contact.csv")
    set_tests_properties(fuelsim_b81_rz_finite_sliding_abaqus_tests PROPERTIES FIXTURES_SETUP b81_rz_production)
    add_executable(fuelsim_rz_friction_history_contract_tests tests/rz_friction_history_contract_tests.cpp)
    target_link_libraries(fuelsim_rz_friction_history_contract_tests PRIVATE fuelsim::io fuelsim_exodus_result_test_support)
    target_include_directories(fuelsim_rz_friction_history_contract_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_enable_warnings(fuelsim_rz_friction_history_contract_tests)
    add_test(NAME fuelsim_b81_rz_friction_history_contract COMMAND fuelsim_rz_friction_history_contract_tests
        "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b81_rz_finite_sliding/verification/fuelsim/transient_b81_rz_finite_sliding.fsi"
        "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b81_rz_finite_sliding/verification/fuelsim/transient_b81_rz_finite_sliding.checkpoint"
        "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b81_rz_finite_sliding/verification/fuelsim/transient_b81_rz_finite_sliding_results.e")
    set_tests_properties(fuelsim_b81_rz_friction_history_contract PROPERTIES
        FIXTURES_REQUIRED b81_rz_production LABELS "unit;contact;checkpoint")

    # B15.1 retains the small-strain capacity/source stages and adds variable conductivity.
    # The constant-conductivity B15.0 small-strain input remains a manual reference.
    foreach(cax4t_thermal_case IN ITEMS b150_cax4t_finite_thermal_operator b151_cax4t_small_conductivity)
        fuelsim_add_production_test("fuelsim_${cax4t_thermal_case}_abaqus_tests" "${cax4t_thermal_case}"
            "transient_${cax4t_thermal_case}" rz-abaqus
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax4t_thermal_case}_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax4t_thermal_case}_points.csv" "probe")
    endforeach()

    # B9.16 retains the B9.0 mechanical state and extends its thermal operator checks;
    # B9.15 separately identifies constant-conductivity modes. B9.0 remains manual.
    foreach(cax4rt_case IN ITEMS b91_cax4rt_finite_probe b92_cax4rt_rectangle_probe
        b93_cax4rt_small_coupled b94_cax4rt_finite_coupled b95_cax4rt_small_creep
        b96_cax4rt_finite_creep b97_cax4rt_small_plastic b98_cax4rt_finite_plastic
        b99_cax4rt_small_thermal b910_cax4rt_finite_thermal b911_cax4rt_small_contact b912_cax4rt_finite_contact)
        string(REGEX REPLACE ".*_" "" cax4rt_mechanism "${cax4rt_case}")
        if(cax4rt_case STREQUAL "b95_cax4rt_small_creep")
            set(cax4rt_mechanism "creep-integration")
        endif()
        if(cax4rt_mechanism STREQUAL "contact")
            set(cax4rt_input "steady_${cax4rt_case}")
        else()
            set(cax4rt_input "transient_${cax4rt_case}")
        endif()
        fuelsim_add_production_test("fuelsim_${cax4rt_case}_abaqus_tests" "${cax4rt_case}"
            "${cax4rt_input}" rz-abaqus
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax4rt_case}_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax4rt_case}_points.csv" "${cax4rt_mechanism}")
    endforeach()

    foreach(cax8t_case IN ITEMS b100_cax8t_small_probe b101_cax8t_finite_probe
        b104_cax8t_small_coupled b105_cax8t_finite_coupled
        b106_cax8t_small_creep b107_cax8t_finite_creep
        b108_cax8t_small_plastic b109_cax8t_finite_plastic
        b110_cax8t_small_thermal b111_cax8t_finite_thermal)
        string(REGEX REPLACE ".*_" "" cax8t_mechanism "${cax8t_case}")
        if(cax8t_case STREQUAL "b106_cax8t_small_creep")
            set(cax8t_mechanism "creep-integration")
        endif()
        fuelsim_add_production_test("fuelsim_${cax8t_case}_abaqus_tests" "${cax8t_case}"
            "transient_${cax8t_case}" rz-abaqus
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax8t_case}_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax8t_case}_points.csv" "${cax8t_mechanism}")
    endforeach()

    foreach(cax8t_case IN ITEMS b102_cax8t_small_contact b103_cax8t_finite_contact)
        fuelsim_add_production_test("fuelsim_${cax8t_case}_abaqus_tests" "${cax8t_case}"
            "steady_${cax8t_case}" rz-abaqus
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax8t_case}_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax8t_case}_points.csv" contact)
    endforeach()

    fuelsim_add_production_test(fuelsim_b112_cax8t_small_friction_abaqus_tests b112_cax8t_small_friction
        transient_b112_cax8t_small_friction rz8-friction-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b112_cax8t_small_friction_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b112_cax8t_small_friction_points.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b112_cax8t_small_friction_contact.csv")
    fuelsim_add_production_test(fuelsim_b113_cax8t_finite_sliding_abaqus_tests b113_cax8t_finite_sliding
        transient_b113_cax8t_finite_sliding rz8-sliding-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b113_cax8t_finite_sliding_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b113_cax8t_finite_sliding_points.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b113_cax8t_finite_sliding_contact.csv")

    fuelsim_add_production_test(fuelsim_b114_cax8t_recovery_abaqus_tests b114_cax8t_recovery
        transient_b114_cax8t_recovery rz8-recovery-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b114_cax8t_recovery_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b114_cax8t_recovery_points.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b114_cax8t_recovery_contact.csv")

    foreach(cax8rt_case IN ITEMS b120_cax8rt_small_probe b121_cax8rt_finite_probe
        b124_cax8rt_small_coupled b125_cax8rt_finite_coupled
        b126_cax8rt_small_creep b127_cax8rt_finite_creep
        b128_cax8rt_small_plastic b129_cax8rt_finite_plastic
        b1210_cax8rt_small_thermal b1211_cax8rt_finite_thermal)
        string(REGEX REPLACE ".*_" "" cax8rt_mechanism "${cax8rt_case}")
        if(cax8rt_case STREQUAL "b126_cax8rt_small_creep")
            set(cax8rt_mechanism "creep-integration")
        endif()
        fuelsim_add_production_test("fuelsim_${cax8rt_case}_abaqus_tests" "${cax8rt_case}"
            "transient_${cax8rt_case}" rz-abaqus
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax8rt_case}_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax8rt_case}_points.csv" "${cax8rt_mechanism}")
    endforeach()

    foreach(cax8rt_case IN ITEMS b122_cax8rt_small_contact b123_cax8rt_finite_contact)
        fuelsim_add_production_test("fuelsim_${cax8rt_case}_abaqus_tests" "${cax8rt_case}"
            "steady_${cax8rt_case}" rz-abaqus
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax8rt_case}_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax8rt_case}_points.csv" contact)
    endforeach()

    fuelsim_add_production_test(fuelsim_b1212_cax8rt_small_friction_abaqus_tests b1212_cax8rt_small_friction
        transient_b1212_cax8rt_small_friction rz8-friction-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b1212_cax8rt_small_friction_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b1212_cax8rt_small_friction_points.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b1212_cax8rt_small_friction_contact.csv")
    fuelsim_add_production_test(fuelsim_b1213_cax8rt_finite_sliding_abaqus_tests b1213_cax8rt_finite_sliding
        transient_b1213_cax8rt_finite_sliding rz8-sliding-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b1213_cax8rt_finite_sliding_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b1213_cax8rt_finite_sliding_points.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b1213_cax8rt_finite_sliding_contact.csv")

    # B10.14 identifies the shared quadratic-edge recovery on prescribed states.
    # B12.0/B12.12/B12.13 retain CAX8RT bulk/contact coverage; B12.14 remains manual.

    set_tests_properties(fuelsim_b105_cax8t_finite_coupled_abaqus_tests PROPERTIES FIXTURES_SETUP b105_cax8t_finite_coupled)
    fuelsim_stage_blackbox_case(rz8_split b105_cax8t_finite_coupled_restart transient_b105_cax8t_finite_coupled_split.fsi)
    fuelsim_stage_blackbox_case(rz8_restart b105_cax8t_finite_coupled_restart transient_b105_cax8t_finite_coupled_restart.fsi)
    get_filename_component(rz8_directory "${rz8_restart}" DIRECTORY)
    set(rz8_result "${rz8_directory}/transient_b105_cax8t_finite_coupled_restart_results.part1.e")
    set(rz8_summary "${rz8_directory}/transient_b105_cax8t_finite_coupled_restart_summary.csv")
    add_test(NAME fuelsim_b105_cax8t_finite_coupled_restart_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DCHECKPOINT_INPUT=${rz8_split}"
            "-DCHECKPOINT_FILE=${rz8_directory}/transient_b105_cax8t_finite_coupled_split.checkpoint"
            "-DINPUT_FILE=${rz8_restart}"
            "-DOUTPUT_FILES=${rz8_result};${rz8_summary}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=rz8-restart;${rz8_result};${rz8_summary};${CMAKE_CURRENT_BINARY_DIR}/blackbox/b105_cax8t_finite_coupled/verification/fuelsim/transient_b105_cax8t_finite_coupled_results.e;3"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_restart.cmake")
    set_tests_properties(fuelsim_b105_cax8t_finite_coupled_restart_tests PROPERTIES
        FIXTURES_REQUIRED b105_cax8t_finite_coupled LABELS "acceptance;production;restart;cax8t" TIMEOUT 120)

    set_tests_properties(fuelsim_b125_cax8rt_finite_coupled_abaqus_tests PROPERTIES FIXTURES_SETUP b125_cax8rt_finite_coupled)
    fuelsim_stage_blackbox_case(rz8_split b125_cax8rt_finite_coupled_restart transient_b125_cax8rt_finite_coupled_split.fsi)
    fuelsim_stage_blackbox_case(rz8_restart b125_cax8rt_finite_coupled_restart transient_b125_cax8rt_finite_coupled_restart.fsi)
    get_filename_component(rz8_directory "${rz8_restart}" DIRECTORY)
    set(rz8_result "${rz8_directory}/transient_b125_cax8rt_finite_coupled_restart_results.part1.e")
    set(rz8_summary "${rz8_directory}/transient_b125_cax8rt_finite_coupled_restart_summary.csv")
    add_test(NAME fuelsim_b125_cax8rt_finite_coupled_restart_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DCHECKPOINT_INPUT=${rz8_split}"
            "-DCHECKPOINT_FILE=${rz8_directory}/transient_b125_cax8rt_finite_coupled_split.checkpoint"
            "-DINPUT_FILE=${rz8_restart}"
            "-DOUTPUT_FILES=${rz8_result};${rz8_summary}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=rz8-restart;${rz8_result};${rz8_summary};${CMAKE_CURRENT_BINARY_DIR}/blackbox/b125_cax8rt_finite_coupled/verification/fuelsim/transient_b125_cax8rt_finite_coupled_results.e;3"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_restart.cmake")
    set_tests_properties(fuelsim_b125_cax8rt_finite_coupled_restart_tests PROPERTIES
        FIXTURES_REQUIRED b125_cax8rt_finite_coupled LABELS "acceptance;production;restart;cax8rt" TIMEOUT 120)

    set_tests_properties(fuelsim_b1213_cax8rt_finite_sliding_abaqus_tests PROPERTIES FIXTURES_SETUP b1213_cax8rt_finite_sliding)
    fuelsim_stage_blackbox_case(rz8_split b1213_cax8rt_finite_sliding_restart transient_b1213_cax8rt_finite_sliding_split.fsi)
    fuelsim_stage_blackbox_case(rz8_restart b1213_cax8rt_finite_sliding_restart transient_b1213_cax8rt_finite_sliding_restart.fsi)
    get_filename_component(rz8_directory "${rz8_restart}" DIRECTORY)
    set(rz8_result "${rz8_directory}/transient_b1213_cax8rt_finite_sliding_restart_results.part1.e")
    set(rz8_summary "${rz8_directory}/transient_b1213_cax8rt_finite_sliding_restart_summary.csv")
    add_test(NAME fuelsim_b1213_cax8rt_finite_sliding_restart_tests
        COMMAND "${CMAKE_COMMAND}"
            "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
            "-DCHECKPOINT_INPUT=${rz8_split}"
            "-DCHECKPOINT_FILE=${rz8_directory}/transient_b1213_cax8rt_finite_sliding_split.checkpoint"
            "-DINPUT_FILE=${rz8_restart}"
            "-DOUTPUT_FILES=${rz8_result};${rz8_summary}"
            "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_result_tests>"
            "-DCHECK_ARGUMENTS=rz8-restart;${rz8_result};${rz8_summary};${CMAKE_CURRENT_BINARY_DIR}/blackbox/b1213_cax8rt_finite_sliding/verification/fuelsim/transient_b1213_cax8rt_finite_sliding_results.e;5"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_restart.cmake")
    set_tests_properties(fuelsim_b1213_cax8rt_finite_sliding_restart_tests PROPERTIES
        FIXTURES_REQUIRED b1213_cax8rt_finite_sliding LABELS "acceptance;production;restart;cax8rt" TIMEOUT 120)

    fuelsim_add_production_test(fuelsim_b915_cax4rt_thermal_identification_abaqus_tests b915_cax4rt_thermal_identification
        transient_b915_cax4rt_thermal_identification rz-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b915_cax4rt_thermal_identification_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b915_cax4rt_thermal_identification_points.csv" probe)

    foreach(cax4rt_thermal_case IN ITEMS b916_cax4rt_variable_conductivity b917_cax4rt_finite_thermal_operators)
        fuelsim_add_production_test(fuelsim_${cax4rt_thermal_case}_abaqus_tests ${cax4rt_thermal_case}
            transient_${cax4rt_thermal_case} rz-abaqus
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax4rt_thermal_case}_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${cax4rt_thermal_case}_points.csv" probe)
    endforeach()

    fuelsim_add_production_test(fuelsim_b913_cax4rt_small_friction_abaqus_tests b913_cax4rt_small_friction
        transient_b913_cax4rt_small_friction rz-friction-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b913_cax4rt_small_friction_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b913_cax4rt_small_friction_points.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b913_cax4rt_small_friction_contact.csv")
    fuelsim_add_production_test(fuelsim_b914_cax4rt_finite_sliding_abaqus_tests b914_cax4rt_finite_sliding
        transient_b914_cax4rt_finite_sliding rz-sliding-abaqus
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b914_cax4rt_finite_sliding_nodes.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b914_cax4rt_finite_sliding_points.csv"
        "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b914_cax4rt_finite_sliding_contact.csv")
    set_tests_properties(fuelsim_b914_cax4rt_finite_sliding_abaqus_tests PROPERTIES FIXTURES_SETUP b914_cax4rt_production)
    add_test(NAME fuelsim_b914_cax4rt_friction_history_contract COMMAND fuelsim_rz_friction_history_contract_tests
        "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b914_cax4rt_finite_sliding/verification/fuelsim/transient_b914_cax4rt_finite_sliding.fsi"
        "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b914_cax4rt_finite_sliding/verification/fuelsim/transient_b914_cax4rt_finite_sliding.checkpoint"
        "${CMAKE_CURRENT_BINARY_DIR}/blackbox/b914_cax4rt_finite_sliding/verification/fuelsim/transient_b914_cax4rt_finite_sliding_results.e")
    set_tests_properties(fuelsim_b914_cax4rt_friction_history_contract PROPERTIES
        FIXTURES_REQUIRED b914_cax4rt_production LABELS "unit;contact;checkpoint")

    foreach(pcmi_group IN ITEMS small finite integrated)
        foreach(pcmi_element IN ITEMS cax4t cax4rt cax8t cax8rt)
            set(pcmi_case "b13_${pcmi_group}_${pcmi_element}")
            fuelsim_stage_blackbox_case(pcmi_input "${pcmi_case}" "transient_${pcmi_case}.fsi")
            get_filename_component(pcmi_directory "${pcmi_input}" DIRECTORY)
            set(pcmi_output "${pcmi_directory}/transient_${pcmi_case}_results.e")
            add_test(NAME "fuelsim_${pcmi_case}_abaqus_tests"
                COMMAND "${CMAKE_COMMAND}"
                    "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
                    "-DINPUT_FILE=${pcmi_input}"
                    "-DOUTPUT_FILE=${pcmi_output}"
                    "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_pcmi_abaqus>"
                    "-DCHECK_ARGUMENTS=${pcmi_output};${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${pcmi_case}_nodes.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${pcmi_case}_points.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${pcmi_case}_contact.csv;${pcmi_group}"
                    -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
            set_tests_properties("fuelsim_${pcmi_case}_abaqus_tests" PROPERTIES
                LABELS "acceptance;production;abaqus;pcmi" TIMEOUT 240)
        endforeach()
    endforeach()
    # Prescribed-state NTS probes need one representative per edge topology.
    # Reduced-integration production contact remains covered by B13.
    foreach(nts_element IN ITEMS cax4t cax8t)
        set(nts_case "b14_nts_${nts_element}")
        fuelsim_stage_blackbox_case(nts_input "${nts_case}" "transient_${nts_case}.fsi")
        get_filename_component(nts_directory "${nts_input}" DIRECTORY)
        set(nts_output "${nts_directory}/transient_${nts_case}_results.e")
        add_test(NAME "fuelsim_${nts_case}_abaqus_tests"
            COMMAND "${CMAKE_COMMAND}"
                "-DFUELSIM_EXECUTABLE=$<TARGET_FILE:fuelsim>"
                "-DINPUT_FILE=${nts_input}"
                "-DOUTPUT_FILE=${nts_output}"
                "-DCHECK_COMMAND=$<TARGET_FILE:fuelsim_production_nts_abaqus>"
                "-DCHECK_ARGUMENTS=${nts_output};${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${nts_case}_nodes.csv;${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/${nts_case}_contact.csv"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fuelsim_case.cmake")
        set_tests_properties("fuelsim_${nts_case}_abaqus_tests" PROPERTIES LABELS "acceptance;production;abaqus;contact" TIMEOUT 120)
    endforeach()

    add_test(NAME fuelsim_transient_jacobian_check
        COMMAND fuelsim -i
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/input_cards/rz_internal.fsi"
            --check-jacobian
    )
    set_tests_properties(fuelsim_transient_jacobian_check
        PROPERTIES LABELS "integration"
    )

    find_program(FUELSIM_BASH_EXECUTABLE bash REQUIRED)
    add_test(NAME fuelsim_reproducible_build_tests
        COMMAND "${FUELSIM_BASH_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/reproducible_build_tests.sh"
            "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    set_tests_properties(fuelsim_reproducible_build_tests PROPERTIES
        LABELS "audit;release"
    )

    add_test(NAME fuelsim_geometry_independence_audit
        COMMAND "${CMAKE_COMMAND}"
            "-DROOT=${CMAKE_CURRENT_SOURCE_DIR}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/geometry_independence_audit.cmake"
    )
    set_tests_properties(fuelsim_geometry_independence_audit PROPERTIES
        LABELS "audit;architecture"
    )

    get_property(FUELSIM_REGISTERED_TESTS DIRECTORY PROPERTY TESTS)
    list(JOIN FUELSIM_REGISTERED_TESTS "\n" FUELSIM_REGISTERED_TEST_TEXT)
    file(GENERATE
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fuelsim_registered_tests.txt"
        CONTENT "${FUELSIM_REGISTERED_TEST_TEXT}\n"
    )
    add_executable(fuelsim_verification_matrix_tests
        tests/verification_matrix_tests.cpp
    )
    fuelsim_enable_warnings(fuelsim_verification_matrix_tests)
    add_test(NAME fuelsim_verification_matrix_tests
        COMMAND fuelsim_verification_matrix_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/verification_matrix.tsv"
            "${CMAKE_CURRENT_SOURCE_DIR}"
            "${CMAKE_CURRENT_BINARY_DIR}/fuelsim_registered_tests.txt"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/c3d8t_validation_contract.tsv"
    )
    set_tests_properties(fuelsim_verification_matrix_tests PROPERTIES
        LABELS "audit;release"
    )
    add_test(NAME fuelsim_moose_reference_sha256
        COMMAND "${CMAKE_COMMAND}"
            "-DROOT=${CMAKE_CURRENT_SOURCE_DIR}/verification/moose"
            "-DMANIFEST=${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/SHA256SUMS"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifySha256.cmake"
    )
    set_tests_properties(fuelsim_moose_reference_sha256 PROPERTIES
        LABELS "audit;release"
    )
    add_test(NAME fuelsim_abaqus_reference_sha256
        COMMAND "${CMAKE_COMMAND}"
            "-DROOT=${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus"
            "-DMANIFEST=${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/SHA256SUMS"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/VerifySha256.cmake"
    )
    set_tests_properties(fuelsim_abaqus_reference_sha256 PROPERTIES
        LABELS "audit;release"
    )

    # The PETSc supplied by the MOOSE environment links OpenMP-enabled HYPRE
    # and OpenBLAS.  CTest otherwise inherits their hardware-sized thread
    # teams for every MPI rank, which makes the small regression problems much
    # slower and their timing dependent on the invoking shell.
    set(FUELSIM_CTEST_THREAD_ENVIRONMENT
        "OMP_NUM_THREADS=1"
        "OPENBLAS_NUM_THREADS=1"
        "MKL_NUM_THREADS=1"
        "NUMEXPR_NUM_THREADS=1"
    )
    get_property(FUELSIM_CTEST_ALL_TESTS DIRECTORY PROPERTY TESTS)
    set_tests_properties(${FUELSIM_CTEST_ALL_TESTS} PROPERTIES
        ENVIRONMENT "${FUELSIM_CTEST_THREAD_ENVIRONMENT}"
    )

endif()
