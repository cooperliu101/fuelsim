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

find_package(adlite 0.2.1 CONFIG REQUIRED)
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

add_library(fuelsim_core
    src/core/rz_assembly.cpp
    src/core/cartesian3d_assembly.cpp
    src/core/cartesian3d_contact.cpp
    src/core/cartesian3d_hex20_contact.cpp
    src/core/contact_search.cpp
    src/core/detail/cartesian3d_mechanics.cpp
    src/core/cartesian3d_kernels.cpp
    src/core/cartesian3d_hex20.cpp
    src/core/rz_kernels.cpp
    src/core/mesh.cpp
    src/core/material.cpp
    src/core/material_functions.cpp
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
target_link_libraries(fuelsim_core PUBLIC adlite::adlite)
fuelsim_enable_warnings(fuelsim_core)

add_library(fuelsim_io
    src/io/case_input.cpp
    src/io/checkpoint.cpp
    src/io/problem_signature.cpp
    src/io/results_io.cpp
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
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")

    add_library(fuelsim_test_mesh_fixture
        tests/support/mesh_fixture.cpp
    )
    target_include_directories(fuelsim_test_mesh_fixture
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    )
    target_link_libraries(fuelsim_test_mesh_fixture PUBLIC fuelsim::core)
    fuelsim_enable_warnings(fuelsim_test_mesh_fixture)

    add_library(fuelsim_moose_test_support
        tests/support/moose_field_comparison.cpp
    )
    target_include_directories(fuelsim_moose_test_support
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    )
    target_link_libraries(fuelsim_moose_test_support
        PUBLIC fuelsim::core
    )
    fuelsim_enable_warnings(fuelsim_moose_test_support)

    add_library(fuelsim_abaqus_hex8_test_support
        tests/support/abaqus_hex8_full_field.cpp
    )
    target_include_directories(fuelsim_abaqus_hex8_test_support
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    )
    target_link_libraries(fuelsim_abaqus_hex8_test_support
        PUBLIC fuelsim_moose_test_support fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_abaqus_hex8_test_support)

    add_executable(fuelsim_m1_single_core_benchmark
        benchmarks/m1_single_core.cpp
    )
    target_link_libraries(fuelsim_m1_single_core_benchmark
        PRIVATE fuelsim_test_mesh_fixture fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m1_single_core_benchmark)

    add_executable(fuelsim_input_tests tests/input_tests.cpp)
    target_link_libraries(fuelsim_input_tests PRIVATE fuelsim::io)
    fuelsim_enable_warnings(fuelsim_input_tests)
    add_test(NAME fuelsim_input_tests
        COMMAND fuelsim_input_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_fuel_cladding.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_fuel_cladding_pcmi.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_finite_strain_pcmi.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_coupled_displacement_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_coupled_traction_moose.fsi"
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
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_fuel_cladding.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_fuel_cladding_pcmi.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_noncoaxial_finite_strain.fsi"
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
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m3_load_boundary_tests)
    add_test(NAME fuelsim_m3_load_boundary_tests
        COMMAND fuelsim_m3_load_boundary_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_heat_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_table_convection_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m31_transient_table_convection_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_fuel_cladding_pcmi.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_pressure_boundary.fsi"
    )
    set_tests_properties(fuelsim_m3_load_boundary_tests PROPERTIES
        LABELS "acceptance"
    )

    add_executable(fuelsim_m3_contact_tests
        tests/m3_contact_tests.cpp
    )
    target_link_libraries(fuelsim_m3_contact_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m3_contact_tests)
    add_test(NAME fuelsim_m3_contact_tests
        COMMAND fuelsim_m3_contact_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_two_pellet_contact_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m33_two_pellet_contact_rz_out_all_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m33_two_pellet_contact_rz_out_lower_surface_0001.csv"
    )
    set_tests_properties(fuelsim_m3_contact_tests PROPERTIES
        LABELS "acceptance"
    )

    add_executable(fuelsim_m33_multi_contact_moose_tests
        tests/m33_multi_contact_moose_tests.cpp
    )
    target_link_libraries(fuelsim_m33_multi_contact_moose_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m33_multi_contact_moose_tests)
    add_test(NAME fuelsim_m33_multi_contact_moose_tests
        COMMAND fuelsim_m33_multi_contact_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_multi_contact_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m33_multi_contact_rz_out_all_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m33_multi_contact_rz_out_pellet_contact_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m33_multi_contact_rz_out_inner_clad_contact_0001.csv"
    )
    set_tests_properties(fuelsim_m33_multi_contact_moose_tests PROPERTIES
        LABELS "acceptance;contact"
    )
    add_test(NAME fuelsim_m33_multi_nonmatching_friction_moose_tests
        COMMAND fuelsim_m33_multi_contact_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_multi_nonmatching_friction_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m33_multi_nonmatching_friction_rz_out_all_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m33_multi_nonmatching_friction_rz_out_pellet_contact_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m33_multi_nonmatching_friction_rz_out_inner_clad_contact_0001.csv"
    )
    set_tests_properties(fuelsim_m33_multi_nonmatching_friction_moose_tests PROPERTIES
        LABELS "acceptance;contact;friction"
    )

    add_executable(fuelsim_m51_coulomb_friction_moose_tests
        tests/m51_coulomb_friction_moose_tests.cpp
    )
    target_link_libraries(fuelsim_m51_coulomb_friction_moose_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m51_coulomb_friction_moose_tests)
    add_test(NAME fuelsim_m51_coulomb_friction_moose_tests
        COMMAND fuelsim_m51_coulomb_friction_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_coulomb_friction_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m51_coulomb_friction_rz_all_nodes_0020.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m51_coulomb_friction_rz_fuel_surface_0020.csv"
    )
    set_tests_properties(fuelsim_m51_coulomb_friction_moose_tests PROPERTIES
        LABELS "acceptance"
    )

    add_executable(fuelsim_m52_large_sliding_contact_moose_tests
        tests/m52_large_sliding_contact_moose_tests.cpp
    )
    target_link_libraries(fuelsim_m52_large_sliding_contact_moose_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m52_large_sliding_contact_moose_tests)
    add_test(NAME fuelsim_m52_large_sliding_contact_moose_tests
        COMMAND fuelsim_m52_large_sliding_contact_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_large_sliding_contact_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m52_large_sliding_contact_rz_all_nodes_0020.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m52_large_sliding_contact_rz_lower_surface_0020.csv"
    )
    set_tests_properties(fuelsim_m52_large_sliding_contact_moose_tests PROPERTIES
        LABELS "acceptance"
    )

    add_executable(fuelsim_m57_integrated_validation_tests
        tests/m57_integrated_validation_tests.cpp
    )
    target_link_libraries(fuelsim_m57_integrated_validation_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m57_integrated_validation_tests)
    add_test(NAME fuelsim_m57_integrated_validation_tests
        COMMAND fuelsim_m57_integrated_validation_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_integrated_fuel_cladding_ctest.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m57_integrated_fuel_cladding_rz_ctest_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m57_integrated_fuel_cladding_rz_ctest_contact_pressure_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m57_integrated_fuel_cladding_rz_ctest_clad_qp_coordinates_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m57_integrated_fuel_cladding_rz_ctest_clad_qp_values_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m57_integrated_fuel_cladding_rz_ctest_out.csv"
    )
    set_tests_properties(fuelsim_m57_integrated_validation_tests PROPERTIES
        LABELS "acceptance"
        TIMEOUT 180
    )

    add_executable(fuelsim_m58_integrated_hex8_benchmark
        benchmarks/m58_integrated_hex8.cpp
    )
    target_link_libraries(fuelsim_m58_integrated_hex8_benchmark
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m58_integrated_hex8_benchmark)

    add_executable(fuelsim_m52_large_sliding_restart_tests
        tests/m52_large_sliding_restart_tests.cpp
    )
    target_link_libraries(fuelsim_m52_large_sliding_restart_tests
        PRIVATE fuelsim::io fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m52_large_sliding_restart_tests)
    add_test(NAME fuelsim_m52_large_sliding_restart_tests
        COMMAND fuelsim_m52_large_sliding_restart_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_large_sliding_contact_restart.fsi"
            "${CMAKE_CURRENT_BINARY_DIR}/m52_checkpoint.bin"
    )
    set_tests_properties(fuelsim_m52_large_sliding_restart_tests PROPERTIES
        LABELS "integration"
    )

    add_executable(fuelsim_m1_exodus_moose_tests
        tests/m1_exodus_moose_tests.cpp
    )
    target_link_libraries(fuelsim_m1_exodus_moose_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m1_exodus_moose_tests)
    add_test(NAME fuelsim_m1_exodus_moose_tests
        COMMAND fuelsim_m1_exodus_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_fuel_cladding.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m1_fuel_cladding_gap_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m1_fuel_surface_final.csv"
    )
    set_tests_properties(fuelsim_m1_exodus_moose_tests PROPERTIES
        LABELS "acceptance"
    )

    add_executable(fuelsim_m54_augmented_contact_moose_tests
        tests/m54_augmented_contact_moose_tests.cpp
    )
    target_link_libraries(fuelsim_m54_augmented_contact_moose_tests
        PRIVATE
            fuelsim::io
            fuelsim::solver
            fuelsim_moose_test_support
    )
    fuelsim_enable_warnings(fuelsim_m54_augmented_contact_moose_tests)
    add_test(NAME fuelsim_m54_augmented_contact_moose_tests
        COMMAND fuelsim_m54_augmented_contact_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_augmented_contact_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m54_augmented_contact_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m54_augmented_contact_rz_fuel_surface_final.csv"
    )
    set_tests_properties(fuelsim_m54_augmented_contact_moose_tests PROPERTIES
        LABELS "acceptance"
    )

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
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_fuel_cladding_pcmi.fsi"
    )
    set_tests_properties(fuelsim_m54_transient_augmented_contact_tests
        PROPERTIES LABELS "integration"
    )

    add_executable(fuelsim_m1_unstructured_moose_tests
        tests/m1_unstructured_moose_tests.cpp
    )
    target_link_libraries(fuelsim_m1_unstructured_moose_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m1_unstructured_moose_tests)
    add_test(NAME fuelsim_m1_unstructured_moose_tests
        COMMAND fuelsim_m1_unstructured_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_fuel_cladding_unstructured.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m1_fuel_cladding_unstructured_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m1_fuel_cladding_unstructured_rz_fuel_surface_final.csv"
    )
    set_tests_properties(fuelsim_m1_unstructured_moose_tests PROPERTIES
        LABELS "acceptance"
    )

    add_executable(fuelsim_m0_moose_tests tests/m0_moose_tests.cpp)
    target_link_libraries(fuelsim_m0_moose_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m0_moose_tests)
    add_test(NAME fuelsim_m0_moose_tests
        COMMAND fuelsim_m0_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_single_fuel_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m0_simple_fuel_rz_all_nodes_final.csv"
    )
    set_tests_properties(fuelsim_m0_moose_tests PROPERTIES
        LABELS "acceptance"
    )

    add_executable(fuelsim_core_tests tests/core_tests.cpp)
    target_link_libraries(fuelsim_core_tests
        PRIVATE fuelsim_test_mesh_fixture fuelsim::core
    )
    fuelsim_enable_warnings(fuelsim_core_tests)
    add_test(NAME fuelsim_core_tests COMMAND fuelsim_core_tests)
    set_tests_properties(fuelsim_core_tests PROPERTIES LABELS "unit")

    add_executable(fuelsim_hex8_tests tests/hex8_tests.cpp)
    target_link_libraries(fuelsim_hex8_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_hex8_tests)
    add_test(NAME fuelsim_hex8_tests COMMAND fuelsim_hex8_tests)
    set_tests_properties(fuelsim_hex8_tests PROPERTIES LABELS "unit;three-dimensional")

    add_executable(fuelsim_hex20_tests tests/hex20_tests.cpp)
    target_link_libraries(fuelsim_hex20_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_hex20_tests)
    add_test(NAME fuelsim_hex20_tests COMMAND fuelsim_hex20_tests)
    set_tests_properties(fuelsim_hex20_tests
        PROPERTIES LABELS "unit;three-dimensional;quadratic;warped-geometry;contact;friction")

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
        PRIVATE fuelsim::io fuelsim_moose_test_support)
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
        PRIVATE fuelsim::io fuelsim_moose_test_support)
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

    add_executable(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests
        tests/b523_hex8_c3d8t_integrated_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests
        PRIVATE fuelsim_abaqus_hex8_test_support fuelsim_moose_test_support fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests)
    add_test(NAME fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests
        COMMAND fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_energy.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus")
    set_tests_properties(fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;finite-strain;plasticity;creep;friction;contact;abaqus")

    set(b524_b525_response_directory "${CMAKE_CURRENT_BINARY_DIR}/b524_b525_responses")
    set(b524_b525_cases
        b524_mesh_coarse
        b524_mesh_medium
        b524_mesh_fine
        b524_time_coarse
        b524_time_fine
        b524_penalty_low
        b524_penalty_high
        b524_friction_low
        b524_friction_high
        b524_slip_low
        b524_slip_high
        b524_thermal_low
        b524_thermal_high
        b525_poisson_030
        b525_poisson_045
        b525_poisson_049
        b525_poisson_0499
        b525_poisson_0499_refined)
    set(b524_b525_functional_cases
        b524_mesh_coarse
        b525_poisson_049
        b525_poisson_0499
        b525_poisson_0499_refined)
    set(b524_b525_response_fixtures)
    foreach(b524_b525_case IN LISTS b524_b525_cases)
        set(b524_b525_test "fuelsim_${b524_b525_case}_abaqus_tests")
        set(b524_b525_fixture "${b524_b525_case}_response")
        add_test(NAME ${b524_b525_test}
            COMMAND fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_nodal.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_integration.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_contact.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_energy.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus"
                "${b524_b525_case}"
                "${b524_b525_response_directory}/${b524_b525_case}.txt")
        set(b524_b525_labels
            acceptance three-dimensional hex8 transient thermal mechanical finite-strain friction contact abaqus split-case)
        list(FIND b524_b525_functional_cases "${b524_b525_case}" b524_b525_functional_index)
        if(b524_b525_functional_index EQUAL -1)
            list(APPEND b524_b525_labels qualification)
        endif()
        set_tests_properties(${b524_b525_test} PROPERTIES
            FIXTURES_SETUP "${b524_b525_fixture}"
            LABELS "${b524_b525_labels}")
        list(APPEND b524_b525_response_fixtures "${b524_b525_fixture}")
    endforeach()
    add_test(NAME fuelsim_b524_b525_response_aggregate_tests
        COMMAND fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_energy.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus"
            --aggregate-responses
            "${b524_b525_response_directory}")
    set_tests_properties(fuelsim_b524_b525_response_aggregate_tests PROPERTIES
        FIXTURES_REQUIRED "${b524_b525_response_fixtures}"
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;finite-strain;friction;contact;abaqus;aggregate;qualification")

    foreach(b526_case IN ITEMS b526_contact_cycle b526_friction_reversal)
        string(REPLACE "b526_" "" b526_suffix "${b526_case}")
        set(b526_test "fuelsim_b526_hex8_c3d8t_${b526_suffix}_abaqus_tests")
        add_test(NAME ${b526_test}
            COMMAND fuelsim_b523_hex8_c3d8t_integrated_abaqus_tests
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_nodal.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_integration.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_contact.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b523_hex8_c3d8t_integrated_path_energy.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus"
                "${b526_case}")
        set_tests_properties(${b526_test} PROPERTIES
            LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;finite-strain;friction;contact;abaqus;transition")
    endforeach()

    add_executable(fuelsim_b527_hex8_engineering_fuel_clad_tests
        tests/b527_hex8_engineering_fuel_clad_tests.cpp)
    target_link_libraries(fuelsim_b527_hex8_engineering_fuel_clad_tests
        PRIVATE fuelsim_abaqus_hex8_test_support fuelsim_moose_test_support fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b527_hex8_engineering_fuel_clad_tests)
    set(b527_response_directory "${CMAKE_CURRENT_BINARY_DIR}/b527_responses")
    set(b527_cases coarse medium fine half_step low_penalty high_penalty)
    set(b527_response_fixtures)
    foreach(b527_case IN LISTS b527_cases)
        set(b527_test "fuelsim_b527_${b527_case}_abaqus_tests")
        set(b527_fixture "b527_${b527_case}_response")
        set(b527_labels
            acceptance three-dimensional hex8 transient thermal mechanical finite-strain friction contact engineering abaqus split-case)
        if(NOT b527_case STREQUAL "medium")
            list(APPEND b527_labels qualification)
        endif()
        add_test(NAME ${b527_test}
            COMMAND fuelsim_b527_hex8_engineering_fuel_clad_tests
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus"
                "${b527_case}"
                "${b527_response_directory}/${b527_case}.txt")
        set_tests_properties(${b527_test} PROPERTIES
            FIXTURES_SETUP "${b527_fixture}"
            LABELS "${b527_labels}")
        list(APPEND b527_response_fixtures "${b527_fixture}")
    endforeach()
    add_test(NAME fuelsim_b527_response_aggregate_tests
        COMMAND fuelsim_b527_hex8_engineering_fuel_clad_tests
            --aggregate-responses
            "${b527_response_directory}")
    set_tests_properties(fuelsim_b527_response_aggregate_tests PROPERTIES
        FIXTURES_REQUIRED "${b527_response_fixtures}"
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;finite-strain;friction;contact;engineering;abaqus;aggregate;qualification")

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

    add_executable(fuelsim_b55_hex8_c3d8t_transient_abaqus_tests
        tests/b55_hex8_c3d8t_transient_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b55_hex8_c3d8t_transient_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b55_hex8_c3d8t_transient_abaqus_tests)
    add_test(NAME fuelsim_b55_hex8_c3d8t_transient_abaqus_tests
        COMMAND fuelsim_b55_hex8_c3d8t_transient_abaqus_tests
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

    add_executable(fuelsim_b58_hex8_c3d8t_multistep_abaqus_tests
        tests/b58_hex8_c3d8t_multistep_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b58_hex8_c3d8t_multistep_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b58_hex8_c3d8t_multistep_abaqus_tests)
    add_test(NAME fuelsim_b58_hex8_c3d8t_multistep_abaqus_tests
        COMMAND fuelsim_b58_hex8_c3d8t_multistep_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b58_hex8_c3d8t_multistep_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b58_hex8_c3d8t_multistep_integration.csv")
    set_tests_properties(fuelsim_b58_hex8_c3d8t_multistep_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;temperature-dependent;abaqus")

    add_executable(fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests
        tests/b59_hex8_c3d8t_multimaterial_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests)
    add_test(NAME fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests
        COMMAND fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b59_hex8_c3d8t_multimaterial_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b59_hex8_c3d8t_multimaterial_integration.csv")
    set_tests_properties(fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;thermal;mechanical;multimaterial;refinement;abaqus")

    add_executable(fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests
        tests/b510_hex8_c3d8t_small_j2_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests)
    add_test(NAME fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests
        COMMAND fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests
            plastic
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b510_hex8_c3d8t_small_j2_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b510_hex8_c3d8t_small_j2_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b510_hex8_c3d8t_small_j2_energy.csv")
    set_tests_properties(fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;abaqus")

    add_test(NAME fuelsim_b512_hex8_c3d8t_small_coupled_abaqus_tests
        COMMAND fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests
            coupled
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b512_hex8_c3d8t_small_coupled_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b512_hex8_c3d8t_small_coupled_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b512_hex8_c3d8t_small_coupled_energy.csv")
    set_tests_properties(fuelsim_b512_hex8_c3d8t_small_coupled_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;creep;abaqus")

    add_test(NAME fuelsim_b513_hex8_c3d8t_small_noncoaxial_abaqus_tests
        COMMAND fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests
            noncoaxial
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b513_hex8_c3d8t_small_noncoaxial_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b513_hex8_c3d8t_small_noncoaxial_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b513_hex8_c3d8t_small_noncoaxial_energy.csv")
    set_tests_properties(fuelsim_b513_hex8_c3d8t_small_noncoaxial_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;creep;noncoaxial;reversal;abaqus")

    add_executable(fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests
        tests/b514_hex8_c3d8t_finite_elastic_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests)
    add_test(NAME fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests
        COMMAND fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b514_hex8_c3d8t_finite_elastic_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b514_hex8_c3d8t_finite_elastic_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b514_hex8_c3d8t_finite_elastic_energy.csv")
    set_tests_properties(fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;elasticity;finite-strain;rotation;abaqus")

    add_test(NAME fuelsim_b515_hex8_c3d8t_finite_j2_abaqus_tests
        COMMAND fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests
            finite_plastic
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b515_hex8_c3d8t_finite_j2_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b515_hex8_c3d8t_finite_j2_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b515_hex8_c3d8t_finite_j2_energy.csv")
    set_tests_properties(fuelsim_b515_hex8_c3d8t_finite_j2_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;finite-strain;reversal;abaqus")

    add_test(NAME fuelsim_b517_hex8_c3d8t_finite_coupled_abaqus_tests
        COMMAND fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests
            finite_coupled
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b517_hex8_c3d8t_finite_coupled_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b517_hex8_c3d8t_finite_coupled_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b517_hex8_c3d8t_finite_coupled_energy.csv")
    set_tests_properties(fuelsim_b517_hex8_c3d8t_finite_coupled_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;creep;finite-strain;abaqus")

    add_test(NAME fuelsim_b518_hex8_c3d8t_finite_noncoaxial_abaqus_tests
        COMMAND fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests
            finite_noncoaxial
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b518_hex8_c3d8t_finite_noncoaxial_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b518_hex8_c3d8t_finite_noncoaxial_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b518_hex8_c3d8t_finite_noncoaxial_energy.csv")
    set_tests_properties(fuelsim_b518_hex8_c3d8t_finite_noncoaxial_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;plasticity;creep;finite-strain;noncoaxial;reversal;abaqus")

    add_executable(fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests
        tests/b519_hex8_c3d8t_finite_selective_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::core)
    fuelsim_enable_warnings(fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests)
    add_test(NAME fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests
        COMMAND fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b519_hex8_c3d8t_finite_selective_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b519_hex8_c3d8t_finite_selective_integration.csv")
    set_tests_properties(fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;mechanical;elasticity;finite-strain;selective-integration;abaqus")

    add_executable(fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests
        tests/b511_hex8_c3d8t_small_norton_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests)
    add_test(NAME fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests
        COMMAND fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests
            small
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b511_hex8_c3d8t_small_norton_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b511_hex8_c3d8t_small_norton_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b511_hex8_c3d8t_small_norton_energy.csv")
    set_tests_properties(fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;creep;abaqus")

    add_test(NAME fuelsim_b516_hex8_c3d8t_finite_norton_abaqus_tests
        COMMAND fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests
            finite
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b516_hex8_c3d8t_finite_norton_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b516_hex8_c3d8t_finite_norton_integration.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b516_hex8_c3d8t_finite_norton_energy.csv")
    set_tests_properties(fuelsim_b516_hex8_c3d8t_finite_norton_abaqus_tests PROPERTIES
        LABELS "acceptance;three-dimensional;hex8;transient;mechanical;creep;finite-strain;abaqus")

    add_executable(fuelsim_b39_hex8_sts_multicase_abaqus_tests
        tests/b39_hex8_sts_multicase_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b39_hex8_sts_multicase_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b39_hex8_sts_multicase_abaqus_tests)
    add_test(NAME fuelsim_b39_hex8_sts_multicase_abaqus_tests
        COMMAND fuelsim_b39_hex8_sts_multicase_abaqus_tests
            matching
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_sts_matching_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_matching_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_matching_force.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_matching_reaction.csv"
            nonmatching
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_sts_nonmatching_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_nonmatching_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_nonmatching_force.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_nonmatching_reaction.csv"
            tilted_gap
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_sts_tilted_gap_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_tilted_gap_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_tilted_gap_force.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b39_hex8_sts_tilted_gap_reaction.csv")
    set_tests_properties(fuelsim_b39_hex8_sts_multicase_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;nonmatching;tilted;abaqus")

    add_executable(fuelsim_b40_hex8_sts_friction_abaqus_tests
        tests/b40_hex8_sts_friction_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b40_hex8_sts_friction_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b40_hex8_sts_friction_abaqus_tests)
    add_test(NAME fuelsim_b40_hex8_sts_friction_abaqus_tests
        COMMAND fuelsim_b40_hex8_sts_friction_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_sts_friction_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b40_hex8_sts_friction_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b40_hex8_sts_friction_reaction.csv"
            "${CMAKE_CURRENT_BINARY_DIR}/b40_hex8_sts_friction_checkpoint.bin")
    set_tests_properties(fuelsim_b40_hex8_sts_friction_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;friction;restart;abaqus")

    add_executable(fuelsim_b41_hex8_finite_sliding_abaqus_tests
        tests/b41_hex8_finite_sliding_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b41_hex8_finite_sliding_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b41_hex8_finite_sliding_abaqus_tests)
    add_test(NAME fuelsim_b41_hex8_finite_sliding_abaqus_tests
        COMMAND fuelsim_b41_hex8_finite_sliding_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b41_hex8_finite_sliding_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b41_hex8_finite_sliding_small.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b41_hex8_finite_sliding_finite.csv")
    set_tests_properties(fuelsim_b41_hex8_finite_sliding_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;friction;finite-sliding;abaqus")
    add_test(NAME fuelsim_b41_hex8_finite_sliding_identification_tests
        COMMAND fuelsim_b38_hex8_sts_identification_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b41_hex8_sts_finite_sliding_operator.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b41_hex8_sts_finite_sliding_contact_summary.txt"
            finite-sliding-probe)
    set_tests_properties(fuelsim_b41_hex8_finite_sliding_identification_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;finite-sliding;abaqus;identification")

    add_executable(fuelsim_b42_hex8_friction_objectivity_abaqus_tests
        tests/b42_hex8_friction_objectivity_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b42_hex8_friction_objectivity_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b42_hex8_friction_objectivity_abaqus_tests)
    add_test(NAME fuelsim_b42_hex8_friction_objectivity_abaqus_tests
        COMMAND fuelsim_b42_hex8_friction_objectivity_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b40_hex8_sts_friction_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b42_hex8_sts_friction_objectivity.csv")
    set_tests_properties(fuelsim_b42_hex8_friction_objectivity_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;friction;finite-sliding;objectivity;abaqus")

    add_executable(fuelsim_b43_hex8_finite_strain_contact_abaqus_tests
        tests/b43_hex8_finite_strain_contact_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b43_hex8_finite_strain_contact_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b43_hex8_finite_strain_contact_abaqus_tests)
    add_test(NAME fuelsim_b43_hex8_finite_strain_contact_abaqus_tests
        COMMAND fuelsim_b43_hex8_finite_strain_contact_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_finite_strain_contact_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b43_hex8_finite_strain_contact_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b43_hex8_finite_strain_contact.csv")
    set_tests_properties(fuelsim_b43_hex8_finite_strain_contact_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;finite-sliding;finite-strain;abaqus")

    add_executable(fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests
        tests/b44_hex8_nonmatching_partial_contact_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b44_hex8_nonmatching_partial_contact_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
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
        PRIVATE fuelsim_moose_test_support fuelsim::io)
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
        PRIVATE fuelsim_moose_test_support fuelsim::io)
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
        PRIVATE fuelsim_moose_test_support fuelsim::io)
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

    add_executable(fuelsim_b48_hex8_multi_contact_path_abaqus_tests
        tests/b48_hex8_multi_contact_path_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b48_hex8_multi_contact_path_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b48_hex8_multi_contact_path_abaqus_tests)
    add_test(NAME fuelsim_b48_hex8_multi_contact_path_abaqus_tests
        COMMAND fuelsim_b48_hex8_multi_contact_path_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_multi_contact_path_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b48_hex8_multi_contact_pair_a.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b48_hex8_multi_contact_pair_b.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b48_hex8_multi_contact_reaction.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b48_hex8_multi_contact_energy.csv"
            "${CMAKE_CURRENT_BINARY_DIR}/b48_hex8_multi_contact_checkpoint.bin"
            "${CMAKE_CURRENT_BINARY_DIR}/b48_hex8_multi_contact_one_rank_state.txt")
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

    add_executable(fuelsim_b6_hex20_moose_tests tests/b6_hex20_moose_tests.cpp)
    target_link_libraries(fuelsim_b6_hex20_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b6_hex20_moose_tests)
    add_test(NAME fuelsim_b6_hex20_moose_tests
        COMMAND fuelsim_b6_hex20_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_u2_t1_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b6_hex20_u2_t1_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b6_hex20_u2_t1_displacement_nodes_0001.csv")
    set_tests_properties(fuelsim_b6_hex20_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic")

    add_executable(fuelsim_hex20_moose_case_tests tests/hex20_moose_case_tests.cpp)
    target_link_libraries(fuelsim_hex20_moose_case_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_hex20_moose_case_tests)
    add_test(NAME fuelsim_h20_02_convection_moose_tests
        COMMAND fuelsim_hex20_moose_case_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_convection_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_02_convection_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_02_convection_displacement_nodes_0001.csv"
            "1e-3" temperature_only)
    set_tests_properties(fuelsim_h20_02_convection_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic")
    add_test(NAME fuelsim_h20_03_pressure_moose_tests
        COMMAND fuelsim_hex20_moose_case_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_pressure_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_03_pressure_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_03_pressure_displacement_nodes_0001.csv"
            "1e-5")
    set_tests_properties(fuelsim_h20_03_pressure_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic")
    add_test(NAME fuelsim_h20_04_traction_moose_tests
        COMMAND fuelsim_hex20_moose_case_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_traction_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_04_traction_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_04_traction_displacement_nodes_0001.csv"
            "1e-5")
    set_tests_properties(fuelsim_h20_04_traction_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic")
    add_test(NAME fuelsim_h20_05_thermo_mechanical_moose_tests
        COMMAND fuelsim_hex20_moose_case_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_thermo_mechanical_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_05_thermo_mechanical_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_05_thermo_mechanical_displacement_nodes_0001.csv"
            "1e-5")
    set_tests_properties(fuelsim_h20_05_thermo_mechanical_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic")

    add_executable(fuelsim_hex20_transient_moose_tests tests/hex20_transient_moose_tests.cpp)
    target_link_libraries(fuelsim_hex20_transient_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_hex20_transient_moose_tests)
    add_test(NAME fuelsim_h20_06_transient_moose_tests
        COMMAND fuelsim_hex20_transient_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_06_transient_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_06_transient_displacement_nodes_0001.csv"
            "1e-3")
    set_tests_properties(fuelsim_h20_06_transient_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic;transient")

    add_executable(fuelsim_hex20_inelastic_moose_tests tests/hex20_inelastic_moose_tests.cpp)
    target_link_libraries(fuelsim_hex20_inelastic_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_hex20_inelastic_moose_tests)
    add_test(NAME fuelsim_h20_07_08_09_inelastic_moose_tests
        COMMAND fuelsim_hex20_inelastic_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_inelastic_plastic_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_plastic_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_plastic_displacement_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_plastic_element_state_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_inelastic_creep_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_creep_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_creep_displacement_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_creep_element_state_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_inelastic_coupled_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_coupled_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_coupled_displacement_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_07_inelastic_coupled_element_state_0001.csv")
    set_tests_properties(fuelsim_h20_07_08_09_inelastic_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic;transient;inelastic")
    add_test(NAME fuelsim_h20_10_finite_strain_moose_tests
        COMMAND fuelsim_hex20_moose_case_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_finite_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_10_finite_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_10_finite_displacement_nodes_0001.csv"
            "1e-3")
    set_tests_properties(fuelsim_h20_10_finite_strain_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic;finite-strain")
    add_executable(fuelsim_h20_24_hex20_nonmatching_contact_tests tests/h20_nonmatching_contact_tests.cpp)
    target_link_libraries(fuelsim_h20_24_hex20_nonmatching_contact_tests PRIVATE fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_24_hex20_nonmatching_contact_tests)
    add_test(NAME fuelsim_h20_24_hex20_nonmatching_contact_tests
        COMMAND fuelsim_h20_24_hex20_nonmatching_contact_tests)
    set_tests_properties(fuelsim_h20_24_hex20_nonmatching_contact_tests
        PROPERTIES LABELS "integration;three-dimensional;quadratic;contact;nonmatching")

    add_executable(fuelsim_b3_hex8_moose_tests tests/b3_hex8_moose_tests.cpp)
    target_link_libraries(fuelsim_b3_hex8_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b3_hex8_moose_tests)
    add_test(NAME fuelsim_b3_hex8_moose_tests
        COMMAND fuelsim_b3_hex8_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_thermoelastic.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_thermoelastic_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_thermoelastic_element_stress_final.csv")
    set_tests_properties(fuelsim_b3_hex8_moose_tests PROPERTIES LABELS "acceptance;three-dimensional")

    add_executable(fuelsim_b35_hex8_shared_nodes_moose_tests tests/b35_hex8_shared_nodes_moose_tests.cpp)
    target_link_libraries(fuelsim_b35_hex8_shared_nodes_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b35_hex8_shared_nodes_moose_tests)
    add_test(NAME fuelsim_b35_hex8_shared_nodes_moose_tests
        COMMAND fuelsim_b35_hex8_shared_nodes_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_shared_meat_clad_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b35_hex8_shared_meat_clad_all_nodes_0001.csv")
    set_tests_properties(fuelsim_b35_hex8_shared_nodes_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional")

    add_executable(fuelsim_b36_plate_meat_clad_moose_tests tests/b36_plate_meat_clad_moose_tests.cpp)
    target_link_libraries(fuelsim_b36_plate_meat_clad_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b36_plate_meat_clad_moose_tests)
    add_test(NAME fuelsim_b36_plate_meat_clad_moose_tests
        COMMAND fuelsim_b36_plate_meat_clad_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_plate_meat_clad_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b36_plate_meat_clad_all_nodes_0040.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b36_plate_meat_clad_element_state_0040.csv")
    set_tests_properties(fuelsim_b36_plate_meat_clad_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;inelastic")

    add_executable(fuelsim_rz_shared_nodes_moose_tests tests/rz_shared_nodes_moose_tests.cpp)
    target_link_libraries(fuelsim_rz_shared_nodes_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_rz_shared_nodes_moose_tests)
    add_test(NAME fuelsim_rz_shared_nodes_moose_tests
        COMMAND fuelsim_rz_shared_nodes_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_rz_shared_meat_clad_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/rz_shared_meat_clad_out_all_nodes_0001.csv")
    set_tests_properties(fuelsim_rz_shared_nodes_moose_tests PROPERTIES LABELS "acceptance;axisymmetric")

    add_executable(fuelsim_b33_hex8_contact_moose_tests tests/b33_hex8_contact_moose_tests.cpp)
    target_link_libraries(fuelsim_b33_hex8_contact_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b33_hex8_contact_moose_tests)
    add_test(NAME fuelsim_b33_hex8_contact_moose_tests
        COMMAND fuelsim_b33_hex8_contact_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_contact_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b33_hex8_thermal_contact_all_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b33_hex8_contact_all_nodes_0010.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b33_hex8_contact_secondary_contact_0010.csv")
    set_tests_properties(fuelsim_b33_hex8_contact_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;contact;friction")

    add_executable(fuelsim_b34_hex8_sliding_contact_abaqus_tests
        tests/b34_hex8_sliding_contact_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b34_hex8_sliding_contact_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b34_hex8_sliding_contact_abaqus_tests)
    add_test(NAME fuelsim_b34_hex8_sliding_contact_abaqus_tests
        COMMAND fuelsim_b34_hex8_sliding_contact_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_sliding_contact_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b34_hex8_sliding_contact_nodes.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b34_hex8_sliding_contact_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b34_hex8_sliding_contact_reaction.csv"
            "${CMAKE_CURRENT_BINARY_DIR}/b34_hex8_sliding_contact_checkpoint.bin")
    set_tests_properties(fuelsim_b34_hex8_sliding_contact_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex8;contact;friction;restart;abaqus")

    add_executable(fuelsim_h20_hex20_contact_reference_tests tests/h20_hex20_contact_reference_tests.cpp)
    target_link_libraries(fuelsim_h20_hex20_contact_reference_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_hex20_contact_reference_tests)
    add_test(NAME fuelsim_h20_hex20_contact_reference_tests
        COMMAND fuelsim_h20_hex20_contact_reference_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_thermal_contact_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_16_hex20_thermal_contact_temperature_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_mechanical_contact_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_17_hex20_mechanical_contact_displacement_nodes_0010.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_17_hex20_mechanical_contact_secondary_contact_0010.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_sliding_contact_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_18_hex20_sliding_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_mortar_area_compare.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_19_hex20_mortar_mechanical_displacement_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_19_hex20_mortar_mechanical.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_23_hex20_surface_friction_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_23_hex20_surface_friction_reaction.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_21_hex20_surface_contact_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_21_hex20_surface_contact_reaction.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_22_hex20_unnormalized_mechanical_displacement_nodes_0009.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_22_hex20_unnormalized_mechanical_secondary_contact_0009.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_22_hex20_unnormalized_mechanical.csv")
    set_tests_properties(fuelsim_h20_hex20_contact_reference_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction")

    add_executable(fuelsim_h20_25_hex20_node_to_face_three_way_tests
        tests/h20_node_to_face_three_way_tests.cpp)
    target_link_libraries(fuelsim_h20_25_hex20_node_to_face_three_way_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_25_hex20_node_to_face_three_way_tests)
    add_test(NAME fuelsim_h20_25_hex20_node_to_face_three_way_tests
        COMMAND fuelsim_h20_25_hex20_node_to_face_three_way_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_mortar_area_compare.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_25_hex20_node_to_face_displacement_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_25_hex20_node_to_face_secondary_contact_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_25_hex20_node_to_face_reaction.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_25_hex20_node_to_surface_contact_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_25_hex20_node_to_surface_contact_pressure.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_25_hex20_node_to_surface_contact_reaction.csv")
    set_tests_properties(fuelsim_h20_25_hex20_node_to_face_three_way_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;abaqus;moose")

    add_executable(fuelsim_h20_26_abaqus_sts_identification_tests
        tests/h20_abaqus_sts_identification_tests.cpp)
    fuelsim_enable_warnings(fuelsim_h20_26_abaqus_sts_identification_tests)
    add_test(NAME fuelsim_h20_26_abaqus_sts_identification_tests
        COMMAND fuelsim_h20_26_abaqus_sts_identification_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_26_hex20_sts_operator.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_26_hex20_sts_history.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_26_hex20_sts_finite_sliding_probe_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_26_hex20_sts_finite_sliding_probe_history.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_26_hex20_sts_tilted_normal_probe_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_26_hex20_sts_tilted_normal_probe_history.csv")
    set_tests_properties(fuelsim_h20_26_abaqus_sts_identification_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;abaqus;identification")

    add_executable(fuelsim_h20_27_abaqus_sts_nonmatching_identification_tests
        tests/h20_abaqus_sts_nonmatching_identification_tests.cpp)
    fuelsim_enable_warnings(fuelsim_h20_27_abaqus_sts_nonmatching_identification_tests)
    add_test(NAME fuelsim_h20_27_abaqus_sts_nonmatching_identification_tests
        COMMAND fuelsim_h20_27_abaqus_sts_nonmatching_identification_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_26_hex20_sts_operator.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_27_hex20_nonmatching_sts_operator.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_27_hex20_nonmatching_sts_history.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_27_hex20_small_sliding_transfer_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_27_hex20_small_sliding_transfer_history.csv")
    set_tests_properties(fuelsim_h20_27_abaqus_sts_nonmatching_identification_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;abaqus;identification")

    add_executable(fuelsim_h20_28_abaqus_sts_smoothing_identification_tests
        tests/h20_abaqus_sts_smoothing_identification_tests.cpp)
    fuelsim_enable_warnings(fuelsim_h20_28_abaqus_sts_smoothing_identification_tests)
    add_test(NAME fuelsim_h20_28_abaqus_sts_smoothing_identification_tests
        COMMAND fuelsim_h20_28_abaqus_sts_smoothing_identification_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_28_hex20_refined_primary_sts_default_operator.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_28_hex20_refined_primary_sts_linear_operator.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_28_hex20_refined_primary_sts_quadratic_operator.csv")
    set_tests_properties(fuelsim_h20_28_abaqus_sts_smoothing_identification_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;abaqus;identification")

    add_executable(fuelsim_h20_29_hex20_sts_multicase_abaqus_tests
        tests/h20_sts_multicase_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_29_hex20_sts_multicase_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_29_hex20_sts_multicase_abaqus_tests)
    foreach(h20_29_case IN ITEMS biaxial tilted graded soft_penalty initial_gap)
        set(h20_29_test "fuelsim_h20_29_hex20_sts_${h20_29_case}_abaqus_tests")
        add_test(NAME ${h20_29_test}
            COMMAND fuelsim_h20_29_hex20_sts_multicase_abaqus_tests
                "${h20_29_case}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_sts_${h20_29_case}_abaqus.fsi"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_29_hex20_sts_${h20_29_case}_displacement.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_29_hex20_sts_${h20_29_case}_force.csv"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_29_hex20_sts_${h20_29_case}_reaction.csv")
        set_tests_properties(${h20_29_test}
            PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;nonmatching;abaqus;split-case")
    endforeach()

    add_executable(fuelsim_h20_30_hex20_curved_sts_abaqus_tests
        tests/h20_curved_sts_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_30_hex20_curved_sts_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_30_hex20_curved_sts_abaqus_tests)
    add_test(NAME fuelsim_h20_30_hex20_curved_sts_abaqus_tests
        COMMAND fuelsim_h20_30_hex20_curved_sts_abaqus_tests
            faceted_cylinder
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_sts_faceted_cylinder_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_faceted_cylinder_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_faceted_cylinder_force.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_faceted_cylinder_resultant.csv"
            quadratic_cylinder
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_sts_quadratic_cylinder_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_quadratic_cylinder_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_quadratic_cylinder_force.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_30_hex20_sts_quadratic_cylinder_resultant.csv")
    set_tests_properties(fuelsim_h20_30_hex20_curved_sts_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;curved;abaqus")

    add_executable(fuelsim_h20_friction_path_abaqus_tests
        tests/h20_friction_path_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_friction_path_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_friction_path_abaqus_tests)
    add_test(NAME fuelsim_h20_friction_path_abaqus_tests
        COMMAND fuelsim_h20_friction_path_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_nonmatching_friction_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_33_hex20_nonmatching_friction_path_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_33_hex20_nonmatching_friction_path_contact.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_33_hex20_nonmatching_friction_path_reaction.csv"
            "${CMAKE_CURRENT_BINARY_DIR}/h20_friction_path_checkpoint.bin")
    set_tests_properties(fuelsim_h20_friction_path_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction;nonmatching;abaqus;restart")

    add_executable(fuelsim_h20_35_curved_friction_abaqus_tests
        tests/h20_curved_friction_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_35_curved_friction_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_35_curved_friction_abaqus_tests)
    add_test(NAME fuelsim_h20_35_curved_friction_abaqus_tests
        COMMAND fuelsim_h20_35_curved_friction_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_sts_quadratic_cylinder_friction_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_35_hex20_sts_quadratic_cylinder_friction_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_35_hex20_sts_quadratic_cylinder_friction_contact.csv")
    set_tests_properties(fuelsim_h20_35_curved_friction_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction;curved;abaqus")

    add_executable(fuelsim_h20_36_curved_friction_path_abaqus_tests
        tests/h20_curved_friction_path_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_36_curved_friction_path_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_36_curved_friction_path_abaqus_tests)
    add_test(NAME fuelsim_h20_36_curved_friction_path_abaqus_tests
        COMMAND fuelsim_h20_36_curved_friction_path_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex20_sts_quadratic_cylinder_friction_path_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_36_hex20_sts_quadratic_cylinder_friction_path_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_36_hex20_sts_quadratic_cylinder_friction_path_contact.csv"
            "${CMAKE_CURRENT_BINARY_DIR}/h20_curved_friction_path_checkpoint.bin")
    set_tests_properties(fuelsim_h20_36_curved_friction_path_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction;curved;abaqus;restart")

    add_executable(fuelsim_h20_38_friction_objectivity_abaqus_tests
        tests/h20_friction_objectivity_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_38_friction_objectivity_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io)
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
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
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
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_40_curved_finite_sliding_abaqus_tests)
    foreach(h20_40_mode IN ITEMS small finite)
        set(h20_40_test "fuelsim_h20_40_curved_finite_sliding_${h20_40_mode}_abaqus_tests")
        add_test(NAME ${h20_40_test}
            COMMAND fuelsim_h20_40_curved_finite_sliding_abaqus_tests
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_40_hex20_curved_finite_sliding_mesh.e"
                "${h20_40_mode}"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_40_hex20_curved_finite_sliding_${h20_40_mode}.csv")
        set_tests_properties(${h20_40_test}
            PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;friction;curved;finite-sliding;abaqus;split-case")
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

    add_executable(fuelsim_h20_24_hex20_nonmatching_abaqus_tests
        tests/h20_nonmatching_contact_abaqus_tests.cpp)
    target_link_libraries(fuelsim_h20_24_hex20_nonmatching_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_h20_24_hex20_nonmatching_abaqus_tests)
    add_test(NAME fuelsim_h20_24_hex20_nonmatching_abaqus_tests
        COMMAND fuelsim_h20_24_hex20_nonmatching_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex20_nonmatching_contact_abaqus.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_24_hex20_nonmatching_surface_contact_displacement.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_24_hex20_nonmatching_surface_contact_reaction.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/h20_24_hex20_nonmatching_surface_contact_pressure.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_24_hex20_nonmatching_mortar_all_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/h20_24_hex20_nonmatching_mortar.csv")
    set_tests_properties(fuelsim_h20_24_hex20_nonmatching_abaqus_tests
        PROPERTIES LABELS "acceptance;three-dimensional;hex20;contact;nonmatching;abaqus;moose")

    add_executable(fuelsim_b37_hex8_multi_contact_moose_tests tests/b37_hex8_multi_contact_moose_tests.cpp)
    target_link_libraries(fuelsim_b37_hex8_multi_contact_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b37_hex8_multi_contact_moose_tests)
    add_test(NAME fuelsim_b37_hex8_multi_contact_moose_tests
        COMMAND fuelsim_b37_hex8_multi_contact_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_multi_contact_friction_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b37_hex8_multi_contact_all_nodes_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b37_hex8_multi_contact_pair_a_contact_0001.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b37_hex8_multi_contact_pair_b_contact_0001.csv")
    set_tests_properties(fuelsim_b37_hex8_multi_contact_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;contact;friction")

    add_executable(fuelsim_b3_hex8_inelastic_moose_tests tests/b3_hex8_inelastic_moose_tests.cpp)
    target_link_libraries(fuelsim_b3_hex8_inelastic_moose_tests
        PRIVATE fuelsim_moose_test_support fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b3_hex8_inelastic_moose_tests)
    add_test(NAME fuelsim_b3_hex8_inelastic_moose_tests
        COMMAND fuelsim_b3_hex8_inelastic_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_plastic_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_plastic_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_plastic_element_state_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_creep_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_creep_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_creep_element_state_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_coupled_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_coupled_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b3_hex8_coupled_element_state_final.csv")
    set_tests_properties(fuelsim_b3_hex8_inelastic_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;inelastic")

    add_test(NAME fuelsim_b32_hex8_finite_inelastic_moose_tests
        COMMAND fuelsim_b3_hex8_inelastic_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_finite_plastic_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b32_hex8_finite_plastic_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b32_hex8_finite_plastic_element_state_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_finite_creep_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b32_hex8_finite_creep_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b32_hex8_finite_creep_element_state_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_finite_coupled_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b32_hex8_finite_coupled_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b32_hex8_finite_coupled_element_state_final.csv")
    set_tests_properties(fuelsim_b32_hex8_finite_inelastic_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;finite-strain;inelastic")

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
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b528_hex8_c3d8rt_holdout_operator_nodal.csv")
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

    add_executable(fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests
        tests/b531_hex8_c3d8rt_transient_abaqus_tests.cpp)
    target_link_libraries(fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests
        PRIVATE fuelsim_moose_test_support fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests)
    add_test(NAME fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests
        COMMAND fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b531_hex8_c3d8rt_transient_nodal.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/abaqus/b531_hex8_c3d8rt_transient_integration.csv")
    set_tests_properties(fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests PROPERTIES
        LABELS "acceptance;abaqus;three-dimensional;reduced-integration;transient;thermo-mechanical")
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
    if(FUELSIM_MPIEXEC_EXECUTABLE)
        add_test(NAME fuelsim_mpi_hex8_solver_tests
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                    "$<TARGET_FILE:fuelsim_hex8_solver_tests>"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex8_steady_results.e"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex8_transient_results.e"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex8_checkpoint.bin"
                    -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps)
        set_tests_properties(fuelsim_mpi_hex8_solver_tests PROPERTIES
            PROCESSORS 2
            RUN_SERIAL TRUE
            LABELS "integration;mpi;three-dimensional")

        add_test(NAME fuelsim_mpi_hex20_solver_tests
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                    "$<TARGET_FILE:fuelsim_hex20_solver_tests>"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex20_mesh.e"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex20_steady_results.e"
                    "${CMAKE_CURRENT_BINARY_DIR}/mpi_hex20_checkpoint.bin"
                    -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps)
        set_tests_properties(fuelsim_mpi_hex20_solver_tests PROPERTIES
            PROCESSORS 2
            RUN_SERIAL TRUE
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
            PRIVATE fuelsim::io fuelsim::solver
        )
        target_include_directories(fuelsim_mpi_equivalence_benchmark PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
        fuelsim_enable_warnings(fuelsim_mpi_equivalence_benchmark)
        add_test(NAME fuelsim_mpi_reference
            COMMAND fuelsim_mpi_equivalence_benchmark
                write
                "${CMAKE_CURRENT_BINARY_DIR}/m34_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_two_pellet_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_mpi_reference PROPERTIES
            FIXTURES_SETUP m34_mpi_reference
            LABELS "integration;parallel"
        )
        add_test(NAME fuelsim_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${CMAKE_CURRENT_BINARY_DIR}/m34_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_two_pellet_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED m34_mpi_reference
            LABELS "acceptance;parallel"
            PROCESSORS 2
        )

        add_test(NAME fuelsim_m52_mpi_reference
            COMMAND fuelsim_mpi_equivalence_benchmark
                write
                "${CMAKE_CURRENT_BINARY_DIR}/m52_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_large_sliding_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_m52_mpi_reference PROPERTIES
            FIXTURES_SETUP m52_mpi_reference
            LABELS "integration;parallel"
        )
        add_test(NAME fuelsim_m52_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${CMAKE_CURRENT_BINARY_DIR}/m52_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_large_sliding_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_m52_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED m52_mpi_reference
            LABELS "acceptance;parallel"
            PROCESSORS 2
        )
        add_test(NAME fuelsim_b33_mpi_reference
            COMMAND fuelsim_mpi_equivalence_benchmark
                write
                "${CMAKE_CURRENT_BINARY_DIR}/b33_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_b33_mpi_reference PROPERTIES
            FIXTURES_SETUP b33_mpi_reference
            LABELS "integration;parallel;three-dimensional;contact"
        )
        add_test(NAME fuelsim_b33_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${CMAKE_CURRENT_BINARY_DIR}/b33_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_b33_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED b33_mpi_reference
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
        add_test(NAME fuelsim_b39_hex8_sts_mpi_reference
            COMMAND fuelsim_mpi_equivalence_benchmark
                write
                "${CMAKE_CURRENT_BINARY_DIR}/b39_hex8_sts_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_sts_nonmatching_abaqus.fsi"
        )
        set_tests_properties(fuelsim_b39_hex8_sts_mpi_reference PROPERTIES
            FIXTURES_SETUP b39_hex8_sts_mpi_reference
            LABELS "integration;parallel;three-dimensional;hex8;contact;surface-to-surface"
        )
        add_test(NAME fuelsim_b39_hex8_sts_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare
                "${CMAKE_CURRENT_BINARY_DIR}/b39_hex8_sts_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_sts_nonmatching_abaqus.fsi"
        )
        set_tests_properties(fuelsim_b39_hex8_sts_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED b39_hex8_sts_mpi_reference
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
        add_test(NAME fuelsim_b48_hex8_multi_contact_path_mpi_equivalence
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_b48_hex8_multi_contact_path_abaqus_tests>
                --compare-mpi-reference
                "${CMAKE_CURRENT_BINARY_DIR}/b48_hex8_multi_contact_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_hex8_multi_contact_path_abaqus.fsi"
        )
        set_tests_properties(fuelsim_b48_hex8_multi_contact_path_mpi_equivalence PROPERTIES
            FIXTURES_REQUIRED b48_hex8_multi_contact_mpi_reference
            LABELS "acceptance;parallel;three-dimensional;hex8;contact;friction;surface-to-surface;finite-sliding;multi-contact"
            PROCESSORS 2
        )
        add_test(NAME fuelsim_mpi_field_split
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare_field_split
                "${CMAKE_CURRENT_BINARY_DIR}/m34_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_two_pellet_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_mpi_field_split PROPERTIES
            FIXTURES_REQUIRED m34_mpi_reference
            LABELS "acceptance;parallel"
            PROCESSORS 2
            TIMEOUT 60
        )
        add_test(NAME fuelsim_mpi_block_jacobi
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare_block_jacobi
                "${CMAKE_CURRENT_BINARY_DIR}/m34_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_two_pellet_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_mpi_block_jacobi PROPERTIES
            FIXTURES_REQUIRED m34_mpi_reference
            LABELS "integration;parallel"
            PROCESSORS 2
            TIMEOUT 60
        )
        add_test(NAME fuelsim_mpi_hypre
            COMMAND "${FUELSIM_MPIEXEC_EXECUTABLE}" -n 2
                $<TARGET_FILE:fuelsim_mpi_equivalence_benchmark>
                compare_hypre
                "${CMAKE_CURRENT_BINARY_DIR}/m34_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_two_pellet_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_mpi_hypre PROPERTIES
            FIXTURES_REQUIRED m34_mpi_reference
            LABELS "integration;parallel"
            PROCESSORS 2
            TIMEOUT 120
        )
        add_test(NAME fuelsim_transient_mpi_reference
            COMMAND fuelsim_mpi_equivalence_benchmark
                write_transient
                "${CMAKE_CURRENT_BINARY_DIR}/m35_transient_one_rank_state.txt"
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_coupled_traction_moose.fsi"
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
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_coupled_traction_moose.fsi"
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
                "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_two_pellet_contact_moose.fsi"
        )
        set_tests_properties(fuelsim_mpi_root_io_failure PROPERTIES
            LABELS "integration;parallel"
            PROCESSORS 2
            TIMEOUT 30
        )
    endif()

    add_executable(fuelsim_m2_solver_tests tests/m2_solver_tests.cpp)
    target_link_libraries(fuelsim_m2_solver_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m2_solver_tests)
    add_test(NAME fuelsim_m2_solver_tests
        COMMAND fuelsim_m2_solver_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_heat_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m21_transient_heat_rz_all_nodes_final.csv"
    )
    set_tests_properties(fuelsim_m2_solver_tests PROPERTIES LABELS "acceptance")

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

    add_executable(fuelsim_m2_inelastic_solver_tests
        tests/m2_inelastic_solver_tests.cpp
    )
    target_link_libraries(fuelsim_m2_inelastic_solver_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m2_inelastic_solver_tests)
    add_test(NAME fuelsim_m2_inelastic_solver_tests
        COMMAND fuelsim_m2_inelastic_solver_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_j2_plastic_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_norton_creep_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_coupled_displacement_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_coupled_traction_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_j2_plastic_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_norton_creep_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_coupled_plastic_creep_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_coupled_plastic_creep_traction_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_j2_unload_reload_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_j2_unload_reload_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_j2_unload_reload_rz_out.csv"
    )
    set_tests_properties(fuelsim_m2_inelastic_solver_tests
        PROPERTIES LABELS "acceptance"
    )

    add_executable(fuelsim_m2_pcmi_solver_tests
        tests/m2_pcmi_solver_tests.cpp
    )
    target_link_libraries(fuelsim_m2_pcmi_solver_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m2_pcmi_solver_tests)
    add_test(NAME fuelsim_m2_pcmi_solver_tests
        COMMAND fuelsim_m2_pcmi_solver_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_fuel_cladding_pcmi.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m23_pcmi_coupled_cladding_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m23_pcmi_coupled_cladding_rz_fuel_surface_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m23_pcmi_coupled_cladding_rz_clad_qp_coordinates_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m23_pcmi_coupled_cladding_rz_clad_qp_values_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m23_pcmi_coupled_cladding_rz_out.csv"
    )
    set_tests_properties(fuelsim_m2_pcmi_solver_tests
        PROPERTIES LABELS "acceptance"
    )

    add_test(NAME fuelsim_m41_finite_strain_pcmi_tests
        COMMAND fuelsim_m2_pcmi_solver_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_finite_strain_pcmi.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m41_finite_strain_pcmi_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m41_finite_strain_pcmi_rz_fuel_surface_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m41_finite_strain_pcmi_rz_clad_qp_coordinates_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m41_finite_strain_pcmi_rz_clad_qp_values_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m41_finite_strain_pcmi_rz_out.csv"
    )
    set_tests_properties(fuelsim_m41_finite_strain_pcmi_tests
        PROPERTIES LABELS "acceptance"
    )

    add_executable(fuelsim_m42_follower_pressure_moose_tests
        tests/m42_follower_pressure_moose_tests.cpp
    )
    target_link_libraries(fuelsim_m42_follower_pressure_moose_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(fuelsim_m42_follower_pressure_moose_tests)
    add_test(NAME fuelsim_m42_follower_pressure_moose_tests
        COMMAND fuelsim_m42_follower_pressure_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_finite_follower_pressure.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m42_follower_pressure_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_finite_left_top_pressure.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m42_left_top_pressure_rz_all_nodes_final.csv"
    )
    set_tests_properties(fuelsim_m42_follower_pressure_moose_tests
        PROPERTIES LABELS "acceptance"
    )

    add_executable(fuelsim_m43_noncoaxial_finite_strain_moose_tests
        tests/m43_noncoaxial_finite_strain_moose_tests.cpp
    )
    target_link_libraries(fuelsim_m43_noncoaxial_finite_strain_moose_tests
        PRIVATE
            fuelsim_moose_test_support
            fuelsim::io
            fuelsim::solver
    )
    fuelsim_enable_warnings(
        fuelsim_m43_noncoaxial_finite_strain_moose_tests
    )
    add_test(NAME fuelsim_m43_noncoaxial_finite_strain_moose_tests
        COMMAND fuelsim_m43_noncoaxial_finite_strain_moose_tests
            production
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_noncoaxial_finite_strain.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_noncoaxial_finite_strain_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_noncoaxial_finite_strain_rz_element_history.csv"
    )
    set_tests_properties(
        fuelsim_m43_noncoaxial_finite_strain_moose_tests
        PROPERTIES LABELS "acceptance"
    )
    add_test(NAME fuelsim_m43_elastic_displacement_moose_tests
        COMMAND fuelsim_m43_noncoaxial_finite_strain_moose_tests
            elastic_displacement
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_noncoaxial_finite_strain_elastic_displacement.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_elastic_displacement_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_elastic_displacement_rz_element_history.csv"
    )
    add_test(NAME fuelsim_m43_plastic_displacement_moose_tests
        COMMAND fuelsim_m43_noncoaxial_finite_strain_moose_tests
            plastic_displacement
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_noncoaxial_finite_strain_plastic_displacement.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_plastic_displacement_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_plastic_displacement_rz_element_history.csv"
    )
    add_test(NAME fuelsim_m43_creep_displacement_moose_tests
        COMMAND fuelsim_m43_noncoaxial_finite_strain_moose_tests
            creep_displacement
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_noncoaxial_finite_strain_creep_displacement.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_creep_displacement_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_creep_displacement_rz_element_history.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_creep_displacement_rz_nodal_history.csv"
    )
    add_test(NAME fuelsim_m43_coupled_displacement_moose_tests
        COMMAND fuelsim_m43_noncoaxial_finite_strain_moose_tests
            coupled_displacement
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_noncoaxial_finite_strain_coupled_displacement.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_coupled_displacement_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_coupled_displacement_rz_element_history.csv"
    )
    add_test(NAME fuelsim_m43_coupled_pressure_moose_tests
        COMMAND fuelsim_m43_noncoaxial_finite_strain_moose_tests
            coupled_pressure
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_noncoaxial_finite_strain_coupled_pressure.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_coupled_pressure_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_coupled_pressure_rz_element_history.csv"
    )
    add_test(NAME fuelsim_m43_material_oracle_moose_tests
        COMMAND fuelsim_m43_noncoaxial_finite_strain_moose_tests
            material_oracle
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_noncoaxial_finite_strain_material_oracle.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_material_oracle_rz_all_nodes_final.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m43_material_oracle_rz_element_history.csv"
    )
    set_tests_properties(
        fuelsim_m43_elastic_displacement_moose_tests
        fuelsim_m43_plastic_displacement_moose_tests
        fuelsim_m43_creep_displacement_moose_tests
        fuelsim_m43_coupled_displacement_moose_tests
        fuelsim_m43_coupled_pressure_moose_tests
        fuelsim_m43_material_oracle_moose_tests
        PROPERTIES LABELS "acceptance"
    )

    add_test(NAME fuelsim_steady_input_card
        COMMAND fuelsim -i
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_fuel_cladding.fsi"
    )
    set_tests_properties(fuelsim_steady_input_card
        PROPERTIES LABELS "acceptance"
    )

    add_test(NAME fuelsim_transient_input_card
        COMMAND fuelsim -i
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_fuel_cladding_pcmi.fsi"
    )
    set_tests_properties(fuelsim_transient_input_card
        PROPERTIES LABELS "acceptance"
    )

    add_test(NAME fuelsim_transient_jacobian_check
        COMMAND fuelsim -i
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_table_convection_moose.fsi"
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

    # Do not run independent two-rank cases concurrently when callers choose
    # ctest -j: each case already owns two MPI processes and is intended to
    # exercise a distributed PETSc path, not process-level oversubscription.
    if(FUELSIM_MPIEXEC_EXECUTABLE)
        set(FUELSIM_CTEST_TWO_RANK_TESTS
            fuelsim_mpi_hex8_solver_tests
            fuelsim_mpi_hex20_solver_tests
            fuelsim_mpi_solver_tests
            fuelsim_mpi_runtime_layout_field_split
            fuelsim_mpi_equivalence
            fuelsim_m52_mpi_equivalence
            fuelsim_b33_mpi_equivalence
            fuelsim_b521_thermal_contact_mpi_equivalence
            fuelsim_b39_hex8_sts_mpi_equivalence
            fuelsim_b41_hex8_finite_sliding_mpi_equivalence
            fuelsim_b48_hex8_multi_contact_path_mpi_equivalence
            fuelsim_mpi_field_split
            fuelsim_mpi_block_jacobi
            fuelsim_mpi_hypre
            fuelsim_transient_mpi_equivalence
            fuelsim_mpi_root_io_failure
        )
        set_tests_properties(${FUELSIM_CTEST_TWO_RANK_TESTS} PROPERTIES
            RESOURCE_LOCK fuelsim_two_rank_mpi
        )
    endif()
endif()
