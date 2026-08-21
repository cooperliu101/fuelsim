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
    src/core/contact_search.cpp
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
    set_tests_properties(fuelsim_hex20_tests PROPERTIES LABELS "unit;three-dimensional;quadratic")

    add_executable(fuelsim_hex8_solver_tests tests/hex8_solver_tests.cpp)
    target_link_libraries(fuelsim_hex8_solver_tests PRIVATE fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_hex8_solver_tests)
    add_test(NAME fuelsim_hex8_solver_tests
        COMMAND fuelsim_hex8_solver_tests
            "${CMAKE_CURRENT_BINARY_DIR}/hex8_steady_results.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex8_transient_results.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex8_checkpoint.bin")
    set_tests_properties(fuelsim_hex8_solver_tests PROPERTIES LABELS "integration;three-dimensional")

    add_executable(fuelsim_hex20_solver_tests tests/hex20_solver_tests.cpp)
    target_link_libraries(fuelsim_hex20_solver_tests PRIVATE fuelsim::io fuelsim::solver)
    fuelsim_enable_warnings(fuelsim_hex20_solver_tests)
    add_test(NAME fuelsim_hex20_solver_tests
        COMMAND fuelsim_hex20_solver_tests
            "${CMAKE_CURRENT_BINARY_DIR}/hex20_mesh.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex20_steady_results.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex20_checkpoint.bin")
    set_tests_properties(fuelsim_hex20_solver_tests PROPERTIES LABELS "integration;three-dimensional;quadratic")

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
    add_test(NAME fuelsim_h20_11_hex20_checkpoint_tests
        COMMAND fuelsim_hex20_solver_tests
            "${CMAKE_CURRENT_BINARY_DIR}/hex20_mesh.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex20_steady_results.e"
            "${CMAKE_CURRENT_BINARY_DIR}/hex20_checkpoint.bin")
    set_tests_properties(fuelsim_h20_11_hex20_checkpoint_tests
        PROPERTIES LABELS "integration;three-dimensional;quadratic;checkpoint")
    add_test(NAME fuelsim_h20_12_multiblock_shared_nodes_tests COMMAND fuelsim_multi_region_tests)
    set_tests_properties(fuelsim_h20_12_multiblock_shared_nodes_tests
        PROPERTIES LABELS "acceptance;three-dimensional;quadratic;multiblock")
    add_test(NAME fuelsim_h20_13_warped_hex20_geometry_tests COMMAND fuelsim_hex20_tests)
    set_tests_properties(fuelsim_h20_13_warped_hex20_geometry_tests
        PROPERTIES LABELS "unit;three-dimensional;quadratic;warped-geometry")

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
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b33_hex8_contact_secondary_contact_0010.csv"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_hex8_sliding_contact_moose.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/b34_hex8_sliding_contact.csv"
            "${CMAKE_CURRENT_BINARY_DIR}/b34_hex8_sliding_contact_checkpoint.bin")
    set_tests_properties(fuelsim_b33_hex8_contact_moose_tests
        PROPERTIES LABELS "acceptance;three-dimensional;contact;friction")

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

    add_test(NAME fuelsim_source_size_audit
        COMMAND "${CMAKE_COMMAND}"
            "-DROOT=${CMAKE_CURRENT_SOURCE_DIR}"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/source_size_audit.cmake"
    )
    set_tests_properties(fuelsim_source_size_audit PROPERTIES
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
