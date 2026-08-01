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

find_package(adlite 0.1.0 CONFIG REQUIRED)
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
    elseif(MSVC)
        target_compile_options("${target}" PRIVATE /W4 /permissive-)
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
    src/dof_map.cpp
    src/inelastic_material.cpp
    src/interface.cpp
    src/steady_problem.cpp
    src/transient_problem.cpp
    src/material.cpp
    src/mesh.cpp
    src/nonlinear_problem.cpp
    src/quad4_rz.cpp
    src/quad4_rz_transient.cpp
)
add_library(fuelsim::core ALIAS fuelsim_core)

target_include_directories(fuelsim_core
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(fuelsim_core PUBLIC adlite::adlite)
fuelsim_enable_warnings(fuelsim_core)

add_library(fuelsim_exodus
    src/exodus_mesh_io.cpp
)
add_library(fuelsim::exodus ALIAS fuelsim_exodus)

target_include_directories(fuelsim_exodus
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(fuelsim_exodus
    PUBLIC fuelsim::core
    PRIVATE SEACASExodus::exodus
)
fuelsim_enable_warnings(fuelsim_exodus)

add_library(fuelsim_input
    src/case_input.cpp
    src/input_file.cpp
)
add_library(fuelsim::input ALIAS fuelsim_input)

target_include_directories(fuelsim_input
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(fuelsim_input PUBLIC fuelsim::core)
fuelsim_enable_warnings(fuelsim_input)

add_library(fuelsim_petsc
    src/problem_solver.cpp
    src/petsc_solver.cpp
)
add_library(fuelsim::petsc ALIAS fuelsim_petsc)

target_include_directories(fuelsim_petsc
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(fuelsim_petsc
    PUBLIC fuelsim::core
)
fuelsim_link_petsc(fuelsim_petsc)
fuelsim_enable_warnings(fuelsim_petsc)

add_executable(fuelsim src/main.cpp)
target_link_libraries(fuelsim
    PRIVATE fuelsim::input fuelsim::exodus fuelsim::petsc
)
fuelsim_enable_warnings(fuelsim)

include(CTest)

if(BUILD_TESTING)
    add_library(fuelsim_test_support
        tests/support/steady_single_region_problem.cpp
        tests/support/steady_fuel_cladding_problem.cpp
        tests/support/steady_fuel_cladding_solver.cpp
        tests/support/transient_fuel_cladding_problem.cpp
        tests/support/transient_fuel_cladding_solver.cpp
    )
    target_include_directories(fuelsim_test_support
        PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    )
    target_link_libraries(fuelsim_test_support
        PUBLIC fuelsim::core fuelsim::petsc
    )
    fuelsim_enable_warnings(fuelsim_test_support)

    add_executable(fuelsim_m1_single_core_benchmark
        benchmarks/m1_single_core.cpp
    )
    target_link_libraries(fuelsim_m1_single_core_benchmark
        PRIVATE fuelsim_test_support
    )
    fuelsim_enable_warnings(fuelsim_m1_single_core_benchmark)

    add_executable(fuelsim_input_tests tests/input_tests.cpp)
    target_link_libraries(fuelsim_input_tests PRIVATE fuelsim::input)
    fuelsim_enable_warnings(fuelsim_input_tests)
    add_test(NAME fuelsim_input_tests
        COMMAND fuelsim_input_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/steady_fuel_cladding.fsi"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/fuelsim/transient_fuel_cladding_pcmi.fsi"
            "${CMAKE_CURRENT_BINARY_DIR}/malformed_input.fsi"
    )
    set_tests_properties(fuelsim_input_tests PROPERTIES LABELS "unit")

    add_executable(fuelsim_exodus_mesh_io_tests
        tests/exodus_mesh_io_tests.cpp
    )
    target_link_libraries(fuelsim_exodus_mesh_io_tests
        PRIVATE fuelsim::exodus
    )
    fuelsim_enable_warnings(fuelsim_exodus_mesh_io_tests)
    add_test(NAME fuelsim_exodus_mesh_io_tests
        COMMAND fuelsim_exodus_mesh_io_tests
            "${CMAKE_CURRENT_BINARY_DIR}/fuelsim_two_quad_smoke.exo"
    )
    set_tests_properties(fuelsim_exodus_mesh_io_tests PROPERTIES
        LABELS "dependency"
    )

    add_executable(fuelsim_m1_exodus_moose_tests
        tests/m1_exodus_moose_tests.cpp
    )
    target_link_libraries(fuelsim_m1_exodus_moose_tests
        PRIVATE fuelsim_test_support fuelsim::exodus
    )
    fuelsim_enable_warnings(fuelsim_m1_exodus_moose_tests)
    add_test(NAME fuelsim_m1_exodus_moose_tests
        COMMAND fuelsim_m1_exodus_moose_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m1_fuel_cladding_gap_rz_mesh.e"
    )
    set_tests_properties(fuelsim_m1_exodus_moose_tests PROPERTIES
        LABELS "acceptance"
    )

    add_executable(fuelsim_core_tests tests/core_tests.cpp)
    target_link_libraries(fuelsim_core_tests PRIVATE fuelsim_test_support)
    fuelsim_enable_warnings(fuelsim_core_tests)
    add_test(NAME fuelsim_core_tests COMMAND fuelsim_core_tests)
    set_tests_properties(fuelsim_core_tests PROPERTIES LABELS "unit")

    add_executable(fuelsim_multi_region_tests tests/multi_region_tests.cpp)
    target_link_libraries(fuelsim_multi_region_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_multi_region_tests)
    add_test(NAME fuelsim_multi_region_tests COMMAND fuelsim_multi_region_tests)
    set_tests_properties(fuelsim_multi_region_tests PROPERTIES LABELS "unit")

    add_executable(fuelsim_m2_core_tests tests/m2_core_tests.cpp)
    target_link_libraries(fuelsim_m2_core_tests PRIVATE fuelsim_test_support)
    fuelsim_enable_warnings(fuelsim_m2_core_tests)
    add_test(NAME fuelsim_m2_core_tests COMMAND fuelsim_m2_core_tests)
    set_tests_properties(fuelsim_m2_core_tests PROPERTIES LABELS "unit")

    add_executable(fuelsim_solver_tests tests/solver_tests.cpp)
    target_link_libraries(fuelsim_solver_tests
        PRIVATE fuelsim_test_support fuelsim::exodus
    )
    fuelsim_enable_warnings(fuelsim_solver_tests)
    add_test(NAME fuelsim_solver_tests
        COMMAND fuelsim_solver_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m0_simple_fuel_rz_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m1_fuel_cladding_gap_rz_mesh.e"
    )
    set_tests_properties(fuelsim_solver_tests PROPERTIES LABELS "acceptance")

    add_executable(fuelsim_m2_solver_tests tests/m2_solver_tests.cpp)
    target_link_libraries(fuelsim_m2_solver_tests
        PRIVATE fuelsim_test_support fuelsim::exodus
    )
    fuelsim_enable_warnings(fuelsim_m2_solver_tests)
    add_test(NAME fuelsim_m2_solver_tests
        COMMAND fuelsim_m2_solver_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m21_transient_heat_rz_mesh.e"
    )
    set_tests_properties(fuelsim_m2_solver_tests PROPERTIES LABELS "acceptance")

    add_executable(fuelsim_m2_inelastic_solver_tests
        tests/m2_inelastic_solver_tests.cpp
    )
    target_link_libraries(fuelsim_m2_inelastic_solver_tests
        PRIVATE fuelsim::exodus fuelsim::petsc
    )
    fuelsim_enable_warnings(fuelsim_m2_inelastic_solver_tests)
    add_test(NAME fuelsim_m2_inelastic_solver_tests
        COMMAND fuelsim_m2_inelastic_solver_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_j2_plastic_rz_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_norton_creep_rz_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_coupled_plastic_creep_rz_mesh.e"
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m22_coupled_plastic_creep_traction_rz_mesh.e"
    )
    set_tests_properties(fuelsim_m2_inelastic_solver_tests
        PROPERTIES LABELS "acceptance"
    )

    add_executable(fuelsim_m2_pcmi_solver_tests
        tests/m2_pcmi_solver_tests.cpp
    )
    target_link_libraries(fuelsim_m2_pcmi_solver_tests
        PRIVATE fuelsim_test_support fuelsim::exodus
    )
    fuelsim_enable_warnings(fuelsim_m2_pcmi_solver_tests)
    add_test(NAME fuelsim_m2_pcmi_solver_tests
        COMMAND fuelsim_m2_pcmi_solver_tests
            "${CMAKE_CURRENT_SOURCE_DIR}/verification/moose/m23_pcmi_coupled_cladding_rz_mesh.e"
    )
    set_tests_properties(fuelsim_m2_pcmi_solver_tests
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
endif()
