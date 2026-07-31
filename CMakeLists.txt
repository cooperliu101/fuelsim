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
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc petsc)

option(FUELSIM_REQUIRE_PETSC_EXODUS
    "Require the selected PETSc installation to provide ExodusII support"
    OFF
)

if(FUELSIM_REQUIRE_PETSC_EXODUS)
    include(CMakePushCheckState)
    include(CheckCXXSourceCompiles)
    cmake_push_check_state(RESET)
    set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES PkgConfig::PETSC)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND PETSC_STATIC_LIBRARY_DIRS)
        list(JOIN PETSC_STATIC_LIBRARY_DIRS ":" petsc_check_link_dirs)
        set(CMAKE_REQUIRED_LINK_OPTIONS
            "LINKER:-rpath-link,${petsc_check_link_dirs}"
        )
    endif()
    check_cxx_source_compiles(
        "#include <petscsys.h>
         #ifndef PETSC_HAVE_EXODUSII
         #error PETSc was built without ExodusII
         #endif
         int main() {
             return PetscGlobalSize < 0;
         }"
        FUELSIM_PETSC_HAS_EXODUS
    )
    cmake_pop_check_state()

    if(NOT FUELSIM_PETSC_HAS_EXODUS)
        message(FATAL_ERROR
            "FUELSIM_REQUIRE_PETSC_EXODUS=ON, but the selected PETSc was "
            "built without ExodusII support or the C++ compiler cannot link "
            "the selected library. Set PKG_CONFIG_PATH to an Exodus-enabled "
            "PETSc prefix."
        )
    endif()
endif()

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
    src/m1_problem.cpp
    src/m2_problem.cpp
    src/material.cpp
    src/mesh.cpp
    src/nonlinear_problem.cpp
    src/problem.cpp
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

add_library(fuelsim_petsc
    src/m1_solver.cpp
    src/m2_solver.cpp
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
target_link_libraries(fuelsim PRIVATE fuelsim::petsc)
fuelsim_enable_warnings(fuelsim)

add_executable(fuelsim_m2 src/m2_main.cpp)
target_link_libraries(fuelsim_m2 PRIVATE fuelsim::petsc)
fuelsim_enable_warnings(fuelsim_m2)

add_executable(fuelsim_m1_single_core_benchmark
    benchmarks/m1_single_core.cpp
)
target_link_libraries(fuelsim_m1_single_core_benchmark PRIVATE fuelsim::petsc)
fuelsim_enable_warnings(fuelsim_m1_single_core_benchmark)

include(CTest)

if(BUILD_TESTING)
    if(FUELSIM_REQUIRE_PETSC_EXODUS)
        add_executable(fuelsim_petsc_exodus_smoke
            tests/petsc_exodus_smoke.cpp
        )
        fuelsim_link_petsc(fuelsim_petsc_exodus_smoke)
        fuelsim_enable_warnings(fuelsim_petsc_exodus_smoke)
        add_test(NAME fuelsim_petsc_exodus_smoke
            COMMAND fuelsim_petsc_exodus_smoke
                "${CMAKE_CURRENT_BINARY_DIR}/fuelsim_two_quad_smoke.exo"
        )
        set_tests_properties(fuelsim_petsc_exodus_smoke PROPERTIES
            LABELS "dependency"
        )
    endif()

    add_executable(fuelsim_core_tests tests/core_tests.cpp)
    target_link_libraries(fuelsim_core_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_core_tests)
    add_test(NAME fuelsim_core_tests COMMAND fuelsim_core_tests)
    set_tests_properties(fuelsim_core_tests PROPERTIES LABELS "unit")

    add_executable(fuelsim_m2_core_tests tests/m2_core_tests.cpp)
    target_link_libraries(fuelsim_m2_core_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_m2_core_tests)
    add_test(NAME fuelsim_m2_core_tests COMMAND fuelsim_m2_core_tests)
    set_tests_properties(fuelsim_m2_core_tests PROPERTIES LABELS "unit")

    add_executable(fuelsim_solver_tests tests/solver_tests.cpp)
    target_link_libraries(fuelsim_solver_tests PRIVATE fuelsim::petsc)
    fuelsim_enable_warnings(fuelsim_solver_tests)
    add_test(NAME fuelsim_solver_tests COMMAND fuelsim_solver_tests)
    set_tests_properties(fuelsim_solver_tests PROPERTIES LABELS "acceptance")

    add_executable(fuelsim_m2_solver_tests tests/m2_solver_tests.cpp)
    target_link_libraries(fuelsim_m2_solver_tests PRIVATE fuelsim::petsc)
    fuelsim_enable_warnings(fuelsim_m2_solver_tests)
    add_test(NAME fuelsim_m2_solver_tests COMMAND fuelsim_m2_solver_tests)
    set_tests_properties(fuelsim_m2_solver_tests PROPERTIES LABELS "acceptance")

    add_executable(fuelsim_m2_inelastic_solver_tests
        tests/m2_inelastic_solver_tests.cpp
    )
    target_link_libraries(fuelsim_m2_inelastic_solver_tests
        PRIVATE fuelsim::petsc
    )
    fuelsim_enable_warnings(fuelsim_m2_inelastic_solver_tests)
    add_test(NAME fuelsim_m2_inelastic_solver_tests
        COMMAND fuelsim_m2_inelastic_solver_tests
    )
    set_tests_properties(fuelsim_m2_inelastic_solver_tests
        PROPERTIES LABELS "acceptance"
    )

    add_executable(fuelsim_m2_pcmi_solver_tests
        tests/m2_pcmi_solver_tests.cpp
    )
    target_link_libraries(fuelsim_m2_pcmi_solver_tests
        PRIVATE fuelsim::petsc
    )
    fuelsim_enable_warnings(fuelsim_m2_pcmi_solver_tests)
    add_test(NAME fuelsim_m2_pcmi_solver_tests
        COMMAND fuelsim_m2_pcmi_solver_tests
    )
    set_tests_properties(fuelsim_m2_pcmi_solver_tests
        PROPERTIES LABELS "acceptance"
    )

    add_test(NAME fuelsim_m2_generic_demonstration COMMAND fuelsim_m2)
    set_tests_properties(fuelsim_m2_generic_demonstration
        PROPERTIES LABELS "acceptance"
    )
endif()
