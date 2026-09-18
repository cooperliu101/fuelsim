cmake_minimum_required(VERSION 3.20)
project(fuelsim_elements VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(adlite 0.2.3 CONFIG REQUIRED)
include(CTest)
option(FUELSIM_ELEMENTS_WARNINGS_AS_ERRORS "Treat element library warnings as errors" ON)

add_library(fuelsim_elements
    src/cpeg8t.cpp
    src/cax2t_gps.cpp
    src/c3d20rt.cpp
    src/c3d20t.cpp
    src/c3d8rt.cpp
    src/c3d8t.cpp
    src/c3d_common.cpp
    src/cax4rt.cpp
    src/cax4t.cpp
    src/cax8rt.cpp
    src/cax8t.cpp
    src/contact_common.cpp
    src/cax_common.cpp
    src/line2_rz.cpp
    src/line3_rz.cpp
    src/line3_plane.cpp
    src/material.cpp
    src/material_functions.cpp
    src/quad4_face.cpp
    src/quad8_face.cpp
    src/ring_gps.cpp
)
add_library(fuelsim::elements ALIAS fuelsim_elements)
target_compile_features(fuelsim_elements PUBLIC cxx_std_17)
set_target_properties(fuelsim_elements PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(fuelsim_elements PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(fuelsim_elements PUBLIC adlite::adlite)
function(fuelsim_element_warnings target)
    if(COMMAND fuelsim_enable_warnings)
        fuelsim_enable_warnings(${target})
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion)
        if(FUELSIM_ELEMENTS_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()
fuelsim_element_warnings(fuelsim_elements)

if(BUILD_TESTING)
    add_test(NAME fuelsim_element_library_boundaries
        COMMAND ${CMAKE_COMMAND}
            -DELEMENTS_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DFUELSIM_SOURCE_DIR=${CMAKE_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/library_boundaries.cmake)
    function(fuelsim_add_element_test name source)
        add_executable(${name} tests/${source})
        target_link_libraries(${name} PRIVATE fuelsim_elements_test_support)
        target_include_directories(${name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
        fuelsim_element_warnings(${name})
        # The parent keeps the same test names and its unified verification matrix.
        if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
            add_test(NAME ${name} COMMAND ${name})
        endif()
    endfunction()
    add_executable(fuelsim_common_math_tests tests/common_math_tests.cpp)
    target_link_libraries(fuelsim_common_math_tests PRIVATE fuelsim::elements)
    target_include_directories(fuelsim_common_math_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    fuelsim_element_warnings(fuelsim_common_math_tests)
    add_test(NAME fuelsim_common_math_tests COMMAND fuelsim_common_math_tests)
    fuelsim_add_element_test(fuelsim_creep_rate_tests creep_rate_tests.cpp)
    add_executable(fuelsim_body_acceleration_tests tests/body_acceleration_tests.cpp)
    target_link_libraries(fuelsim_body_acceleration_tests PRIVATE fuelsim::elements)
    target_include_directories(fuelsim_body_acceleration_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_element_warnings(fuelsim_body_acceleration_tests)
    if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
        add_test(NAME fuelsim_body_acceleration_tests COMMAND fuelsim_body_acceleration_tests)
    endif()
    # This internal contract test deliberately exercises private failure branches.
    add_executable(fuelsim_contact_common_tests tests/contact_common_tests.cpp)
    target_link_libraries(fuelsim_contact_common_tests PRIVATE fuelsim::elements)
    target_include_directories(fuelsim_contact_common_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    fuelsim_element_warnings(fuelsim_contact_common_tests)
    add_test(NAME fuelsim_contact_common_tests COMMAND fuelsim_contact_common_tests)
    add_library(fuelsim_elements_test_support STATIC tests/support/cax8_common_tests.cpp)
    target_link_libraries(fuelsim_elements_test_support PUBLIC fuelsim::elements)
    target_include_directories(fuelsim_elements_test_support PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests")
    fuelsim_element_warnings(fuelsim_elements_test_support)
    fuelsim_add_element_test(fuelsim_cax4t_tests cax4t_tests.cpp)
    fuelsim_add_element_test(fuelsim_cpeg8t_tests cpeg8t_tests.cpp)
    fuelsim_add_element_test(fuelsim_cax2t_gps_tests cax2t_gps_tests.cpp)
    fuelsim_add_element_test(fuelsim_ring_gps_tests ring_gps_tests.cpp)
    fuelsim_add_element_test(fuelsim_cax4rt_tests cax4rt_tests.cpp)
    fuelsim_add_element_test(fuelsim_cax8t_tests cax8t_tests.cpp)
    fuelsim_add_element_test(fuelsim_cax8rt_tests cax8rt_tests.cpp)
    fuelsim_add_element_test(fuelsim_c3d8t_tests c3d8t_tests.cpp)
    fuelsim_add_element_test(fuelsim_c3d8rt_tests c3d8rt_tests.cpp)
    fuelsim_add_element_test(fuelsim_c3d20t_tests c3d20t_tests.cpp)
    fuelsim_add_element_test(fuelsim_c3d20rt_tests c3d20rt_tests.cpp)
    fuelsim_add_element_test(fuelsim_line2_rz_tests line2_rz_tests.cpp)
    fuelsim_add_element_test(fuelsim_line3_rz_tests line3_rz_tests.cpp)
    fuelsim_add_element_test(fuelsim_quad4_face_tests quad4_face_tests.cpp)
    fuelsim_add_element_test(fuelsim_quad8_face_tests quad8_face_tests.cpp)
endif()
