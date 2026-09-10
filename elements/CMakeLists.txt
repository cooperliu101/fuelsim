cmake_minimum_required(VERSION 3.20)
project(fuelsim_elements VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(adlite 0.2.3 CONFIG REQUIRED)
include(CTest)
option(FUELSIM_ELEMENTS_WARNINGS_AS_ERRORS "Treat element library warnings as errors" ON)

add_library(fuelsim_elements
    src/rz_kernels.cpp
    src/rz_quad8.cpp
    src/cax4rt.cpp
    src/cartesian3d_kernels.cpp
    src/cartesian3d_hex20.cpp
    src/cartesian3d_mechanics.cpp
    src/cartesian3d_contact.cpp
    src/cartesian3d_hex20_contact.cpp
    src/line3_rz_contact.cpp
    src/line3_rz_boundary.cpp
    src/cax4t.cpp
    src/rz_geometry.cpp
    src/material.cpp
    src/material_functions.cpp
)
add_library(fuelsim::elements ALIAS fuelsim_elements)
target_compile_features(fuelsim_elements PUBLIC cxx_std_17)
set_target_properties(fuelsim_elements PROPERTIES CXX_EXTENSIONS OFF)
target_include_directories(fuelsim_elements PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
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
    function(fuelsim_add_element_test name source)
        add_executable(${name} tests/${source})
        target_link_libraries(${name} PRIVATE fuelsim::elements)
        fuelsim_element_warnings(${name})
        # The parent keeps the same test names and its unified verification matrix.
        if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
            add_test(NAME ${name} COMMAND ${name})
        endif()
    endfunction()
    fuelsim_add_element_test(fuelsim_cax4t_element_tests cax4t_tests.cpp)
    fuelsim_add_element_test(fuelsim_rz_element_tests rz_quad4_tests.cpp)
    fuelsim_add_element_test(fuelsim_quad8_rz_tests quad8_rz_tests.cpp)
    fuelsim_add_element_test(fuelsim_hex8_tests hex8_tests.cpp)
    fuelsim_add_element_test(fuelsim_hex20_tests hex20_tests.cpp)
endif()
