cmake_minimum_required(VERSION 3.20)
project(fuelsim_elements VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(adlite 0.2.3 CONFIG REQUIRED)
include(CTest)
option(FUELSIM_ELEMENTS_WARNINGS_AS_ERRORS "Treat element library warnings as errors" ON)

add_library(fuelsim_elements
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
if(COMMAND fuelsim_enable_warnings)
    fuelsim_enable_warnings(fuelsim_elements)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fuelsim_elements PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion)
    if(FUELSIM_ELEMENTS_WARNINGS_AS_ERRORS)
        target_compile_options(fuelsim_elements PRIVATE -Werror)
    endif()
endif()

if(BUILD_TESTING)
    add_executable(fuelsim_cax4t_element_tests tests/cax4t_tests.cpp)
    target_link_libraries(fuelsim_cax4t_element_tests PRIVATE fuelsim::elements)
    if(COMMAND fuelsim_enable_warnings)
        fuelsim_enable_warnings(fuelsim_cax4t_element_tests)
    endif()
    # The parent registers this same test in its unified verification matrix.
    if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
        add_test(NAME fuelsim_cax4t_element_tests COMMAND fuelsim_cax4t_element_tests)
    endif()
endif()
