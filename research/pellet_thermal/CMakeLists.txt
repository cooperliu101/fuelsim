target_sources(fuelsim_core PRIVATE src/pellet_thermal.cpp src/pellet_assembly.cpp src/pellet_mlp.cpp)
target_include_directories(fuelsim_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

if(BUILD_TESTING)
    add_executable(fuelsim_pellet_response_tests tests/response.cpp)
    target_link_libraries(fuelsim_pellet_response_tests PRIVATE fuelsim::io)
    target_include_directories(fuelsim_pellet_response_tests PRIVATE "${PROJECT_SOURCE_DIR}/src")
    fuelsim_enable_warnings(fuelsim_pellet_response_tests)
    add_executable(fuelsim_pellet_exact_tests tests/exact.cpp)
    target_link_libraries(fuelsim_pellet_exact_tests PRIVATE fuelsim::core)
    fuelsim_enable_warnings(fuelsim_pellet_exact_tests)
    add_test(NAME fuelsim_pellet_exact COMMAND fuelsim_pellet_exact_tests)
endif()

add_executable(fuelsim_pellet_sample src/sample.cpp)
target_link_libraries(fuelsim_pellet_sample PRIVATE fuelsim::solver fuelsim::io)
target_include_directories(fuelsim_pellet_sample PRIVATE "${PROJECT_SOURCE_DIR}/src")
fuelsim_enable_warnings(fuelsim_pellet_sample)

add_executable(fuelsim_pellet_infer src/infer.cpp)
target_link_libraries(fuelsim_pellet_infer PRIVATE fuelsim::core)
fuelsim_enable_warnings(fuelsim_pellet_infer)

if(BUILD_TESTING)
    add_test(NAME fuelsim_pellet_python_cpp COMMAND "${FUELSIM_SECTION_PYTHON}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/parity.py"
        --infer "$<TARGET_FILE:fuelsim_pellet_infer>"
        --model "${CMAKE_CURRENT_SOURCE_DIR}/models/pellet_mlp.txt"
        --work "${CMAKE_CURRENT_BINARY_DIR}/parity")
    foreach(study single coupled)
        add_test(NAME fuelsim_pellet_${study} COMMAND "${FUELSIM_SECTION_PYTHON}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/production.py"
            --fuelsim "$<TARGET_FILE:fuelsim>" --source "${CMAKE_CURRENT_SOURCE_DIR}"
            --response-check "$<TARGET_FILE:fuelsim_pellet_response_tests>"
            --work "${CMAKE_CURRENT_BINARY_DIR}/${study}" --study ${study})
    endforeach()
    set_tests_properties(fuelsim_pellet_exact fuelsim_pellet_python_cpp fuelsim_pellet_single fuelsim_pellet_coupled
        PROPERTIES LABELS "research;pellet_thermal" TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
endif()
