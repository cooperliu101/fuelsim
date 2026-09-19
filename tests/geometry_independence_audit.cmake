if(NOT DEFINED ROOT)
    message(FATAL_ERROR "geometry-independence audit requires ROOT")
endif()

set(common_files
    "${ROOT}/include/core/nonlinear_problem.hpp"
    "${ROOT}/include/solver/petsc_solver.hpp"
    "${ROOT}/include/solver/solve_workflows.hpp"
    "${ROOT}/src/solver/solver.cpp"
    "${ROOT}/src/solver/solve_workflows.cpp"
)

foreach(path IN LISTS common_files)
    file(READ "${path}" contents)
    foreach(forbidden
            "cax4_local_dof_count"
            "Cax4LocalDofs"
            "Cax4LocalValues"
            "Cax4LocalResidual"
            "LocalSystem"
            "std::array<double, 3>"
            "dof_count() / 3"
            "dof_count() / 3U"
            "global_count / 3"
            "global_count / 3U"
            "[T, ur, uz]"
            "rz_problem_access.hpp"
            "rz::ProblemAccess")
        string(FIND "${contents}" "${forbidden}" location)
        if(NOT location EQUAL -1)
            message(FATAL_ERROR "${path} reintroduced geometry-specific solver text: ${forbidden}")
        endif()
    endforeach()
endforeach()

set(common_problem_headers
    "${ROOT}/include/core/steady_problem.hpp"
    "${ROOT}/include/core/transient_problem.hpp"
)

foreach(path IN LISTS common_problem_headers)
    file(READ "${path}" contents)
    foreach(forbidden
            "DofMap"
            "RegionMesh"
            "Quad4Rz"
            "Cax4LocalDofs"
            "Cax4LocalValues"
            "Cax4LocalResidual"
            "LocalSystem"
            "SpatialAssembly"
            "ContactPointHistory"
            "AxisymmetricStressValues"
            "Quad4MaterialHistory"
            "TransientCommittedState")
        string(FIND "${contents}" "${forbidden}" location)
        if(NOT location EQUAL -1)
            message(FATAL_ERROR "${path} reintroduced an RZ query or state layout: ${forbidden}")
        endif()
    endforeach()
endforeach()

file(READ "${ROOT}/src/solver/solver.cpp" common_solver)
file(READ "${ROOT}/src/solver/solve_workflows.cpp" common_workflows)
string(APPEND common_solver "\n${common_workflows}")
foreach(forbidden_state_layout
        "TransientCommittedState"
        "Quad4MaterialHistory"
        "AxisymmetricStressValues"
        "material_histories"
        "ContactPointHistory")
    string(FIND "${common_solver}" "${forbidden_state_layout}" location)
    if(NOT location EQUAL -1)
        message(FATAL_ERROR
            "The common steady/time solver reintroduced backend state layout: ${forbidden_state_layout}")
    endif()
endforeach()

file(READ "${ROOT}/include/core/nonlinear_problem.hpp" problem_port)
foreach(required "ContributionWorkspace" "FieldDescriptor" "field_layout" "contribution_dofs"
                 "discretization_identity")
    string(FIND "${problem_port}" "${required}" location)
    if(location EQUAL -1)
        message(FATAL_ERROR "NonlinearProblem runtime port is missing: ${required}")
    endif()
endforeach()

file(READ "${ROOT}/tests/solver_tests.cpp" runtime_test)
foreach(required "RuntimeLayoutProblem" "return 32" "_narrow_dofs" "displacement_z"
                 "Preconditioner::field_split")
    string(FIND "${runtime_test}" "${required}" location)
    if(location EQUAL -1)
        message(FATAL_ERROR "Runtime-layout behavior test is missing: ${required}")
    endif()
endforeach()

file(READ "${ROOT}/src/core/rz_assembly.hpp" rz_assembly)
string(FIND "${rz_assembly}" "namespace fuelsim::rz" rz_namespace)
if(rz_namespace EQUAL -1)
    message(FATAL_ERROR "The RZ spatial assembly is no longer isolated in namespace fuelsim::rz")
endif()

file(READ "${ROOT}/src/core/problem_backend_access.hpp" backend_access)
foreach(required "namespace rz" "struct TransientCommittedState" "Quad4MaterialHistory"
                 "class BackendAccess" "static const cartesian::SpatialAssembly& cartesian_spatial")
    string(FIND "${backend_access}" "${required}" location)
    if(location EQUAL -1)
        message(FATAL_ERROR "The private concrete backend access layer is missing: ${required}")
    endif()
endforeach()

foreach(path "${ROOT}/tests/support/rz_problem_access.hpp" "${ROOT}/tests/support/cartesian3d_problem_access.hpp")
    file(READ "${path}" test_access)
    string(FIND "${test_access}" "class ProblemAccess" location)
    if(location EQUAL -1)
        message(FATAL_ERROR "The test-only concrete problem access layer is missing: ${path}")
    endif()
endforeach()
string(FIND "${backend_access}" "Hex8TransientKernel" obsolete_kernel)
if(NOT obsolete_kernel EQUAL -1)
    message(FATAL_ERROR "The obsolete separate HEX8 transient kernel was reintroduced")
endif()

file(READ "${ROOT}/src/core/spatial_backend.hpp" spatial_backend)
foreach(required "class SpatialBackend" "required_state_dofs" "contribution_metadata_is_fixed"
                 "sparsity_contribution_jacobian_pattern" "prepare_time_step" "publish_material_history")
    string(FIND "${spatial_backend}" "${required}" location)
    if(location EQUAL -1)
        message(FATAL_ERROR "The private spatial backend protocol is missing: ${required}")
    endif()
endforeach()
foreach(forbidden "Quad4MaterialHistory" "Quad8MaterialHistory" "CartesianMaterialHistory"
                  "Cax2tGpsMaterialHistory" "thermal_assembly.hpp" "cartesian3d_assembly.hpp"
                  "petsc" "template<" "template <")
    string(FIND "${spatial_backend}" "${forbidden}" location)
    if(NOT location EQUAL -1)
        message(FATAL_ERROR "The backend interface acquired a concrete runtime dependency: ${forbidden}")
    endif()
endforeach()

file(READ "${ROOT}/src/core/spatial_problem.cpp" spatial_problem)
string(FIND "${spatial_problem}" "std::unique_ptr<SpatialBackend> _backend" location)
if(location EQUAL -1)
    message(FATAL_ERROR "SpatialProblemStorage must own a single private backend")
endif()
foreach(forbidden "std::unique_ptr<thermal::SpatialAssembly>" "std::unique_ptr<rz::SpatialAssembly>"
                  "std::unique_ptr<cartesian::SpatialAssembly>" "evaluate_cax4(" "compute_cax8("
                  "elements::evaluate_cpeg8t(" "transient_update(" "body_point_work(")
    string(FIND "${spatial_problem}" "${forbidden}" location)
    if(NOT location EQUAL -1)
        message(FATAL_ERROR "The public problem lifecycle reacquired backend implementation: ${forbidden}")
    endif()
endforeach()

message(STATUS "Geometry-independent solver, problem port and spatial backend source audit passed")
