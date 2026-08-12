if(NOT DEFINED ROOT)
    message(FATAL_ERROR "geometry-independence audit requires ROOT")
endif()

set(common_files
    "${ROOT}/include/fuelsim/nonlinear_problem.hpp"
    "${ROOT}/include/fuelsim/diagnostics.hpp"
    "${ROOT}/include/fuelsim/petsc_solver.hpp"
    "${ROOT}/include/fuelsim/problem_solver.hpp"
    "${ROOT}/src/solver.cpp"
)

foreach(path IN LISTS common_files)
    file(READ "${path}" contents)
    foreach(forbidden
            "local_dof_count"
            "LocalDofs"
            "LocalValues"
            "LocalResidual"
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
    "${ROOT}/include/fuelsim/steady_problem.hpp"
    "${ROOT}/include/fuelsim/transient_problem.hpp"
)

foreach(path IN LISTS common_problem_headers)
    file(READ "${path}" contents)
    foreach(forbidden
            "DofMap"
            "RegionMesh"
            "Quad4Rz"
            "LocalDofs"
            "LocalValues"
            "LocalResidual"
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

file(READ "${ROOT}/src/solver.cpp" common_solver)
foreach(forbidden_state_layout
        "TransientCommittedState"
        "Quad4MaterialHistory"
        "AxisymmetricStressValues"
        "material_histories"
        "material_stresses"
        "ContactPointHistory")
    string(FIND "${common_solver}" "${forbidden_state_layout}" location)
    if(NOT location EQUAL -1)
        message(FATAL_ERROR
            "The common steady/time solver reintroduced backend state layout: ${forbidden_state_layout}")
    endif()
endforeach()

file(READ "${ROOT}/include/fuelsim/nonlinear_problem.hpp" problem_port)
foreach(required "ContributionWorkspace" "FieldDescriptor" "field_layout" "contribution_dof_count"
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

file(READ "${ROOT}/src/assembly.hpp" rz_assembly)
string(FIND "${rz_assembly}" "namespace rz" rz_namespace)
if(rz_namespace EQUAL -1)
    message(FATAL_ERROR "The RZ spatial assembly is no longer isolated in namespace fuelsim::rz")
endif()

file(READ "${ROOT}/include/fuelsim/rz_problem_access.hpp" rz_access)
foreach(required "class ProblemAccess" "struct TransientCommittedState" "Quad4RzGeometry"
                 "ContactPointHistory" "contribution_state")
    string(FIND "${rz_access}" "${required}" location)
    if(location EQUAL -1)
        message(FATAL_ERROR "The explicit RZ problem access layer is missing: ${required}")
    endif()
endforeach()

message(STATUS "Geometry-independent solver and problem-port source audit passed")
