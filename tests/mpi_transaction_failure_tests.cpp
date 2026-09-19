#include "core/problem_backend_access.hpp"
#include "core/spatial_backend.hpp"
#include "io/case_input.hpp"
#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include "solver/solve_workflows.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <petscsys.h>
#include <stdexcept>

namespace {
using namespace fuelsim;

void require(bool value, const char* message) {
    if (!value)
        throw std::runtime_error(message);
}

void sum(std::vector<double>& values) {
    require(MPIU_Allreduce(MPI_IN_PLACE,
                values.data(),
                static_cast<int>(values.size()),
                MPI_DOUBLE,
                MPI_SUM,
                PETSC_COMM_WORLD)
                == MPI_SUCCESS,
        "Commit reduction failed");
}

void commit(TransientProblem& problem, const SolveResult& result) {
    problem.commit_time_step(result.state, result.local_contribution_begin, result.local_contribution_end, sum);
}

std::vector<char> checkpoint(TransientProblem& mirror, const TransientProblem& problem, const std::string& path) {
    // The writer deliberately rejects active steps. Export only committed data
    // into a separate inactive problem to inspect it before rollback as well.
    BackendAccess::restore_committed_state(mirror, BackendAccess::committed_state(problem));
    write_transient_checkpoint(path, mirror, 0.02);
    std::ifstream file(path, std::ios::binary);
    require(static_cast<bool>(file), "Cannot read transaction checkpoint");
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

void check_mpi(const FuelSimCaseDefinition& input, const std::string& path, int rank) {
    const auto mesh = read_exodus_hex8(input.mesh_file);
    TransientProblem problem(input.spatial, mesh), mirror(input.spatial, mesh);
    problem.track_previous_committed_solution(true);
    PetscSolver solver;
    for (int step = 1; step <= 10; ++step) {
        problem.begin_time_step({step * 0.02, 1.0});
        auto result = solver.solve(problem, problem.committed_solution(), input.solver);
        require(result.converged, "Initial contact/material loading did not converge");
        commit(problem, result);
    }
    double inelastic = 0.0, friction = 0.0;
    const auto initial = BackendAccess::committed_state(problem);
    for (const auto& region : initial.cartesian_material_histories)
        for (const auto& element : region)
            for (const auto& point : element)
                inelastic += point.equivalent_creep_strain + point.equivalent_plastic_strain;
    for (const auto& contact : initial.contact_histories)
        for (const auto& point : contact)
            for (double component : point.cartesian_total_tangential_slip)
                friction += std::abs(component);
    require(inelastic > 0.0 && friction > 0.0, "Failure test requires nonzero inelastic and friction history");
    const auto baseline = problem.capture_state();
    const auto before = checkpoint(mirror, problem, path);
    problem.begin_time_step({0.22, 1.0});
    const auto valid = solver.solve(problem, problem.committed_solution(), input.solver);
    require(valid.converged, "Reference continuation did not converge");
    commit(problem, valid);
    const auto expected = checkpoint(mirror, problem, path);
    problem.restore_state(baseline);
    const CommitPreparationStage stages[] = {CommitPreparationStage::nodal_state,
        CommitPreparationStage::material,
        CommitPreparationStage::contact,
        CommitPreparationStage::exchange_buffer,
        CommitPreparationStage::final_diagnostics};
    for (std::size_t scenario = 0; scenario < 6; ++scenario) {
        problem.begin_time_step({0.22, 1.0});
        bool material_prepared = false, contact_prepared = false;
        if (scenario < 5)
            BackendAccess::set_commit_test_hook(problem, [&](CommitPreparationStage stage) {
                material_prepared = material_prepared || stage == CommitPreparationStage::material;
                contact_prepared = contact_prepared || stage == CommitPreparationStage::contact;
                if (rank == 1 && stage == stages[scenario]) {
                    if (scenario == 1)
                        throw std::domain_error("Injected staged material failure");
                    if (scenario == 2)
                        throw std::runtime_error("Injected staged contact failure");
                    throw std::bad_alloc();
                }
            });
        auto candidate = valid.state;
        if (scenario == 5 && rank == 1)
            candidate.pop_back();
        bool rejected = false;
        unsigned collectives = 0;
        try {
            problem.commit_time_step(candidate,
                valid.local_contribution_begin,
                valid.local_contribution_end,
                [&](std::vector<double>& values) {
                    ++collectives;
                    sum(values);
                });
        } catch (const std::exception&) {
            rejected = true;
        }
        BackendAccess::set_commit_test_hook(problem, {});
        int failed_ranks = rejected ? 1 : 0;
        MPIU_Allreduce(MPI_IN_PLACE, &failed_ranks, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        require(failed_ranks == 2, "Both ranks must reject the failed commit");
        require(collectives == (scenario == 4 ? 3U : 1U), "Failure entered an unexpected collective sequence");
        if (scenario == 1)
            require(material_prepared, "Material failure hook ran before material preparation");
        if (scenario == 2)
            require(material_prepared && contact_prepared, "Contact failure hook ran before real staging");
        require(problem.time_step_active(), "Failed commit must allow rollback");
        require(before == checkpoint(mirror, problem, path), "Failure changed committed checkpoint before rollback");
        problem.rollback_time_step();
        require(before == checkpoint(mirror, problem, path), "Rollback changed committed checkpoint");
        problem.begin_time_step({0.22, 1.0});
        auto retry = solver.solve(problem, problem.committed_solution(), input.solver);
        require(retry.converged, "Solve after failed transaction did not converge");
        commit(problem, retry);
        require(expected == checkpoint(mirror, problem, path), "Retry differs from failure-free continuation");
        problem.restore_state(baseline);
        if (rank == 0)
            std::cout << "scenario=" << scenario << " rejected_ranks=" << failed_ranks
                      << " checkpoint_equal=1 retry_equal=1 collectives=" << collectives << '\n';
    }
    if (rank == 0)
        std::cout << "inelastic_history_sum=" << inelastic << " friction_history_sum=" << friction << '\n';
}
} // namespace

int main(int argc, char** argv) {
    fuelsim::PetscSession session(argc, argv, "Distributed transaction failure contracts");
    try {
        require(argc == 3 && session.size() == 2, "Expected input card and scratch path, with exactly two MPI ranks");
        check_mpi(fuelsim::read_case_input(argv[1]),
            std::string(argv[2]) + std::to_string(session.rank()),
            session.rank());
    } catch (const std::exception& error) {
        std::cerr << "rank " << session.rank() << ": " << error.what() << '\n';
        (void)PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,
            __LINE__,
            __func__,
            __FILE__,
            PETSC_ERR_LIB,
            PETSC_ERROR_INITIAL,
            error.what(),
            nullptr);
        return 1;
    }
}
