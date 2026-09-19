#include "core/problem_backend_access.hpp"
#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
void append_point(std::vector<double>& values, const fuelsim::MaterialPointState& point) {
    for (const auto* tensor : {&point.elastic_strain, &point.plastic_strain, &point.creep_strain})
        values.insert(values.end(), tensor->begin(), tensor->end());
    values.insert(values.end(),
        {point.equivalent_plastic_strain,
            point.equivalent_creep_strain,
            point.stress.rr,
            point.stress.zz,
            point.stress.hoop,
            point.stress.rz});
}

std::vector<double> committed_values(const fuelsim::TransientProblem& problem) {
    const auto state = fuelsim::BackendAccess::committed_state(problem);
    std::vector<double> values{state.time, state.load_factor, state.previous_time};
    for (const auto* field :
        {&state.solution, &state.previous_solution, &state.raw_residual, &state.external_load_residual})
        values.insert(values.end(), field->begin(), field->end());
    for (const auto& field : fuelsim::transient_conservation_fields)
        values.push_back(state.conservation.*field.member);
    for (const auto& region : state.material_histories)
        for (const auto& history : region)
            for (const auto& point : history)
                append_point(values, point);
    for (const auto& region : state.quad8_material_histories)
        for (const auto& history : region)
            for (const auto& point : history)
                append_point(values, point);
    for (const auto& interface : state.contact_histories)
        for (const auto& point : interface)
            values.insert(values.end(),
                {point.elastic_tangential_slip,
                    point.sliding ? 1.0 : 0.0,
                    point.normal_multiplier,
                    point.total_tangential_slip});
    return values;
}

void require_close(const std::vector<double>& actual, const std::vector<double>& expected, double tolerance) {
    if (actual.size() != expected.size())
        throw std::runtime_error("Partitioned commit changed state dimensions");
    for (std::size_t i = 0; i < actual.size(); ++i)
        if (!std::isfinite(actual[i]) || std::abs(actual[i] - expected[i]) > tolerance * (1.0 + std::abs(expected[i])))
            throw std::runtime_error("Partitioned commit changed state entry " + std::to_string(i));
}

void check(const std::string& path) {
    const auto input = fuelsim::read_case_input(path);
    const bool quadratic =
        input.spatial.regions.front().rz_element_formulation == fuelsim::RzElementFormulation::cax8t
        || input.spatial.regions.front().rz_element_formulation == fuelsim::RzElementFormulation::cax8rt;
    std::unique_ptr<fuelsim::TransientProblem> instance;
    if (quadratic)
        instance =
            std::make_unique<fuelsim::TransientProblem>(input.spatial, fuelsim::read_exodus_quad8(input.mesh_file));
    else
        instance =
            std::make_unique<fuelsim::TransientProblem>(input.spatial, fuelsim::read_exodus_quad4(input.mesh_file));
    auto& problem = *instance;
    problem.track_previous_committed_solution(true);
    for (int step = 1; step <= 2; ++step) {
        const auto before = problem.capture_state();
        const auto before_values = committed_values(problem);
        auto candidate = problem.committed_solution();
        for (const auto& field : problem.field_layout())
            for (std::size_t dof = field.begin; dof < field.end; ++dof)
                candidate[dof] +=
                    field.category == fuelsim::FieldCategory::thermal ? 1.0 : 1e-8 * static_cast<double>(1 + dof % 7);
        const fuelsim::TransientStepInput increment{step * 0.0625, 1.0};
        problem.begin_time_step(increment);
        problem.commit_time_step(candidate);
        const auto expected = committed_values(problem);
        problem.restore_state(before);

        // Capture two disjoint local buffers without committing either one.
        std::vector<double> sum;
        for (std::size_t rank = 0; rank < 2; ++rank) {
            problem.begin_time_step(increment);
            const auto interval = problem.contribution_partition(rank, 2);
            bool rejected = false;
            try {
                problem.commit_time_step(candidate,
                    interval.first,
                    interval.second,
                    [&sum](std::vector<double>& values) {
                        if (values.size() == 2)
                            return;
                        if (sum.empty())
                            sum.assign(values.size(), 0.0);
                        for (std::size_t i = 0; i < values.size(); ++i)
                            sum[i] += values[i];
                        values[0] = 1.0;
                    });
            } catch (const std::domain_error&) {
                rejected = true;
            }
            if (!rejected)
                throw std::runtime_error("Remote partition failure was not rejected");
            require_close(committed_values(problem), before_values, 0.0);
            problem.rollback_time_step();
            require_close(committed_values(problem), before_values, 0.0);
        }
        if (sum[0] != 0.0 || sum[1] != 0.0)
            throw std::runtime_error("Local material update unexpectedly failed");
        // An empty owner must still receive every history and diagnostic.
        problem.begin_time_step(increment);
        problem.commit_time_step(candidate, 0, 0, [&sum](std::vector<double>& values) {
            if (values.size() != 2)
                values = sum;
        });
        require_close(committed_values(problem), expected, 2e-12);
    }
    std::cout << "[PASS] Partition ownership, remote failure, empty owner, history and diagnostics: " << path << '\n';
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5)
            throw std::runtime_error("Expected four axisymmetric input cards");
        for (int i = 1; i < argc; ++i)
            check(argv[i]);
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
