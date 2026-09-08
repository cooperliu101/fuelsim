#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace fuelsim;

struct Report {
    void value(const std::string& name, double value) { std::cout << name << "=" << value << "\n"; }
};

bool write_jacobian_check(const NonlinearProblem& problem,
    const std::vector<double>& state,
    Report& output,
    double step) {
    problem.validate_discretization();
    problem.validate_state(state);
    std::vector<double> direction(problem.dof_count()), analytic(problem.dof_count()), plus = state, minus = state;
    for (const FieldDescriptor& field : problem.field_layout())
        for (std::size_t dof = field.begin; dof < field.end; ++dof) {
            const double index = static_cast<double>((dof - field.begin) % 7);
            direction[dof] =
                field.category == FieldCategory::thermal ? 0.25 + 0.05 * index : 1.0e-6 * (0.4 + 0.1 * index);
            plus[dof] += step * direction[dof];
            minus[dof] -= step * direction[dof];
        }
    ContributionWorkspace workspace;
    for (std::size_t entry = 0; entry < problem.contribution_count(); ++entry) {
        problem.evaluate_contribution(entry, state, workspace, true);
        const std::size_t local_count = workspace.dofs.size();
        for (std::size_t row = 0; row < local_count; ++row)
            for (std::size_t column = 0; column < local_count; ++column)
                analytic[workspace.dofs[row]] +=
                    workspace.jacobian[row * local_count + column] * direction[workspace.dofs[column]];
    }
    const auto residual = [&](const std::vector<double>& values) {
        std::vector<double> result(problem.dof_count());
        ContributionWorkspace local_workspace;
        problem.validate_state(values);
        for (std::size_t entry = 0; entry < problem.contribution_count(); ++entry) {
            problem.evaluate_contribution(entry, values, local_workspace, false);
            for (std::size_t local = 0; local < local_workspace.dofs.size(); ++local)
                result[local_workspace.dofs[local]] += local_workspace.residual[local];
        }
        for (const DirichletCondition& condition : problem.dirichlet_conditions())
            result[condition.dof] = values[condition.dof] - condition.value;
        return result;
    };
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        analytic[condition.dof] = direction[condition.dof];
    const std::vector<double> current_residual = residual(state), plus_residual = residual(plus),
                              minus_residual = residual(minus);
    bool passed = true;
    for (const FieldDescriptor& field : problem.field_layout()) {
        double residual_l2 = 0.0, analytic_l2 = 0.0, reference = 0.0, difference = 0.0, reference_maximum = 0.0,
               difference_maximum = 0.0;
        for (std::size_t dof = field.begin; dof < field.end; ++dof) {
            const double finite_difference = (plus_residual[dof] - minus_residual[dof]) / (2.0 * step),
                         error = analytic[dof] - finite_difference;
            residual_l2 = std::hypot(residual_l2, current_residual[dof]);
            analytic_l2 = std::hypot(analytic_l2, analytic[dof]);
            reference = std::hypot(reference, finite_difference);
            difference = std::hypot(difference, error);
            reference_maximum = std::max(reference_maximum, std::abs(finite_difference));
            difference_maximum = std::max(difference_maximum, std::abs(error));
        }
        const std::string prefix = "jacobian." + field.name + ".";
        const double relative = reference > 0.0 ? difference / reference
                                                : (difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
        output.value(prefix + "residual_l2", residual_l2);
        output.value(prefix + "analytic_l2", analytic_l2);
        output.value(prefix + "finite_difference_l2", reference);
        output.value(prefix + "difference_l2", difference);
        output.value(prefix + "relative_l2", relative);
        output.value(prefix + "maximum_absolute_difference", difference_maximum);
        passed = difference <= 1.0e-6 * (1.0 + reference) && difference_maximum <= 1.0e-6 * (1.0 + reference_maximum)
                 && passed;
    }
    output.value("jacobian.check_passed", passed);
    return passed;
}

int main(int argc, char** argv) {
    if (argc < 2 || argc > 4)
        throw std::invalid_argument(
            "Expected restart input card, optional contribution step and optional captured state");
    const auto input = read_case_input(argv[1]);
    if (input.restart_file.empty())
        throw std::invalid_argument("Contact tangent probe requires a restart checkpoint");
    const auto mesh = read_exodus_hex20(input.mesh_file);
    TransientProblem problem(input.spatial, mesh);
    const double next_step = restore_transient_checkpoint(input.restart_file, problem);
    if (!(input.transient_execution.end_time > problem.committed_time()))
        throw std::invalid_argument("Contact tangent probe requires a remaining time step");
    double end = std::min(input.transient_execution.end_time, problem.committed_time() + next_step);
    for (const double event : problem.time_events())
        if (event > problem.committed_time() && event < end)
            end = event;
    const double factor = input.transient_execution.load_ramp_time == 0.0
                              ? 1.0
                              : std::min(end / input.transient_execution.load_ramp_time, 1.0);
    problem.begin_time_step({end, factor, input.transient_execution.include_thermal_time_term});
    std::cout << "committed_time=" << problem.committed_time() << "\nattempted_end_time=" << end << "\n";
    auto state = problem.committed_solution();
    for (const auto& bc : problem.dirichlet_conditions())
        state[bc.dof] = bc.value;
    if (argc == 4) {
        std::ifstream captured(argv[3]);
        for (double& value : state)
            if (!(captured >> value) || !std::isfinite(value))
                throw std::runtime_error("Invalid captured state");
        std::string extra;
        if (captured >> extra)
            throw std::runtime_error("Captured state has extra entries");
    }
    Report report;
    for (double step : {1e-4, 1e-3, 1e-2, 1e-1}) {
        std::cout << "directional_step=" << step << "\n";
        write_jacobian_check(problem, state, report, step);
    }
    const std::size_t size = problem.dof_count(), count = problem.contribution_count();
    std::vector<double> direction(size);
    for (const auto& field : problem.field_layout())
        for (std::size_t i = field.begin; i < field.end; ++i) {
            const double v = static_cast<double>((i - field.begin) % 7);
            direction[i] = field.category == FieldCategory::thermal ? 0.25 + 0.05 * v : 1e-6 * (0.4 + 0.1 * v);
        }
    const double h = argc >= 3 ? std::stod(argv[2]) : 0.01;
    if (!(h > 0.0) || !std::isfinite(h))
        throw std::invalid_argument("Contribution step must be finite and positive");
    std::cout << "contribution_step=" << h << "\n";
    std::vector<std::vector<double>> values(4 * count, std::vector<double>(size, 0.0));
    for (std::size_t phase = 0; phase < 3; ++phase) {
        auto trial = state;
        for (std::size_t i = 0; i < size; ++i)
            trial[i] += (phase == 1 ? h : (phase == 2 ? -h : 0.0)) * direction[i];
        problem.validate_state(trial);
        const auto contact = cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, trial);
        const auto source = cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
        for (std::size_t i = 0; i < contact.size(); ++i) {
            const auto& node = contact[i];
            std::cout << "contact_phase=" << phase << " node=" << source[i] + 1 << " primary_face=" << node.primary_face
                      << " sliding=" << node.sliding << " gap=" << node.gap
                      << " constraint_pressure=" << node.constraint_pressure << " elastic_slip="
                      << std::hypot(node.elastic_tangential_slip[0],
                             node.elastic_tangential_slip[1],
                             node.elastic_tangential_slip[2])
                      << "\n";
        }
        ContributionWorkspace work;
        for (std::size_t e = 0; e < count; ++e) {
            problem.evaluate_contribution(e, trial, work, phase == 0);
            for (std::size_t r = 0; r < work.dofs.size(); ++r) {
                double value = work.residual[r];
                if (phase == 0)
                    values[3 * count + e][work.dofs[r]] += value;
                if (phase == 0) {
                    value = 0;
                    for (std::size_t c = 0; c < work.dofs.size(); ++c)
                        value += work.jacobian[r * work.dofs.size() + c] * direction[work.dofs[c]];
                }
                values[phase * count + e][work.dofs[r]] += value;
            }
        }
    }
    for (std::size_t e = 0; e < count; ++e) {
        double error = 0, scale = 0, forward_error = 0, backward_error = 0;
        for (std::size_t i = 0; i < size; ++i) {
            double v = (values[count + e][i] - values[2 * count + e][i]) / (2 * h);
            error = std::hypot(error, values[e][i] - v);
            scale = std::hypot(scale, v);
            forward_error =
                std::hypot(forward_error, values[e][i] - (values[count + e][i] - values[3 * count + e][i]) / h);
            backward_error =
                std::hypot(backward_error, values[e][i] - (values[3 * count + e][i] - values[2 * count + e][i]) / h);
        }
        std::cout << "contribution=" << e << " error=" << error << " scale=" << scale
                  << " relative=" << (scale > 0 ? error / scale : 0) << " forward_error=" << forward_error
                  << " backward_error=" << backward_error << "\n";
    }
    problem.rollback_time_step();
}
