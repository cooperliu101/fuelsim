#pragma once
#include "core/nonlinear_problem.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fuelsim::test {
struct FieldNorms final {
    std::vector<double> l2, maximum_absolute;
};

struct DirectionalJacobianCheck final {
    FieldNorms residual, analytic_directional_derivative;
    FieldNorms finite_difference_directional_derivative, difference;
};

inline std::vector<double> assembled_residual(const NonlinearProblem& problem, const std::vector<double>& state) {
    std::vector<double> result(problem.dof_count());
    ContributionWorkspace workspace;
    problem.validate_state(state);
    for (std::size_t entry = 0; entry < problem.contribution_count(); ++entry) {
        problem.evaluate_contribution(entry, state, workspace, false);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            result[workspace.dofs[local]] += workspace.residual[local];
    }
    return result;
}

inline std::vector<double> constrained_residual(const NonlinearProblem& problem, const std::vector<double>& state) {
    std::vector<double> result = assembled_residual(problem, state);
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        result.at(condition.dof) = state.at(condition.dof) - condition.value;
    return result;
}

inline FieldNorms field_norms(const NonlinearProblem& problem, const std::vector<double>& values) {
    FieldNorms result;
    result.l2.resize(problem.field_layout().size());
    result.maximum_absolute.resize(problem.field_layout().size());
    for (std::size_t field = 0; field < problem.field_layout().size(); ++field) {
        const FieldDescriptor& descriptor = problem.field_layout()[field];
        for (std::size_t dof = descriptor.begin; dof < descriptor.end; ++dof) {
            result.l2[field] = std::hypot(result.l2[field], values.at(dof));
            result.maximum_absolute[field] = std::max(result.maximum_absolute[field], std::abs(values.at(dof)));
        }
    }
    return result;
}

inline DirectionalJacobianCheck check_directional_jacobian(const NonlinearProblem& problem,
    const std::vector<double>& state,
    const std::vector<double>& direction,
    double step) {
    if (state.size() != problem.dof_count() || direction.size() != problem.dof_count() || !std::isfinite(step)
        || !(step > 0.0))
        throw std::invalid_argument("Directional Jacobian check inputs do not match the problem");
    problem.validate_discretization();
    problem.validate_state(state);
    std::vector<double> analytic(problem.dof_count(), 0.0), plus = state, minus = state;
    ContributionWorkspace workspace;
    for (std::size_t entry = 0; entry < problem.contribution_count(); ++entry) {
        problem.evaluate_contribution(entry, state, workspace, true);
        const std::size_t local_count = workspace.dofs.size();
        for (std::size_t row = 0; row < local_count; ++row)
            for (std::size_t column = 0; column < local_count; ++column)
                analytic.at(workspace.dofs[row]) +=
                    workspace.jacobian[row * local_count + column] * direction.at(workspace.dofs[column]);
    }
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        analytic.at(condition.dof) = direction.at(condition.dof);
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const std::vector<double> residual = constrained_residual(problem, state);
    const std::vector<double> plus_residual = constrained_residual(problem, plus);
    const std::vector<double> minus_residual = constrained_residual(problem, minus);
    std::vector<double> finite_difference(problem.dof_count()), difference(problem.dof_count());
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        finite_difference[dof] = (plus_residual[dof] - minus_residual[dof]) / (2.0 * step);
        difference[dof] = analytic[dof] - finite_difference[dof];
    }
    return {field_norms(problem, residual),
        field_norms(problem, analytic),
        field_norms(problem, finite_difference),
        field_norms(problem, difference)};
}
} // namespace fuelsim::test
