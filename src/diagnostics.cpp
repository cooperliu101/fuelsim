#include "fuelsim/diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

std::vector<double>
analytic_directional_derivative(const NonlinearProblem& problem,
                                const std::vector<double>& state,
                                const std::vector<double>& direction) {
    problem.validate_state(state);
    std::vector<double> result(problem.dof_count(), 0.0);
    for (std::size_t contribution = 0;
         contribution < problem.contribution_count(); ++contribution) {
        const LocalDofs dofs = problem.contribution_dofs(contribution);
        const LocalValues local_state =
            problem.contribution_state(contribution, state);
        const LocalSystem local =
            problem.linearize_contribution(contribution, local_state);
        for (std::size_t row = 0; row < local_dof_count; ++row) {
            double value = 0.0;
            for (std::size_t column = 0; column < local_dof_count; ++column)
                value += local.jacobian[row * local_dof_count + column] *
                         direction.at(dofs[column]);
            result.at(dofs[row]) += value;
        }
    }
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        result.at(condition.dof) = direction.at(condition.dof);
    return result;
}

} // namespace

std::vector<double> constrained_residual(const NonlinearProblem& problem,
                                         const std::vector<double>& state) {
    std::vector<double> result;
    problem.assemble_residual(state, result);
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        result.at(condition.dof) = state.at(condition.dof) - condition.value;
    return result;
}

FieldNorms field_norms(const DofMap& dof_map,
                       const std::vector<double>& values) {
    if (values.size() != dof_map.dof_count())
        throw std::invalid_argument(
            "Field norm vector size does not match the DOF map");
    FieldNorms result;
    for (std::size_t node = 0; node < dof_map.node_count(); ++node) {
        const std::array<std::size_t, 3> dofs = {
            dof_map.temperature(node), dof_map.radial_displacement(node),
            dof_map.axial_displacement(node)};
        for (std::size_t field = 0; field < dofs.size(); ++field) {
            const double value = values[dofs[field]];
            result.l2[field] = std::hypot(result.l2[field], value);
            result.maximum_absolute[field] =
                std::max(result.maximum_absolute[field], std::abs(value));
        }
    }
    return result;
}

DirectionalJacobianCheck
check_directional_jacobian(const NonlinearProblem& problem,
                           const DofMap& dof_map,
                           const std::vector<double>& state,
                           const std::vector<double>& direction, double step) {
    if (state.size() != problem.dof_count() ||
        direction.size() != problem.dof_count())
        throw std::invalid_argument(
            "Directional Jacobian vectors do not match the problem");
    if (!std::isfinite(step) || !(step > 0.0))
        throw std::invalid_argument(
            "Directional Jacobian step must be finite and positive");
    std::vector<double> plus = state;
    std::vector<double> minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const std::vector<double> residual = constrained_residual(problem, state);
    const std::vector<double> plus_residual =
        constrained_residual(problem, plus);
    const std::vector<double> minus_residual =
        constrained_residual(problem, minus);
    const std::vector<double> analytic =
        analytic_directional_derivative(problem, state, direction);
    std::vector<double> finite_difference(problem.dof_count(), 0.0);
    std::vector<double> difference(problem.dof_count(), 0.0);
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        finite_difference[dof] =
            (plus_residual[dof] - minus_residual[dof]) / (2.0 * step);
        difference[dof] = analytic[dof] - finite_difference[dof];
    }
    return {field_norms(dof_map, residual), field_norms(dof_map, analytic),
            field_norms(dof_map, finite_difference),
            field_norms(dof_map, difference)};
}

} // namespace fuelsim
