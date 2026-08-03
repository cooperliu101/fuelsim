#include "fuelsim/nonlinear_problem.hpp"

#include <stdexcept>

namespace fuelsim {

void NonlinearProblem::validate_state(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "NonlinearProblem validation state size does not match problem");
}

LocalValues NonlinearProblem::contribution_state(
    std::size_t contribution_index,
    const std::vector<double>& global_state) const {
    if (global_state.size() != dof_count())
        throw std::invalid_argument(
            "NonlinearProblem global state size does not match DOF count");

    const LocalDofs dofs = contribution_dofs(contribution_index);
    LocalValues local_state{};
    for (std::size_t local = 0; local < dofs.size(); ++local) {
        if (dofs[local] >= global_state.size())
            throw std::out_of_range(
                "NonlinearProblem contribution DOF is out of range");
        local_state[local] = global_state[dofs[local]];
    }
    return local_state;
}

void NonlinearProblem::assemble_residual(const std::vector<double>& state,
                                         std::vector<double>& residual) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "NonlinearProblem state size does not match DOF count");

    residual.assign(dof_count(), 0.0);
    for (std::size_t contribution = 0; contribution < contribution_count();
         ++contribution) {
        const LocalDofs dofs = contribution_dofs(contribution);
        const LocalValues local_state = contribution_state(contribution, state);
        const LocalResidual local_residual =
            contribution_residual(contribution, local_state);

        for (std::size_t local = 0; local < dofs.size(); ++local) {
            if (dofs[local] >= residual.size())
                throw std::out_of_range(
                    "NonlinearProblem contribution DOF is out of range");
            residual[dofs[local]] += local_residual[local];
        }
    }

    add_state_independent_residual(residual);
    if (residual.size() != dof_count())
        throw std::logic_error(
            "NonlinearProblem residual hook changed the residual size");
}

void NonlinearProblem::assemble_state_independent_residual(
    std::vector<double>& residual) const {
    residual.assign(dof_count(), 0.0);
    add_state_independent_residual(residual);
    if (residual.size() != dof_count())
        throw std::logic_error(
            "NonlinearProblem residual hook changed the residual size");
}

} // namespace fuelsim
