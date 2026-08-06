#include "fuelsim/nonlinear_problem.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace fuelsim {

GlobalStateView::GlobalStateView(const std::vector<double>& dense_values)
    : _global_size(dense_values.size()), _dense_values(&dense_values),
      _global_dofs(nullptr), _sparse_values(nullptr) {}

GlobalStateView::GlobalStateView(
    std::size_t global_size, const std::vector<std::uint32_t>& global_dofs,
    const std::vector<double>& values)
    : _global_size(global_size), _dense_values(nullptr),
      _global_dofs(&global_dofs), _sparse_values(&values) {
    if (global_dofs.size() != values.size())
        throw std::invalid_argument(
            "GlobalStateView sparse index and value sizes differ");
    if (!std::is_sorted(global_dofs.begin(), global_dofs.end()) ||
        std::adjacent_find(global_dofs.begin(), global_dofs.end()) !=
            global_dofs.end())
        throw std::invalid_argument(
            "GlobalStateView sparse DOFs must be sorted and unique");
    if (!global_dofs.empty() && global_dofs.back() >= global_size)
        throw std::out_of_range(
            "GlobalStateView sparse DOF exceeds the global size");
}

std::size_t GlobalStateView::global_size() const noexcept {
    return _global_size;
}

std::size_t GlobalStateView::local_size() const noexcept {
    return _dense_values != nullptr ? _dense_values->size()
                                    : _sparse_values->size();
}

bool GlobalStateView::contains(std::size_t global_dof) const {
    if (global_dof >= _global_size)
        return false;
    if (_dense_values != nullptr)
        return true;
    if (global_dof >
        static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
        return false;
    return std::binary_search(_global_dofs->begin(), _global_dofs->end(),
                              static_cast<std::uint32_t>(global_dof));
}

double GlobalStateView::value(std::size_t global_dof) const {
    if (global_dof >= _global_size)
        throw std::out_of_range(
            "GlobalStateView requested DOF exceeds the global size");
    if (_dense_values != nullptr)
        return _dense_values->at(global_dof);
    if (global_dof >
        static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
        throw std::out_of_range(
            "GlobalStateView requested DOF exceeds sparse index range");
    const auto found =
        std::lower_bound(_global_dofs->begin(), _global_dofs->end(),
                         static_cast<std::uint32_t>(global_dof));
    if (found == _global_dofs->end() ||
        static_cast<std::size_t>(*found) != global_dof)
        throw std::out_of_range(
            "GlobalStateView requested DOF is absent from the shadow state");
    return _sparse_values->at(
        static_cast<std::size_t>(found - _global_dofs->begin()));
}

void NonlinearProblem::validate_state(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "NonlinearProblem validation state size does not match problem");
}

std::vector<std::size_t> NonlinearProblem::required_state_dofs(
    std::size_t contribution_begin, std::size_t contribution_end) const {
    if (contribution_begin > contribution_end ||
        contribution_end > contribution_count())
        throw std::out_of_range(
            "NonlinearProblem contribution range is invalid");
    std::vector<std::size_t> result;
    result.reserve((contribution_end - contribution_begin) * local_dof_count);
    for (std::size_t contribution = contribution_begin;
         contribution < contribution_end; ++contribution) {
        const LocalDofs dofs = contribution_dofs(contribution);
        for (const std::size_t dof : dofs) {
            if (dof >= dof_count())
                throw std::out_of_range(
                    "NonlinearProblem contribution DOF is out of range");
            result.push_back(dof);
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

void NonlinearProblem::validate_local_state(
    std::size_t contribution_begin, std::size_t contribution_end,
    const GlobalStateView& state) const {
    if (contribution_begin > contribution_end ||
        contribution_end > contribution_count())
        throw std::out_of_range(
            "NonlinearProblem local contribution range is invalid");
    if (state.global_size() != dof_count())
        throw std::invalid_argument(
            "NonlinearProblem shadow state size does not match problem");
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

LocalValues NonlinearProblem::contribution_state(
    std::size_t contribution_index,
    const GlobalStateView& global_state) const {
    if (global_state.global_size() != dof_count())
        throw std::invalid_argument(
            "NonlinearProblem shadow state size does not match DOF count");

    const LocalDofs dofs = contribution_dofs(contribution_index);
    LocalValues local_state{};
    for (std::size_t local = 0; local < dofs.size(); ++local)
        local_state[local] = global_state.value(dofs[local]);
    return local_state;
}

void NonlinearProblem::assemble_residual(const std::vector<double>& state,
                                         std::vector<double>& residual) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "NonlinearProblem state size does not match DOF count");
    validate_state(state);

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
