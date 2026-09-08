#include "fuelsim/core/nonlinear_problem.hpp"
#include "fuelsim/core/spatial_definition.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace fuelsim {
void ContributionWorkspace::reserve(std::size_t maximum_dof_count) {
    if (maximum_dof_count > 0 && maximum_dof_count > std::numeric_limits<std::size_t>::max() / maximum_dof_count)
        throw std::length_error("Local contribution Jacobian size overflows");
    dofs.reserve(maximum_dof_count);
    state.reserve(maximum_dof_count);
    residual.reserve(maximum_dof_count);
    jacobian.reserve(maximum_dof_count * maximum_dof_count);
}

void ContributionWorkspace::resize(std::size_t dof_count_value, bool include_jacobian) {
    if (dof_count_value == 0)
        throw std::invalid_argument("Local contribution must contain at least one DOF");
    if (dof_count_value > std::numeric_limits<std::size_t>::max() / dof_count_value)
        throw std::length_error("Local contribution Jacobian size overflows");
    dofs.resize(dof_count_value);
    state.resize(dof_count_value);
    residual.assign(dof_count_value, 0.0);
    if (include_jacobian)
        jacobian.assign(dof_count_value * dof_count_value, 0.0);
    else
        jacobian.clear();
}

void NonlinearProblem::validate_state(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("NonlinearProblem validation state size does not match problem");
}

bool NonlinearProblem::uses_augmented_contact() const noexcept {
    return false;
}

std::size_t NonlinearProblem::sparsity_contribution_count() const noexcept {
    return contribution_count();
}

bool NonlinearProblem::jacobian_sparsity_is_state_dependent() const noexcept {
    return false;
}

bool NonlinearProblem::contribution_metadata_is_fixed() const noexcept {
    return false;
}

std::pair<std::size_t, std::size_t> NonlinearProblem::contribution_partition(std::size_t partition,
    std::size_t partition_count) const {
    if (partition_count == 0 || partition >= partition_count)
        throw std::out_of_range("NonlinearProblem contribution partition is out of range");
    return {contribution_count() * partition / partition_count,
        contribution_count() * (partition + 1U) / partition_count};
}

void NonlinearProblem::contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
    std::vector<std::size_t> dofs;
    contribution_dofs(index, dofs);
    pattern.assign(dofs.size() * dofs.size(), 1U);
}

void NonlinearProblem::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    contribution_dofs(index, dofs);
}

void NonlinearProblem::sparsity_contribution_jacobian_pattern(std::size_t index,
    std::vector<unsigned char>& pattern) const {
    contribution_jacobian_pattern(index, pattern);
}

AugmentedContactUpdate NonlinearProblem::update_augmented_contact_multipliers(const std::vector<double>&, std::size_t) {
    throw std::logic_error("NonlinearProblem does not support augmented contact");
}

void NonlinearProblem::validate_discretization() const {
    const std::vector<FieldDescriptor>& fields = field_layout();
    if (dof_count() == 0 || fields.empty())
        throw std::invalid_argument("NonlinearProblem requires DOFs and field metadata");
    std::size_t expected_begin = 0;
    for (std::size_t field = 0; field < fields.size(); ++field) {
        const FieldDescriptor& descriptor = fields[field];
        if (descriptor.name.empty() || descriptor.begin != expected_begin || descriptor.end <= descriptor.begin
            || descriptor.end > dof_count())
            throw std::invalid_argument("NonlinearProblem fields must be named, nonempty, contiguous ranges");
        for (std::size_t previous = 0; previous < field; ++previous)
            if (fields[previous].name == descriptor.name)
                throw std::invalid_argument("NonlinearProblem field names must be unique");
        expected_begin = descriptor.end;
    }
    if (expected_begin != dof_count())
        throw std::invalid_argument("NonlinearProblem field ranges must cover every DOF");
    std::vector<std::size_t> dofs;
    for (std::size_t entry = 0; entry < contribution_count(); ++entry) {
        contribution_dofs(entry, dofs);
        const std::size_t local_count = dofs.size();
        if (local_count == 0 || local_count > std::numeric_limits<std::size_t>::max() / local_count)
            throw std::invalid_argument("NonlinearProblem contribution size is invalid");
    }
    for (std::size_t entry = 0; entry < sparsity_contribution_count(); ++entry) {
        sparsity_contribution_dofs(entry, dofs);
        if (dofs.empty())
            throw std::invalid_argument("NonlinearProblem sparsity contribution is empty");
        for (const std::size_t dof : dofs)
            if (dof >= dof_count())
                throw std::out_of_range("NonlinearProblem sparsity contribution DOF is out of range");
    }
}

std::size_t NonlinearProblem::field_index(std::size_t dof) const {
    if (dof >= dof_count())
        throw std::out_of_range("NonlinearProblem field lookup DOF is out of range");
    const std::vector<FieldDescriptor>& fields = field_layout();
    const auto found =
        std::find_if(fields.begin(), fields.end(), [dof](const FieldDescriptor& field) { return dof < field.end; });
    if (found == fields.end() || dof < found->begin)
        throw std::logic_error("NonlinearProblem field metadata does not cover a DOF");
    return static_cast<std::size_t>(found - fields.begin());
}

std::vector<std::size_t> NonlinearProblem::required_state_dofs(std::size_t first, std::size_t last) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("NonlinearProblem contribution range is invalid");
    std::vector<std::size_t> result;
    std::vector<std::size_t> dofs;
    for (std::size_t entry = first; entry < last; ++entry) {
        contribution_dofs(entry, dofs);
        for (const std::size_t dof : dofs) {
            if (dof >= dof_count())
                throw std::out_of_range("NonlinearProblem contribution DOF is out of range");
            result.push_back(dof);
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

void NonlinearProblem::validate_local_state(std::size_t first,
    std::size_t last,
    const std::vector<double>& state) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("NonlinearProblem local contribution range is invalid");
    if (state.size() != dof_count())
        throw std::invalid_argument("NonlinearProblem shadow state size does not match problem");
}

void NonlinearProblem::evaluate_contribution(std::size_t index,
    const std::vector<double>& global_state,
    ContributionWorkspace& workspace,
    bool linearize) const {
    if (global_state.size() != dof_count())
        throw std::invalid_argument("NonlinearProblem shadow state size does not match DOF count");
    if (index >= contribution_count())
        throw std::out_of_range("NonlinearProblem contribution index is out of range");
    contribution_dofs(index, workspace.dofs);
    const std::size_t local_count = workspace.dofs.size();
    workspace.resize(local_count, linearize);
    for (std::size_t local = 0; local < local_count; ++local) {
        if (workspace.dofs[local] >= dof_count())
            throw std::out_of_range("NonlinearProblem contribution DOF is out of range");
        workspace.state[local] = global_state.at(workspace.dofs[local]);
    }
    compute_contribution(index, workspace.state, workspace.residual, linearize ? &workspace.jacobian : nullptr);
    if (workspace.residual.size() != local_count
        || (linearize && workspace.jacobian.size() != local_count * local_count))
        throw std::logic_error("NonlinearProblem contribution output has the wrong size");
}

} // namespace fuelsim
