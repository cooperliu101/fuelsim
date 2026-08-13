#include "assembly.hpp"
#include "cartesian3d_assembly.hpp"
#include "fuelsim/hex8.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/spatial_definition.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"
#include "problem_backend_access.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
namespace fuelsim {
namespace {
void finalize_conservation(const NonlinearProblem& problem, const std::vector<double>& current,
    const std::vector<double>& old, const std::vector<double>& raw_residual, TransientConservationSummary& result) {
    std::vector<bool> constrained(problem.dof_count(), false);
    for (const DirichletCondition& condition : problem.dirichlet_conditions()) {
        constrained[condition.dof] = true;
        if (problem.field_layout()[problem.field_index(condition.dof)].category == FieldCategory::thermal)
            result.dirichlet_heat_input_rate += raw_residual[condition.dof];
        else
            result.dirichlet_reaction_work_increment +=
                raw_residual[condition.dof] * (current[condition.dof] - old[condition.dof]);
    }
    double thermal_residual_squared = 0.0, mechanical_residual_squared = 0.0;
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        if (constrained[dof]) continue;
        double& norm_squared = problem.field_layout()[problem.field_index(dof)].category == FieldCategory::thermal
                                   ? thermal_residual_squared
                                   : mechanical_residual_squared;
        norm_squared += raw_residual[dof] * raw_residual[dof];
    }
    result.unconstrained_thermal_residual_l2 = std::sqrt(thermal_residual_squared);
    result.unconstrained_mechanical_residual_l2 = std::sqrt(mechanical_residual_squared);
    result.global_thermal_balance = result.stored_heat_rate + result.convection_heat_rate +
                                    result.interface_heat_imbalance - result.generated_heat_rate -
                                    result.dirichlet_heat_input_rate;
    const double thermal_scale = std::abs(result.stored_heat_rate) + std::abs(result.convection_heat_rate) +
                                 std::abs(result.interface_heat_imbalance) + std::abs(result.generated_heat_rate) +
                                 std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0 ? std::abs(result.global_thermal_balance) / thermal_scale : 0.0;
    result.mechanical_work_balance = result.internal_mechanical_work_increment + result.contact_work_increment -
                                     result.pressure_traction_work_increment - result.dirichlet_reaction_work_increment;
    const double mechanical_scale =
        std::abs(result.internal_mechanical_work_increment) + std::abs(result.contact_work_increment) +
        std::abs(result.pressure_traction_work_increment) + std::abs(result.dirichlet_reaction_work_increment);
    result.relative_mechanical_work_balance =
        mechanical_scale > 0.0 ? std::abs(result.mechanical_work_balance) / mechanical_scale : 0.0;
}
} // namespace
namespace rz {
class TransientConservationCalculator final {
  public:
    static TransientConservationSummary summarize(const TransientProblem& problem,
        const std::vector<double>& converged_solution,
        const std::vector<std::vector<Quad4MaterialHistory>>& staged_histories,
        const std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>& staged_stresses);
};
} // namespace rz
GlobalStateView::GlobalStateView(const std::vector<double>& dense_values)
    : _global_size(dense_values.size()), _dense_values(&dense_values), _global_dofs(nullptr), _sparse_values(nullptr) {}
GlobalStateView::GlobalStateView(
    std::size_t global_size, const std::vector<std::uint32_t>& global_dofs, const std::vector<double>& values)
    : _global_size(global_size), _dense_values(nullptr), _global_dofs(&global_dofs), _sparse_values(&values) {
    if (global_dofs.size() != values.size())
        throw std::invalid_argument("GlobalStateView sparse index and value sizes differ");
    if (!std::is_sorted(global_dofs.begin(), global_dofs.end()) ||
        std::adjacent_find(global_dofs.begin(), global_dofs.end()) != global_dofs.end())
        throw std::invalid_argument("GlobalStateView sparse DOFs must be sorted and unique");
    if (!global_dofs.empty() && global_dofs.back() >= global_size)
        throw std::out_of_range("GlobalStateView sparse DOF exceeds the global size");
}
std::size_t GlobalStateView::global_size() const noexcept { return _global_size; }
double GlobalStateView::value(std::size_t global_dof) const {
    if (global_dof >= _global_size) throw std::out_of_range("GlobalStateView requested DOF exceeds the global size");
    if (_dense_values != nullptr) return _dense_values->at(global_dof);
    if (global_dof > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
        throw std::out_of_range("GlobalStateView requested DOF exceeds sparse index range");
    const auto found =
        std::lower_bound(_global_dofs->begin(), _global_dofs->end(), static_cast<std::uint32_t>(global_dof));
    if (found == _global_dofs->end() || static_cast<std::size_t>(*found) != global_dof)
        throw std::out_of_range("GlobalStateView requested DOF is absent from the shadow state");
    return _sparse_values->at(static_cast<std::size_t>(found - _global_dofs->begin()));
}
void ContributionWorkspace::reserve(std::size_t maximum_dof_count) {
    if (maximum_dof_count > 0 && maximum_dof_count > std::numeric_limits<std::size_t>::max() / maximum_dof_count)
        throw std::length_error("Local contribution Jacobian size overflows");
    dofs.reserve(maximum_dof_count);
    state.reserve(maximum_dof_count);
    residual.reserve(maximum_dof_count);
    jacobian.reserve(maximum_dof_count * maximum_dof_count);
}
void ContributionWorkspace::resize(std::size_t dof_count_value, bool include_jacobian) {
    if (dof_count_value == 0) throw std::invalid_argument("Local contribution must contain at least one DOF");
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
bool NonlinearProblem::uses_augmented_contact() const noexcept { return false; }
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
        if (descriptor.name.empty() || descriptor.begin != expected_begin || descriptor.end <= descriptor.begin ||
            descriptor.end > dof_count())
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
}
std::size_t NonlinearProblem::field_index(std::size_t dof) const {
    if (dof >= dof_count()) throw std::out_of_range("NonlinearProblem field lookup DOF is out of range");
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
            if (dof >= dof_count()) throw std::out_of_range("NonlinearProblem contribution DOF is out of range");
            result.push_back(dof);
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}
void NonlinearProblem::validate_local_state(std::size_t first, std::size_t last, const GlobalStateView& state) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("NonlinearProblem local contribution range is invalid");
    if (state.global_size() != dof_count())
        throw std::invalid_argument("NonlinearProblem shadow state size does not match problem");
}
void NonlinearProblem::gather_contribution_state(std::size_t index, const GlobalStateView& global_state,
    ContributionWorkspace& workspace, bool include_jacobian) const {
    if (global_state.global_size() != dof_count())
        throw std::invalid_argument("NonlinearProblem shadow state size does not match DOF count");
    if (index >= contribution_count()) throw std::out_of_range("NonlinearProblem contribution index is out of range");
    contribution_dofs(index, workspace.dofs);
    const std::size_t local_count = workspace.dofs.size();
    workspace.resize(local_count, include_jacobian);
    for (std::size_t local = 0; local < local_count; ++local) {
        if (workspace.dofs[local] >= dof_count())
            throw std::out_of_range("NonlinearProblem contribution DOF is out of range");
        workspace.state[local] = global_state.value(workspace.dofs[local]);
    }
}
void NonlinearProblem::evaluate_contribution_residual(
    std::size_t index, const GlobalStateView& global_state, ContributionWorkspace& workspace) const {
    gather_contribution_state(index, global_state, workspace, false);
    compute_contribution_residual(index, workspace.state, workspace.residual);
    if (workspace.residual.size() != workspace.dofs.size())
        throw std::logic_error("NonlinearProblem contribution residual has the wrong size");
}
void NonlinearProblem::evaluate_contribution_system(
    std::size_t index, const GlobalStateView& global_state, ContributionWorkspace& workspace) const {
    gather_contribution_state(index, global_state, workspace, true);
    compute_contribution_system(index, workspace.state, workspace.residual, workspace.jacobian);
    if (workspace.residual.size() != workspace.dofs.size() ||
        workspace.jacobian.size() != workspace.dofs.size() * workspace.dofs.size())
        throw std::logic_error("NonlinearProblem contribution system has the wrong size");
}
void NonlinearProblem::assemble_residual(const std::vector<double>& state, std::vector<double>& residual) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("NonlinearProblem state size does not match DOF count");
    validate_discretization();
    validate_state(state);
    residual.assign(dof_count(), 0.0);
    const GlobalStateView global_state(state);
    ContributionWorkspace workspace;
    for (std::size_t entry = 0; entry < contribution_count(); ++entry) {
        evaluate_contribution_residual(entry, global_state, workspace);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            residual[workspace.dofs[local]] += workspace.residual[local];
    }
}
namespace {
LocalValues rz_local_values(const std::vector<double>& values) {
    if (values.size() != local_dof_count) throw std::invalid_argument("RZ contribution state must contain 12 DOFs");
    LocalValues result{};
    std::copy(values.begin(), values.end(), result.begin());
    return result;
}
LocalValues gather_rz_state(
    const rz::SpatialAssembly& spatial, std::size_t index, const GlobalStateView& global_state) {
    const LocalDofs dofs = spatial.contribution_dofs(index);
    LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local) result[local] = global_state.value(dofs[local]);
    return result;
}
LocalValues gather_rz_state(
    const rz::SpatialAssembly& spatial, std::size_t index, const std::vector<double>& global_state) {
    return gather_rz_state(spatial, index, GlobalStateView(global_state));
}
LocalResidual steady_rz_residual(const rz::SpatialAssembly& spatial,
    const std::vector<Quad4RzThermoelasticData>& kernel_data, std::size_t index, const LocalValues& state) {
    if (index >= spatial.volume_contribution_count()) return spatial.contribution_residual(index, state);
    const auto location = spatial.element_location(index);
    return compute_quad4_rz_thermoelastic_residual(
        kernel_data[location.first], spatial.region_element_geometry(location.first, location.second), state);
}
LocalSystem steady_rz_system(const rz::SpatialAssembly& spatial,
    const std::vector<Quad4RzThermoelasticData>& kernel_data, std::size_t index, const LocalValues& state) {
    if (index >= spatial.volume_contribution_count()) return spatial.linearize_contribution(index, state);
    const auto location = spatial.element_location(index);
    return compute_quad4_rz_thermoelastic_system(
        kernel_data[location.first], spatial.region_element_geometry(location.first, location.second), state);
}
LocalResidual transient_rz_residual(
    const rz::TransientBackendView& backend, std::size_t index, const LocalValues& state) {
    if (index >= backend.spatial.volume_contribution_count())
        return backend.spatial.contribution_residual(index, state);
    const auto location = backend.spatial.element_location(index);
    return compute_quad4_rz_transient_residual(backend.kernel_data[location.first],
        backend.spatial.region_element_geometry(location.first, location.second), state,
        gather_rz_state(backend.spatial, index, backend.committed_solution),
        backend.histories[location.first][location.second], backend.active_time_step);
}
LocalSystem transient_rz_system(const rz::TransientBackendView& backend, std::size_t index, const LocalValues& state) {
    if (index >= backend.spatial.volume_contribution_count())
        return backend.spatial.linearize_contribution(index, state);
    const auto location = backend.spatial.element_location(index);
    return compute_quad4_rz_transient_system(backend.kernel_data[location.first],
        backend.spatial.region_element_geometry(location.first, location.second), state,
        gather_rz_state(backend.spatial, index, backend.committed_solution),
        backend.histories[location.first][location.second], backend.active_time_step);
}
void copy_rz_residual(const LocalResidual& source, std::vector<double>& destination) {
    destination.assign(source.begin(), source.end());
}
void copy_rz_system(const LocalSystem& source, std::vector<double>& residual, std::vector<double>& jacobian) {
    residual.assign(source.residual.begin(), source.residual.end());
    jacobian.assign(source.jacobian.begin(), source.jacobian.end());
}
class SpatialProblemBackend {
  public:
    SpatialProblemBackend(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
        : rz(std::make_unique<rz::SpatialAssembly>(std::move(definition), source_mesh)) {}
    SpatialProblemBackend(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
        : cartesian(std::make_unique<cartesian::SpatialAssembly>(std::move(definition), source_mesh)) {}
    bool is_cartesian() const noexcept { return cartesian != nullptr; }
    const spatial_detail::SpatialLayout& layout() const noexcept {
        return is_cartesian() ? static_cast<const spatial_detail::SpatialLayout&>(*cartesian) : *rz;
    }
    const SpatialDefinition& definition() const noexcept { return layout().definition(); }
    std::size_t contribution_count() const noexcept {
        return is_cartesian() ? cartesian->contribution_count() : rz->contribution_count();
    }
    SpatialContributionType contribution_type(std::size_t index) const {
        return is_cartesian() ? cartesian->contribution_type(index) : rz->contribution_type(index);
    }
    void set_load_factor(double value) {
        if (is_cartesian())
            cartesian->set_load_factor(value);
        else
            rz->set_load_factor(value);
    }
    void set_time(double value) {
        if (is_cartesian())
            cartesian->set_time(value);
        else
            rz->set_time(value);
    }
    void validate_state(const std::vector<double>& state) const {
        if (is_cartesian())
            cartesian->validate_state(state);
        else
            rz->validate_state(state);
    }
    std::vector<std::size_t> required_state_dofs(
        const NonlinearProblem& problem, std::size_t first, std::size_t last) const {
        return is_cartesian() ? problem.NonlinearProblem::required_state_dofs(first, last)
                              : rz->required_state_dofs(first, last);
    }
    void validate_local_state(
        const NonlinearProblem& problem, std::size_t first, std::size_t last, const GlobalStateView& state) const {
        if (is_cartesian())
            problem.NonlinearProblem::validate_local_state(first, last, state);
        else
            rz->validate_local_state(first, last, state);
    }
    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
        if (is_cartesian()) {
            cartesian->contribution_dofs(index, dofs);
            return;
        }
        const LocalDofs fixed = rz->contribution_dofs(index);
        dofs.assign(fixed.begin(), fixed.end());
    }
    std::unique_ptr<rz::SpatialAssembly> rz;
    std::unique_ptr<cartesian::SpatialAssembly> cartesian;
};
} // namespace
class SteadyProblem::Implementation final : public SpatialProblemBackend {
  public:
    Implementation(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
        : SpatialProblemBackend(std::move(definition), source_mesh) {
        rz_kernel_data.reserve(rz->region_count());
        for (std::size_t region = 0; region < rz->region_count(); ++region) {
            const RegionDefinition& value = rz->region(region);
            rz_kernel_data.push_back({IsotropicThermoelasticMaterial(value.material), rz->region_heat_source(region),
                0.0, value.strain_formulation});
        }
    }
    Implementation(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
        : SpatialProblemBackend(std::move(definition), source_mesh) {}
    std::vector<Quad4RzThermoelasticData> rz_kernel_data;
};
SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : _impl(std::make_unique<Implementation>(std::move(definition), source_mesh)) {}
SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : _impl(std::make_unique<Implementation>(std::move(definition), source_mesh)) {}
SteadyProblem::~SteadyProblem() = default;
const cartesian::SpatialAssembly& cartesian::BackendAccess::spatial(const SteadyProblem& problem) noexcept {
    return *problem._impl->cartesian;
}
rz::SteadyBackendView rz::BackendAccess::steady(const SteadyProblem& problem) noexcept {
    return {*problem._impl->rz, problem._impl->rz_kernel_data};
}
bool SteadyProblem::uses_augmented_contact() const noexcept {
    return !_impl->is_cartesian() && _impl->rz->uses_augmented_contact();
}
AugmentedContactUpdate SteadyProblem::update_augmented_contact_multipliers(
    const std::vector<double>& state, std::size_t completed_updates) {
    if (_impl->is_cartesian())
        throw std::logic_error("Cartesian three-dimensional stage B does not support augmented contact");
    return _impl->rz->update_augmented_contact_multipliers(state, completed_updates);
}
void SteadyProblem::set_load_factor(double value) {
    _impl->set_load_factor(value);
    refresh_region_heat_sources();
}
void SteadyProblem::refresh_region_heat_sources() {
    if (_impl->is_cartesian()) return;
    for (std::size_t region = 0; region < _impl->rz->region_count(); ++region)
        _impl->rz_kernel_data[region].volumetric_heat_source = _impl->rz->region_heat_source(region);
}
double SteadyProblem::load_factor() const noexcept {
    return _impl->is_cartesian() ? _impl->cartesian->load_factor() : _impl->rz->load_factor();
}
void SteadyProblem::set_time(double value) {
    if (_impl->is_cartesian()) {
        _impl->cartesian->set_time(value);
        return;
    }
    _impl->set_time(value);
    for (Quad4RzThermoelasticData& data : _impl->rz_kernel_data) data.time = value;
    refresh_region_heat_sources();
}
std::vector<double> SteadyProblem::initial_state() const { return _impl->layout().initial_state(); }
struct SteadyStateSnapshot::Storage final {
    Storage(std::shared_ptr<const void> owner_value, std::vector<std::vector<ContactPointHistory>> histories_value)
        : owner(std::move(owner_value)), histories(std::move(histories_value)) {}
    std::shared_ptr<const void> owner;
    std::vector<std::vector<ContactPointHistory>> histories;
};
SteadyStateSnapshot::SteadyStateSnapshot(std::shared_ptr<const Storage> storage) : _storage(std::move(storage)) {}
SteadyStateSnapshot SteadyProblem::capture_internal_state() const {
    const auto histories = _impl->is_cartesian() ? std::vector<std::vector<ContactPointHistory>>{}
                                                 : _impl->rz->committed_contact_histories();
    return SteadyStateSnapshot(std::make_shared<SteadyStateSnapshot::Storage>(discretization_identity(), histories));
}
void SteadyProblem::restore_internal_state(const SteadyStateSnapshot& snapshot, const std::vector<double>& state) {
    if (snapshot.empty()) throw std::invalid_argument("SteadyProblem cannot restore an empty internal-state snapshot");
    if (snapshot._storage->owner != discretization_identity())
        throw std::invalid_argument("SteadyProblem cannot restore a snapshot from another problem");
    if (_impl->is_cartesian()) {
        if (state.size() != dof_count()) throw std::invalid_argument("SteadyProblem restore state size mismatch");
        return;
    }
    _impl->rz->restore_contact_state(state, snapshot._storage->histories);
}
void SteadyProblem::commit_internal_state(const std::vector<double>& state) {
    if (_impl->is_cartesian()) {
        if (state.size() != dof_count()) throw std::invalid_argument("SteadyProblem commit state size mismatch");
        return;
    }
    _impl->rz->commit_contact_state(state);
}
std::size_t SteadyProblem::dof_count() const noexcept { return _impl->layout().dof_count(); }
std::size_t SteadyProblem::contribution_count() const noexcept { return _impl->contribution_count(); }
const std::vector<FieldDescriptor>& SteadyProblem::field_layout() const noexcept {
    return _impl->layout().dof_map().field_layout();
}
const std::vector<DirichletCondition>& SteadyProblem::dirichlet_conditions() const noexcept {
    return _impl->layout().dirichlet_conditions();
}
void SteadyProblem::validate_state(const std::vector<double>& state) const { _impl->validate_state(state); }
std::vector<std::size_t> SteadyProblem::required_state_dofs(std::size_t first, std::size_t last) const {
    return _impl->required_state_dofs(*this, first, last);
}
void SteadyProblem::validate_local_state(std::size_t first, std::size_t last, const GlobalStateView& state) const {
    _impl->validate_local_state(*this, first, last, state);
}
void SteadyProblem::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->contribution_dofs(index, dofs);
}
void SteadyProblem::compute_contribution_residual(
    std::size_t index, const std::vector<double>& state, std::vector<double>& residual) const {
    if (_impl->is_cartesian()) {
        _impl->cartesian->contribution_residual(index, state, nullptr, 0.0, residual);
        return;
    }
    copy_rz_residual(steady_rz_residual(*_impl->rz, _impl->rz_kernel_data, index, rz_local_values(state)), residual);
}
void SteadyProblem::compute_contribution_system(std::size_t index, const std::vector<double>& state,
    std::vector<double>& residual, std::vector<double>& jacobian) const {
    if (_impl->is_cartesian()) {
        _impl->cartesian->contribution_system(index, state, nullptr, 0.0, residual, jacobian);
        return;
    }
    copy_rz_system(
        steady_rz_system(*_impl->rz, _impl->rz_kernel_data, index, rz_local_values(state)), residual, jacobian);
}
class TransientProblem::Implementation final : public SpatialProblemBackend {
  public:
    Implementation(SpatialDefinition definition_value, const UnstructuredQuad4Mesh& source_mesh)
        : SpatialProblemBackend(std::move(definition_value), source_mesh) {}
    Implementation(SpatialDefinition definition_value, const UnstructuredHex8Mesh& source_mesh)
        : SpatialProblemBackend(std::move(definition_value), source_mesh) {}
    std::vector<Quad4RzTransientData> rz_kernel_data;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> material_stresses;
    TransientConservationSummary last_conservation_summary;
    std::vector<double> committed_solution;
    double committed_time = 0.0, committed_load_factor = 0.0, active_time_step = 0.0, active_end_time = 0.0,
           active_load_factor = 0.0;
    std::vector<std::vector<ContactPointHistory>> active_contact_histories;
    bool time_step_active = false;
};
std::vector<double> TransientProblem::accumulate_contribution_conservation(
    const std::vector<double>& solution, TransientConservationSummary& summary) const {
    std::vector<double> raw_residual(dof_count(), 0.0);
    const GlobalStateView state(solution);
    ContributionWorkspace workspace;
    for (std::size_t entry = 0; entry < contribution_count(); ++entry) {
        evaluate_contribution_residual(entry, state, workspace);
        const SpatialContributionType type = _impl->contribution_type(entry);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local) {
            const std::size_t dof = workspace.dofs[local];
            const double residual = workspace.residual[local];
            raw_residual[dof] += residual;
            if (field_layout()[field_index(dof)].category == FieldCategory::thermal) {
                if (type == SpatialContributionType::thermal_contact)
                    summary.interface_heat_imbalance += residual;
                else if (type == SpatialContributionType::convection)
                    summary.convection_heat_rate += residual;
                continue;
            }
            const double work = residual * (solution[dof] - _impl->committed_solution[dof]);
            if (type == SpatialContributionType::volume)
                summary.internal_mechanical_work_increment += work;
            else if (type == SpatialContributionType::mechanical_contact)
                summary.contact_work_increment += work;
            else if (type == SpatialContributionType::pressure || type == SpatialContributionType::traction)
                summary.pressure_traction_work_increment -= work;
        }
    }
    return raw_residual;
}
namespace rz {
namespace {
double stress_strain_inner_product(
    const AxisymmetricStressValues& stress, const std::array<double, 4>& strain) noexcept {
    return stress.rr * strain[0] + stress.zz * strain[1] + stress.hoop * strain[2] + 2.0 * stress.rz * strain[3];
}
std::array<double, 4> strain_difference(const std::array<double, 4>& current, const std::array<double, 4>& old) {
    std::array<double, 4> result{};
    for (std::size_t component = 0; component < result.size(); ++component)
        result[component] = current[component] - old[component];
    return result;
}
} // namespace
TransientConservationSummary TransientConservationCalculator::summarize(const TransientProblem& problem,
    const std::vector<double>& converged_solution,
    const std::vector<std::vector<Quad4MaterialHistory>>& staged_histories,
    const std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>& staged_stresses) {
    problem.require_active_time_step();
    const TransientBackendView backend = BackendAccess::transient(problem);
    TransientConservationSummary result;
    const std::vector<double> raw_residual = problem.accumulate_contribution_conservation(converged_solution, result);
    for (std::size_t region = 0; region < backend.spatial.region_count(); ++region) {
        const Quad4RzTransientData& kernel_data = backend.kernel_data[region];
        const std::size_t offset = backend.spatial.region_element_offset(region);
        for (std::size_t element = 0; element < staged_histories[region].size(); ++element) {
            const LocalValues current = gather_rz_state(backend.spatial, offset + element, converged_solution);
            const LocalValues old = gather_rz_state(backend.spatial, offset + element, backend.committed_solution);
            const Quad4RzGeometry& geometry = backend.spatial.region_element_geometry(region, element);
            for (std::size_t q = 0; q < geometry.points.size(); ++q) {
                const RzQuadraturePoint& point = geometry.points[q];
                double current_temperature = 0.0, old_temperature = 0.0;
                for (std::size_t node = 0; node < quad4_node_count; ++node) {
                    current_temperature += point.shape[node] * current[node];
                    old_temperature += point.shape[node] * old[node];
                }
                const double heat_capacity = compute_quad4_rz_transient_heat_capacity(
                    kernel_data, current_temperature, point.radius, point.axial_coordinate);
                result.stored_heat_rate += point.weighted_measure * heat_capacity *
                                           (current_temperature - old_temperature) / backend.active_time_step;
                result.generated_heat_rate += point.weighted_measure * kernel_data.volumetric_heat_source;
                const MaterialPointState &old_history = backend.histories[region][element][q],
                                         &new_history = staged_histories[region][element][q];
                const AxisymmetricStressValues &old_stress = backend.stresses[region][element][q],
                                               &new_stress = staged_stresses[region][element][q];
                result.elastic_energy_change +=
                    0.5 * point.weighted_measure *
                    (stress_strain_inner_product(new_stress, new_history.elastic_strain) -
                        stress_strain_inner_product(old_stress, old_history.elastic_strain));
                result.plastic_dissipation_increment +=
                    point.weighted_measure *
                    stress_strain_inner_product(
                        new_stress, strain_difference(new_history.plastic_strain, old_history.plastic_strain));
                result.creep_dissipation_increment +=
                    point.weighted_measure * stress_strain_inner_product(new_stress,
                                                 strain_difference(new_history.creep_strain, old_history.creep_strain));
            }
        }
    }
    finalize_conservation(problem, converged_solution, backend.committed_solution, raw_residual, result);
    return result;
}
} // namespace rz
namespace {
bool finite_stress(const AxisymmetricStressValues& stress) {
    return std::isfinite(stress.rr) && std::isfinite(stress.zz) && std::isfinite(stress.hoop) &&
           std::isfinite(stress.rz);
}
bool valid_material_state(const MaterialPointState& state) {
    for (std::size_t component = 0; component < 4; ++component)
        if (!std::isfinite(state.elastic_strain[component]) || !std::isfinite(state.plastic_strain[component]) ||
            !std::isfinite(state.creep_strain[component]))
            return false;
    return std::isfinite(state.equivalent_plastic_strain) && state.equivalent_plastic_strain >= 0.0 &&
           std::isfinite(state.equivalent_creep_strain) && state.equivalent_creep_strain >= 0.0;
}
} // namespace
TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : _impl(std::make_unique<Implementation>(std::move(definition), source_mesh)) {
    const std::size_t regions = _impl->definition().regions.size();
    _impl->rz_kernel_data.reserve(regions);
    _impl->material_histories.resize(regions);
    _impl->material_stresses.resize(regions);
    for (std::size_t region = 0; region < regions; ++region) {
        _impl->rz_kernel_data.push_back({IsotropicInelasticMaterial(_impl->definition().regions[region].material), 0.0,
            0.0, _impl->definition().regions[region].strain_formulation});
        _impl->material_histories[region].resize(_impl->rz->region_element_count(region));
        _impl->material_stresses[region].resize(_impl->rz->region_element_count(region));
    }
    apply_spatial_controls(0.0, 0.0);
    _impl->committed_solution = _impl->rz->initial_state();
    _impl->rz->restore_contact_state(_impl->committed_solution, _impl->rz->committed_contact_histories());
}
TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : _impl(std::make_unique<Implementation>(std::move(definition), source_mesh)) {
    apply_spatial_controls(0.0, 0.0);
    _impl->committed_solution = _impl->cartesian->initial_state();
}
TransientProblem::~TransientProblem() = default;
bool TransientProblem::is_cartesian_3d() const noexcept { return _impl->is_cartesian(); }
const cartesian::SpatialAssembly& cartesian::BackendAccess::spatial(const TransientProblem& problem) noexcept {
    return *problem._impl->cartesian;
}
TransientCommittedState cartesian::BackendAccess::committed_state(const TransientProblem& problem) {
    return {problem._impl->committed_solution, {}, {}, {}, problem._impl->last_conservation_summary,
        problem._impl->committed_time, problem._impl->committed_load_factor};
}
void cartesian::BackendAccess::restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
    if (problem._impl->time_step_active)
        throw std::logic_error("TransientProblem cannot restore during an active time step");
    if (state.solution.size() != problem.dof_count() || !std::isfinite(state.time) || state.time < 0.0 ||
        !std::isfinite(state.load_factor) || state.load_factor < 0.0)
        throw std::invalid_argument("Cartesian committed state layout is invalid");
    problem._impl->cartesian->validate_state(state.solution);
    problem._impl->committed_solution = std::move(state.solution);
    problem._impl->last_conservation_summary = state.conservation;
    problem._impl->committed_time = state.time;
    problem._impl->committed_load_factor = state.load_factor;
    problem.clear_active_time_step();
    problem.apply_spatial_controls(state.time, state.load_factor);
}
rz::TransientBackendView rz::BackendAccess::transient(const TransientProblem& problem) noexcept {
    return {*problem._impl->rz, problem._impl->rz_kernel_data, problem._impl->material_histories,
        problem._impl->material_stresses, problem._impl->committed_solution, problem._impl->active_time_step,
        problem._impl->time_step_active};
}
const SpatialDefinition& TransientProblem::definition() const noexcept { return _impl->definition(); }
const std::vector<double>& TransientProblem::committed_solution() const noexcept { return _impl->committed_solution; }
double TransientProblem::committed_time() const noexcept { return _impl->committed_time; }
double TransientProblem::committed_load_factor() const noexcept { return _impl->committed_load_factor; }
bool TransientProblem::time_step_active() const noexcept { return _impl->time_step_active; }
std::vector<double> TransientProblem::time_events() const {
    std::vector<double> result;
    for (const PiecewiseLinearTimeTable& table : _impl->definition().time_tables)
        result.insert(result.end(), table.times().begin(), table.times().end());
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}
RegionStateSummary TransientProblem::summarize_region(std::size_t region) const {
    if (region >= _impl->definition().regions.size())
        throw std::out_of_range("TransientProblem region summary index is out of range");
    const std::size_t node_offset = _impl->layout().region_node_offset(region);
    const std::size_t node_count = _impl->is_cartesian() ? _impl->cartesian->region_mesh(region).nodes().size()
                                                         : _impl->rz->region_mesh(region).nodes().size();
    const auto temperature = std::find_if(field_layout().begin(), field_layout().end(),
        [](const FieldDescriptor& field) { return field.category == FieldCategory::thermal; });
    RegionStateSummary result{-std::numeric_limits<double>::infinity(), 0.0, 0.0};
    for (std::size_t node = 0; node < node_count; ++node)
        result.maximum_temperature =
            std::max(result.maximum_temperature, _impl->committed_solution.at(temperature->begin + node_offset + node));
    if (_impl->is_cartesian()) return result;
    for (const Quad4MaterialHistory& element : _impl->material_histories[region])
        for (const MaterialPointState& point : element) {
            result.maximum_equivalent_plastic_strain =
                std::max(result.maximum_equivalent_plastic_strain, point.equivalent_plastic_strain);
            result.maximum_equivalent_creep_strain =
                std::max(result.maximum_equivalent_creep_strain, point.equivalent_creep_strain);
        }
    return result;
}
struct TransientStateSnapshot::Storage final {
    Storage(std::shared_ptr<const void> owner_value, TransientCommittedState state_value)
        : owner(std::move(owner_value)), state(std::move(state_value)) {}
    std::shared_ptr<const void> owner;
    TransientCommittedState state;
};
TransientStateSnapshot::TransientStateSnapshot(std::shared_ptr<const Storage> storage) : _storage(std::move(storage)) {}
TransientCommittedState rz::BackendAccess::committed_state(const TransientProblem& problem) {
    return {problem._impl->committed_solution, problem._impl->material_histories, problem._impl->material_stresses,
        problem._impl->rz->committed_contact_histories(), problem._impl->last_conservation_summary,
        problem._impl->committed_time, problem._impl->committed_load_factor};
}
TransientStateSnapshot TransientProblem::capture_state() const {
    if (_impl->time_step_active) throw std::logic_error("TransientProblem cannot capture an active time step");
    TransientCommittedState state = _impl->is_cartesian() ? cartesian::BackendAccess::committed_state(*this)
                                                          : rz::BackendAccess::committed_state(*this);
    return TransientStateSnapshot(
        std::make_shared<TransientStateSnapshot::Storage>(discretization_identity(), std::move(state)));
}
void TransientProblem::restore_state(const TransientStateSnapshot& snapshot) {
    if (snapshot.empty()) throw std::invalid_argument("TransientProblem cannot restore an empty state snapshot");
    if (snapshot._storage->owner != discretization_identity())
        throw std::invalid_argument("TransientProblem cannot restore a snapshot from another problem");
    if (_impl->is_cartesian())
        cartesian::BackendAccess::restore_committed_state(*this, snapshot._storage->state);
    else
        rz::BackendAccess::restore_committed_state(*this, snapshot._storage->state);
}
void rz::BackendAccess::restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
    if (problem._impl->time_step_active)
        throw std::logic_error("TransientProblem cannot restore during an active time step");
    if (state.solution.size() != problem.dof_count() ||
        state.material_histories.size() != problem._impl->rz->region_count() ||
        state.material_stresses.size() != problem._impl->rz->region_count() ||
        state.contact_histories.size() != problem._impl->definition().contacts.size())
        throw std::invalid_argument("Transient committed state layout does not match the problem");
    if (!std::isfinite(state.time) || state.time < 0.0 || !std::isfinite(state.load_factor) || state.load_factor < 0.0)
        throw std::invalid_argument("Transient committed time and load factor must be valid");
    const DofMap& dofs = problem._impl->rz->dof_map();
    for (std::size_t node = 0; node < dofs.node_count(); ++node) {
        const double temperature = state.solution.at(dofs.temperature(node));
        if (!std::isfinite(temperature) || !(temperature > 0.0) ||
            !std::isfinite(state.solution.at(dofs.radial_displacement(node))) ||
            !std::isfinite(state.solution.at(dofs.axial_displacement(node))))
            throw std::invalid_argument("Transient committed nodal state must be finite with positive temperatures");
    }
    for (std::size_t region = 0; region < problem._impl->rz->region_count(); ++region) {
        const std::size_t elements = problem._impl->rz->region_element_count(region);
        if (state.material_histories[region].size() != elements || state.material_stresses[region].size() != elements)
            throw std::invalid_argument("Transient committed element state layout does not match");
        for (std::size_t element = 0; element < elements; ++element) {
            for (std::size_t q = 0; q < 4; ++q)
                if (!valid_material_state(state.material_histories[region][element][q]) ||
                    !finite_stress(state.material_stresses[region][element][q]))
                    throw std::invalid_argument("Transient committed integration-point state is invalid");
        }
    }
    problem._impl->rz->restore_contact_state(state.solution, state.contact_histories);
    problem._impl->committed_solution = std::move(state.solution);
    problem._impl->material_histories = std::move(state.material_histories);
    problem._impl->material_stresses = std::move(state.material_stresses);
    problem._impl->last_conservation_summary = state.conservation;
    problem._impl->committed_time = state.time;
    problem._impl->committed_load_factor = state.load_factor;
    problem.clear_active_time_step();
    problem.apply_spatial_controls(problem._impl->committed_time, problem._impl->committed_load_factor);
}
namespace rz {
namespace {
struct TimeErrorAccumulator final {
    double difference_squared = 0.0, solution_squared = 0.0;
    std::size_t count = 0;
};
void accumulate_time_error(TimeErrorAccumulator& accumulator, double full_step, double two_half_steps) {
    const double difference = two_half_steps - full_step;
    accumulator.difference_squared += difference * difference;
    accumulator.solution_squared += two_half_steps * two_half_steps;
    ++accumulator.count;
}
double normalized_time_error(
    const TimeErrorAccumulator& accumulator, double absolute_tolerance, double relative_tolerance) {
    if (accumulator.count == 0) return 0.0;
    const double denominator = absolute_tolerance * std::sqrt(static_cast<double>(accumulator.count)) +
                               relative_tolerance * std::sqrt(accumulator.solution_squared);
    return std::sqrt(accumulator.difference_squared) / denominator;
}
TransientTimeErrorEstimate nodal_time_error(const TransientCommittedState& full, const TransientCommittedState& half,
    const std::vector<FieldDescriptor>& fields, const TransientTimeOptions& options) {
    if (full.solution.size() != half.solution.size()) throw std::logic_error("step-doubling nodal layouts differ");
    TransientTimeErrorEstimate result;
    for (const FieldDescriptor& field : fields) {
        TimeErrorAccumulator error;
        for (std::size_t dof = field.begin; dof < field.end; ++dof)
            accumulate_time_error(error, full.solution.at(dof), half.solution.at(dof));
        const double tolerance = field.category == FieldCategory::thermal
                                     ? options.temperature_time_absolute_tolerance
                                     : options.displacement_time_absolute_tolerance;
        const std::string name =
            field.name == "radial" || field.name == "axial" ? field.name + "_displacement" : field.name;
        result.nodal_fields.push_back(
            {name, normalized_time_error(error, tolerance, options.time_error_relative_tolerance)});
        result.maximum = std::max(result.maximum, result.nodal_fields.back().value);
    }
    return result;
}
TransientConservationSummary combine_rz_half_step_conservation(
    const TransientConservationSummary& first, const TransientConservationSummary& second) {
    TransientConservationSummary result;
    const auto average = [](double left, double right) { return 0.5 * (left + right); };
    result.generated_heat_rate = average(first.generated_heat_rate, second.generated_heat_rate);
    result.stored_heat_rate = average(first.stored_heat_rate, second.stored_heat_rate);
    result.convection_heat_rate = average(first.convection_heat_rate, second.convection_heat_rate);
    result.interface_heat_imbalance = average(first.interface_heat_imbalance, second.interface_heat_imbalance);
    result.dirichlet_heat_input_rate = average(first.dirichlet_heat_input_rate, second.dirichlet_heat_input_rate);
    result.global_thermal_balance = result.stored_heat_rate + result.convection_heat_rate +
                                    result.interface_heat_imbalance - result.generated_heat_rate -
                                    result.dirichlet_heat_input_rate;
    const double thermal_scale = std::abs(result.generated_heat_rate) + std::abs(result.stored_heat_rate) +
                                 std::abs(result.convection_heat_rate) + std::abs(result.interface_heat_imbalance) +
                                 std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0 ? std::abs(result.global_thermal_balance) / thermal_scale : 0.0;
    result.unconstrained_thermal_residual_l2 =
        std::max(first.unconstrained_thermal_residual_l2, second.unconstrained_thermal_residual_l2);
    result.internal_mechanical_work_increment =
        first.internal_mechanical_work_increment + second.internal_mechanical_work_increment;
    result.pressure_traction_work_increment =
        first.pressure_traction_work_increment + second.pressure_traction_work_increment;
    result.dirichlet_reaction_work_increment =
        first.dirichlet_reaction_work_increment + second.dirichlet_reaction_work_increment;
    result.contact_work_increment = first.contact_work_increment + second.contact_work_increment;
    result.mechanical_work_balance = result.internal_mechanical_work_increment + result.contact_work_increment -
                                     result.pressure_traction_work_increment - result.dirichlet_reaction_work_increment;
    const double mechanical_scale =
        std::abs(result.internal_mechanical_work_increment) + std::abs(result.pressure_traction_work_increment) +
        std::abs(result.dirichlet_reaction_work_increment) + std::abs(result.contact_work_increment);
    result.relative_mechanical_work_balance =
        mechanical_scale > 0.0 ? std::abs(result.mechanical_work_balance) / mechanical_scale : 0.0;
    result.unconstrained_mechanical_residual_l2 =
        std::max(first.unconstrained_mechanical_residual_l2, second.unconstrained_mechanical_residual_l2);
    result.elastic_energy_change = first.elastic_energy_change + second.elastic_energy_change;
    result.plastic_dissipation_increment = first.plastic_dissipation_increment + second.plastic_dissipation_increment;
    result.creep_dissipation_increment = first.creep_dissipation_increment + second.creep_dissipation_increment;
    return result;
}
TransientTimeErrorEstimate compare_step_doubling_states(const TransientCommittedState& full_step,
    const TransientCommittedState& two_half_steps, const std::vector<FieldDescriptor>& fields,
    std::size_t expected_dof_count, const TransientTimeOptions& options) {
    if (full_step.solution.size() != two_half_steps.solution.size() || full_step.solution.size() != expected_dof_count)
        throw std::logic_error("step-doubling nodal-state layouts differ");
    if (full_step.material_histories.size() != two_half_steps.material_histories.size() ||
        full_step.material_stresses.size() != two_half_steps.material_stresses.size())
        throw std::logic_error("step-doubling material-state region layouts differ");
    TimeErrorAccumulator elastic, plastic;
    TimeErrorAccumulator creep, equivalent_plastic;
    TimeErrorAccumulator equivalent_creep, stress;
    TimeErrorAccumulator contact_friction, contact_normal_multiplier;
    bool contact_state_mismatch = false;
    for (std::size_t region = 0; region < full_step.material_histories.size(); ++region) {
        const auto &full_history = full_step.material_histories[region],
                   &half_history = two_half_steps.material_histories[region];
        const auto &full_stress = full_step.material_stresses[region],
                   &half_stress = two_half_steps.material_stresses[region];
        if (full_history.size() != half_history.size() || full_stress.size() != half_stress.size() ||
            full_history.size() != full_stress.size())
            throw std::logic_error("step-doubling material-state element layouts differ");
        for (std::size_t element = 0; element < full_history.size(); ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const MaterialPointState &full_point = full_history[element][q], &half_point = half_history[element][q];
                for (std::size_t component = 0; component < 4; ++component) {
                    accumulate_time_error(
                        elastic, full_point.elastic_strain[component], half_point.elastic_strain[component]);
                    accumulate_time_error(
                        plastic, full_point.plastic_strain[component], half_point.plastic_strain[component]);
                    accumulate_time_error(
                        creep, full_point.creep_strain[component], half_point.creep_strain[component]);
                }
                accumulate_time_error(
                    equivalent_plastic, full_point.equivalent_plastic_strain, half_point.equivalent_plastic_strain);
                accumulate_time_error(
                    equivalent_creep, full_point.equivalent_creep_strain, half_point.equivalent_creep_strain);
                const AxisymmetricStressValues &full_value = full_stress[element][q],
                                               &half_value = half_stress[element][q];
                accumulate_time_error(stress, full_value.rr, half_value.rr);
                accumulate_time_error(stress, full_value.zz, half_value.zz);
                accumulate_time_error(stress, full_value.hoop, half_value.hoop);
                accumulate_time_error(stress, full_value.rz, half_value.rz);
            }
        }
    }
    if (full_step.contact_histories.size() != two_half_steps.contact_histories.size())
        throw std::logic_error("step-doubling contact-history layouts differ");
    for (std::size_t contact = 0; contact < full_step.contact_histories.size(); ++contact) {
        if (full_step.contact_histories[contact].size() != two_half_steps.contact_histories[contact].size())
            throw std::logic_error("step-doubling contact-node history layouts differ");
        for (std::size_t node = 0; node < full_step.contact_histories[contact].size(); ++node) {
            const ContactPointHistory &full = full_step.contact_histories[contact][node],
                                      &half = two_half_steps.contact_histories[contact][node];
            accumulate_time_error(contact_friction, full.elastic_tangential_slip, half.elastic_tangential_slip);
            accumulate_time_error(contact_normal_multiplier, full.normal_multiplier, half.normal_multiplier);
            contact_state_mismatch = contact_state_mismatch || full.sliding != half.sliding;
        }
    }
    TransientTimeErrorEstimate result = nodal_time_error(full_step, two_half_steps, fields, options);
    result.elastic_strain = normalized_time_error(
        elastic, options.strain_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.plastic_strain = normalized_time_error(
        plastic, options.strain_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.creep_strain = normalized_time_error(
        creep, options.strain_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.equivalent_plastic_strain = normalized_time_error(
        equivalent_plastic, options.strain_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.equivalent_creep_strain = normalized_time_error(
        equivalent_creep, options.strain_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.stress = normalized_time_error(
        stress, options.stress_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.contact_friction =
        contact_state_mismatch ? std::numeric_limits<double>::infinity()
                               : normalized_time_error(contact_friction, options.displacement_time_absolute_tolerance,
                                     options.time_error_relative_tolerance);
    result.contact_normal_multiplier = normalized_time_error(contact_normal_multiplier,
        options.stress_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.maximum =
        std::max({result.elastic_strain, result.plastic_strain, result.creep_strain, result.equivalent_plastic_strain,
            result.equivalent_creep_strain, result.stress, result.contact_friction, result.contact_normal_multiplier});
    for (const TransientFieldTimeError& field : result.nodal_fields)
        result.maximum = std::max(result.maximum, field.value);
    return result;
}
} // namespace
} // namespace rz
TransientTimeErrorEstimate TransientProblem::step_doubling_error(const TransientStateSnapshot& full_snapshot,
    const TransientStateSnapshot& half_snapshot, const TransientTimeOptions& options) const {
    if (full_snapshot.empty() || half_snapshot.empty())
        throw std::invalid_argument("step-doubling requires two complete state snapshots");
    if (full_snapshot._storage->owner != discretization_identity() ||
        half_snapshot._storage->owner != discretization_identity())
        throw std::invalid_argument("step-doubling snapshots belong to another problem");
    if (_impl->is_cartesian()) {
        const TransientCommittedState &full = full_snapshot._storage->state, &half = half_snapshot._storage->state;
        if (full.solution.size() != dof_count() || half.solution.size() != dof_count())
            throw std::logic_error("Cartesian step-doubling snapshot layouts differ");
        TransientTimeErrorEstimate result = rz::nodal_time_error(full, half, field_layout(), options);
        double stress_difference_squared = 0.0, stress_solution_squared = 0.0;
        std::size_t stress_count = 0;
        for (std::size_t region = 0; region < _impl->cartesian->region_count(); ++region) {
            for (std::size_t element = 0; element < _impl->cartesian->region_element_count(region); ++element) {
                const auto full_stresses = _impl->cartesian->stress(region, element, full.solution);
                const auto half_stresses = _impl->cartesian->stress(region, element, half.solution);
                for (std::size_t q = 0; q < 8; ++q) {
                    const SymmetricTensor3Values &full_value = full_stresses[q], &half_value = half_stresses[q];
                    const std::array<double, 6> full_components = {
                        full_value.xx, full_value.yy, full_value.zz, full_value.xy, full_value.yz, full_value.xz};
                    const std::array<double, 6> half_components = {
                        half_value.xx, half_value.yy, half_value.zz, half_value.xy, half_value.yz, half_value.xz};
                    for (std::size_t component = 0; component < 6; ++component) {
                        const double difference = half_components[component] - full_components[component];
                        stress_difference_squared += difference * difference;
                        stress_solution_squared += half_components[component] * half_components[component];
                        ++stress_count;
                    }
                }
            }
        }
        result.stress = rz::normalized_time_error({stress_difference_squared, stress_solution_squared, stress_count},
            options.stress_history_time_absolute_tolerance, options.time_error_relative_tolerance);
        result.maximum = std::max(result.maximum, result.stress);
        return result;
    }
    return rz::compare_step_doubling_states(
        full_snapshot._storage->state, half_snapshot._storage->state, field_layout(), dof_count(), options);
}
void TransientProblem::combine_last_half_step_conservation(const TransientConservationSummary& first_half) {
    if (_impl->time_step_active)
        throw std::logic_error("TransientProblem cannot combine conservation during an active time step");
    _impl->last_conservation_summary =
        rz::combine_rz_half_step_conservation(first_half, _impl->last_conservation_summary);
}
void TransientProblem::begin_time_step(const TransientStepInput& input) {
    if (_impl->time_step_active) throw std::logic_error("TransientProblem already has an active time step");
    if (!std::isfinite(input.end_time) || input.end_time <= _impl->committed_time)
        throw std::invalid_argument("TransientProblem end time must exceed committed time");
    if (!std::isfinite(input.load_factor) || input.load_factor < 0.0)
        throw std::invalid_argument("TransientProblem load factor must be finite and nonnegative");
    _impl->active_time_step = input.end_time - _impl->committed_time;
    _impl->active_end_time = input.end_time;
    _impl->active_load_factor = input.load_factor;
    if (!_impl->is_cartesian() && _impl->rz->uses_augmented_contact())
        _impl->active_contact_histories = _impl->rz->committed_contact_histories();
    try {
        apply_spatial_controls(input.end_time, input.load_factor);
    } catch (...) {
        if (!_impl->is_cartesian() && !_impl->active_contact_histories.empty())
            _impl->rz->restore_contact_state(_impl->committed_solution, std::move(_impl->active_contact_histories));
        apply_spatial_controls(_impl->committed_time, _impl->committed_load_factor);
        clear_active_time_step();
        throw;
    }
    _impl->time_step_active = true;
}
void TransientProblem::commit_time_step(const std::vector<double>& converged_solution) {
    require_active_time_step();
    if (converged_solution.size() != dof_count())
        throw std::invalid_argument("TransientProblem committed solution size mismatch");
    if (!std::all_of(
            converged_solution.begin(), converged_solution.end(), [](double value) { return std::isfinite(value); }))
        throw std::domain_error("TransientProblem committed solution must be finite");
    for (const FieldDescriptor& field : field_layout())
        if (field.category == FieldCategory::thermal)
            for (std::size_t dof = field.begin; dof < field.end; ++dof)
                if (!(converged_solution[dof] > 0.0))
                    throw std::domain_error("TransientProblem committed temperatures must be positive");
    if (_impl->is_cartesian()) {
        TransientConservationSummary conservation;
        const std::vector<double> raw_residual = accumulate_contribution_conservation(converged_solution, conservation);
        for (std::size_t region = 0; region < _impl->cartesian->region_count(); ++region) {
            const std::size_t offset = _impl->cartesian->region_element_offset(region);
            for (std::size_t element = 0; element < _impl->cartesian->region_element_count(region); ++element) {
                const Hex8LocalValues current = _impl->cartesian->volume_state(offset + element, converged_solution);
                const Hex8LocalValues old = _impl->cartesian->volume_state(offset + element, _impl->committed_solution);
                const Hex8Geometry& geometry = _impl->cartesian->region_element_geometry(region, element);
                for (const Hex8QuadraturePoint& point : geometry.points) {
                    double current_temperature = 0.0, old_temperature = 0.0;
                    for (std::size_t node = 0; node < 8; ++node) {
                        current_temperature += point.shape[node] * current[node];
                        old_temperature += point.shape[node] * old[node];
                    }
                    conservation.stored_heat_rate +=
                        point.weighted_measure *
                        _impl->cartesian->heat_capacity(region, current_temperature, point.position) *
                        (current_temperature - old_temperature) / _impl->active_time_step;
                    conservation.generated_heat_rate +=
                        point.weighted_measure * _impl->cartesian->region_heat_source(region);
                }
            }
        }
        finalize_conservation(*this, converged_solution, _impl->committed_solution, raw_residual, conservation);
        _impl->last_conservation_summary = conservation;
    } else {
        const std::size_t regions = _impl->rz->region_count();
        std::vector<std::vector<Quad4MaterialHistory>> staged(regions);
        std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> staged_stresses(regions);
        for (std::size_t region = 0; region < regions; ++region) {
            staged[region].resize(_impl->rz->region_element_count(region));
            staged_stresses[region].resize(_impl->rz->region_element_count(region));
            const std::size_t offset = _impl->rz->region_element_offset(region);
            for (std::size_t element = 0; element < staged[region].size(); ++element) {
                const LocalValues state = gather_rz_state(*_impl->rz, offset + element, converged_solution);
                const LocalValues committed_state =
                    gather_rz_state(*_impl->rz, offset + element, _impl->committed_solution);
                staged[region][element] = compute_quad4_rz_transient_trial_state(_impl->rz_kernel_data[region],
                    _impl->rz->region_element_geometry(region, element), state, committed_state,
                    _impl->material_histories[region][element], _impl->active_time_step);
                staged_stresses[region][element] = compute_quad4_rz_transient_stress(_impl->rz_kernel_data[region],
                    _impl->rz->region_element_geometry(region, element), state, committed_state,
                    _impl->material_histories[region][element], _impl->active_time_step);
            }
        }
        _impl->last_conservation_summary =
            rz::TransientConservationCalculator::summarize(*this, converged_solution, staged, staged_stresses);
        _impl->rz->commit_contact_state(converged_solution);
        _impl->material_histories.swap(staged);
        _impl->material_stresses.swap(staged_stresses);
    }
    _impl->committed_solution = converged_solution;
    _impl->committed_time = _impl->active_end_time;
    _impl->committed_load_factor = _impl->active_load_factor;
    clear_active_time_step();
}
void TransientProblem::rollback_time_step() noexcept {
    if (!_impl->time_step_active) return;
    apply_spatial_controls(_impl->committed_time, _impl->committed_load_factor);
    if (!_impl->is_cartesian() && !_impl->active_contact_histories.empty())
        _impl->rz->restore_contact_state(_impl->committed_solution, std::move(_impl->active_contact_histories));
    clear_active_time_step();
}
void TransientProblem::apply_spatial_controls(double time, double load_factor) {
    _impl->set_time(time);
    _impl->set_load_factor(load_factor);
    if (_impl->is_cartesian()) return;
    for (Quad4RzTransientData& kernel_data : _impl->rz_kernel_data) kernel_data.time = time;
    refresh_region_heat_sources();
}
void TransientProblem::clear_active_time_step() noexcept {
    _impl->active_time_step = 0.0;
    _impl->active_end_time = _impl->committed_time;
    _impl->active_load_factor = _impl->committed_load_factor;
    _impl->active_contact_histories.clear();
    _impl->time_step_active = false;
}
void TransientProblem::refresh_region_heat_sources() {
    if (_impl->is_cartesian()) return;
    for (std::size_t region = 0; region < _impl->rz->region_count(); ++region)
        _impl->rz_kernel_data[region].volumetric_heat_source = _impl->rz->region_heat_source(region);
}
bool TransientProblem::uses_augmented_contact() const noexcept {
    return !_impl->is_cartesian() && _impl->rz->uses_augmented_contact();
}
AugmentedContactUpdate TransientProblem::update_augmented_contact_multipliers(
    const std::vector<double>& state, std::size_t completed_updates) {
    require_active_time_step();
    if (_impl->is_cartesian())
        throw std::logic_error("Cartesian three-dimensional stage B does not support augmented contact");
    return _impl->rz->update_augmented_contact_multipliers(state, completed_updates);
}
const TransientConservationSummary& TransientProblem::last_conservation_summary() const noexcept {
    return _impl->last_conservation_summary;
}
std::size_t TransientProblem::dof_count() const noexcept { return _impl->layout().dof_count(); }
std::size_t TransientProblem::contribution_count() const noexcept { return _impl->contribution_count(); }
const std::vector<FieldDescriptor>& TransientProblem::field_layout() const noexcept {
    return _impl->layout().dof_map().field_layout();
}
const std::vector<DirichletCondition>& TransientProblem::dirichlet_conditions() const noexcept {
    return _impl->layout().dirichlet_conditions();
}
void TransientProblem::validate_state(const std::vector<double>& state) const {
    require_active_time_step();
    _impl->validate_state(state);
}
std::vector<std::size_t> TransientProblem::required_state_dofs(std::size_t first, std::size_t last) const {
    return _impl->required_state_dofs(*this, first, last);
}
void TransientProblem::validate_local_state(std::size_t first, std::size_t last, const GlobalStateView& state) const {
    require_active_time_step();
    _impl->validate_local_state(*this, first, last, state);
}
void TransientProblem::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->contribution_dofs(index, dofs);
}
void TransientProblem::compute_contribution_residual(
    std::size_t index, const std::vector<double>& state, std::vector<double>& residual) const {
    if (_impl->is_cartesian()) {
        require_active_time_step();
        _impl->cartesian->contribution_residual(
            index, state, &_impl->committed_solution, _impl->active_time_step, residual);
        return;
    }
    require_active_time_step();
    copy_rz_residual(
        transient_rz_residual(rz::BackendAccess::transient(*this), index, rz_local_values(state)), residual);
}
void TransientProblem::compute_contribution_system(std::size_t index, const std::vector<double>& state,
    std::vector<double>& residual, std::vector<double>& jacobian) const {
    if (_impl->is_cartesian()) {
        require_active_time_step();
        _impl->cartesian->contribution_system(
            index, state, &_impl->committed_solution, _impl->active_time_step, residual, jacobian);
        return;
    }
    require_active_time_step();
    copy_rz_system(
        transient_rz_system(rz::BackendAccess::transient(*this), index, rz_local_values(state)), residual, jacobian);
}
void TransientProblem::require_active_time_step() const {
    if (!_impl->time_step_active)
        throw std::logic_error("TransientProblem residual evaluation requires an active time step");
}
namespace {
std::vector<double> analytic_directional_derivative(
    const NonlinearProblem& problem, const std::vector<double>& state, const std::vector<double>& direction) {
    problem.validate_discretization();
    problem.validate_state(state);
    std::vector<double> result(problem.dof_count(), 0.0);
    const GlobalStateView global_state(state);
    ContributionWorkspace workspace;
    for (std::size_t entry = 0; entry < problem.contribution_count(); ++entry) {
        problem.evaluate_contribution_system(entry, global_state, workspace);
        const std::size_t local_count = workspace.dofs.size();
        for (std::size_t row = 0; row < local_count; ++row) {
            double value = 0.0;
            for (std::size_t column = 0; column < local_count; ++column)
                value += workspace.jacobian[row * local_count + column] * direction.at(workspace.dofs[column]);
            result.at(workspace.dofs[row]) += value;
        }
    }
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        result.at(condition.dof) = direction.at(condition.dof);
    return result;
}
} // namespace
std::vector<double> constrained_residual(const NonlinearProblem& problem, const std::vector<double>& state) {
    std::vector<double> result;
    problem.assemble_residual(state, result);
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        result.at(condition.dof) = state.at(condition.dof) - condition.value;
    return result;
}
FieldNorms field_norms(const NonlinearProblem& problem, const std::vector<double>& values) {
    if (values.size() != problem.dof_count())
        throw std::invalid_argument("Field norm vector size does not match the problem");
    problem.validate_discretization();
    FieldNorms result;
    result.l2.resize(problem.field_layout().size());
    result.maximum_absolute.resize(problem.field_layout().size());
    for (std::size_t field = 0; field < problem.field_layout().size(); ++field) {
        const FieldDescriptor& descriptor = problem.field_layout()[field];
        for (std::size_t dof = descriptor.begin; dof < descriptor.end; ++dof) {
            const double value = values[dof];
            result.l2[field] = std::hypot(result.l2[field], value);
            result.maximum_absolute[field] = std::max(result.maximum_absolute[field], std::abs(value));
        }
    }
    return result;
}
DirectionalJacobianCheck check_directional_jacobian(const NonlinearProblem& problem, const std::vector<double>& state,
    const std::vector<double>& direction, double step) {
    if (state.size() != problem.dof_count() || direction.size() != problem.dof_count())
        throw std::invalid_argument("Directional Jacobian vectors do not match the problem");
    if (!std::isfinite(step) || !(step > 0.0))
        throw std::invalid_argument("Directional Jacobian step must be finite and positive");
    std::vector<double> plus = state;
    std::vector<double> minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const std::vector<double> residual = constrained_residual(problem, state);
    const std::vector<double> plus_residual = constrained_residual(problem, plus);
    const std::vector<double> minus_residual = constrained_residual(problem, minus);
    const std::vector<double> analytic = analytic_directional_derivative(problem, state, direction);
    std::vector<double> finite_difference(problem.dof_count(), 0.0);
    std::vector<double> difference(problem.dof_count(), 0.0);
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        finite_difference[dof] = (plus_residual[dof] - minus_residual[dof]) / (2.0 * step);
        difference[dof] = analytic[dof] - finite_difference[dof];
    }
    return {field_norms(problem, residual), field_norms(problem, analytic), field_norms(problem, finite_difference),
        field_norms(problem, difference)};
}
PiecewiseLinearTimeTable::PiecewiseLinearTimeTable(
    std::string name, std::vector<double> times, std::vector<double> values)
    : _name(std::move(name)), _times(std::move(times)), _values(std::move(values)) {
    if (_name.empty()) throw std::invalid_argument("Time-table name must not be empty");
    if (_times.size() < 2 || _times.size() != _values.size())
        throw std::invalid_argument("Time table requires at least two time/value pairs: " + _name);
    for (std::size_t index = 0; index < _times.size(); ++index) {
        if (!std::isfinite(_times[index]) || _times[index] < 0.0 || !std::isfinite(_values[index]))
            throw std::invalid_argument("Time-table entries must be finite with nonnegative times: " + _name);
        if (index > 0 && !(_times[index] > _times[index - 1]))
            throw std::invalid_argument("Time-table times must be strictly increasing: " + _name);
    }
}
double PiecewiseLinearTimeTable::value(double time) const {
    if (!std::isfinite(time) || time < 0.0)
        throw std::invalid_argument("Time-table evaluation time must be finite and nonnegative");
    if (time <= _times.front()) return _values.front();
    if (time >= _times.back()) return _values.back();
    const auto upper = std::upper_bound(_times.begin(), _times.end(), time);
    const std::size_t right = static_cast<std::size_t>(upper - _times.begin()), left = right - 1;
    const double fraction = (time - _times[left]) / (_times[right] - _times[left]);
    return (1.0 - fraction) * _values[left] + fraction * _values[right];
}
} // namespace fuelsim
