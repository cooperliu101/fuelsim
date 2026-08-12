#include "assembly.hpp"
#include "fuelsim/diagnostics.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/time_table.hpp"
#include "fuelsim/transient_problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fuelsim {

namespace rz {

class TransientConservationCalculator final {
  public:
    static TransientConservationSummary
    summarize(const TransientProblem& problem, const std::vector<double>& converged_solution,
              const std::vector<std::vector<Quad4MaterialHistory>>& staged_histories,
              const std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>& staged_stresses);
};

} // namespace rz

GlobalStateView::GlobalStateView(const std::vector<double>& dense_values)
    : _global_size(dense_values.size()), _dense_values(&dense_values), _global_dofs(nullptr), _sparse_values(nullptr) {}

GlobalStateView::GlobalStateView(std::size_t global_size, const std::vector<std::uint32_t>& global_dofs,
                                 const std::vector<double>& values)
    : _global_size(global_size), _dense_values(nullptr), _global_dofs(&global_dofs), _sparse_values(&values) {
    if (global_dofs.size() != values.size())
        throw std::invalid_argument("GlobalStateView sparse index and value sizes differ");
    if (!std::is_sorted(global_dofs.begin(), global_dofs.end()) ||
        std::adjacent_find(global_dofs.begin(), global_dofs.end()) != global_dofs.end())
        throw std::invalid_argument("GlobalStateView sparse DOFs must be sorted and unique");
    if (!global_dofs.empty() && global_dofs.back() >= global_size)
        throw std::out_of_range("GlobalStateView sparse DOF exceeds the global size");
}

std::size_t GlobalStateView::global_size() const noexcept {
    return _global_size;
}

std::size_t GlobalStateView::local_size() const noexcept {
    return _dense_values != nullptr ? _dense_values->size() : _sparse_values->size();
}

bool GlobalStateView::contains(std::size_t global_dof) const {
    if (global_dof >= _global_size)
        return false;
    if (_dense_values != nullptr)
        return true;
    if (global_dof > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))
        return false;
    return std::binary_search(_global_dofs->begin(), _global_dofs->end(), static_cast<std::uint32_t>(global_dof));
}

double GlobalStateView::value(std::size_t global_dof) const {
    if (global_dof >= _global_size)
        throw std::out_of_range("GlobalStateView requested DOF exceeds the global size");
    if (_dense_values != nullptr)
        return _dense_values->at(global_dof);
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
    if (dof_count_value == 0)
        throw std::invalid_argument("Local contribution must contain at least one DOF");
    if (dof_count_value > std::numeric_limits<std::size_t>::max() / dof_count_value)
        throw std::length_error("Local contribution Jacobian size overflows");
    dofs.assign(dof_count_value, std::numeric_limits<std::size_t>::max());
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
        for (std::size_t previous = 0; previous < field; ++previous) {
            if (fields[previous].name == descriptor.name)
                throw std::invalid_argument("NonlinearProblem field names must be unique");
        }
        expected_begin = descriptor.end;
    }
    if (expected_begin != dof_count())
        throw std::invalid_argument("NonlinearProblem field ranges must cover every DOF");
    for (std::size_t contribution = 0; contribution < contribution_count(); ++contribution) {
        const std::size_t local_count = contribution_dof_count(contribution);
        if (local_count == 0 || local_count > std::numeric_limits<std::size_t>::max() / local_count)
            throw std::invalid_argument("NonlinearProblem contribution size is invalid");
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

std::vector<std::size_t> NonlinearProblem::required_state_dofs(std::size_t contribution_begin,
                                                               std::size_t contribution_end) const {
    if (contribution_begin > contribution_end || contribution_end > contribution_count())
        throw std::out_of_range("NonlinearProblem contribution range is invalid");
    std::vector<std::size_t> result;
    ContributionWorkspace workspace;
    for (std::size_t contribution = contribution_begin; contribution < contribution_end; ++contribution) {
        workspace.resize(contribution_dof_count(contribution), false);
        fill_contribution_dofs(contribution, workspace.dofs);
        if (workspace.dofs.size() != contribution_dof_count(contribution))
            throw std::logic_error("NonlinearProblem contribution DOF buffer has the wrong size");
        for (const std::size_t dof : workspace.dofs) {
            if (dof >= dof_count())
                throw std::out_of_range("NonlinearProblem contribution DOF is out of range");
            result.push_back(dof);
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

void NonlinearProblem::validate_local_state(std::size_t contribution_begin, std::size_t contribution_end,
                                            const GlobalStateView& state) const {
    if (contribution_begin > contribution_end || contribution_end > contribution_count())
        throw std::out_of_range("NonlinearProblem local contribution range is invalid");
    if (state.global_size() != dof_count())
        throw std::invalid_argument("NonlinearProblem shadow state size does not match problem");
}

void NonlinearProblem::gather_contribution_state(std::size_t contribution_index, const GlobalStateView& global_state,
                                                 ContributionWorkspace& workspace, bool include_jacobian) const {
    if (global_state.global_size() != dof_count())
        throw std::invalid_argument("NonlinearProblem shadow state size does not match DOF count");
    if (contribution_index >= contribution_count())
        throw std::out_of_range("NonlinearProblem contribution index is out of range");
    const std::size_t local_count = contribution_dof_count(contribution_index);
    workspace.resize(local_count, include_jacobian);
    fill_contribution_dofs(contribution_index, workspace.dofs);
    if (workspace.dofs.size() != local_count)
        throw std::logic_error("NonlinearProblem contribution DOF buffer has the wrong size");
    for (std::size_t local = 0; local < local_count; ++local) {
        if (workspace.dofs[local] >= dof_count())
            throw std::out_of_range("NonlinearProblem contribution DOF is out of range");
        workspace.state[local] = global_state.value(workspace.dofs[local]);
    }
}

void NonlinearProblem::evaluate_contribution_residual(std::size_t contribution_index,
                                                      const GlobalStateView& global_state,
                                                      ContributionWorkspace& workspace) const {
    gather_contribution_state(contribution_index, global_state, workspace, false);
    compute_contribution_residual(contribution_index, workspace.state, workspace.residual);
    if (workspace.residual.size() != workspace.dofs.size())
        throw std::logic_error("NonlinearProblem contribution residual has the wrong size");
}

void NonlinearProblem::evaluate_contribution_system(std::size_t contribution_index, const GlobalStateView& global_state,
                                                    ContributionWorkspace& workspace) const {
    gather_contribution_state(contribution_index, global_state, workspace, true);
    compute_contribution_system(contribution_index, workspace.state, workspace.residual, workspace.jacobian);
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
    for (std::size_t contribution = 0; contribution < contribution_count(); ++contribution) {
        evaluate_contribution_residual(contribution, global_state, workspace);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            residual[workspace.dofs[local]] += workspace.residual[local];
    }
}

namespace {

LocalValues rz_local_values(const std::vector<double>& values) {
    if (values.size() != local_dof_count)
        throw std::invalid_argument("RZ contribution state must contain 12 DOFs");
    LocalValues result{};
    std::copy(values.begin(), values.end(), result.begin());
    return result;
}

void copy_rz_residual(const LocalResidual& source, std::vector<double>& destination) {
    destination.assign(source.begin(), source.end());
}

void copy_rz_system(const LocalSystem& source, std::vector<double>& residual, std::vector<double>& jacobian) {
    residual.assign(source.residual.begin(), source.residual.end());
    jacobian.assign(source.jacobian.begin(), source.jacobian.end());
}

} // namespace

// Steady nonlinear problem.

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : _spatial(std::make_unique<rz::SpatialAssembly>(std::move(definition), source_mesh)) {
    _region_kernels.reserve(_spatial->region_count());
    for (std::size_t region_value = 0; region_value < _spatial->region_count(); ++region_value) {
        const RegionDefinition& value = _spatial->region(region_value);
        _region_kernels.emplace_back(IsotropicThermoelasticMaterial(value.material),
                                     _spatial->region_heat_source(region_value), value.strain_formulation);
    }
}

SteadyProblem::~SteadyProblem() = default;

const SpatialDefinition& SteadyProblem::definition() const noexcept {
    return _spatial->definition();
}

const DofMap& SteadyProblem::dof_map() const noexcept {
    return _spatial->dof_map();
}

std::size_t SteadyProblem::region_count() const noexcept {
    return _spatial->region_count();
}

std::size_t SteadyProblem::region_index(const std::string& name) const {
    return _spatial->region_index(name);
}

const RegionDefinition& SteadyProblem::region(std::size_t index) const {
    return _spatial->region(index);
}

const RegionMesh& SteadyProblem::region_mesh(std::size_t index) const {
    return _spatial->region_mesh(index);
}

const Quad4RzThermoelasticKernel& SteadyProblem::region_kernel(std::size_t index) const {
    return _region_kernels.at(index);
}

std::size_t SteadyProblem::region_node_offset(std::size_t index) const {
    return _spatial->region_node_offset(index);
}

std::size_t SteadyProblem::region_element_count(std::size_t index) const {
    return _spatial->region_element_count(index);
}

std::size_t SteadyProblem::region_element_offset(std::size_t index) const {
    return _spatial->region_element_offset(index);
}

std::size_t SteadyProblem::volume_contribution_count() const noexcept {
    return _spatial->volume_contribution_count();
}

SpatialContributionType SteadyProblem::contribution_type(std::size_t contribution_index) const {
    return _spatial->contribution_type(contribution_index);
}

const Quad4RzGeometry& SteadyProblem::region_element_geometry(std::size_t region_value,
                                                              std::size_t element_index) const {
    return _spatial->region_element_geometry(region_value, element_index);
}

std::size_t SteadyProblem::contact_count() const noexcept {
    return _spatial->contact_count();
}

const ContactDefinition& SteadyProblem::contact(std::size_t index) const {
    return _spatial->contact(index);
}

const std::vector<std::vector<ContactPointHistory>>& SteadyProblem::committed_contact_histories() const noexcept {
    return _spatial->committed_contact_histories();
}

void SteadyProblem::commit_contact_state(const std::vector<double>& state) {
    _spatial->commit_contact_state(state);
}

bool SteadyProblem::uses_augmented_contact() const noexcept {
    return _spatial->uses_augmented_contact();
}

AugmentedContactUpdate SteadyProblem::update_augmented_contact_multipliers(const std::vector<double>& state,
                                                                           std::size_t completed_updates) {
    return _spatial->update_augmented_contact_multipliers(state, completed_updates);
}

void SteadyProblem::restore_contact_state(const std::vector<double>& state,
                                          std::vector<std::vector<ContactPointHistory>> histories) {
    _spatial->restore_contact_state(state, std::move(histories));
}

void SteadyProblem::set_load_factor(double value) {
    _spatial->set_load_factor(value);
    refresh_region_heat_sources();
}

void SteadyProblem::refresh_region_heat_sources() {
    for (std::size_t region_value = 0; region_value < region_count(); ++region_value)
        _region_kernels[region_value].set_volumetric_heat_source(_spatial->region_heat_source(region_value));
}

double SteadyProblem::load_factor() const noexcept {
    return _spatial->load_factor();
}

void SteadyProblem::set_time(double value) {
    _spatial->set_time(value);
    for (Quad4RzThermoelasticKernel& kernel : _region_kernels)
        kernel.set_time(value);
    refresh_region_heat_sources();
}

std::vector<double> SteadyProblem::initial_state() const {
    return _spatial->initial_state();
}

std::vector<ContactNodeSummary> SteadyProblem::summarize_contact_nodes(std::size_t contact_index,
                                                                       const std::vector<double>& state) const {
    return _spatial->summarize_contact_nodes(contact_index, state);
}

std::vector<std::size_t> SteadyProblem::contact_secondary_source_nodes(std::size_t contact_index) const {
    return _spatial->contact_secondary_source_nodes(contact_index);
}

InterfaceSummary SteadyProblem::summarize_interface(std::size_t contact_index, const std::vector<double>& state) const {
    return _spatial->summarize_interface(contact_index, state);
}

namespace rz {

struct SteadyStateStorage final {
    explicit SteadyStateStorage(std::vector<std::vector<ContactPointHistory>> value)
        : contact_histories(std::move(value)) {}

    std::vector<std::vector<ContactPointHistory>> contact_histories;
};

} // namespace rz

struct SteadyStateSnapshot::Storage final {
    Storage(std::shared_ptr<const void> owner_value, std::shared_ptr<const rz::SteadyStateStorage> value)
        : owner(std::move(owner_value)), rz_state(std::move(value)) {}

    std::shared_ptr<const void> owner;
    std::shared_ptr<const rz::SteadyStateStorage> rz_state;
};

SteadyStateSnapshot::SteadyStateSnapshot() = default;
SteadyStateSnapshot::~SteadyStateSnapshot() = default;
SteadyStateSnapshot::SteadyStateSnapshot(const SteadyStateSnapshot& other) = default;
SteadyStateSnapshot& SteadyStateSnapshot::operator=(const SteadyStateSnapshot& other) = default;
SteadyStateSnapshot::SteadyStateSnapshot(SteadyStateSnapshot&& other) noexcept = default;
SteadyStateSnapshot& SteadyStateSnapshot::operator=(SteadyStateSnapshot&& other) noexcept = default;

SteadyStateSnapshot::SteadyStateSnapshot(std::shared_ptr<const Storage> storage) : _storage(std::move(storage)) {}

bool SteadyStateSnapshot::empty() const noexcept {
    return _storage == nullptr;
}

std::shared_ptr<const void> SteadyStateSnapshot::snapshot_owner() const noexcept {
    return _storage != nullptr ? _storage->owner : nullptr;
}

SteadyStateSnapshot SteadyProblem::capture_internal_state() const {
    const auto state = std::make_shared<rz::SteadyStateStorage>(_spatial->committed_contact_histories());
    return SteadyStateSnapshot(
        std::make_shared<SteadyStateSnapshot::Storage>(discretization_identity(), std::move(state)));
}

void SteadyProblem::restore_internal_state(const SteadyStateSnapshot& snapshot, const std::vector<double>& state) {
    if (snapshot.empty())
        throw std::invalid_argument("SteadyProblem cannot restore an empty internal-state snapshot");
    if (snapshot.snapshot_owner() != discretization_identity())
        throw std::invalid_argument("SteadyProblem cannot restore a snapshot from another problem");
    _spatial->restore_contact_state(state, snapshot._storage->rz_state->contact_histories);
}

void SteadyProblem::commit_internal_state(const std::vector<double>& state) {
    _spatial->commit_contact_state(state);
}

std::size_t SteadyProblem::dof_count() const noexcept {
    return _spatial->dof_count();
}

std::size_t SteadyProblem::contribution_count() const noexcept {
    return _spatial->contribution_count();
}

const std::vector<FieldDescriptor>& SteadyProblem::field_layout() const noexcept {
    return _spatial->dof_map().field_layout();
}

const std::vector<DirichletCondition>& SteadyProblem::dirichlet_conditions() const noexcept {
    return _spatial->dirichlet_conditions();
}

void SteadyProblem::validate_state(const std::vector<double>& state) const {
    _spatial->validate_state(state);
}

std::vector<std::size_t> SteadyProblem::required_state_dofs(std::size_t contribution_begin,
                                                            std::size_t contribution_end) const {
    return _spatial->required_state_dofs(contribution_begin, contribution_end);
}

void SteadyProblem::validate_local_state(std::size_t contribution_begin, std::size_t contribution_end,
                                         const GlobalStateView& state) const {
    _spatial->validate_local_state(contribution_begin, contribution_end, state);
}

std::size_t SteadyProblem::contribution_dof_count(std::size_t contribution_index) const {
    if (contribution_index >= contribution_count())
        throw std::out_of_range("SteadyProblem contribution index is out of range");
    return local_dof_count;
}

void SteadyProblem::fill_contribution_dofs(std::size_t contribution_index, std::vector<std::size_t>& dofs) const {
    const LocalDofs fixed = contribution_dofs(contribution_index);
    dofs.assign(fixed.begin(), fixed.end());
}

void SteadyProblem::compute_contribution_residual(std::size_t contribution_index, const std::vector<double>& state,
                                                  std::vector<double>& residual) const {
    copy_rz_residual(contribution_residual(contribution_index, rz_local_values(state)), residual);
}

void SteadyProblem::compute_contribution_system(std::size_t contribution_index, const std::vector<double>& state,
                                                std::vector<double>& residual, std::vector<double>& jacobian) const {
    copy_rz_system(linearize_contribution(contribution_index, rz_local_values(state)), residual, jacobian);
}

LocalDofs SteadyProblem::contribution_dofs(std::size_t contribution_index) const {
    return _spatial->contribution_dofs(contribution_index);
}

LocalValues SteadyProblem::contribution_state(std::size_t contribution_index,
                                              const std::vector<double>& global_state) const {
    return contribution_state(contribution_index, GlobalStateView(global_state));
}

LocalValues SteadyProblem::contribution_state(std::size_t contribution_index,
                                              const GlobalStateView& global_state) const {
    if (global_state.global_size() != dof_count())
        throw std::invalid_argument("SteadyProblem contribution state has the wrong global size");
    const LocalDofs dofs = contribution_dofs(contribution_index);
    LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local)
        result[local] = global_state.value(dofs[local]);
    return result;
}

LocalResidual SteadyProblem::contribution_residual(std::size_t contribution_index, const LocalValues& state) const {
    if (contribution_index < volume_contribution_count()) {
        const auto location = _spatial->element_location(contribution_index);
        return _region_kernels[location.first].residual(region_element_geometry(location.first, location.second),
                                                        state);
    }
    return _spatial->contribution_residual(contribution_index, state);
}

LocalSystem SteadyProblem::linearize_contribution(std::size_t contribution_index, const LocalValues& state) const {
    if (contribution_index < volume_contribution_count()) {
        const auto location = _spatial->element_location(contribution_index);
        return _region_kernels[location.first].linearize(region_element_geometry(location.first, location.second),
                                                         state);
    }
    return _spatial->linearize_contribution(contribution_index, state);
}

// Transient conservation diagnostics.
namespace rz {
namespace {

double stress_strain_inner_product(const AxisymmetricStressValues& stress,
                                   const std::array<double, 4>& strain) noexcept {
    return stress.rr * strain[0] + stress.zz * strain[1] + stress.hoop * strain[2] + 2.0 * stress.rz * strain[3];
}

std::array<double, 4> strain_difference(const std::array<double, 4>& current, const std::array<double, 4>& old) {
    std::array<double, 4> result{};
    for (std::size_t component = 0; component < result.size(); ++component)
        result[component] = current[component] - old[component];
    return result;
}

} // namespace

TransientConservationSummary TransientConservationCalculator::summarize(
    const TransientProblem& problem, const std::vector<double>& converged_solution,
    const std::vector<std::vector<Quad4MaterialHistory>>& staged_histories,
    const std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>& staged_stresses) {
    problem.require_active_time_step();
    TransientConservationSummary result;
    std::vector<double> raw_residual(problem.dof_count(), 0.0);

    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        const LocalDofs dofs = problem.contribution_dofs(contribution);
        const LocalValues state = problem.contribution_state(contribution, converged_solution);
        const LocalResidual residual = problem.contribution_residual(contribution, state);
        const SpatialContributionType type = problem._spatial->contribution_type(contribution);
        for (std::size_t local = 0; local < local_dof_count; ++local) {
            raw_residual[dofs[local]] += residual[local];
            const double increment = converged_solution[dofs[local]] - problem._committed_solution[dofs[local]];
            if (local < 4) {
                if (type == SpatialContributionType::thermal_contact)
                    result.interface_heat_imbalance += residual[local];
                else if (type == SpatialContributionType::convection)
                    result.convection_heat_rate += residual[local];
                continue;
            }
            const double work = residual[local] * increment;
            if (type == SpatialContributionType::volume)
                result.internal_mechanical_work_increment += work;
            else if (type == SpatialContributionType::mechanical_contact)
                result.contact_work_increment += work;
            else if (type == SpatialContributionType::pressure || type == SpatialContributionType::traction)
                result.pressure_traction_work_increment -= work;
        }
    }
    for (std::size_t region_value = 0; region_value < problem.region_count(); ++region_value) {
        const Quad4RzTransientKernel& kernel = problem._region_kernels[region_value];
        const std::size_t offset = problem._spatial->region_element_offset(region_value);
        for (std::size_t element = 0; element < staged_histories[region_value].size(); ++element) {
            const LocalValues current = problem.contribution_state(offset + element, converged_solution);
            const LocalValues old = problem.contribution_state(offset + element, problem._committed_solution);
            const Quad4RzGeometry& geometry = problem.region_element_geometry(region_value, element);
            for (std::size_t q = 0; q < geometry.points.size(); ++q) {
                const RzQuadraturePoint& point = geometry.points[q];
                double current_temperature = 0.0;
                double old_temperature = 0.0;
                for (std::size_t node = 0; node < quad4_node_count; ++node) {
                    current_temperature += point.shape[node] * current[node];
                    old_temperature += point.shape[node] * old[node];
                }
                const double heat_capacity =
                    kernel.heat_capacity(current_temperature, point.radius, point.axial_coordinate);
                result.stored_heat_rate += point.weighted_measure * heat_capacity *
                                           (current_temperature - old_temperature) / problem._active_time_step;
                result.generated_heat_rate += point.weighted_measure * kernel.volumetric_heat_source();

                const MaterialPointState& old_history = problem._material_histories[region_value][element][q];
                const MaterialPointState& new_history = staged_histories[region_value][element][q];
                const AxisymmetricStressValues& old_stress = problem._material_stresses[region_value][element][q];
                const AxisymmetricStressValues& new_stress = staged_stresses[region_value][element][q];
                result.elastic_energy_change += 0.5 * point.weighted_measure *
                                                (stress_strain_inner_product(new_stress, new_history.elastic_strain) -
                                                 stress_strain_inner_product(old_stress, old_history.elastic_strain));
                result.plastic_dissipation_increment +=
                    point.weighted_measure *
                    stress_strain_inner_product(
                        new_stress, strain_difference(new_history.plastic_strain, old_history.plastic_strain));
                result.creep_dissipation_increment +=
                    point.weighted_measure *
                    stress_strain_inner_product(new_stress,
                                                strain_difference(new_history.creep_strain, old_history.creep_strain));
            }
        }
    }

    std::vector<bool> constrained(problem.dof_count(), false);
    for (const DirichletCondition& condition : problem.dirichlet_conditions()) {
        constrained[condition.dof] = true;
        const double increment = converged_solution[condition.dof] - problem._committed_solution[condition.dof];
        if (condition.dof < problem.dof_map().node_count())
            result.dirichlet_heat_input_rate += raw_residual[condition.dof];
        else
            result.dirichlet_reaction_work_increment += raw_residual[condition.dof] * increment;
    }
    double thermal_residual_squared = 0.0;
    double mechanical_residual_squared = 0.0;
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        if (constrained[dof])
            continue;
        const double value = raw_residual[dof];
        if (dof < problem.dof_map().node_count())
            thermal_residual_squared += value * value;
        else
            mechanical_residual_squared += value * value;
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
    return result;
}

} // namespace rz

// Transient state transactions and nonlinear problem.
namespace {

bool finite_stress(const AxisymmetricStressValues& stress) {
    return std::isfinite(stress.rr) && std::isfinite(stress.zz) && std::isfinite(stress.hoop) &&
           std::isfinite(stress.rz);
}

bool valid_material_state(const MaterialPointState& state) {
    for (std::size_t component = 0; component < 4; ++component) {
        if (!std::isfinite(state.elastic_strain[component]) || !std::isfinite(state.plastic_strain[component]) ||
            !std::isfinite(state.creep_strain[component]))
            return false;
    }
    return std::isfinite(state.equivalent_plastic_strain) && state.equivalent_plastic_strain >= 0.0 &&
           std::isfinite(state.equivalent_creep_strain) && state.equivalent_creep_strain >= 0.0;
}

RegionInelasticSummary summarize_history(const std::vector<Quad4MaterialHistory>& history) noexcept {
    RegionInelasticSummary summary{0.0, 0.0};
    for (const Quad4MaterialHistory& element : history) {
        for (const MaterialPointState& point : element) {
            summary.maximum_equivalent_plastic_strain =
                std::max(summary.maximum_equivalent_plastic_strain, point.equivalent_plastic_strain);
            summary.maximum_equivalent_creep_strain =
                std::max(summary.maximum_equivalent_creep_strain, point.equivalent_creep_strain);
        }
    }
    return summary;
}

void validate_definition(const TransientProblemDefinition& definition) {
    if (definition.regions.size() != definition.spatial.regions.size())
        throw std::invalid_argument("TransientProblem requires one transient material per region");
    for (std::size_t region = 0; region < definition.regions.size(); ++region) {
        if (definition.regions[region].region != definition.spatial.regions[region].name)
            throw std::invalid_argument("TransientProblem region material order does not match the "
                                        "spatial regions");
    }
}

} // namespace

TransientProblem::TransientProblem(TransientProblemDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : _definition(std::move(definition)),
      _spatial(std::make_unique<rz::SpatialAssembly>(_definition.spatial, source_mesh)), _committed_time(0.0),
      _committed_load_factor(0.0), _active_time_step(0.0), _active_end_time(0.0), _active_load_factor(0.0),
      _time_step_active(false) {
    validate_definition(_definition);
    _region_kernels.reserve(region_count());
    _material_histories.resize(region_count());
    _material_stresses.resize(region_count());
    for (std::size_t region_value = 0; region_value < region_count(); ++region_value) {
        _region_kernels.emplace_back(IsotropicInelasticMaterial(_definition.spatial.regions[region_value].material,
                                                                _definition.regions[region_value].material),
                                     0.0, _definition.spatial.regions[region_value].strain_formulation);
        _material_histories[region_value].resize(_spatial->region_element_count(region_value));
        _material_stresses[region_value].resize(_spatial->region_element_count(region_value));
    }
    apply_spatial_controls(0.0, 0.0);
    _committed_solution = _spatial->initial_state();
    _spatial->restore_contact_state(_committed_solution, _spatial->committed_contact_histories());
}

TransientProblem::~TransientProblem() = default;

const TransientProblemDefinition& TransientProblem::definition() const noexcept {
    return _definition;
}

const DofMap& TransientProblem::dof_map() const noexcept {
    return _spatial->dof_map();
}

std::size_t TransientProblem::region_count() const noexcept {
    return _definition.spatial.regions.size();
}

std::size_t TransientProblem::region_index(const std::string& name) const {
    return _spatial->region_index(name);
}

std::size_t TransientProblem::region_node_offset(std::size_t region_value) const {
    return _spatial->region_node_offset(region_value);
}

const RegionDefinition& TransientProblem::region(std::size_t index) const {
    return _spatial->region(index);
}

const RegionMesh& TransientProblem::region_mesh(std::size_t index) const {
    return _spatial->region_mesh(index);
}

const Quad4RzTransientKernel& TransientProblem::region_kernel(std::size_t index) const {
    return _region_kernels.at(index);
}

const Quad4RzGeometry& TransientProblem::region_element_geometry(std::size_t region_value,
                                                                 std::size_t element_index) const {
    return _spatial->region_element_geometry(region_value, element_index);
}

const std::vector<double>& TransientProblem::committed_solution() const noexcept {
    return _committed_solution;
}

double TransientProblem::committed_time() const noexcept {
    return _committed_time;
}

double TransientProblem::committed_load_factor() const noexcept {
    return _committed_load_factor;
}

bool TransientProblem::time_step_active() const noexcept {
    return _time_step_active;
}

std::vector<double> TransientProblem::time_events() const {
    std::vector<double> result;
    for (const PiecewiseLinearTimeTable& table : _definition.spatial.time_tables)
        result.insert(result.end(), table.times().begin(), table.times().end());
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

namespace rz {

struct TransientStateStorage final {
    explicit TransientStateStorage(TransientCommittedState value) : state(std::move(value)) {}

    TransientCommittedState state;
};

} // namespace rz

struct TransientStateSnapshot::Storage final {
    Storage(std::shared_ptr<const void> owner_value, std::shared_ptr<const rz::TransientStateStorage> value)
        : owner(std::move(owner_value)), rz_state(std::move(value)) {}

    std::shared_ptr<const void> owner;
    std::shared_ptr<const rz::TransientStateStorage> rz_state;
};

TransientStateSnapshot::TransientStateSnapshot() = default;
TransientStateSnapshot::~TransientStateSnapshot() = default;
TransientStateSnapshot::TransientStateSnapshot(const TransientStateSnapshot& other) = default;
TransientStateSnapshot& TransientStateSnapshot::operator=(const TransientStateSnapshot& other) = default;
TransientStateSnapshot::TransientStateSnapshot(TransientStateSnapshot&& other) noexcept = default;
TransientStateSnapshot& TransientStateSnapshot::operator=(TransientStateSnapshot&& other) noexcept = default;

TransientStateSnapshot::TransientStateSnapshot(std::shared_ptr<const Storage> storage) : _storage(std::move(storage)) {}

bool TransientStateSnapshot::empty() const noexcept {
    return _storage == nullptr;
}

std::shared_ptr<const void> TransientStateSnapshot::snapshot_owner() const noexcept {
    return _storage != nullptr ? _storage->owner : nullptr;
}

TransientCommittedState TransientProblem::committed_state() const {
    return {
        _committed_solution,        _material_histories, _material_stresses,    _spatial->committed_contact_histories(),
        _last_conservation_summary, _committed_time,     _committed_load_factor};
}

TransientStateSnapshot TransientProblem::capture_state() const {
    if (_time_step_active)
        throw std::logic_error("TransientProblem cannot capture an active time step");
    const auto state = std::make_shared<rz::TransientStateStorage>(committed_state());
    return TransientStateSnapshot(
        std::make_shared<TransientStateSnapshot::Storage>(discretization_identity(), std::move(state)));
}

void TransientProblem::restore_state(const TransientStateSnapshot& snapshot) {
    if (snapshot.empty())
        throw std::invalid_argument("TransientProblem cannot restore an empty state snapshot");
    if (snapshot.snapshot_owner() != discretization_identity())
        throw std::invalid_argument("TransientProblem cannot restore a snapshot from another problem");
    restore_committed_state(snapshot._storage->rz_state->state);
}

void TransientProblem::restore_committed_state(TransientCommittedState state) {
    if (_time_step_active)
        throw std::logic_error("TransientProblem cannot restore during an active time step");
    if (state.solution.size() != dof_count() || state.material_histories.size() != region_count() ||
        state.material_stresses.size() != region_count() ||
        state.contact_histories.size() != _definition.spatial.contacts.size())
        throw std::invalid_argument("Transient committed state layout does not match the problem");
    if (!std::isfinite(state.time) || state.time < 0.0 || !std::isfinite(state.load_factor) || state.load_factor < 0.0)
        throw std::invalid_argument("Transient committed time and load factor must be valid");
    for (std::size_t node = 0; node < dof_map().node_count(); ++node) {
        const double temperature = state.solution.at(dof_map().temperature(node));
        if (!std::isfinite(temperature) || !(temperature > 0.0) ||
            !std::isfinite(state.solution.at(dof_map().radial_displacement(node))) ||
            !std::isfinite(state.solution.at(dof_map().axial_displacement(node))))
            throw std::invalid_argument("Transient committed nodal state must be finite with "
                                        "positive temperatures");
    }
    for (std::size_t region_value = 0; region_value < region_count(); ++region_value) {
        const std::size_t elements = region_mesh(region_value).elements().size();
        if (state.material_histories[region_value].size() != elements ||
            state.material_stresses[region_value].size() != elements)
            throw std::invalid_argument("Transient committed element state layout does not match");
        for (std::size_t element = 0; element < elements; ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                if (!valid_material_state(state.material_histories[region_value][element][q]) ||
                    !finite_stress(state.material_stresses[region_value][element][q]))
                    throw std::invalid_argument("Transient committed integration-point state is "
                                                "invalid");
            }
        }
    }

    _spatial->restore_contact_state(state.solution, state.contact_histories);
    _committed_solution = std::move(state.solution);
    _material_histories = std::move(state.material_histories);
    _material_stresses = std::move(state.material_stresses);
    _last_conservation_summary = state.conservation;
    _committed_time = state.time;
    _committed_load_factor = state.load_factor;
    clear_active_time_step();
    apply_spatial_controls(_committed_time, _committed_load_factor);
}

namespace rz {
namespace {

struct TimeErrorAccumulator final {
    double difference_squared = 0.0;
    double solution_squared = 0.0;
    std::size_t count = 0;
};

void accumulate_time_error(TimeErrorAccumulator& accumulator, double full_step, double two_half_steps) {
    const double difference = two_half_steps - full_step;
    accumulator.difference_squared += difference * difference;
    accumulator.solution_squared += two_half_steps * two_half_steps;
    ++accumulator.count;
}

double normalized_time_error(const TimeErrorAccumulator& accumulator, double absolute_tolerance,
                             double relative_tolerance) {
    if (accumulator.count == 0)
        return 0.0;
    const double denominator = absolute_tolerance * std::sqrt(static_cast<double>(accumulator.count)) +
                               relative_tolerance * std::sqrt(accumulator.solution_squared);
    return std::sqrt(accumulator.difference_squared) / denominator;
}

TransientConservationSummary combine_rz_half_step_conservation(const TransientConservationSummary& first,
                                                               const TransientConservationSummary& second) {
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
                                                        const TransientCommittedState& two_half_steps,
                                                        const std::vector<FieldDescriptor>& fields,
                                                        std::size_t expected_dof_count,
                                                        const TransientTimeOptions& options) {
    if (full_step.solution.size() != two_half_steps.solution.size() || full_step.solution.size() != expected_dof_count)
        throw std::logic_error("step-doubling nodal-state layouts differ");
    if (full_step.material_histories.size() != two_half_steps.material_histories.size() ||
        full_step.material_stresses.size() != two_half_steps.material_stresses.size())
        throw std::logic_error("step-doubling material-state region layouts differ");
    if (fields.size() != 3 || fields[0].name != "temperature" || fields[1].name != "radial" ||
        fields[2].name != "axial")
        throw std::logic_error("RZ step-doubling field layout is invalid");

    std::vector<TimeErrorAccumulator> nodal(fields.size());
    for (std::size_t field = 0; field < fields.size(); ++field) {
        const FieldDescriptor& descriptor = fields[field];
        for (std::size_t dof = descriptor.begin; dof < descriptor.end; ++dof)
            accumulate_time_error(nodal[field], full_step.solution[dof], two_half_steps.solution[dof]);
    }

    TimeErrorAccumulator elastic;
    TimeErrorAccumulator plastic;
    TimeErrorAccumulator creep;
    TimeErrorAccumulator equivalent_plastic;
    TimeErrorAccumulator equivalent_creep;
    TimeErrorAccumulator stress;
    TimeErrorAccumulator contact_friction;
    TimeErrorAccumulator contact_normal_multiplier;
    bool contact_state_mismatch = false;
    for (std::size_t region = 0; region < full_step.material_histories.size(); ++region) {
        const auto& full_history = full_step.material_histories[region];
        const auto& half_history = two_half_steps.material_histories[region];
        const auto& full_stress = full_step.material_stresses[region];
        const auto& half_stress = two_half_steps.material_stresses[region];
        if (full_history.size() != half_history.size() || full_stress.size() != half_stress.size() ||
            full_history.size() != full_stress.size())
            throw std::logic_error("step-doubling material-state element layouts differ");
        for (std::size_t element = 0; element < full_history.size(); ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const MaterialPointState& full_point = full_history[element][q];
                const MaterialPointState& half_point = half_history[element][q];
                for (std::size_t component = 0; component < 4; ++component) {
                    accumulate_time_error(elastic, full_point.elastic_strain[component],
                                          half_point.elastic_strain[component]);
                    accumulate_time_error(plastic, full_point.plastic_strain[component],
                                          half_point.plastic_strain[component]);
                    accumulate_time_error(creep, full_point.creep_strain[component],
                                          half_point.creep_strain[component]);
                }
                accumulate_time_error(equivalent_plastic, full_point.equivalent_plastic_strain,
                                      half_point.equivalent_plastic_strain);
                accumulate_time_error(equivalent_creep, full_point.equivalent_creep_strain,
                                      half_point.equivalent_creep_strain);

                const AxisymmetricStressValues& full_value = full_stress[element][q];
                const AxisymmetricStressValues& half_value = half_stress[element][q];
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
            const ContactPointHistory& full = full_step.contact_histories[contact][node];
            const ContactPointHistory& half = two_half_steps.contact_histories[contact][node];
            accumulate_time_error(contact_friction, full.elastic_tangential_slip, half.elastic_tangential_slip);
            accumulate_time_error(contact_normal_multiplier, full.normal_multiplier, half.normal_multiplier);
            contact_state_mismatch = contact_state_mismatch || full.sliding != half.sliding;
        }
    }

    TransientTimeErrorEstimate result;
    result.nodal_fields = {
        {"temperature", normalized_time_error(nodal[0], options.temperature_time_absolute_tolerance,
                                              options.time_error_relative_tolerance)},
        {"radial_displacement", normalized_time_error(nodal[1], options.displacement_time_absolute_tolerance,
                                                      options.time_error_relative_tolerance)},
        {"axial_displacement", normalized_time_error(nodal[2], options.displacement_time_absolute_tolerance,
                                                     options.time_error_relative_tolerance)},
    };
    result.elastic_strain = normalized_time_error(elastic, options.strain_history_time_absolute_tolerance,
                                                  options.time_error_relative_tolerance);
    result.plastic_strain = normalized_time_error(plastic, options.strain_history_time_absolute_tolerance,
                                                  options.time_error_relative_tolerance);
    result.creep_strain = normalized_time_error(creep, options.strain_history_time_absolute_tolerance,
                                                options.time_error_relative_tolerance);
    result.equivalent_plastic_strain = normalized_time_error(
        equivalent_plastic, options.strain_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.equivalent_creep_strain = normalized_time_error(
        equivalent_creep, options.strain_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.stress = normalized_time_error(stress, options.stress_history_time_absolute_tolerance,
                                          options.time_error_relative_tolerance);
    result.contact_friction =
        contact_state_mismatch ? std::numeric_limits<double>::infinity()
                               : normalized_time_error(contact_friction, options.displacement_time_absolute_tolerance,
                                                       options.time_error_relative_tolerance);
    result.contact_normal_multiplier =
        normalized_time_error(contact_normal_multiplier, options.stress_history_time_absolute_tolerance,
                              options.time_error_relative_tolerance);
    result.maximum = std::max({result.elastic_strain, result.plastic_strain, result.creep_strain,
                               result.equivalent_plastic_strain, result.equivalent_creep_strain, result.stress,
                               result.contact_friction, result.contact_normal_multiplier});
    for (const TransientFieldTimeError& field : result.nodal_fields)
        result.maximum = std::max(result.maximum, field.value);
    return result;
}

} // namespace
} // namespace rz

TransientTimeErrorEstimate TransientProblem::step_doubling_error(const TransientStateSnapshot& full_snapshot,
                                                                 const TransientStateSnapshot& half_snapshot,
                                                                 const TransientTimeOptions& options) const {
    if (full_snapshot.empty() || half_snapshot.empty())
        throw std::invalid_argument("step-doubling requires two complete state snapshots");
    if (full_snapshot.snapshot_owner() != discretization_identity() ||
        half_snapshot.snapshot_owner() != discretization_identity())
        throw std::invalid_argument("step-doubling snapshots belong to another problem");
    return rz::compare_step_doubling_states(full_snapshot._storage->rz_state->state,
                                            half_snapshot._storage->rz_state->state, field_layout(), dof_count(),
                                            options);
}

void TransientProblem::combine_last_half_step_conservation(const TransientConservationSummary& first_half) {
    if (_time_step_active)
        throw std::logic_error("TransientProblem cannot combine conservation during an active time step");
    _last_conservation_summary = rz::combine_rz_half_step_conservation(first_half, _last_conservation_summary);
}

void TransientProblem::begin_time_step(const TransientStepInput& input) {
    if (_time_step_active)
        throw std::logic_error("TransientProblem already has an active time step");
    if (!std::isfinite(input.end_time) || input.end_time <= _committed_time)
        throw std::invalid_argument("TransientProblem end time must exceed committed time");
    if (!std::isfinite(input.load_factor) || input.load_factor < 0.0)
        throw std::invalid_argument("TransientProblem load factor must be finite and nonnegative");
    _active_time_step = input.end_time - _committed_time;
    _active_end_time = input.end_time;
    _active_load_factor = input.load_factor;
    if (_spatial->uses_augmented_contact())
        _active_contact_histories = _spatial->committed_contact_histories();
    try {
        apply_spatial_controls(input.end_time, input.load_factor);
    } catch (...) {
        if (!_active_contact_histories.empty())
            _spatial->restore_contact_state(_committed_solution, std::move(_active_contact_histories));
        apply_spatial_controls(_committed_time, _committed_load_factor);
        clear_active_time_step();
        throw;
    }
    _time_step_active = true;
}

void TransientProblem::commit_time_step(const std::vector<double>& converged_solution) {
    require_active_time_step();
    if (converged_solution.size() != dof_count())
        throw std::invalid_argument("TransientProblem committed solution size mismatch");
    if (!std::all_of(converged_solution.begin(), converged_solution.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::domain_error("TransientProblem committed solution must be finite");
    for (std::size_t node = 0; node < dof_map().node_count(); ++node) {
        if (!(converged_solution[dof_map().temperature(node)] > 0.0))
            throw std::domain_error("TransientProblem committed temperatures must be positive");
    }

    std::vector<std::vector<Quad4MaterialHistory>> staged(region_count());
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> staged_stresses(region_count());
    for (std::size_t region_value = 0; region_value < region_count(); ++region_value) {
        staged[region_value].resize(_spatial->region_element_count(region_value));
        staged_stresses[region_value].resize(_spatial->region_element_count(region_value));
        const std::size_t offset = _spatial->region_element_offset(region_value);
        for (std::size_t element = 0; element < staged[region_value].size(); ++element) {
            const LocalValues state = contribution_state(offset + element, converged_solution);
            const LocalValues committed_state = contribution_state(offset + element, _committed_solution);
            staged[region_value][element] = _region_kernels[region_value].trial_state_values(
                region_element_geometry(region_value, element), state, committed_state,
                _material_histories[region_value][element], _active_time_step);
            staged_stresses[region_value][element] = _region_kernels[region_value].stress_values(
                region_element_geometry(region_value, element), state, committed_state,
                _material_histories[region_value][element], _active_time_step);
        }
    }

    const TransientConservationSummary conservation =
        rz::TransientConservationCalculator::summarize(*this, converged_solution, staged, staged_stresses);
    _spatial->commit_contact_state(converged_solution);
    _material_histories.swap(staged);
    _material_stresses.swap(staged_stresses);
    _last_conservation_summary = conservation;
    _committed_solution = converged_solution;
    _committed_time = _active_end_time;
    _committed_load_factor = _active_load_factor;
    clear_active_time_step();
}

void TransientProblem::rollback_time_step() noexcept {
    if (!_time_step_active)
        return;
    apply_spatial_controls(_committed_time, _committed_load_factor);
    if (!_active_contact_histories.empty())
        _spatial->restore_contact_state(_committed_solution, std::move(_active_contact_histories));
    clear_active_time_step();
}

void TransientProblem::apply_spatial_controls(double time, double load_factor) {
    _spatial->set_time(time);
    _spatial->set_load_factor(load_factor);
    for (Quad4RzTransientKernel& kernel : _region_kernels)
        kernel.set_time(time);
    refresh_region_heat_sources();
}

void TransientProblem::clear_active_time_step() noexcept {
    _active_time_step = 0.0;
    _active_end_time = _committed_time;
    _active_load_factor = _committed_load_factor;
    _active_contact_histories.clear();
    _time_step_active = false;
}

void TransientProblem::refresh_region_heat_sources() {
    for (std::size_t region_value = 0; region_value < region_count(); ++region_value)
        _region_kernels[region_value].set_volumetric_heat_source(_spatial->region_heat_source(region_value));
}

bool TransientProblem::uses_augmented_contact() const noexcept {
    return _spatial->uses_augmented_contact();
}

AugmentedContactUpdate TransientProblem::update_augmented_contact_multipliers(const std::vector<double>& state,
                                                                              std::size_t completed_updates) {
    require_active_time_step();
    return _spatial->update_augmented_contact_multipliers(state, completed_updates);
}

const Quad4MaterialHistory& TransientProblem::material_history(std::size_t region_value,
                                                               std::size_t element_index) const {
    return _material_histories.at(region_value).at(element_index);
}

const std::array<AxisymmetricStressValues, 4>& TransientProblem::material_stress(std::size_t region_value,
                                                                                 std::size_t element_index) const {
    return _material_stresses.at(region_value).at(element_index);
}

RegionInelasticSummary TransientProblem::summarize_region_history(std::size_t region_value) const {
    return summarize_history(_material_histories.at(region_value));
}

const TransientConservationSummary& TransientProblem::last_conservation_summary() const noexcept {
    return _last_conservation_summary;
}

InterfaceSummary TransientProblem::summarize_interface(std::size_t contact_index,
                                                       const std::vector<double>& state) const {
    return _spatial->summarize_interface(contact_index, state);
}

std::vector<ContactNodeSummary> TransientProblem::summarize_contact_nodes(std::size_t contact_index,
                                                                          const std::vector<double>& state) const {
    return _spatial->summarize_contact_nodes(contact_index, state);
}

std::vector<std::size_t> TransientProblem::contact_secondary_source_nodes(std::size_t contact_index) const {
    return _spatial->contact_secondary_source_nodes(contact_index);
}

std::size_t TransientProblem::dof_count() const noexcept {
    return _spatial->dof_count();
}

std::size_t TransientProblem::contribution_count() const noexcept {
    return _spatial->contribution_count();
}

const std::vector<FieldDescriptor>& TransientProblem::field_layout() const noexcept {
    return _spatial->dof_map().field_layout();
}

const std::vector<DirichletCondition>& TransientProblem::dirichlet_conditions() const noexcept {
    return _spatial->dirichlet_conditions();
}

void TransientProblem::validate_state(const std::vector<double>& state) const {
    require_active_time_step();
    _spatial->validate_state(state);
}

std::vector<std::size_t> TransientProblem::required_state_dofs(std::size_t contribution_begin,
                                                               std::size_t contribution_end) const {
    return _spatial->required_state_dofs(contribution_begin, contribution_end);
}

void TransientProblem::validate_local_state(std::size_t contribution_begin, std::size_t contribution_end,
                                            const GlobalStateView& state) const {
    require_active_time_step();
    _spatial->validate_local_state(contribution_begin, contribution_end, state);
}

std::size_t TransientProblem::contribution_dof_count(std::size_t contribution_index) const {
    if (contribution_index >= contribution_count())
        throw std::out_of_range("TransientProblem contribution index is out of range");
    return local_dof_count;
}

void TransientProblem::fill_contribution_dofs(std::size_t contribution_index, std::vector<std::size_t>& dofs) const {
    const LocalDofs fixed = contribution_dofs(contribution_index);
    dofs.assign(fixed.begin(), fixed.end());
}

void TransientProblem::compute_contribution_residual(std::size_t contribution_index, const std::vector<double>& state,
                                                     std::vector<double>& residual) const {
    copy_rz_residual(contribution_residual(contribution_index, rz_local_values(state)), residual);
}

void TransientProblem::compute_contribution_system(std::size_t contribution_index, const std::vector<double>& state,
                                                   std::vector<double>& residual, std::vector<double>& jacobian) const {
    copy_rz_system(linearize_contribution(contribution_index, rz_local_values(state)), residual, jacobian);
}

LocalDofs TransientProblem::contribution_dofs(std::size_t contribution_index) const {
    return _spatial->contribution_dofs(contribution_index);
}

LocalValues TransientProblem::contribution_state(std::size_t contribution_index,
                                                 const std::vector<double>& global_state) const {
    return contribution_state(contribution_index, GlobalStateView(global_state));
}

LocalValues TransientProblem::contribution_state(std::size_t contribution_index,
                                                 const GlobalStateView& global_state) const {
    if (global_state.global_size() != dof_count())
        throw std::invalid_argument("TransientProblem contribution state has the wrong global size");
    const LocalDofs dofs = contribution_dofs(contribution_index);
    LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local)
        result[local] = global_state.value(dofs[local]);
    return result;
}

LocalResidual TransientProblem::contribution_residual(std::size_t contribution_index, const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _spatial->volume_contribution_count()) {
        const auto location = _spatial->element_location(contribution_index);
        return _region_kernels[location.first].residual(region_element_geometry(location.first, location.second), state,
                                                        contribution_state(contribution_index, _committed_solution),
                                                        _material_histories[location.first][location.second],
                                                        _active_time_step);
    }
    return _spatial->contribution_residual(contribution_index, state);
}

LocalSystem TransientProblem::linearize_contribution(std::size_t contribution_index, const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _spatial->volume_contribution_count()) {
        const auto location = _spatial->element_location(contribution_index);
        return _region_kernels[location.first].linearize(
            region_element_geometry(location.first, location.second), state,
            contribution_state(contribution_index, _committed_solution),
            _material_histories[location.first][location.second], _active_time_step);
    }
    return _spatial->linearize_contribution(contribution_index, state);
}

void TransientProblem::require_active_time_step() const {
    if (!_time_step_active)
        throw std::logic_error("TransientProblem residual evaluation requires "
                               "an active time step");
}

// Directional Jacobian diagnostics.
namespace {

std::vector<double> analytic_directional_derivative(const NonlinearProblem& problem, const std::vector<double>& state,
                                                    const std::vector<double>& direction) {
    problem.validate_discretization();
    problem.validate_state(state);
    std::vector<double> result(problem.dof_count(), 0.0);
    const GlobalStateView global_state(state);
    ContributionWorkspace workspace;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        problem.evaluate_contribution_system(contribution, global_state, workspace);
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

// Piecewise-linear time functions.

PiecewiseLinearTimeTable::PiecewiseLinearTimeTable(std::string name, std::vector<double> times,
                                                   std::vector<double> values)
    : _name(std::move(name)), _times(std::move(times)), _values(std::move(values)) {
    if (_name.empty())
        throw std::invalid_argument("Time-table name must not be empty");
    if (_times.size() < 2 || _times.size() != _values.size())
        throw std::invalid_argument("Time table requires at least two time/value pairs: " + _name);
    for (std::size_t index = 0; index < _times.size(); ++index) {
        if (!std::isfinite(_times[index]) || _times[index] < 0.0 || !std::isfinite(_values[index]))
            throw std::invalid_argument("Time-table entries must be finite with nonnegative times: " + _name);
        if (index > 0 && !(_times[index] > _times[index - 1]))
            throw std::invalid_argument("Time-table times must be strictly increasing: " + _name);
    }
}

const std::string& PiecewiseLinearTimeTable::name() const noexcept {
    return _name;
}

const std::vector<double>& PiecewiseLinearTimeTable::times() const noexcept {
    return _times;
}

const std::vector<double>& PiecewiseLinearTimeTable::values() const noexcept {
    return _values;
}

double PiecewiseLinearTimeTable::value(double time) const {
    if (!std::isfinite(time) || time < 0.0)
        throw std::invalid_argument("Time-table evaluation time must be finite and nonnegative");
    if (time <= _times.front())
        return _values.front();
    if (time >= _times.back())
        return _values.back();
    const auto upper = std::upper_bound(_times.begin(), _times.end(), time);
    const std::size_t right = static_cast<std::size_t>(upper - _times.begin());
    const std::size_t left = right - 1;
    const double fraction = (time - _times[left]) / (_times[right] - _times[left]);
    return (1.0 - fraction) * _values[left] + fraction * _values[right];
}

} // namespace fuelsim
