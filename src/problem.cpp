#include "assembly.hpp"
#include "cartesian3d_assembly.hpp"
#include "fuelsim/cartesian3d_problem_access.hpp"
#include "fuelsim/diagnostics.hpp"
#include "fuelsim/hex8_thermoelastic.hpp"
#include "fuelsim/nonlinear_problem.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/rz_problem_access.hpp"
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

std::size_t GlobalStateView::global_size() const noexcept { return _global_size; }

std::size_t GlobalStateView::local_size() const noexcept {
    return _dense_values != nullptr ? _dense_values->size() : _sparse_values->size();
}

bool GlobalStateView::contains(std::size_t global_dof) const {
    if (global_dof >= _global_size) return false;
    if (_dense_values != nullptr) return true;
    if (global_dof > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) return false;
    return std::binary_search(_global_dofs->begin(), _global_dofs->end(), static_cast<std::uint32_t>(global_dof));
}

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
    if (dof >= dof_count()) throw std::out_of_range("NonlinearProblem field lookup DOF is out of range");
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
            if (dof >= dof_count()) throw std::out_of_range("NonlinearProblem contribution DOF is out of range");
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
    if (values.size() != local_dof_count) throw std::invalid_argument("RZ contribution state must contain 12 DOFs");
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

Hex8LocalValues hex8_local_values(const std::vector<double>& values) {
    if (values.size() != hex8_local_dof_count)
        throw std::invalid_argument("HEX8 contribution state must contain 32 DOFs");
    Hex8LocalValues result{};
    std::copy(values.begin(), values.end(), result.begin());
    return result;
}

Quad4FaceLocalValues face_local_values(const std::vector<double>& values) {
    if (values.size() != quad4_face_local_dof_count)
        throw std::invalid_argument("Three-dimensional face state must contain 16 DOFs");
    Quad4FaceLocalValues result{};
    std::copy(values.begin(), values.end(), result.begin());
    return result;
}

void copy_hex8_residual(const Hex8LocalResidual& source, std::vector<double>& destination) {
    destination.assign(source.begin(), source.end());
}

void copy_hex8_system(const Hex8LocalSystem& source, std::vector<double>& residual, std::vector<double>& jacobian) {
    residual.assign(source.residual.begin(), source.residual.end());
    jacobian.assign(source.jacobian.begin(), source.jacobian.end());
}

void copy_face_residual(const Quad4FaceLocalResidual& source, std::vector<double>& destination) {
    destination.assign(source.begin(), source.end());
}

void copy_face_system(const Quad4FaceLocalSystem& source, std::vector<double>& residual,
                      std::vector<double>& jacobian) {
    residual.assign(source.residual.begin(), source.residual.end());
    jacobian.assign(source.jacobian.begin(), source.jacobian.end());
}

Hex8LocalValues gather_hex8_state(const cartesian3d::SpatialAssembly& spatial, std::size_t contribution_index,
                                  const std::vector<double>& global_state) {
    std::vector<std::size_t> dofs;
    spatial.fill_contribution_dofs(contribution_index, dofs);
    if (dofs.size() != hex8_local_dof_count)
        throw std::logic_error("HEX8 volume contribution has an invalid DOF layout");
    Hex8LocalValues result{};
    for (std::size_t local = 0; local < result.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}

} // namespace

// Steady nonlinear problem.

class SteadyProblem::Implementation final {
  public:
    Implementation(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
        : spatial(std::make_unique<rz::SpatialAssembly>(std::move(definition), source_mesh)) {
        region_kernels.reserve(spatial->region_count());
        for (std::size_t region_value = 0; region_value < spatial->region_count(); ++region_value) {
            const RegionDefinition& value = spatial->region(region_value);
            region_kernels.emplace_back(IsotropicThermoelasticMaterial(value.material),
                                        spatial->region_heat_source(region_value), value.strain_formulation);
        }
    }

    Implementation(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
        : cartesian_spatial(std::make_unique<cartesian3d::SpatialAssembly>(std::move(definition), source_mesh)) {
        cartesian_region_kernels.reserve(cartesian_spatial->region_count());
        for (std::size_t region_value = 0; region_value < cartesian_spatial->region_count(); ++region_value) {
            const RegionDefinition& value = cartesian_spatial->region(region_value);
            cartesian_region_kernels.emplace_back(IsotropicThermoelasticMaterial(value.material),
                                                  cartesian_spatial->region_heat_source(region_value));
        }
    }

    bool is_cartesian() const noexcept { return cartesian_spatial != nullptr; }

    std::unique_ptr<rz::SpatialAssembly> spatial;
    std::vector<Quad4RzThermoelasticKernel> region_kernels;
    std::unique_ptr<cartesian3d::SpatialAssembly> cartesian_spatial;
    std::vector<Hex8ThermoelasticKernel> cartesian_region_kernels;
};

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : _implementation(std::make_unique<Implementation>(std::move(definition), source_mesh)) {}

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : _implementation(std::make_unique<Implementation>(std::move(definition), source_mesh)) {}

SteadyProblem::~SteadyProblem() = default;

bool SteadyProblem::is_cartesian_3d() const noexcept { return _implementation->is_cartesian(); }

const Hex8DofMap& cartesian3d::ProblemAccess::dof_map(const SteadyProblem& problem) noexcept {
    return problem._implementation->cartesian_spatial->dof_map();
}
std::size_t cartesian3d::ProblemAccess::region_count(const SteadyProblem& problem) noexcept {
    return problem._implementation->cartesian_spatial->region_count();
}
const RegionDefinition& cartesian3d::ProblemAccess::region(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->cartesian_spatial->region(index);
}
const Hex8RegionMesh& cartesian3d::ProblemAccess::region_mesh(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->cartesian_spatial->region_mesh(index);
}
const Hex8ThermoelasticKernel& cartesian3d::ProblemAccess::region_kernel(const SteadyProblem& problem,
                                                                         std::size_t index) {
    return problem._implementation->cartesian_region_kernels.at(index);
}
std::size_t cartesian3d::ProblemAccess::region_node_offset(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->cartesian_spatial->region_node_offset(index);
}
const Hex8Geometry& cartesian3d::ProblemAccess::region_element_geometry(const SteadyProblem& problem,
                                                                        std::size_t region_value,
                                                                        std::size_t element_index) {
    return problem._implementation->cartesian_spatial->region_element_geometry(region_value, element_index);
}

const SpatialDefinition& rz::ProblemAccess::definition(const SteadyProblem& problem) noexcept {
    return problem._implementation->spatial->definition();
}

const DofMap& rz::ProblemAccess::dof_map(const SteadyProblem& problem) noexcept {
    return problem._implementation->spatial->dof_map();
}

std::size_t rz::ProblemAccess::region_count(const SteadyProblem& problem) noexcept {
    return problem._implementation->spatial->region_count();
}

std::size_t rz::ProblemAccess::region_index(const SteadyProblem& problem, const std::string& name) {
    return problem._implementation->spatial->region_index(name);
}

const RegionDefinition& rz::ProblemAccess::region(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->spatial->region(index);
}

const RegionMesh& rz::ProblemAccess::region_mesh(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->spatial->region_mesh(index);
}

const Quad4RzThermoelasticKernel& rz::ProblemAccess::region_kernel(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->region_kernels.at(index);
}

std::size_t rz::ProblemAccess::region_node_offset(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->spatial->region_node_offset(index);
}

std::size_t rz::ProblemAccess::region_element_count(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->spatial->region_element_count(index);
}

std::size_t rz::ProblemAccess::region_element_offset(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->spatial->region_element_offset(index);
}

std::size_t rz::ProblemAccess::volume_contribution_count(const SteadyProblem& problem) noexcept {
    return problem._implementation->spatial->volume_contribution_count();
}

SpatialContributionType rz::ProblemAccess::contribution_type(const SteadyProblem& problem,
                                                             std::size_t contribution_index) {
    return problem._implementation->spatial->contribution_type(contribution_index);
}

const Quad4RzGeometry& rz::ProblemAccess::region_element_geometry(const SteadyProblem& problem,
                                                                  std::size_t region_value, std::size_t element_index) {
    return problem._implementation->spatial->region_element_geometry(region_value, element_index);
}

std::size_t rz::ProblemAccess::contact_count(const SteadyProblem& problem) noexcept {
    return problem._implementation->spatial->contact_count();
}

const ContactDefinition& rz::ProblemAccess::contact(const SteadyProblem& problem, std::size_t index) {
    return problem._implementation->spatial->contact(index);
}

const std::vector<std::vector<ContactPointHistory>>&
rz::ProblemAccess::committed_contact_histories(const SteadyProblem& problem) noexcept {
    return problem._implementation->spatial->committed_contact_histories();
}

void rz::ProblemAccess::commit_contact_state(SteadyProblem& problem, const std::vector<double>& state) {
    problem._implementation->spatial->commit_contact_state(state);
}

bool SteadyProblem::uses_augmented_contact() const noexcept {
    return !_implementation->is_cartesian() && _implementation->spatial->uses_augmented_contact();
}

AugmentedContactUpdate SteadyProblem::update_augmented_contact_multipliers(const std::vector<double>& state,
                                                                           std::size_t completed_updates) {
    if (_implementation->is_cartesian())
        throw std::logic_error("Cartesian three-dimensional stage B does not support augmented contact");
    return _implementation->spatial->update_augmented_contact_multipliers(state, completed_updates);
}

void rz::ProblemAccess::restore_contact_state(SteadyProblem& problem, const std::vector<double>& state,
                                              std::vector<std::vector<ContactPointHistory>> histories) {
    problem._implementation->spatial->restore_contact_state(state, std::move(histories));
}

void SteadyProblem::set_load_factor(double value) {
    if (_implementation->is_cartesian())
        _implementation->cartesian_spatial->set_load_factor(value);
    else
        _implementation->spatial->set_load_factor(value);
    refresh_region_heat_sources();
}

void SteadyProblem::refresh_region_heat_sources() {
    if (_implementation->is_cartesian()) {
        for (std::size_t region_value = 0; region_value < _implementation->cartesian_spatial->region_count();
             ++region_value)
            _implementation->cartesian_region_kernels[region_value].set_volumetric_heat_source(
                _implementation->cartesian_spatial->region_heat_source(region_value));
        return;
    }
    for (std::size_t region_value = 0; region_value < _implementation->spatial->region_count(); ++region_value)
        _implementation->region_kernels[region_value].set_volumetric_heat_source(
            _implementation->spatial->region_heat_source(region_value));
}

double SteadyProblem::load_factor() const noexcept {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->load_factor()
                                           : _implementation->spatial->load_factor();
}

void SteadyProblem::set_time(double value) {
    if (_implementation->is_cartesian()) {
        _implementation->cartesian_spatial->set_time(value);
        for (Hex8ThermoelasticKernel& kernel : _implementation->cartesian_region_kernels) kernel.set_time(value);
        refresh_region_heat_sources();
        return;
    }
    _implementation->spatial->set_time(value);
    for (Quad4RzThermoelasticKernel& kernel : _implementation->region_kernels) kernel.set_time(value);
    refresh_region_heat_sources();
}

std::vector<double> SteadyProblem::initial_state() const {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->initial_state()
                                           : _implementation->spatial->initial_state();
}

std::vector<ContactNodeSummary> rz::ProblemAccess::summarize_contact_nodes(const SteadyProblem& problem,
                                                                           std::size_t contact_index,
                                                                           const std::vector<double>& state) {
    return problem._implementation->spatial->summarize_contact_nodes(contact_index, state);
}

std::vector<std::size_t> rz::ProblemAccess::contact_secondary_source_nodes(const SteadyProblem& problem,
                                                                           std::size_t contact_index) {
    return problem._implementation->spatial->contact_secondary_source_nodes(contact_index);
}

InterfaceSummary rz::ProblemAccess::summarize_interface(const SteadyProblem& problem, std::size_t contact_index,
                                                        const std::vector<double>& state) {
    return problem._implementation->spatial->summarize_interface(contact_index, state);
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

bool SteadyStateSnapshot::empty() const noexcept { return _storage == nullptr; }

std::shared_ptr<const void> SteadyStateSnapshot::snapshot_owner() const noexcept {
    return _storage != nullptr ? _storage->owner : nullptr;
}

SteadyStateSnapshot SteadyProblem::capture_internal_state() const {
    std::shared_ptr<const rz::SteadyStateStorage> state;
    if (!_implementation->is_cartesian())
        state = std::make_shared<rz::SteadyStateStorage>(_implementation->spatial->committed_contact_histories());
    return SteadyStateSnapshot(
        std::make_shared<SteadyStateSnapshot::Storage>(discretization_identity(), std::move(state)));
}

void SteadyProblem::restore_internal_state(const SteadyStateSnapshot& snapshot, const std::vector<double>& state) {
    if (snapshot.empty()) throw std::invalid_argument("SteadyProblem cannot restore an empty internal-state snapshot");
    if (snapshot.snapshot_owner() != discretization_identity())
        throw std::invalid_argument("SteadyProblem cannot restore a snapshot from another problem");
    if (_implementation->is_cartesian()) {
        if (state.size() != dof_count()) throw std::invalid_argument("SteadyProblem restore state size mismatch");
        return;
    }
    _implementation->spatial->restore_contact_state(state, snapshot._storage->rz_state->contact_histories);
}

void SteadyProblem::commit_internal_state(const std::vector<double>& state) {
    if (_implementation->is_cartesian()) {
        if (state.size() != dof_count()) throw std::invalid_argument("SteadyProblem commit state size mismatch");
        return;
    }
    _implementation->spatial->commit_contact_state(state);
}

std::size_t SteadyProblem::dof_count() const noexcept {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->dof_count()
                                           : _implementation->spatial->dof_count();
}

std::size_t SteadyProblem::contribution_count() const noexcept {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->contribution_count()
                                           : _implementation->spatial->contribution_count();
}

const std::vector<FieldDescriptor>& SteadyProblem::field_layout() const noexcept {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->dof_map().field_layout()
                                           : _implementation->spatial->dof_map().field_layout();
}

const std::vector<DirichletCondition>& SteadyProblem::dirichlet_conditions() const noexcept {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->dirichlet_conditions()
                                           : _implementation->spatial->dirichlet_conditions();
}

void SteadyProblem::validate_state(const std::vector<double>& state) const {
    if (_implementation->is_cartesian())
        _implementation->cartesian_spatial->validate_state(state);
    else
        _implementation->spatial->validate_state(state);
}

std::vector<std::size_t> SteadyProblem::required_state_dofs(std::size_t contribution_begin,
                                                            std::size_t contribution_end) const {
    if (_implementation->is_cartesian())
        return NonlinearProblem::required_state_dofs(contribution_begin, contribution_end);
    return _implementation->spatial->required_state_dofs(contribution_begin, contribution_end);
}

void SteadyProblem::validate_local_state(std::size_t contribution_begin, std::size_t contribution_end,
                                         const GlobalStateView& state) const {
    if (_implementation->is_cartesian()) {
        NonlinearProblem::validate_local_state(contribution_begin, contribution_end, state);
        return;
    }
    _implementation->spatial->validate_local_state(contribution_begin, contribution_end, state);
}

std::size_t SteadyProblem::contribution_dof_count(std::size_t contribution_index) const {
    if (_implementation->is_cartesian())
        return _implementation->cartesian_spatial->contribution_dof_count(contribution_index);
    if (contribution_index >= contribution_count())
        throw std::out_of_range("SteadyProblem contribution index is out of range");
    return local_dof_count;
}

void SteadyProblem::fill_contribution_dofs(std::size_t contribution_index, std::vector<std::size_t>& dofs) const {
    if (_implementation->is_cartesian()) {
        _implementation->cartesian_spatial->fill_contribution_dofs(contribution_index, dofs);
        return;
    }
    const LocalDofs fixed = rz::ProblemAccess::contribution_dofs(*this, contribution_index);
    dofs.assign(fixed.begin(), fixed.end());
}

void SteadyProblem::compute_contribution_residual(std::size_t contribution_index, const std::vector<double>& state,
                                                  std::vector<double>& residual) const {
    if (_implementation->is_cartesian()) {
        if (contribution_index < _implementation->cartesian_spatial->volume_contribution_count()) {
            const auto location = _implementation->cartesian_spatial->element_location(contribution_index);
            copy_hex8_residual(
                _implementation->cartesian_region_kernels[location.first].residual(
                    _implementation->cartesian_spatial->region_element_geometry(location.first, location.second),
                    hex8_local_values(state)),
                residual);
        } else {
            copy_face_residual(
                _implementation->cartesian_spatial->boundary_residual(contribution_index, face_local_values(state)),
                residual);
        }
        return;
    }
    copy_rz_residual(rz::ProblemAccess::contribution_residual(*this, contribution_index, rz_local_values(state)),
                     residual);
}

void SteadyProblem::compute_contribution_system(std::size_t contribution_index, const std::vector<double>& state,
                                                std::vector<double>& residual, std::vector<double>& jacobian) const {
    if (_implementation->is_cartesian()) {
        if (contribution_index < _implementation->cartesian_spatial->volume_contribution_count()) {
            const auto location = _implementation->cartesian_spatial->element_location(contribution_index);
            copy_hex8_system(
                _implementation->cartesian_region_kernels[location.first].linearize(
                    _implementation->cartesian_spatial->region_element_geometry(location.first, location.second),
                    hex8_local_values(state)),
                residual, jacobian);
        } else {
            copy_face_system(
                _implementation->cartesian_spatial->boundary_system(contribution_index, face_local_values(state)),
                residual, jacobian);
        }
        return;
    }
    copy_rz_system(rz::ProblemAccess::linearize_contribution(*this, contribution_index, rz_local_values(state)),
                   residual, jacobian);
}

LocalDofs rz::ProblemAccess::contribution_dofs(const SteadyProblem& problem, std::size_t contribution_index) {
    return problem._implementation->spatial->contribution_dofs(contribution_index);
}

LocalValues rz::ProblemAccess::contribution_state(const SteadyProblem& problem, std::size_t contribution_index,
                                                  const std::vector<double>& global_state) {
    return contribution_state(problem, contribution_index, GlobalStateView(global_state));
}

LocalValues rz::ProblemAccess::contribution_state(const SteadyProblem& problem, std::size_t contribution_index,
                                                  const GlobalStateView& global_state) {
    if (global_state.global_size() != problem.dof_count())
        throw std::invalid_argument("SteadyProblem contribution state has the wrong global size");
    const LocalDofs dofs = contribution_dofs(problem, contribution_index);
    LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local) result[local] = global_state.value(dofs[local]);
    return result;
}

LocalResidual rz::ProblemAccess::contribution_residual(const SteadyProblem& problem, std::size_t contribution_index,
                                                       const LocalValues& state) {
    if (contribution_index < volume_contribution_count(problem)) {
        const auto location = problem._implementation->spatial->element_location(contribution_index);
        return problem._implementation->region_kernels[location.first].residual(
            region_element_geometry(problem, location.first, location.second), state);
    }
    return problem._implementation->spatial->contribution_residual(contribution_index, state);
}

LocalSystem rz::ProblemAccess::linearize_contribution(const SteadyProblem& problem, std::size_t contribution_index,
                                                      const LocalValues& state) {
    if (contribution_index < volume_contribution_count(problem)) {
        const auto location = problem._implementation->spatial->element_location(contribution_index);
        return problem._implementation->region_kernels[location.first].linearize(
            region_element_geometry(problem, location.first, location.second), state);
    }
    return problem._implementation->spatial->linearize_contribution(contribution_index, state);
}

class TransientProblem::Implementation final {
  public:
    Implementation(TransientProblemDefinition definition_value, const UnstructuredQuad4Mesh& source_mesh)
        : definition(std::move(definition_value)),
          spatial(std::make_unique<rz::SpatialAssembly>(definition.spatial, source_mesh)), committed_time(0.0),
          committed_load_factor(0.0), active_time_step(0.0), active_end_time(0.0), active_load_factor(0.0),
          time_step_active(false) {}

    Implementation(TransientProblemDefinition definition_value, const UnstructuredHex8Mesh& source_mesh)
        : definition(std::move(definition_value)),
          cartesian_spatial(std::make_unique<cartesian3d::SpatialAssembly>(definition.spatial, source_mesh)),
          committed_time(0.0), committed_load_factor(0.0), active_time_step(0.0), active_end_time(0.0),
          active_load_factor(0.0), time_step_active(false) {}

    bool is_cartesian() const noexcept { return cartesian_spatial != nullptr; }

    TransientProblemDefinition definition;
    std::unique_ptr<rz::SpatialAssembly> spatial;
    std::vector<Quad4RzTransientKernel> region_kernels;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> material_stresses;
    std::unique_ptr<cartesian3d::SpatialAssembly> cartesian_spatial;
    std::vector<Hex8ThermoelasticKernel> cartesian_region_kernels;
    std::vector<std::vector<std::array<SymmetricTensor3Values, 8>>> cartesian_stresses;
    TransientConservationSummary last_conservation_summary;
    std::vector<double> committed_solution;
    double committed_time;
    double committed_load_factor;
    double active_time_step;
    double active_end_time;
    double active_load_factor;
    std::vector<std::vector<ContactPointHistory>> active_contact_histories;
    bool time_step_active;
};

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
        const LocalDofs dofs = ProblemAccess::contribution_dofs(problem, contribution);
        const LocalValues state = ProblemAccess::contribution_state(problem, contribution, converged_solution);
        const LocalResidual residual = ProblemAccess::contribution_residual(problem, contribution, state);
        const SpatialContributionType type = problem._implementation->spatial->contribution_type(contribution);
        for (std::size_t local = 0; local < local_dof_count; ++local) {
            raw_residual[dofs[local]] += residual[local];
            const double increment =
                converged_solution[dofs[local]] - problem._implementation->committed_solution[dofs[local]];
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
    for (std::size_t region_value = 0; region_value < ProblemAccess::region_count(problem); ++region_value) {
        const Quad4RzTransientKernel& kernel = problem._implementation->region_kernels[region_value];
        const std::size_t offset = problem._implementation->spatial->region_element_offset(region_value);
        for (std::size_t element = 0; element < staged_histories[region_value].size(); ++element) {
            const LocalValues current =
                ProblemAccess::contribution_state(problem, offset + element, converged_solution);
            const LocalValues old = ProblemAccess::contribution_state(problem, offset + element,
                                                                      problem._implementation->committed_solution);
            const Quad4RzGeometry& geometry = ProblemAccess::region_element_geometry(problem, region_value, element);
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
                                           (current_temperature - old_temperature) /
                                           problem._implementation->active_time_step;
                result.generated_heat_rate += point.weighted_measure * kernel.volumetric_heat_source();

                const MaterialPointState& old_history =
                    problem._implementation->material_histories[region_value][element][q];
                const MaterialPointState& new_history = staged_histories[region_value][element][q];
                const AxisymmetricStressValues& old_stress =
                    problem._implementation->material_stresses[region_value][element][q];
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
        const double increment =
            converged_solution[condition.dof] - problem._implementation->committed_solution[condition.dof];
        if (condition.dof < ProblemAccess::dof_map(problem).node_count())
            result.dirichlet_heat_input_rate += raw_residual[condition.dof];
        else
            result.dirichlet_reaction_work_increment += raw_residual[condition.dof] * increment;
    }
    double thermal_residual_squared = 0.0;
    double mechanical_residual_squared = 0.0;
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        if (constrained[dof]) continue;
        const double value = raw_residual[dof];
        if (dof < ProblemAccess::dof_map(problem).node_count())
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
    : _implementation(std::make_unique<Implementation>(std::move(definition), source_mesh)) {
    validate_definition(_implementation->definition);
    const std::size_t regions = _implementation->definition.spatial.regions.size();
    _implementation->region_kernels.reserve(regions);
    _implementation->material_histories.resize(regions);
    _implementation->material_stresses.resize(regions);
    for (std::size_t region_value = 0; region_value < regions; ++region_value) {
        _implementation->region_kernels.emplace_back(
            IsotropicInelasticMaterial(_implementation->definition.spatial.regions[region_value].material,
                                       _implementation->definition.regions[region_value].material),
            0.0, _implementation->definition.spatial.regions[region_value].strain_formulation);
        _implementation->material_histories[region_value].resize(
            _implementation->spatial->region_element_count(region_value));
        _implementation->material_stresses[region_value].resize(
            _implementation->spatial->region_element_count(region_value));
    }
    apply_spatial_controls(0.0, 0.0);
    _implementation->committed_solution = _implementation->spatial->initial_state();
    _implementation->spatial->restore_contact_state(_implementation->committed_solution,
                                                    _implementation->spatial->committed_contact_histories());
}

TransientProblem::TransientProblem(TransientProblemDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : _implementation(std::make_unique<Implementation>(std::move(definition), source_mesh)) {
    validate_definition(_implementation->definition);
    const std::size_t regions = _implementation->definition.spatial.regions.size();
    _implementation->cartesian_region_kernels.reserve(regions);
    _implementation->cartesian_stresses.resize(regions);
    for (std::size_t region_value = 0; region_value < regions; ++region_value) {
        const TransientInelasticProperties& transient = _implementation->definition.regions[region_value].material;
        if (transient.behavior != InelasticBehavior::elastic)
            throw std::invalid_argument("Cartesian three-dimensional stage B supports only elastic materials");
        _implementation->cartesian_region_kernels.emplace_back(
            IsotropicThermoelasticMaterial(_implementation->definition.spatial.regions[region_value].material),
            transient.density * transient.specific_heat, 0.0);
        _implementation->cartesian_stresses[region_value].resize(
            _implementation->cartesian_spatial->region_element_count(region_value));
    }
    apply_spatial_controls(0.0, 0.0);
    _implementation->committed_solution = _implementation->cartesian_spatial->initial_state();
}

TransientProblem::~TransientProblem() = default;

bool TransientProblem::is_cartesian_3d() const noexcept { return _implementation->is_cartesian(); }

const TransientProblemDefinition& cartesian3d::ProblemAccess::definition(const TransientProblem& problem) noexcept {
    return problem._implementation->definition;
}
const Hex8DofMap& cartesian3d::ProblemAccess::dof_map(const TransientProblem& problem) noexcept {
    return problem._implementation->cartesian_spatial->dof_map();
}
std::size_t cartesian3d::ProblemAccess::region_count(const TransientProblem& problem) noexcept {
    return problem._implementation->cartesian_spatial->region_count();
}
const RegionDefinition& cartesian3d::ProblemAccess::region(const TransientProblem& problem, std::size_t index) {
    return problem._implementation->cartesian_spatial->region(index);
}
const Hex8RegionMesh& cartesian3d::ProblemAccess::region_mesh(const TransientProblem& problem, std::size_t index) {
    return problem._implementation->cartesian_spatial->region_mesh(index);
}
const Hex8ThermoelasticKernel& cartesian3d::ProblemAccess::region_kernel(const TransientProblem& problem,
                                                                         std::size_t index) {
    return problem._implementation->cartesian_region_kernels.at(index);
}
std::size_t cartesian3d::ProblemAccess::region_node_offset(const TransientProblem& problem, std::size_t index) {
    return problem._implementation->cartesian_spatial->region_node_offset(index);
}
const Hex8Geometry& cartesian3d::ProblemAccess::region_element_geometry(const TransientProblem& problem,
                                                                        std::size_t region_value,
                                                                        std::size_t element_index) {
    return problem._implementation->cartesian_spatial->region_element_geometry(region_value, element_index);
}
const std::array<SymmetricTensor3Values, 8>& cartesian3d::ProblemAccess::stress(const TransientProblem& problem,
                                                                                std::size_t region_value,
                                                                                std::size_t element_index) {
    return problem._implementation->cartesian_stresses.at(region_value).at(element_index);
}
cartesian3d::TransientCommittedState cartesian3d::ProblemAccess::committed_state(const TransientProblem& problem) {
    return {problem._implementation->committed_solution, problem._implementation->cartesian_stresses,
            problem._implementation->last_conservation_summary, problem._implementation->committed_time,
            problem._implementation->committed_load_factor};
}
void cartesian3d::ProblemAccess::restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
    if (problem._implementation->time_step_active)
        throw std::logic_error("TransientProblem cannot restore during an active time step");
    if (state.solution.size() != problem.dof_count() || state.stresses.size() != region_count(problem) ||
        !std::isfinite(state.time) || state.time < 0.0 || !std::isfinite(state.load_factor) || state.load_factor < 0.0)
        throw std::invalid_argument("Cartesian committed state layout is invalid");
    problem._implementation->cartesian_spatial->validate_state(state.solution);
    for (std::size_t region_value = 0; region_value < state.stresses.size(); ++region_value) {
        if (state.stresses[region_value].size() !=
            problem._implementation->cartesian_spatial->region_element_count(region_value))
            throw std::invalid_argument("Cartesian committed stress layout is invalid");
        for (const auto& element : state.stresses[region_value]) {
            for (const SymmetricTensor3Values& stress : element) {
                if (!std::isfinite(stress.xx) || !std::isfinite(stress.yy) || !std::isfinite(stress.zz) ||
                    !std::isfinite(stress.xy) || !std::isfinite(stress.yz) || !std::isfinite(stress.xz))
                    throw std::invalid_argument("Cartesian committed stress values must be finite");
            }
        }
    }
    problem._implementation->committed_solution = std::move(state.solution);
    problem._implementation->cartesian_stresses = std::move(state.stresses);
    problem._implementation->last_conservation_summary = state.conservation;
    problem._implementation->committed_time = state.time;
    problem._implementation->committed_load_factor = state.load_factor;
    problem.clear_active_time_step();
    problem.apply_spatial_controls(state.time, state.load_factor);
}

const TransientProblemDefinition& rz::ProblemAccess::definition(const TransientProblem& problem) noexcept {
    return problem._implementation->definition;
}

const DofMap& rz::ProblemAccess::dof_map(const TransientProblem& problem) noexcept {
    return problem._implementation->spatial->dof_map();
}

std::size_t rz::ProblemAccess::region_count(const TransientProblem& problem) noexcept {
    return problem._implementation->definition.spatial.regions.size();
}

std::size_t rz::ProblemAccess::region_index(const TransientProblem& problem, const std::string& name) {
    return problem._implementation->spatial->region_index(name);
}

std::size_t rz::ProblemAccess::region_node_offset(const TransientProblem& problem, std::size_t region_value) {
    return problem._implementation->spatial->region_node_offset(region_value);
}

const RegionDefinition& rz::ProblemAccess::region(const TransientProblem& problem, std::size_t index) {
    return problem._implementation->spatial->region(index);
}

const RegionMesh& rz::ProblemAccess::region_mesh(const TransientProblem& problem, std::size_t index) {
    return problem._implementation->spatial->region_mesh(index);
}

const Quad4RzTransientKernel& rz::ProblemAccess::region_kernel(const TransientProblem& problem, std::size_t index) {
    return problem._implementation->region_kernels.at(index);
}

const Quad4RzGeometry& rz::ProblemAccess::region_element_geometry(const TransientProblem& problem,
                                                                  std::size_t region_value, std::size_t element_index) {
    return problem._implementation->spatial->region_element_geometry(region_value, element_index);
}

const std::vector<double>& TransientProblem::committed_solution() const noexcept {
    return _implementation->committed_solution;
}

double TransientProblem::committed_time() const noexcept { return _implementation->committed_time; }

double TransientProblem::committed_load_factor() const noexcept { return _implementation->committed_load_factor; }

bool TransientProblem::time_step_active() const noexcept { return _implementation->time_step_active; }

std::vector<double> TransientProblem::time_events() const {
    std::vector<double> result;
    for (const PiecewiseLinearTimeTable& table : _implementation->definition.spatial.time_tables)
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

namespace cartesian3d {

struct TransientStateStorage final {
    explicit TransientStateStorage(TransientCommittedState value) : state(std::move(value)) {}

    TransientCommittedState state;
};

} // namespace cartesian3d

struct TransientStateSnapshot::Storage final {
    Storage(std::shared_ptr<const void> owner_value, std::shared_ptr<const rz::TransientStateStorage> rz_value,
            std::shared_ptr<const cartesian3d::TransientStateStorage> cartesian_value)
        : owner(std::move(owner_value)), rz_state(std::move(rz_value)), cartesian_state(std::move(cartesian_value)) {}

    std::shared_ptr<const void> owner;
    std::shared_ptr<const rz::TransientStateStorage> rz_state;
    std::shared_ptr<const cartesian3d::TransientStateStorage> cartesian_state;
};

TransientStateSnapshot::TransientStateSnapshot() = default;
TransientStateSnapshot::~TransientStateSnapshot() = default;
TransientStateSnapshot::TransientStateSnapshot(const TransientStateSnapshot& other) = default;
TransientStateSnapshot& TransientStateSnapshot::operator=(const TransientStateSnapshot& other) = default;
TransientStateSnapshot::TransientStateSnapshot(TransientStateSnapshot&& other) noexcept = default;
TransientStateSnapshot& TransientStateSnapshot::operator=(TransientStateSnapshot&& other) noexcept = default;

TransientStateSnapshot::TransientStateSnapshot(std::shared_ptr<const Storage> storage) : _storage(std::move(storage)) {}

bool TransientStateSnapshot::empty() const noexcept { return _storage == nullptr; }

std::shared_ptr<const void> TransientStateSnapshot::snapshot_owner() const noexcept {
    return _storage != nullptr ? _storage->owner : nullptr;
}

rz::TransientCommittedState rz::ProblemAccess::committed_state(const TransientProblem& problem) {
    return {problem._implementation->committed_solution,
            problem._implementation->material_histories,
            problem._implementation->material_stresses,
            problem._implementation->spatial->committed_contact_histories(),
            problem._implementation->last_conservation_summary,
            problem._implementation->committed_time,
            problem._implementation->committed_load_factor};
}

TransientStateSnapshot TransientProblem::capture_state() const {
    if (_implementation->time_step_active)
        throw std::logic_error("TransientProblem cannot capture an active time step");
    if (_implementation->is_cartesian()) {
        const auto state = std::make_shared<cartesian3d::TransientStateStorage>(cartesian3d::TransientCommittedState{
            _implementation->committed_solution, _implementation->cartesian_stresses,
            _implementation->last_conservation_summary, _implementation->committed_time,
            _implementation->committed_load_factor});
        return TransientStateSnapshot(
            std::make_shared<TransientStateSnapshot::Storage>(discretization_identity(), nullptr, std::move(state)));
    }
    const auto state = std::make_shared<rz::TransientStateStorage>(rz::ProblemAccess::committed_state(*this));
    return TransientStateSnapshot(
        std::make_shared<TransientStateSnapshot::Storage>(discretization_identity(), std::move(state), nullptr));
}

void TransientProblem::restore_state(const TransientStateSnapshot& snapshot) {
    if (snapshot.empty()) throw std::invalid_argument("TransientProblem cannot restore an empty state snapshot");
    if (snapshot.snapshot_owner() != discretization_identity())
        throw std::invalid_argument("TransientProblem cannot restore a snapshot from another problem");
    if (_implementation->is_cartesian()) {
        const cartesian3d::TransientCommittedState& state = snapshot._storage->cartesian_state->state;
        if (state.solution.size() != dof_count() ||
            state.stresses.size() != _implementation->cartesian_stresses.size() || !std::isfinite(state.time) ||
            state.time < 0.0 || !std::isfinite(state.load_factor) || state.load_factor < 0.0)
            throw std::invalid_argument("Cartesian transient snapshot layout is invalid");
        for (std::size_t region_value = 0; region_value < state.stresses.size(); ++region_value) {
            if (state.stresses[region_value].size() !=
                _implementation->cartesian_spatial->region_element_count(region_value))
                throw std::invalid_argument("Cartesian transient snapshot element layout is invalid");
        }
        _implementation->cartesian_spatial->validate_state(state.solution);
        _implementation->committed_solution = state.solution;
        _implementation->cartesian_stresses = state.stresses;
        _implementation->last_conservation_summary = state.conservation;
        _implementation->committed_time = state.time;
        _implementation->committed_load_factor = state.load_factor;
        clear_active_time_step();
        apply_spatial_controls(state.time, state.load_factor);
        return;
    }
    rz::ProblemAccess::restore_committed_state(*this, snapshot._storage->rz_state->state);
}

void rz::ProblemAccess::restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
    if (problem._implementation->time_step_active)
        throw std::logic_error("TransientProblem cannot restore during an active time step");
    if (state.solution.size() != problem.dof_count() || state.material_histories.size() != region_count(problem) ||
        state.material_stresses.size() != region_count(problem) ||
        state.contact_histories.size() != problem._implementation->definition.spatial.contacts.size())
        throw std::invalid_argument("Transient committed state layout does not match the problem");
    if (!std::isfinite(state.time) || state.time < 0.0 || !std::isfinite(state.load_factor) || state.load_factor < 0.0)
        throw std::invalid_argument("Transient committed time and load factor must be valid");
    const DofMap& dofs = dof_map(problem);
    for (std::size_t node = 0; node < dofs.node_count(); ++node) {
        const double temperature = state.solution.at(dofs.temperature(node));
        if (!std::isfinite(temperature) || !(temperature > 0.0) ||
            !std::isfinite(state.solution.at(dofs.radial_displacement(node))) ||
            !std::isfinite(state.solution.at(dofs.axial_displacement(node))))
            throw std::invalid_argument("Transient committed nodal state must be finite with "
                                        "positive temperatures");
    }
    for (std::size_t region_value = 0; region_value < region_count(problem); ++region_value) {
        const std::size_t elements = region_mesh(problem, region_value).elements().size();
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

    problem._implementation->spatial->restore_contact_state(state.solution, state.contact_histories);
    problem._implementation->committed_solution = std::move(state.solution);
    problem._implementation->material_histories = std::move(state.material_histories);
    problem._implementation->material_stresses = std::move(state.material_stresses);
    problem._implementation->last_conservation_summary = state.conservation;
    problem._implementation->committed_time = state.time;
    problem._implementation->committed_load_factor = state.load_factor;
    problem.clear_active_time_step();
    problem.apply_spatial_controls(problem._implementation->committed_time,
                                   problem._implementation->committed_load_factor);
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
    if (accumulator.count == 0) return 0.0;
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
    if (_implementation->is_cartesian()) {
        const cartesian3d::TransientCommittedState& full = full_snapshot._storage->cartesian_state->state;
        const cartesian3d::TransientCommittedState& half = half_snapshot._storage->cartesian_state->state;
        if (full.solution.size() != dof_count() || half.solution.size() != dof_count() ||
            full.stresses.size() != half.stresses.size())
            throw std::logic_error("Cartesian step-doubling snapshot layouts differ");
        const auto normalized = [](double difference_squared, double solution_squared, std::size_t count,
                                   double absolute_tolerance, double relative_tolerance) {
            if (count == 0) return 0.0;
            const double denominator = absolute_tolerance * std::sqrt(static_cast<double>(count)) +
                                       relative_tolerance * std::sqrt(solution_squared);
            return std::sqrt(difference_squared) / denominator;
        };
        TransientTimeErrorEstimate result;
        for (const FieldDescriptor& field : field_layout()) {
            double difference_squared = 0.0;
            double solution_squared = 0.0;
            for (std::size_t dof = field.begin; dof < field.end; ++dof) {
                const double difference = half.solution[dof] - full.solution[dof];
                difference_squared += difference * difference;
                solution_squared += half.solution[dof] * half.solution[dof];
            }
            const double absolute_tolerance = field.category == FieldCategory::thermal
                                                  ? options.temperature_time_absolute_tolerance
                                                  : options.displacement_time_absolute_tolerance;
            const double value = normalized(difference_squared, solution_squared, field.end - field.begin,
                                            absolute_tolerance, options.time_error_relative_tolerance);
            result.nodal_fields.push_back({field.name, value});
            result.maximum = std::max(result.maximum, value);
        }
        double stress_difference_squared = 0.0;
        double stress_solution_squared = 0.0;
        std::size_t stress_count = 0;
        for (std::size_t region_value = 0; region_value < full.stresses.size(); ++region_value) {
            if (full.stresses[region_value].size() != half.stresses[region_value].size())
                throw std::logic_error("Cartesian step-doubling stress layouts differ");
            for (std::size_t element = 0; element < full.stresses[region_value].size(); ++element) {
                for (std::size_t q = 0; q < 8; ++q) {
                    const SymmetricTensor3Values& full_value = full.stresses[region_value][element][q];
                    const SymmetricTensor3Values& half_value = half.stresses[region_value][element][q];
                    const std::array<double, 6> full_components = {full_value.xx, full_value.yy, full_value.zz,
                                                                   full_value.xy, full_value.yz, full_value.xz};
                    const std::array<double, 6> half_components = {half_value.xx, half_value.yy, half_value.zz,
                                                                   half_value.xy, half_value.yz, half_value.xz};
                    for (std::size_t component = 0; component < 6; ++component) {
                        const double difference = half_components[component] - full_components[component];
                        stress_difference_squared += difference * difference;
                        stress_solution_squared += half_components[component] * half_components[component];
                        ++stress_count;
                    }
                }
            }
        }
        result.stress =
            normalized(stress_difference_squared, stress_solution_squared, stress_count,
                       options.stress_history_time_absolute_tolerance, options.time_error_relative_tolerance);
        result.maximum = std::max(result.maximum, result.stress);
        return result;
    }
    return rz::compare_step_doubling_states(full_snapshot._storage->rz_state->state,
                                            half_snapshot._storage->rz_state->state, field_layout(), dof_count(),
                                            options);
}

void TransientProblem::combine_last_half_step_conservation(const TransientConservationSummary& first_half) {
    if (_implementation->time_step_active)
        throw std::logic_error("TransientProblem cannot combine conservation during an active time step");
    _implementation->last_conservation_summary =
        rz::combine_rz_half_step_conservation(first_half, _implementation->last_conservation_summary);
}

void TransientProblem::begin_time_step(const TransientStepInput& input) {
    if (_implementation->time_step_active) throw std::logic_error("TransientProblem already has an active time step");
    if (!std::isfinite(input.end_time) || input.end_time <= _implementation->committed_time)
        throw std::invalid_argument("TransientProblem end time must exceed committed time");
    if (!std::isfinite(input.load_factor) || input.load_factor < 0.0)
        throw std::invalid_argument("TransientProblem load factor must be finite and nonnegative");
    _implementation->active_time_step = input.end_time - _implementation->committed_time;
    _implementation->active_end_time = input.end_time;
    _implementation->active_load_factor = input.load_factor;
    if (!_implementation->is_cartesian() && _implementation->spatial->uses_augmented_contact())
        _implementation->active_contact_histories = _implementation->spatial->committed_contact_histories();
    try {
        apply_spatial_controls(input.end_time, input.load_factor);
    } catch (...) {
        if (!_implementation->is_cartesian() && !_implementation->active_contact_histories.empty())
            _implementation->spatial->restore_contact_state(_implementation->committed_solution,
                                                            std::move(_implementation->active_contact_histories));
        apply_spatial_controls(_implementation->committed_time, _implementation->committed_load_factor);
        clear_active_time_step();
        throw;
    }
    _implementation->time_step_active = true;
}

void TransientProblem::commit_time_step(const std::vector<double>& converged_solution) {
    require_active_time_step();
    if (converged_solution.size() != dof_count())
        throw std::invalid_argument("TransientProblem committed solution size mismatch");
    if (!std::all_of(converged_solution.begin(), converged_solution.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::domain_error("TransientProblem committed solution must be finite");
    if (_implementation->is_cartesian()) {
        const Hex8DofMap& dofs = _implementation->cartesian_spatial->dof_map();
        for (std::size_t node = 0; node < dofs.node_count(); ++node) {
            if (!(converged_solution[dofs.temperature(node)] > 0.0))
                throw std::domain_error("TransientProblem committed temperatures must be positive");
        }
        std::vector<std::vector<std::array<SymmetricTensor3Values, 8>>> staged(
            _implementation->cartesian_spatial->region_count());
        for (std::size_t region_value = 0; region_value < staged.size(); ++region_value) {
            staged[region_value].resize(_implementation->cartesian_spatial->region_element_count(region_value));
            const std::size_t offset = _implementation->cartesian_spatial->region_element_offset(region_value);
            for (std::size_t element = 0; element < staged[region_value].size(); ++element) {
                staged[region_value][element] = _implementation->cartesian_region_kernels[region_value].stress_values(
                    _implementation->cartesian_spatial->region_element_geometry(region_value, element),
                    gather_hex8_state(*_implementation->cartesian_spatial, offset + element, converged_solution));
            }
        }

        TransientConservationSummary conservation;
        std::vector<double> raw_residual(dof_count(), 0.0);
        std::vector<std::size_t> local_dofs;
        for (std::size_t contribution = 0; contribution < contribution_count(); ++contribution) {
            local_dofs.clear();
            fill_contribution_dofs(contribution, local_dofs);
            std::vector<double> local_state(local_dofs.size());
            for (std::size_t local = 0; local < local_dofs.size(); ++local)
                local_state[local] = converged_solution[local_dofs[local]];
            std::vector<double> local_residual;
            compute_contribution_residual(contribution, local_state, local_residual);
            const SpatialContributionType type = _implementation->cartesian_spatial->contribution_type(contribution);
            const std::size_t thermal_count = type == SpatialContributionType::volume ? 8 : 4;
            for (std::size_t local = 0; local < local_dofs.size(); ++local) {
                raw_residual[local_dofs[local]] += local_residual[local];
                if (local < thermal_count) {
                    if (type == SpatialContributionType::convection)
                        conservation.convection_heat_rate += local_residual[local];
                    continue;
                }
                const double increment =
                    converged_solution[local_dofs[local]] - _implementation->committed_solution[local_dofs[local]];
                if (type == SpatialContributionType::volume)
                    conservation.internal_mechanical_work_increment += local_residual[local] * increment;
                else if (type == SpatialContributionType::pressure || type == SpatialContributionType::traction)
                    conservation.pressure_traction_work_increment -= local_residual[local] * increment;
            }
        }
        for (std::size_t region_value = 0; region_value < _implementation->cartesian_spatial->region_count();
             ++region_value) {
            const Hex8ThermoelasticKernel& kernel = _implementation->cartesian_region_kernels[region_value];
            const std::size_t offset = _implementation->cartesian_spatial->region_element_offset(region_value);
            for (std::size_t element = 0;
                 element < _implementation->cartesian_spatial->region_element_count(region_value); ++element) {
                const Hex8LocalValues current =
                    gather_hex8_state(*_implementation->cartesian_spatial, offset + element, converged_solution);
                const Hex8LocalValues old = gather_hex8_state(*_implementation->cartesian_spatial, offset + element,
                                                              _implementation->committed_solution);
                const Hex8Geometry& geometry =
                    _implementation->cartesian_spatial->region_element_geometry(region_value, element);
                for (const Hex8QuadraturePoint& point : geometry.points) {
                    double current_temperature = 0.0;
                    double old_temperature = 0.0;
                    for (std::size_t node = 0; node < 8; ++node) {
                        current_temperature += point.shape[node] * current[node];
                        old_temperature += point.shape[node] * old[node];
                    }
                    conservation.stored_heat_rate += point.weighted_measure *
                                                     kernel.heat_capacity(current_temperature, point.position.x,
                                                                          point.position.y, point.position.z) *
                                                     (current_temperature - old_temperature) /
                                                     _implementation->active_time_step;
                    conservation.generated_heat_rate += point.weighted_measure * kernel.volumetric_heat_source();
                }
            }
        }
        std::vector<bool> constrained(dof_count(), false);
        for (const DirichletCondition& condition : dirichlet_conditions()) {
            constrained[condition.dof] = true;
            const double increment =
                converged_solution[condition.dof] - _implementation->committed_solution[condition.dof];
            if (field_layout()[field_index(condition.dof)].category == FieldCategory::thermal)
                conservation.dirichlet_heat_input_rate += raw_residual[condition.dof];
            else
                conservation.dirichlet_reaction_work_increment += raw_residual[condition.dof] * increment;
        }
        double thermal_residual_squared = 0.0;
        double mechanical_residual_squared = 0.0;
        for (std::size_t dof = 0; dof < dof_count(); ++dof) {
            if (constrained[dof]) continue;
            if (field_layout()[field_index(dof)].category == FieldCategory::thermal)
                thermal_residual_squared += raw_residual[dof] * raw_residual[dof];
            else
                mechanical_residual_squared += raw_residual[dof] * raw_residual[dof];
        }
        conservation.unconstrained_thermal_residual_l2 = std::sqrt(thermal_residual_squared);
        conservation.unconstrained_mechanical_residual_l2 = std::sqrt(mechanical_residual_squared);
        conservation.global_thermal_balance = conservation.stored_heat_rate + conservation.convection_heat_rate -
                                              conservation.generated_heat_rate - conservation.dirichlet_heat_input_rate;
        const double thermal_scale =
            std::abs(conservation.stored_heat_rate) + std::abs(conservation.convection_heat_rate) +
            std::abs(conservation.generated_heat_rate) + std::abs(conservation.dirichlet_heat_input_rate);
        conservation.relative_thermal_balance =
            thermal_scale > 0.0 ? std::abs(conservation.global_thermal_balance) / thermal_scale : 0.0;
        conservation.mechanical_work_balance = conservation.internal_mechanical_work_increment -
                                               conservation.pressure_traction_work_increment -
                                               conservation.dirichlet_reaction_work_increment;
        const double mechanical_scale = std::abs(conservation.internal_mechanical_work_increment) +
                                        std::abs(conservation.pressure_traction_work_increment) +
                                        std::abs(conservation.dirichlet_reaction_work_increment);
        conservation.relative_mechanical_work_balance =
            mechanical_scale > 0.0 ? std::abs(conservation.mechanical_work_balance) / mechanical_scale : 0.0;

        _implementation->cartesian_stresses.swap(staged);
        _implementation->last_conservation_summary = conservation;
        _implementation->committed_solution = converged_solution;
        _implementation->committed_time = _implementation->active_end_time;
        _implementation->committed_load_factor = _implementation->active_load_factor;
        clear_active_time_step();
        return;
    }
    const DofMap& dofs = _implementation->spatial->dof_map();
    for (std::size_t node = 0; node < dofs.node_count(); ++node) {
        if (!(converged_solution[dofs.temperature(node)] > 0.0))
            throw std::domain_error("TransientProblem committed temperatures must be positive");
    }

    const std::size_t regions = _implementation->spatial->region_count();
    std::vector<std::vector<Quad4MaterialHistory>> staged(regions);
    std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>> staged_stresses(regions);
    for (std::size_t region_value = 0; region_value < regions; ++region_value) {
        staged[region_value].resize(_implementation->spatial->region_element_count(region_value));
        staged_stresses[region_value].resize(_implementation->spatial->region_element_count(region_value));
        const std::size_t offset = _implementation->spatial->region_element_offset(region_value);
        for (std::size_t element = 0; element < staged[region_value].size(); ++element) {
            const LocalValues state =
                rz::ProblemAccess::contribution_state(*this, offset + element, converged_solution);
            const LocalValues committed_state =
                rz::ProblemAccess::contribution_state(*this, offset + element, _implementation->committed_solution);
            staged[region_value][element] = _implementation->region_kernels[region_value].trial_state_values(
                rz::ProblemAccess::region_element_geometry(*this, region_value, element), state, committed_state,
                _implementation->material_histories[region_value][element], _implementation->active_time_step);
            staged_stresses[region_value][element] = _implementation->region_kernels[region_value].stress_values(
                rz::ProblemAccess::region_element_geometry(*this, region_value, element), state, committed_state,
                _implementation->material_histories[region_value][element], _implementation->active_time_step);
        }
    }

    const TransientConservationSummary conservation =
        rz::TransientConservationCalculator::summarize(*this, converged_solution, staged, staged_stresses);
    _implementation->spatial->commit_contact_state(converged_solution);
    _implementation->material_histories.swap(staged);
    _implementation->material_stresses.swap(staged_stresses);
    _implementation->last_conservation_summary = conservation;
    _implementation->committed_solution = converged_solution;
    _implementation->committed_time = _implementation->active_end_time;
    _implementation->committed_load_factor = _implementation->active_load_factor;
    clear_active_time_step();
}

void TransientProblem::rollback_time_step() noexcept {
    if (!_implementation->time_step_active) return;
    apply_spatial_controls(_implementation->committed_time, _implementation->committed_load_factor);
    if (!_implementation->is_cartesian() && !_implementation->active_contact_histories.empty())
        _implementation->spatial->restore_contact_state(_implementation->committed_solution,
                                                        std::move(_implementation->active_contact_histories));
    clear_active_time_step();
}

void TransientProblem::apply_spatial_controls(double time, double load_factor) {
    if (_implementation->is_cartesian()) {
        _implementation->cartesian_spatial->set_time(time);
        _implementation->cartesian_spatial->set_load_factor(load_factor);
        for (Hex8ThermoelasticKernel& kernel : _implementation->cartesian_region_kernels) kernel.set_time(time);
        refresh_region_heat_sources();
        return;
    }
    _implementation->spatial->set_time(time);
    _implementation->spatial->set_load_factor(load_factor);
    for (Quad4RzTransientKernel& kernel : _implementation->region_kernels) kernel.set_time(time);
    refresh_region_heat_sources();
}

void TransientProblem::clear_active_time_step() noexcept {
    _implementation->active_time_step = 0.0;
    _implementation->active_end_time = _implementation->committed_time;
    _implementation->active_load_factor = _implementation->committed_load_factor;
    _implementation->active_contact_histories.clear();
    _implementation->time_step_active = false;
}

void TransientProblem::refresh_region_heat_sources() {
    if (_implementation->is_cartesian()) {
        for (std::size_t region_value = 0; region_value < _implementation->cartesian_spatial->region_count();
             ++region_value)
            _implementation->cartesian_region_kernels[region_value].set_volumetric_heat_source(
                _implementation->cartesian_spatial->region_heat_source(region_value));
        return;
    }
    for (std::size_t region_value = 0; region_value < _implementation->spatial->region_count(); ++region_value)
        _implementation->region_kernels[region_value].set_volumetric_heat_source(
            _implementation->spatial->region_heat_source(region_value));
}

bool TransientProblem::uses_augmented_contact() const noexcept {
    return !_implementation->is_cartesian() && _implementation->spatial->uses_augmented_contact();
}

AugmentedContactUpdate TransientProblem::update_augmented_contact_multipliers(const std::vector<double>& state,
                                                                              std::size_t completed_updates) {
    require_active_time_step();
    if (_implementation->is_cartesian())
        throw std::logic_error("Cartesian three-dimensional stage B does not support augmented contact");
    return _implementation->spatial->update_augmented_contact_multipliers(state, completed_updates);
}

const Quad4MaterialHistory& rz::ProblemAccess::material_history(const TransientProblem& problem,
                                                                std::size_t region_value, std::size_t element_index) {
    return problem._implementation->material_histories.at(region_value).at(element_index);
}

const std::array<AxisymmetricStressValues, 4>& rz::ProblemAccess::material_stress(const TransientProblem& problem,
                                                                                  std::size_t region_value,
                                                                                  std::size_t element_index) {
    return problem._implementation->material_stresses.at(region_value).at(element_index);
}

RegionInelasticSummary rz::ProblemAccess::summarize_region_history(const TransientProblem& problem,
                                                                   std::size_t region_value) {
    return summarize_history(problem._implementation->material_histories.at(region_value));
}

const TransientConservationSummary& TransientProblem::last_conservation_summary() const noexcept {
    return _implementation->last_conservation_summary;
}

InterfaceSummary rz::ProblemAccess::summarize_interface(const TransientProblem& problem, std::size_t contact_index,
                                                        const std::vector<double>& state) {
    return problem._implementation->spatial->summarize_interface(contact_index, state);
}

std::vector<ContactNodeSummary> rz::ProblemAccess::summarize_contact_nodes(const TransientProblem& problem,
                                                                           std::size_t contact_index,
                                                                           const std::vector<double>& state) {
    return problem._implementation->spatial->summarize_contact_nodes(contact_index, state);
}

std::vector<std::size_t> rz::ProblemAccess::contact_secondary_source_nodes(const TransientProblem& problem,
                                                                           std::size_t contact_index) {
    return problem._implementation->spatial->contact_secondary_source_nodes(contact_index);
}

std::size_t TransientProblem::dof_count() const noexcept {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->dof_count()
                                           : _implementation->spatial->dof_count();
}

std::size_t TransientProblem::contribution_count() const noexcept {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->contribution_count()
                                           : _implementation->spatial->contribution_count();
}

const std::vector<FieldDescriptor>& TransientProblem::field_layout() const noexcept {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->dof_map().field_layout()
                                           : _implementation->spatial->dof_map().field_layout();
}

const std::vector<DirichletCondition>& TransientProblem::dirichlet_conditions() const noexcept {
    return _implementation->is_cartesian() ? _implementation->cartesian_spatial->dirichlet_conditions()
                                           : _implementation->spatial->dirichlet_conditions();
}

void TransientProblem::validate_state(const std::vector<double>& state) const {
    require_active_time_step();
    if (_implementation->is_cartesian())
        _implementation->cartesian_spatial->validate_state(state);
    else
        _implementation->spatial->validate_state(state);
}

std::vector<std::size_t> TransientProblem::required_state_dofs(std::size_t contribution_begin,
                                                               std::size_t contribution_end) const {
    if (_implementation->is_cartesian())
        return NonlinearProblem::required_state_dofs(contribution_begin, contribution_end);
    return _implementation->spatial->required_state_dofs(contribution_begin, contribution_end);
}

void TransientProblem::validate_local_state(std::size_t contribution_begin, std::size_t contribution_end,
                                            const GlobalStateView& state) const {
    require_active_time_step();
    if (_implementation->is_cartesian()) {
        NonlinearProblem::validate_local_state(contribution_begin, contribution_end, state);
        return;
    }
    _implementation->spatial->validate_local_state(contribution_begin, contribution_end, state);
}

std::size_t TransientProblem::contribution_dof_count(std::size_t contribution_index) const {
    if (_implementation->is_cartesian())
        return _implementation->cartesian_spatial->contribution_dof_count(contribution_index);
    if (contribution_index >= contribution_count())
        throw std::out_of_range("TransientProblem contribution index is out of range");
    return local_dof_count;
}

void TransientProblem::fill_contribution_dofs(std::size_t contribution_index, std::vector<std::size_t>& dofs) const {
    if (_implementation->is_cartesian()) {
        _implementation->cartesian_spatial->fill_contribution_dofs(contribution_index, dofs);
        return;
    }
    const LocalDofs fixed = rz::ProblemAccess::contribution_dofs(*this, contribution_index);
    dofs.assign(fixed.begin(), fixed.end());
}

void TransientProblem::compute_contribution_residual(std::size_t contribution_index, const std::vector<double>& state,
                                                     std::vector<double>& residual) const {
    if (_implementation->is_cartesian()) {
        require_active_time_step();
        if (contribution_index < _implementation->cartesian_spatial->volume_contribution_count()) {
            const auto location = _implementation->cartesian_spatial->element_location(contribution_index);
            copy_hex8_residual(
                _implementation->cartesian_region_kernels[location.first].residual(
                    _implementation->cartesian_spatial->region_element_geometry(location.first, location.second),
                    hex8_local_values(state),
                    gather_hex8_state(*_implementation->cartesian_spatial, contribution_index,
                                      _implementation->committed_solution),
                    _implementation->active_time_step),
                residual);
        } else {
            copy_face_residual(
                _implementation->cartesian_spatial->boundary_residual(contribution_index, face_local_values(state)),
                residual);
        }
        return;
    }
    copy_rz_residual(rz::ProblemAccess::contribution_residual(*this, contribution_index, rz_local_values(state)),
                     residual);
}

void TransientProblem::compute_contribution_system(std::size_t contribution_index, const std::vector<double>& state,
                                                   std::vector<double>& residual, std::vector<double>& jacobian) const {
    if (_implementation->is_cartesian()) {
        require_active_time_step();
        if (contribution_index < _implementation->cartesian_spatial->volume_contribution_count()) {
            const auto location = _implementation->cartesian_spatial->element_location(contribution_index);
            copy_hex8_system(
                _implementation->cartesian_region_kernels[location.first].linearize(
                    _implementation->cartesian_spatial->region_element_geometry(location.first, location.second),
                    hex8_local_values(state),
                    gather_hex8_state(*_implementation->cartesian_spatial, contribution_index,
                                      _implementation->committed_solution),
                    _implementation->active_time_step),
                residual, jacobian);
        } else {
            copy_face_system(
                _implementation->cartesian_spatial->boundary_system(contribution_index, face_local_values(state)),
                residual, jacobian);
        }
        return;
    }
    copy_rz_system(rz::ProblemAccess::linearize_contribution(*this, contribution_index, rz_local_values(state)),
                   residual, jacobian);
}

LocalDofs rz::ProblemAccess::contribution_dofs(const TransientProblem& problem, std::size_t contribution_index) {
    return problem._implementation->spatial->contribution_dofs(contribution_index);
}

LocalValues rz::ProblemAccess::contribution_state(const TransientProblem& problem, std::size_t contribution_index,
                                                  const std::vector<double>& global_state) {
    return contribution_state(problem, contribution_index, GlobalStateView(global_state));
}

LocalValues rz::ProblemAccess::contribution_state(const TransientProblem& problem, std::size_t contribution_index,
                                                  const GlobalStateView& global_state) {
    if (global_state.global_size() != problem.dof_count())
        throw std::invalid_argument("TransientProblem contribution state has the wrong global size");
    const LocalDofs dofs = contribution_dofs(problem, contribution_index);
    LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local) result[local] = global_state.value(dofs[local]);
    return result;
}

LocalResidual rz::ProblemAccess::contribution_residual(const TransientProblem& problem, std::size_t contribution_index,
                                                       const LocalValues& state) {
    problem.require_active_time_step();
    if (contribution_index < problem._implementation->spatial->volume_contribution_count()) {
        const auto location = problem._implementation->spatial->element_location(contribution_index);
        return problem._implementation->region_kernels[location.first].residual(
            region_element_geometry(problem, location.first, location.second), state,
            contribution_state(problem, contribution_index, problem._implementation->committed_solution),
            problem._implementation->material_histories[location.first][location.second],
            problem._implementation->active_time_step);
    }
    return problem._implementation->spatial->contribution_residual(contribution_index, state);
}

LocalSystem rz::ProblemAccess::linearize_contribution(const TransientProblem& problem, std::size_t contribution_index,
                                                      const LocalValues& state) {
    problem.require_active_time_step();
    if (contribution_index < problem._implementation->spatial->volume_contribution_count()) {
        const auto location = problem._implementation->spatial->element_location(contribution_index);
        return problem._implementation->region_kernels[location.first].linearize(
            region_element_geometry(problem, location.first, location.second), state,
            contribution_state(problem, contribution_index, problem._implementation->committed_solution),
            problem._implementation->material_histories[location.first][location.second],
            problem._implementation->active_time_step);
    }
    return problem._implementation->spatial->linearize_contribution(contribution_index, state);
}

void TransientProblem::require_active_time_step() const {
    if (!_implementation->time_step_active)
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

const std::string& PiecewiseLinearTimeTable::name() const noexcept { return _name; }

const std::vector<double>& PiecewiseLinearTimeTable::times() const noexcept { return _times; }

const std::vector<double>& PiecewiseLinearTimeTable::values() const noexcept { return _values; }

double PiecewiseLinearTimeTable::value(double time) const {
    if (!std::isfinite(time) || time < 0.0)
        throw std::invalid_argument("Time-table evaluation time must be finite and nonnegative");
    if (time <= _times.front()) return _values.front();
    if (time >= _times.back()) return _values.back();
    const auto upper = std::upper_bound(_times.begin(), _times.end(), time);
    const std::size_t right = static_cast<std::size_t>(upper - _times.begin());
    const std::size_t left = right - 1;
    const double fraction = (time - _times[left]) / (_times[right] - _times[left]);
    return (1.0 - fraction) * _values[left] + fraction * _values[right];
}

} // namespace fuelsim
