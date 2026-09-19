#include "core/steady_problem.hpp"
#include "core/transient_problem.hpp"
#include "solver/solve_workflows.hpp"
#include "spatial_backend_common.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim {
class SpatialProblemStorage final {
  public:
    SpatialProblemStorage(const SpatialProblemStorage&) = delete;
    SpatialProblemStorage& operator=(const SpatialProblemStorage&) = delete;
    SpatialProblemStorage(SpatialProblemStorage&&) = delete;
    SpatialProblemStorage& operator=(SpatialProblemStorage&&) = delete;

    SpatialProblemStorage(SpatialDefinition definition,
        const UnstructuredPlaneQuad8Mesh& mesh,
        bool transient = false) {
        _backend = std::make_unique<PlaneBackend>(_state, std::move(definition), mesh, transient);
    }

    SpatialProblemStorage(SpatialDefinition definition, const UnstructuredBar2Mesh& mesh, bool transient = false) {
        _backend = std::make_unique<RadialBackend>(_state, std::move(definition), mesh, transient);
    }

    SpatialProblemStorage(SpatialDefinition definition, const UnstructuredQuad4Mesh& mesh, bool transient = false) {
        if (definition.physics == Physics::thermal)
            _backend = std::make_unique<ThermalBackend>(_state, std::move(definition), mesh, transient);
        else
            _backend = std::make_unique<Rz4Backend>(_state, std::move(definition), mesh, transient);
    }

    SpatialProblemStorage(SpatialDefinition definition, const UnstructuredQuad8Mesh& mesh, bool transient = false) {
        if (definition.physics == Physics::thermal)
            _backend = std::make_unique<ThermalBackend>(_state, std::move(definition), mesh, transient);
        else
            _backend = std::make_unique<Rz8Backend>(_state, std::move(definition), mesh, transient);
    }

    SpatialProblemStorage(SpatialDefinition definition, const UnstructuredHex8Mesh& mesh, bool transient = false) {
        if (definition.physics == Physics::thermal)
            _backend = std::make_unique<ThermalBackend>(_state, std::move(definition), mesh, transient);
        else
            _backend = std::make_unique<CartesianBackend>(_state, std::move(definition), mesh, transient);
    }

    SpatialProblemStorage(SpatialDefinition definition, const UnstructuredHex20Mesh& mesh, bool transient = false) {
        if (definition.physics == Physics::thermal)
            _backend = std::make_unique<ThermalBackend>(_state, std::move(definition), mesh, transient);
        else
            _backend = std::make_unique<CartesianBackend>(_state, std::move(definition), mesh, transient);
    }

    const spatial_detail::SpatialLayout& layout() const noexcept { return _backend->layout(); }

    void initialize_steady_strain_formulations() {
        _steady_strain_formulations.reserve(layout().region_count());
        for (const auto& region : layout().definition().regions)
            _steady_strain_formulations.push_back(region.strain_formulation);
    }

    void set_small_strain_predictor_active(bool active) {
        if (layout().definition().physics == Physics::thermal) {
            if (active)
                throw std::invalid_argument("Thermal physics does not use a mechanical strain predictor");
            return;
        }
        if (!_steady_strain_formulations.empty() && active && !layout().definition().contacts.empty())
            throw std::invalid_argument("small-strain steady predictor does not support contact");
        const bool has_finite_strain = std::any_of(_steady_strain_formulations.begin(),
            _steady_strain_formulations.end(),
            [](StrainFormulation formulation) { return formulation == StrainFormulation::finite; });
        if (active && !has_finite_strain)
            throw std::invalid_argument("small-strain steady predictor requires at least one finite-strain region");
        if (active && _steady_sparsity_patterns.empty()) {
            _steady_sparsity_patterns.resize(_backend->sparsity_contribution_count());
            for (std::size_t index = 0; index < _steady_sparsity_patterns.size(); ++index)
                _backend->sparsity_contribution_jacobian_pattern(index, _steady_sparsity_patterns[index]);
        }
        for (std::size_t region = 0; region < _steady_strain_formulations.size(); ++region)
            _backend->set_region_strain_formulation(region,
                active ? StrainFormulation::small : _steady_strain_formulations[region]);
    }

    void sparsity_contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
        if (!_steady_sparsity_patterns.empty())
            pattern = _steady_sparsity_patterns.at(index);
        else
            _backend->sparsity_contribution_jacobian_pattern(index, pattern);
    }

    SpatialTimeState _state;
    std::unique_ptr<SpatialBackend> _backend;
    std::vector<StrainFormulation> _steady_strain_formulations;
    std::vector<std::vector<unsigned char>> _steady_sparsity_patterns;
};

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredPlaneQuad8Mesh& source)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source)) {
    for (const auto& region : _impl->layout().definition().regions)
        if (region.material.functions->has_creep() || region.material.functions->has_plasticity()
            || region.strain_formulation != StrainFormulation::small)
            throw std::invalid_argument("Steady CPEG8T currently supports small-strain elasticity");
    _impl->initialize_steady_strain_formulations();
}

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredBar2Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh)) {
    for (const auto& region : _impl->layout().definition().regions)
        if (region.material.functions->has_creep() || region.material.functions->has_plasticity())
            throw std::invalid_argument(
                "Steady radial GPS supports elasticity; use transient execution for inelastic materials");
    _impl->initialize_steady_strain_formulations();
}

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh)) {
    _impl->initialize_steady_strain_formulations();
}

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, false)) {
    _impl->initialize_steady_strain_formulations();
}

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredHex20Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, false)) {
    _impl->initialize_steady_strain_formulations();
}

SteadyProblem::~SteadyProblem() = default;

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredQuad8Mesh& source)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source)) {
    for (const auto& region : _impl->layout().definition().regions)
        if (region.material.functions->has_creep() || region.material.functions->has_plasticity())
            throw std::invalid_argument(
                "Steady CAX8T supports elasticity; use transient execution for inelastic materials");
    _impl->initialize_steady_strain_formulations();
}

bool BackendAccess::uses_radial_gps(const SteadyProblem& problem) {
    return dynamic_cast<const RadialBackend*>(problem._impl->_backend.get()) != nullptr;
}

bool BackendAccess::uses_radial_gps(const TransientProblem& problem) {
    return dynamic_cast<const RadialBackend*>(problem._impl->_backend.get()) != nullptr;
}

const radial::SpatialAssembly& BackendAccess::radial_spatial(const SteadyProblem& problem) {
    return dynamic_cast<const RadialBackend&>(*problem._impl->_backend).spatial();
}

const radial::SpatialAssembly& BackendAccess::radial_spatial(const TransientProblem& problem) {
    return dynamic_cast<const RadialBackend&>(*problem._impl->_backend).spatial();
}

const std::vector<std::vector<Cax2tGpsMaterialHistory>>& BackendAccess::radial_material_histories(
    const TransientProblem& problem) {
    return dynamic_cast<const RadialBackend&>(*problem._impl->_backend).histories();
}

const plane::SpatialAssembly& BackendAccess::plane_spatial(const SteadyProblem& problem) {
    return dynamic_cast<const PlaneBackend&>(*problem._impl->_backend).spatial();
}

const rz8::SpatialAssembly& BackendAccess::quad8_spatial(const SteadyProblem& problem) {
    return dynamic_cast<const Rz8Backend&>(*problem._impl->_backend).spatial();
}

const rz8::SpatialAssembly& BackendAccess::quad8_spatial(const TransientProblem& problem) {
    return dynamic_cast<const Rz8Backend&>(*problem._impl->_backend).spatial();
}

bool BackendAccess::uses_quad8(const SteadyProblem& problem) {
    return dynamic_cast<const Rz8Backend*>(problem._impl->_backend.get()) != nullptr;
}

const std::vector<AxisymmetricRegionData>& BackendAccess::quad8_kernel_data(const SteadyProblem& problem) {
    return dynamic_cast<const Rz8Backend&>(*problem._impl->_backend).kernel_data();
}

const std::vector<std::vector<Quad8MaterialHistory>>& BackendAccess::quad8_material_histories(
    const TransientProblem& problem) {
    return dynamic_cast<const Rz8Backend&>(*problem._impl->_backend).histories();
}

const cartesian::SpatialAssembly& BackendAccess::cartesian_spatial(const SteadyProblem& problem) {
    return dynamic_cast<const CartesianBackend&>(*problem._impl->_backend).spatial();
}

rz::SteadyBackendView BackendAccess::steady(const SteadyProblem& problem) {
    const auto& backend = dynamic_cast<const Rz4Backend&>(*problem._impl->_backend);
    return {backend.spatial(), backend.kernel_data()};
}

bool SteadyProblem::uses_augmented_contact() const noexcept {
    return _impl->_backend->uses_augmented_contact();
}

AugmentedContactUpdate SteadyProblem::update_augmented_contact_multipliers(const std::vector<double>& state,
    std::size_t completed_updates) {
    return _impl->_backend->update_augmented_contact_multipliers(state, completed_updates);
}

void SteadyProblem::set_load_factor(double value) {
    _impl->_backend->set_load_factor(value);
}

double SteadyProblem::load_factor() const noexcept {
    return _impl->layout().load_factor();
}

void SteadyProblem::set_time(double value) {
    _impl->_backend->set_time(value);
}

void SteadyProblem::set_small_strain_predictor_active(bool active) {
    _impl->set_small_strain_predictor_active(active);
}

std::vector<double> SteadyProblem::initial_state() const {
    return _impl->layout().initial_state();
}

ProblemStateSnapshot SteadyProblem::capture_internal_state() const {
    auto histories = std::make_shared<const std::vector<std::vector<ContactPointHistory>>>(
        _impl->_backend->committed_contact_histories());
    return ProblemStateSnapshot(discretization_identity(), histories);
}

void SteadyProblem::restore_internal_state(const ProblemStateSnapshot& snapshot, const std::vector<double>& state) {
    if (snapshot.empty())
        throw std::invalid_argument("SteadyProblem cannot restore an empty internal-state snapshot");
    if (snapshot._owner != discretization_identity())
        throw std::invalid_argument("SteadyProblem cannot restore a snapshot from another problem");
    const auto histories =
        std::static_pointer_cast<const std::vector<std::vector<ContactPointHistory>>>(snapshot._state);
    _impl->_backend->restore_contact_state(state, *histories);
}

void SteadyProblem::commit_internal_state(const std::vector<double>& state) {
    _impl->_backend->commit_contact_state(state);
}

std::size_t SteadyProblem::dof_count() const noexcept {
    return _impl->layout().dof_count();
}

std::size_t SteadyProblem::contribution_count() const noexcept {
    return _impl->_backend->contribution_count();
}

std::size_t SteadyProblem::sparsity_contribution_count() const noexcept {
    return _impl->_backend->sparsity_contribution_count();
}

bool SteadyProblem::jacobian_sparsity_is_state_dependent() const noexcept {
    return _impl->_backend->jacobian_sparsity_is_state_dependent();
}

bool SteadyProblem::contribution_metadata_is_fixed() const noexcept {
    return _impl->_backend->contribution_metadata_is_fixed();
}

std::pair<std::size_t, std::size_t> SteadyProblem::contribution_partition(std::size_t partition,
    std::size_t partition_count) const {
    return _impl->_backend->contribution_partition(partition, partition_count);
}

const std::vector<FieldDescriptor>& SteadyProblem::field_layout() const noexcept {
    return _impl->layout().field_layout();
}

const std::vector<DirichletCondition>& SteadyProblem::dirichlet_conditions() const noexcept {
    return _impl->layout().dirichlet_conditions();
}

void SteadyProblem::validate_state(const std::vector<double>& state) const {
    _impl->_backend->validate_state(state);
}

std::vector<std::size_t> SteadyProblem::required_state_dofs(std::size_t first, std::size_t last) const {
    return _impl->_backend->required_state_dofs(first, last);
}

void SteadyProblem::validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const {
    _impl->_backend->validate_local_state(first, last, state);
}

void SteadyProblem::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->_backend->contribution_dofs(index, dofs);
}

void SteadyProblem::contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
    _impl->_backend->contribution_jacobian_pattern(index, pattern);
}

void SteadyProblem::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->_backend->sparsity_contribution_dofs(index, dofs);
}

void SteadyProblem::sparsity_contribution_jacobian_pattern(std::size_t index,
    std::vector<unsigned char>& pattern) const {
    _impl->sparsity_contribution_jacobian_pattern(index, pattern);
}

void SteadyProblem::compute_contribution(std::size_t index,
    const std::vector<double>& state,
    std::vector<double>& residual,
    std::vector<double>* jacobian) const {
    _impl->_backend->compute_contribution(index, state, residual, jacobian, false);
}

std::vector<double> TransientProblem::accumulate_contribution_conservation(const std::vector<double>& solution,
    TransientConservationSummary& summary,
    std::vector<double>* external_load_residual,
    std::size_t first,
    std::size_t last) const {
    std::vector<double> raw_residual(dof_count(), 0.0);
    if (external_load_residual != nullptr)
        external_load_residual->assign(dof_count(), 0.0);
    ContributionWorkspace workspace;
    for (std::size_t entry = first; entry < last; ++entry) {
        evaluate_contribution(entry, solution, workspace, false);
        const SpatialContributionType type = _impl->_backend->contribution_type(entry);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local) {
            const std::size_t dof = workspace.dofs[local];
            const double residual = workspace.residual[local];
            raw_residual[dof] += residual;
            if (field_layout()[field_index(dof)].category == FieldCategory::thermal) {
                if (type == SpatialContributionType::thermal_contact
                    || (_impl->_backend->combined_contact_contribution()
                        && type == SpatialContributionType::mechanical_contact))
                    summary.interface_heat_imbalance += residual;
                else if (type == SpatialContributionType::convection)
                    summary.convection_heat_rate += residual;
                else if (type == SpatialContributionType::heat_flux)
                    summary.surface_heat_input_rate -= residual;
                continue;
            }
            const double work = residual * (solution[dof] - _impl->_state.committed_solution[dof]);
            if (type == SpatialContributionType::volume)
                summary.internal_mechanical_work_increment += work;
            else if (type == SpatialContributionType::mechanical_contact)
                summary.contact_work_increment += work;
            else if (type == SpatialContributionType::pressure || type == SpatialContributionType::traction)
                summary.pressure_traction_work_increment -= work;
            if (external_load_residual != nullptr
                && (type == SpatialContributionType::pressure || type == SpatialContributionType::traction))
                (*external_load_residual)[dof] += residual;
        }
    }
    return raw_residual;
}

TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredPlaneQuad8Mesh& source)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source, true)) {
    initialize_committed_state();
}

bool TransientProblem::uses_plane_quad8() const noexcept {
    return dynamic_cast<const PlaneBackend*>(_impl->_backend.get()) != nullptr;
}

const plane::SpatialAssembly& BackendAccess::plane_spatial(const TransientProblem& problem) {
    return dynamic_cast<const PlaneBackend&>(*problem._impl->_backend).spatial();
}

TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredBar2Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, true)) {
    initialize_committed_state();
}

TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, true)) {
    initialize_committed_state();
}

TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, true)) {
    initialize_committed_state();
}

TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredHex20Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, true)) {
    initialize_committed_state();
}

TransientProblem::~TransientProblem() = default;

TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredQuad8Mesh& source)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source, true)) {
    initialize_committed_state();
}

void TransientProblem::initialize_committed_state() {
    apply_spatial_controls(0.0, 0.0);
    _impl->_state.committed_solution = _impl->layout().initial_state();
    _impl->_state.committed_raw_residual.assign(dof_count(), 0.0);
    _impl->_state.committed_external_load_residual.assign(dof_count(), 0.0);
    _impl->_backend->restore_contact_state(_impl->_state.committed_solution,
        _impl->_backend->committed_contact_histories());
}

bool TransientProblem::uses_quad8() const noexcept {
    return dynamic_cast<const Rz8Backend*>(_impl->_backend.get()) != nullptr;
}

bool TransientProblem::uses_radial_gps() const noexcept {
    return dynamic_cast<const RadialBackend*>(_impl->_backend.get()) != nullptr;
}

bool TransientProblem::is_cartesian_3d() const noexcept {
    return dynamic_cast<const CartesianBackend*>(_impl->_backend.get()) != nullptr;
}

const cartesian::SpatialAssembly& BackendAccess::cartesian_spatial(const TransientProblem& problem) {
    return dynamic_cast<const CartesianBackend&>(*problem._impl->_backend).spatial();
}

const std::vector<std::vector<CartesianMaterialHistory>>& BackendAccess::cartesian_material_histories(
    const TransientProblem& problem) {
    if (const auto* plane = dynamic_cast<const PlaneBackend*>(problem._impl->_backend.get()))
        return plane->histories();
    return dynamic_cast<const CartesianBackend&>(*problem._impl->_backend).histories();
}

rz::TransientBackendView BackendAccess::transient(const TransientProblem& problem) {
    const auto& backend = dynamic_cast<const Rz4Backend&>(*problem._impl->_backend);
    const auto& runtime = problem._impl->_state;
    return {backend.spatial(),
        backend.kernel_data(),
        backend.histories(),
        runtime.committed_solution,
        runtime.active_time_step,
        runtime.time_step_active,
        runtime.include_thermal_time_term};
}

const SpatialDefinition& TransientProblem::definition() const noexcept {
    return _impl->layout().definition();
}

std::vector<double> TransientProblem::initial_solution() const {
    return _impl->layout().initial_state();
}

const std::vector<double>& TransientProblem::committed_solution() const noexcept {
    return _impl->_state.committed_solution;
}

bool TransientProblem::has_previous_committed_solution() const noexcept {
    return !_impl->_state.previous_committed_solution.empty();
}

const std::vector<double>& TransientProblem::previous_committed_solution() const noexcept {
    return _impl->_state.previous_committed_solution;
}

double TransientProblem::previous_committed_time() const noexcept {
    return _impl->_state.previous_committed_time;
}

void TransientProblem::track_previous_committed_solution(bool enabled) {
    if (_impl->_state.time_step_active)
        throw std::logic_error("TransientProblem cannot change predictor tracking during an active time step");
    _impl->_state.track_previous_committed_solution = enabled;
    if (!enabled) {
        _impl->_state.previous_committed_solution.clear();
        _impl->_state.previous_committed_time = 0.0;
    }
}

double TransientProblem::committed_time() const noexcept {
    return _impl->_state.committed_time;
}

double TransientProblem::committed_load_factor() const noexcept {
    return _impl->_state.committed_load_factor;
}

bool TransientProblem::time_step_active() const noexcept {
    return _impl->_state.time_step_active;
}

std::vector<double> TransientProblem::time_events() const {
    std::vector<double> result;
    for (const PiecewiseLinearTimeTable& table : _impl->layout().definition().time_tables)
        result.insert(result.end(), table.times().begin(), table.times().end());
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

RegionStateSummary TransientProblem::summarize_region(std::size_t region) const {
    return _impl->_backend->summarize_region(region);
}

const std::vector<double>& BackendAccess::committed_raw_residual(const TransientProblem& problem) {
    return problem._impl->_state.committed_raw_residual;
}

TransientCommittedState BackendAccess::committed_state(const TransientProblem& problem) {
    const auto& runtime = problem._impl->_state;
    TransientCommittedState state;
    state.solution = runtime.committed_solution;
    state.previous_solution = runtime.previous_committed_solution;
    state.raw_residual = runtime.committed_raw_residual;
    state.external_load_residual = runtime.committed_external_load_residual;
    state.conservation = runtime.last_conservation_summary;
    state.time = runtime.committed_time;
    state.load_factor = runtime.committed_load_factor;
    state.previous_time = runtime.previous_committed_time;
    problem._impl->_backend->export_histories(state);
    return state;
}

void BackendAccess::restore_committed_state(TransientProblem& problem, TransientCommittedState state) {

    SpatialProblemStorage& storage = *problem._impl;
    if (storage._state.time_step_active)
        throw std::logic_error("TransientProblem cannot restore during an active time step");
    const bool previous_valid = state.previous_solution.empty()
                                    ? state.previous_time == 0.0
                                    : state.previous_solution.size() == problem.dof_count()
                                          && std::isfinite(state.previous_time) && state.previous_time >= 0.0
                                          && state.previous_time < state.time
                                          && std::all_of(
                                              state.previous_solution.begin(),
                                              state.previous_solution.end(),
                                              [](double value) { return std::isfinite(value); });
    if (state.solution.size() != problem.dof_count() || !previous_valid || !std::isfinite(state.time)
        || state.time < 0.0 || !std::isfinite(state.load_factor) || state.load_factor < 0.0
        || state.raw_residual.size() != problem.dof_count()
        || state.external_load_residual.size() != problem.dof_count()
        || !std::all_of(
            state.raw_residual.begin(),
            state.raw_residual.end(),
            [](double value) { return std::isfinite(value); })
        || !std::all_of(
            state.external_load_residual.begin(),
            state.external_load_residual.end(),
            [](double value) { return std::isfinite(value); }))
        throw std::invalid_argument("Transient committed state layout does not match the problem");
    storage._backend->restore_histories(state);
    storage._state.committed_solution = std::move(state.solution);
    storage._state.previous_committed_solution = std::move(state.previous_solution);
    storage._state.committed_raw_residual = std::move(state.raw_residual);
    storage._state.committed_external_load_residual = std::move(state.external_load_residual);
    storage._state.last_conservation_summary = state.conservation;
    storage._state.committed_time = state.time;
    storage._state.committed_load_factor = state.load_factor;
    storage._state.previous_committed_time = state.previous_time;
    problem.clear_active_time_step();
    problem.apply_spatial_controls(storage._state.committed_time, storage._state.committed_load_factor);
}

ProblemStateSnapshot TransientProblem::capture_state() const {
    if (_impl->_state.time_step_active)
        throw std::logic_error("TransientProblem cannot capture an active time step");
    auto state = std::make_shared<const TransientCommittedState>(BackendAccess::committed_state(*this));
    return ProblemStateSnapshot(discretization_identity(), state);
}

void TransientProblem::restore_state(const ProblemStateSnapshot& snapshot) {
    if (snapshot.empty())
        throw std::invalid_argument("TransientProblem cannot restore an empty state snapshot");
    if (snapshot._owner != discretization_identity())
        throw std::invalid_argument("TransientProblem cannot restore a snapshot from another problem");
    BackendAccess::restore_committed_state(*this,
        *std::static_pointer_cast<const TransientCommittedState>(snapshot._state));
}

namespace {
struct TimeErrorAccumulator final {
    double difference_squared = 0.0, solution_squared = 0.0;
    std::size_t count = 0;
};

struct MaterialTimeErrors final {
    TimeErrorAccumulator elastic, plastic, creep, equivalent_plastic, equivalent_creep, stress;
};

void accumulate_time_error(TimeErrorAccumulator& accumulator, double full_step, double two_half_steps) {
    const double difference = two_half_steps - full_step;
    accumulator.difference_squared += difference * difference;
    accumulator.solution_squared += two_half_steps * two_half_steps;
    ++accumulator.count;
}

double
normalized_time_error(const TimeErrorAccumulator& accumulator, double absolute_tolerance, double relative_tolerance) {
    if (accumulator.count == 0)
        return 0.0;
    const double denominator = absolute_tolerance * std::sqrt(static_cast<double>(accumulator.count))
                               + relative_tolerance * std::sqrt(accumulator.solution_squared);
    return std::sqrt(accumulator.difference_squared) / denominator;
}

void accumulate_time_error(TimeErrorAccumulator& accumulator,
    const double* full_step,
    const double* two_half_steps,
    std::size_t count) {
    for (std::size_t component = 0; component < count; ++component)
        accumulate_time_error(accumulator, full_step[component], two_half_steps[component]);
}

void accumulate_material_time_error(MaterialTimeErrors& errors,
    const MaterialPointState& full,
    const MaterialPointState& half) {
    accumulate_time_error(errors.elastic, full.elastic_strain.data(), half.elastic_strain.data(), 4);
    accumulate_time_error(errors.plastic, full.plastic_strain.data(), half.plastic_strain.data(), 4);
    accumulate_time_error(errors.creep, full.creep_strain.data(), half.creep_strain.data(), 4);
    const double full_stress[] = {full.stress.rr, full.stress.zz, full.stress.hoop, full.stress.rz};
    const double half_stress[] = {half.stress.rr, half.stress.zz, half.stress.hoop, half.stress.rz};
    accumulate_time_error(errors.stress, full_stress, half_stress, 4);
    accumulate_time_error(errors.equivalent_plastic, full.equivalent_plastic_strain, half.equivalent_plastic_strain);
    accumulate_time_error(errors.equivalent_creep, full.equivalent_creep_strain, half.equivalent_creep_strain);
}

void accumulate_material_time_error(MaterialTimeErrors& errors,
    const CartesianMaterialPointState& full,
    const CartesianMaterialPointState& half) {
    accumulate_time_error(errors.elastic, full.elastic_strain.data(), half.elastic_strain.data(), 6);
    accumulate_time_error(errors.plastic, full.plastic_strain.data(), half.plastic_strain.data(), 6);
    accumulate_time_error(errors.creep, full.creep_strain.data(), half.creep_strain.data(), 6);
    const auto full_stress = cartesian::components(full.stress), half_stress = cartesian::components(half.stress);
    accumulate_time_error(errors.stress, full_stress.data(), half_stress.data(), 6);
    accumulate_time_error(errors.equivalent_plastic, full.equivalent_plastic_strain, half.equivalent_plastic_strain);
    accumulate_time_error(errors.equivalent_creep, full.equivalent_creep_strain, half.equivalent_creep_strain);
}

void assign_material_time_errors(TransientTimeErrorEstimate& result,
    const MaterialTimeErrors& errors,
    const TransientTimeOptions& options) {
    const double strain = options.strain_history_time_absolute_tolerance,
                 relative = options.time_error_relative_tolerance;
    result.elastic_strain = normalized_time_error(errors.elastic, strain, relative);
    result.plastic_strain = normalized_time_error(errors.plastic, strain, relative);
    result.creep_strain = normalized_time_error(errors.creep, strain, relative);
    result.equivalent_plastic_strain = normalized_time_error(errors.equivalent_plastic, strain, relative);
    result.equivalent_creep_strain = normalized_time_error(errors.equivalent_creep, strain, relative);
    result.stress = normalized_time_error(errors.stress, options.stress_history_time_absolute_tolerance, relative);
    result.maximum = std::max({result.maximum,
        result.elastic_strain,
        result.plastic_strain,
        result.creep_strain,
        result.equivalent_plastic_strain,
        result.equivalent_creep_strain,
        result.stress});
}

TransientTimeErrorEstimate nodal_time_error(const TransientCommittedState& full,
    const TransientCommittedState& half,
    const std::vector<FieldDescriptor>& fields,
    const TransientTimeOptions& options) {
    if (full.solution.size() != half.solution.size())
        throw std::logic_error("step-doubling nodal layouts differ");
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

TransientConservationSummary combine_half_step_conservation(const TransientConservationSummary& first,
    const TransientConservationSummary& second) {
    TransientConservationSummary result;
    for (std::size_t index = 0; index < 6; ++index) {
        double TransientConservationSummary::* member = transient_conservation_fields[index].member;
        result.*member = 0.5 * (first.*member + second.*member);
    }
    result.global_thermal_balance = result.stored_heat_rate + result.convection_heat_rate
                                    + result.interface_heat_imbalance - result.generated_heat_rate
                                    - result.surface_heat_input_rate - result.dirichlet_heat_input_rate;
    const double thermal_scale = std::abs(result.generated_heat_rate) + std::abs(result.stored_heat_rate)
                                 + std::abs(result.convection_heat_rate) + std::abs(result.interface_heat_imbalance)
                                 + std::abs(result.surface_heat_input_rate)
                                 + std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0 ? std::abs(result.global_thermal_balance) / thermal_scale : 0.0;
    result.unconstrained_thermal_residual_l2 =
        std::max(first.unconstrained_thermal_residual_l2, second.unconstrained_thermal_residual_l2);
    for (std::size_t index = 9; index < 13; ++index) {
        double TransientConservationSummary::* member = transient_conservation_fields[index].member;
        result.*member = first.*member + second.*member;
    }
    result.body_force_work_increment = first.body_force_work_increment + second.body_force_work_increment;
    result.mechanical_work_balance = result.internal_mechanical_work_increment + result.contact_work_increment
                                     - result.pressure_traction_work_increment
                                     - result.dirichlet_reaction_work_increment - result.body_force_work_increment;
    const double mechanical_scale =
        std::abs(result.body_force_work_increment) + std::abs(result.internal_mechanical_work_increment)
        + std::abs(result.pressure_traction_work_increment) + std::abs(result.dirichlet_reaction_work_increment)
        + std::abs(result.contact_work_increment);
    result.relative_mechanical_work_balance =
        mechanical_scale > 0.0 ? std::abs(result.mechanical_work_balance) / mechanical_scale : 0.0;
    result.unconstrained_mechanical_residual_l2 =
        std::max(first.unconstrained_mechanical_residual_l2, second.unconstrained_mechanical_residual_l2);
    for (std::size_t index = 16; index < 22; ++index) {
        double TransientConservationSummary::* member = transient_conservation_fields[index].member;
        result.*member = first.*member + second.*member;
    }
    result.mechanical_hourglass_energy = second.mechanical_hourglass_energy;
    result.mechanical_hourglass_energy_change =
        first.mechanical_hourglass_energy_change + second.mechanical_hourglass_energy_change;
    return result;
}

TransientTimeErrorEstimate compare_rz_step_doubling_states(const TransientCommittedState& full_step,
    const TransientCommittedState& two_half_steps,
    const std::vector<FieldDescriptor>& fields,
    std::size_t expected_dof_count,
    const TransientTimeOptions& options,
    const SpatialDefinition& definition,
    const rz8::SpatialAssembly* quad8) {
    if (full_step.solution.size() != two_half_steps.solution.size() || full_step.solution.size() != expected_dof_count)
        throw std::logic_error("step-doubling nodal-state layouts differ");
    if (full_step.material_histories.size() != two_half_steps.material_histories.size())
        throw std::logic_error("step-doubling material-state region layouts differ");
    MaterialTimeErrors material;
    TimeErrorAccumulator contact_friction, contact_normal_multiplier;
    bool contact_state_mismatch = false;
    if (full_step.radial_material_histories.size() != two_half_steps.radial_material_histories.size())
        throw std::logic_error("Radial GPS step-doubling material region layouts differ");
    for (std::size_t region = 0; region < full_step.radial_material_histories.size(); ++region) {
        const auto& full_region = full_step.radial_material_histories[region];
        const auto& half_region = two_half_steps.radial_material_histories[region];
        if (full_region.size() != half_region.size())
            throw std::logic_error("Radial GPS step-doubling material element layouts differ");
        for (std::size_t element = 0; element < full_region.size(); ++element)
            for (std::size_t q = 0; q < cax2t_gps_material_point_count; ++q)
                accumulate_material_time_error(material, full_region[element][q], half_region[element][q]);
    }
    if (full_step.quad8_material_histories.size() != two_half_steps.quad8_material_histories.size())
        throw std::logic_error("CAX8T step-doubling material region layouts differ");
    for (std::size_t r = 0; r < full_step.quad8_material_histories.size(); ++r) {
        const auto &first = full_step.quad8_material_histories[r], &second = two_half_steps.quad8_material_histories[r];
        if (first.size() != second.size())
            throw std::logic_error("CAX8T step-doubling element layouts differ");
        for (std::size_t e = 0; e < first.size(); ++e)
            for (std::size_t q = 0; q < quad8->region_element_geometry(r, e).point_count; ++q)
                accumulate_material_time_error(material, first[e][q], second[e][q]);
    }
    for (std::size_t region = 0; region < full_step.material_histories.size(); ++region) {
        const auto &full_history = full_step.material_histories[region],
                   &half_history = two_half_steps.material_histories[region];
        if (full_history.size() != half_history.size())
            throw std::logic_error("step-doubling material-state element layouts differ");
        const std::size_t point_count =
            definition.regions.at(region).rz_element_formulation == RzElementFormulation::cax4rt ? 1 : 4;
        for (std::size_t element = 0; element < full_history.size(); ++element)
            for (std::size_t q = 0; q < point_count; ++q)
                accumulate_material_time_error(material, full_history[element][q], half_history[element][q]);
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
            if (!full_step.radial_material_histories.empty())
                accumulate_time_error(contact_friction, full.total_tangential_slip, half.total_tangential_slip);
            for (std::size_t component = 0; component < full.cartesian_elastic_tangential_slip.size(); ++component)
                accumulate_time_error(contact_friction,
                    full.cartesian_elastic_tangential_slip[component],
                    half.cartesian_elastic_tangential_slip[component]);
            for (std::size_t component = 0; component < full.cartesian_total_tangential_slip.size(); ++component)
                accumulate_time_error(contact_friction,
                    full.cartesian_total_tangential_slip[component],
                    half.cartesian_total_tangential_slip[component]);
            for (std::size_t component = 0; component < full.cartesian_contact_normal.size(); ++component) {
                accumulate_time_error(contact_friction,
                    full.cartesian_contact_normal[component],
                    half.cartesian_contact_normal[component]);
                accumulate_time_error(contact_friction,
                    full.cartesian_contact_tangent_first[component],
                    half.cartesian_contact_tangent_first[component]);
            }
            accumulate_time_error(contact_normal_multiplier, full.normal_multiplier, half.normal_multiplier);
            contact_state_mismatch =
                contact_state_mismatch || full.sliding != half.sliding
                || full.cartesian_tangent_basis_initialized != half.cartesian_tangent_basis_initialized;
        }
    }
    TransientTimeErrorEstimate result = nodal_time_error(full_step, two_half_steps, fields, options);
    assign_material_time_errors(result, material, options);
    result.contact_friction = contact_state_mismatch ? std::numeric_limits<double>::infinity()
                                                     : normalized_time_error(contact_friction,
                                                           options.displacement_time_absolute_tolerance,
                                                           options.time_error_relative_tolerance);
    result.contact_normal_multiplier = normalized_time_error(contact_normal_multiplier,
        options.stress_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.maximum = std::max({result.maximum, result.contact_friction, result.contact_normal_multiplier});
    for (const TransientFieldTimeError& field : result.nodal_fields)
        result.maximum = std::max(result.maximum, field.value);
    return result;
}
} // namespace

std::vector<double> TransientProblem::committed_creep_rates() const {
    if (time_step_active())
        throw std::logic_error("Creep rates require an accepted state");
    return _impl->_backend->committed_creep_rates();
}

TransientTimeErrorEstimate TransientProblem::step_doubling_error(const ProblemStateSnapshot& full_snapshot,
    const ProblemStateSnapshot& half_snapshot,
    const TransientTimeOptions& options) const {
    if (full_snapshot.empty() || half_snapshot.empty())
        throw std::invalid_argument("step-doubling requires two complete state snapshots");
    if (full_snapshot._owner != discretization_identity() || half_snapshot._owner != discretization_identity())
        throw std::invalid_argument("step-doubling snapshots belong to another problem");
    const TransientCommittedState& full =
        *std::static_pointer_cast<const TransientCommittedState>(full_snapshot._state);
    const TransientCommittedState& half =
        *std::static_pointer_cast<const TransientCommittedState>(half_snapshot._state);
    if (BackendAccess::uses_thermal(*this))
        return nodal_time_error(full, half, field_layout(), options);
    if (is_cartesian_3d() || uses_plane_quad8()) {
        if (full.solution.size() != dof_count() || half.solution.size() != dof_count())
            throw std::logic_error("Cartesian step-doubling snapshot layouts differ");
        TransientTimeErrorEstimate result = nodal_time_error(full, half, field_layout(), options);
        if (full.cartesian_material_histories.size() != half.cartesian_material_histories.size())
            throw std::logic_error("Cartesian step-doubling material region layouts differ");
        MaterialTimeErrors material;
        for (std::size_t region = 0; region < full.cartesian_material_histories.size(); ++region) {
            const auto &full_region = full.cartesian_material_histories[region],
                       &half_region = half.cartesian_material_histories[region];
            if (full_region.size() != half_region.size())
                throw std::logic_error("Cartesian step-doubling material element layouts differ");
            for (std::size_t element = 0; element < full_region.size(); ++element) {
                if (full_region[element].size() != half_region[element].size())
                    throw std::logic_error("Cartesian step-doubling integration-point layouts differ");
                for (std::size_t q = 0; q < full_region[element].size(); ++q)
                    accumulate_material_time_error(material, full_region[element][q], half_region[element][q]);
            }
        }
        assign_material_time_errors(result, material, options);
        return result;
    }
    return compare_rz_step_doubling_states(full,
        half,
        field_layout(),
        dof_count(),
        options,
        definition(),
        uses_quad8() ? &BackendAccess::quad8_spatial(*this) : nullptr);
}

void TransientProblem::combine_last_half_step_conservation(const TransientConservationSummary& first_half) {
    if (_impl->_state.time_step_active)
        throw std::logic_error("TransientProblem cannot combine conservation during an active time step");
    _impl->_state.last_conservation_summary =
        combine_half_step_conservation(first_half, _impl->_state.last_conservation_summary);
}

void TransientProblem::begin_time_step(const TransientStepInput& input) {
    if (_impl->_state.time_step_active)
        throw std::logic_error("TransientProblem already has an active time step");
    if (!std::isfinite(input.end_time) || input.end_time <= _impl->_state.committed_time)
        throw std::invalid_argument("TransientProblem end time must exceed committed time");
    if (!std::isfinite(input.load_factor) || input.load_factor < 0.0)
        throw std::invalid_argument("TransientProblem load factor must be finite and nonnegative");
    _impl->_state.active_time_step = input.end_time - _impl->_state.committed_time;
    _impl->_state.active_end_time = input.end_time;
    _impl->_state.active_load_factor = input.load_factor;
    _impl->_state.include_thermal_time_term = input.include_thermal_time_term;
    _impl->_backend->capture_active_contact_state();
    try {
        apply_spatial_controls(input.end_time, input.load_factor);
        _impl->_backend->set_heat_source_interval(_impl->_state.committed_time, input.end_time);
    } catch (...) {
        _impl->_backend->restore_active_contact_state();
        apply_spatial_controls(_impl->_state.committed_time, _impl->_state.committed_load_factor);
        clear_active_time_step();
        throw;
    }
    _impl->_state.time_step_active = true;
}

void TransientProblem::commit_time_step(const std::vector<double>& converged_solution) {
    commit_time_step(converged_solution, 0, contribution_count(), {});
}

void TransientProblem::commit_time_step(const std::vector<double>& converged_solution,
    std::size_t first_contribution,
    std::size_t last_contribution,
    const std::function<void(std::vector<double>&)>& sum_partitions) {

    require_active_time_step();
    if (first_contribution > last_contribution || last_contribution > contribution_count())
        throw std::invalid_argument("Transient commit contribution range is invalid");
    if (!sum_partitions && (first_contribution != 0 || last_contribution != contribution_count()))
        throw std::invalid_argument("Partial transient commit requires a partition sum");
    if (sum_partitions && uses_radial_gps())
        throw std::invalid_argument("Partitioned transient commit does not support radial GPS geometry");
    if (converged_solution.size() != dof_count())
        throw std::invalid_argument("TransientProblem committed solution size mismatch");
    if (!std::all_of(converged_solution.begin(), converged_solution.end(), [](double value) {
            return std::isfinite(value);
        }))
        throw std::domain_error("TransientProblem committed solution must be finite");
    for (const FieldDescriptor& field : field_layout())
        if (field.category == FieldCategory::thermal)
            for (std::size_t dof = field.begin; dof < field.end; ++dof)
                if (!(converged_solution[dof] > 0.0))
                    throw std::domain_error("TransientProblem committed temperatures must be positive");
    // Solver callbacks validate only the contributions owned by this MPI rank.
    // Committed diagnostics and histories are replicated and traverse all
    // contributions, including contact candidates outside that local partition.
    // Refresh their projections from the complete converged state before
    // assembling reactions or committing any history.
    TransientConservationSummary conservation;
    std::vector<double> external_load_residual(dof_count(), 0.0), raw_residual(dof_count(), 0.0);
    std::exception_ptr partition_failure;
    try {
        _impl->_backend->validate_state(converged_solution);
        raw_residual = accumulate_contribution_conservation(converged_solution,
            conservation,
            &external_load_residual,
            first_contribution,
            last_contribution);
    } catch (...) {
        if (!sum_partitions)
            throw;
        partition_failure = std::current_exception();
    }
    // Allocate the next nodal state before changing any committed backend state.
    auto accepted_solution = converged_solution;
    auto previous_solution =
        _impl->_state.track_previous_committed_solution ? _impl->_state.committed_solution : std::vector<double>{};
    _impl->_backend->prepare_time_step(*this,
        converged_solution,
        first_contribution,
        last_contribution,
        sum_partitions,
        conservation,
        raw_residual,
        external_load_residual,
        partition_failure);
    _impl->_backend->complete_contact_diagnostics(conservation);
    _impl->_backend->publish_contact_history();
    _impl->_backend->publish_material_history();
    _impl->_state.last_conservation_summary = conservation;
    _impl->_state.committed_raw_residual.swap(raw_residual);
    _impl->_state.committed_external_load_residual.swap(external_load_residual);
    if (_impl->_state.track_previous_committed_solution) {
        _impl->_state.previous_committed_solution.swap(previous_solution);
        _impl->_state.previous_committed_time = _impl->_state.committed_time;
    }
    _impl->_state.committed_solution.swap(accepted_solution);
    _impl->_state.committed_time = _impl->_state.active_end_time;
    _impl->_state.committed_load_factor = _impl->_state.active_load_factor;
    clear_active_time_step();
}

void TransientProblem::rollback_time_step() noexcept {
    if (!_impl->_state.time_step_active)
        return;
    apply_spatial_controls(_impl->_state.committed_time, _impl->_state.committed_load_factor);
    _impl->_backend->restore_active_contact_state();
    clear_active_time_step();
}

void TransientProblem::apply_spatial_controls(double time, double load_factor) {
    _impl->_backend->set_time(time);
    _impl->_backend->set_load_factor(load_factor);
}

void TransientProblem::clear_active_time_step() noexcept {
    _impl->_state.active_time_step = 0.0;
    _impl->_state.active_end_time = _impl->_state.committed_time;
    _impl->_state.active_load_factor = _impl->_state.committed_load_factor;
    _impl->_backend->clear_active_contact_state();
    _impl->_state.include_thermal_time_term = true;
    _impl->_state.time_step_active = false;
}

bool TransientProblem::uses_augmented_contact() const noexcept {
    return _impl->_backend->uses_augmented_contact();
}

AugmentedContactUpdate TransientProblem::update_augmented_contact_multipliers(const std::vector<double>& state,
    std::size_t completed_updates) {
    require_active_time_step();
    return _impl->_backend->update_augmented_contact_multipliers(state, completed_updates);
}

const TransientConservationSummary& TransientProblem::last_conservation_summary() const noexcept {
    return _impl->_state.last_conservation_summary;
}

std::size_t TransientProblem::dof_count() const noexcept {
    return _impl->layout().dof_count();
}

std::size_t TransientProblem::contribution_count() const noexcept {
    return _impl->_backend->contribution_count();
}

std::size_t TransientProblem::sparsity_contribution_count() const noexcept {
    return _impl->_backend->sparsity_contribution_count();
}

bool TransientProblem::jacobian_sparsity_is_state_dependent() const noexcept {
    return _impl->_backend->jacobian_sparsity_is_state_dependent();
}

bool TransientProblem::contribution_metadata_is_fixed() const noexcept {
    return _impl->_backend->contribution_metadata_is_fixed();
}

std::pair<std::size_t, std::size_t> TransientProblem::contribution_partition(std::size_t partition,
    std::size_t partition_count) const {
    return _impl->_backend->contribution_partition(partition, partition_count);
}

const std::vector<FieldDescriptor>& TransientProblem::field_layout() const noexcept {
    return _impl->layout().field_layout();
}

const std::vector<DirichletCondition>& TransientProblem::dirichlet_conditions() const noexcept {
    return _impl->layout().dirichlet_conditions();
}

void TransientProblem::validate_state(const std::vector<double>& state) const {
    require_active_time_step();
    _impl->_backend->validate_state(state);
}

std::vector<std::size_t> TransientProblem::required_state_dofs(std::size_t first, std::size_t last) const {
    return _impl->_backend->required_state_dofs(first, last);
}

void TransientProblem::validate_local_state(std::size_t first,
    std::size_t last,
    const std::vector<double>& state) const {
    require_active_time_step();
    _impl->_backend->validate_local_state(first, last, state);
}

void TransientProblem::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->_backend->contribution_dofs(index, dofs);
}

void TransientProblem::contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
    _impl->_backend->contribution_jacobian_pattern(index, pattern);
}

void TransientProblem::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->_backend->sparsity_contribution_dofs(index, dofs);
}

void TransientProblem::sparsity_contribution_jacobian_pattern(std::size_t index,
    std::vector<unsigned char>& pattern) const {
    _impl->sparsity_contribution_jacobian_pattern(index, pattern);
}

void TransientProblem::compute_contribution(std::size_t index,
    const std::vector<double>& state,
    std::vector<double>& residual,
    std::vector<double>* jacobian) const {
    require_active_time_step();
    _impl->_backend->compute_contribution(index, state, residual, jacobian, true);
}

void TransientProblem::require_active_time_step() const {
    if (!_impl->_state.time_step_active)
        throw std::logic_error("TransientProblem residual evaluation requires an active time step");
}

PiecewiseLinearTimeTable::PiecewiseLinearTimeTable(std::string name,
    std::vector<double> times,
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

double PiecewiseLinearTimeTable::value(double time) const {
    if (!std::isfinite(time) || time < 0.0)
        throw std::invalid_argument("Time-table evaluation time must be finite and nonnegative");
    if (time <= _times.front())
        return _values.front();
    if (time >= _times.back())
        return _values.back();
    const auto upper = std::upper_bound(_times.begin(), _times.end(), time);
    const std::size_t right = static_cast<std::size_t>(upper - _times.begin()), left = right - 1;
    const double fraction = (time - _times[left]) / (_times[right] - _times[left]);
    return (1.0 - fraction) * _values[left] + fraction * _values[right];
}

double PiecewiseLinearTimeTable::average_value(double begin_time, double end_time) const {
    if (!std::isfinite(begin_time) || !std::isfinite(end_time) || begin_time < 0.0 || !(end_time > begin_time))
        throw std::invalid_argument("Time-table averaging interval must be finite, nonnegative, and increasing");
    double integral = 0.0, left = begin_time;
    while (left < end_time) {
        const auto upper = std::upper_bound(_times.begin(), _times.end(), left);
        const double right = upper == _times.end() ? end_time : std::min(end_time, *upper);
        integral += 0.5 * (value(left) + value(right)) * (right - left);
        left = right;
    }
    return integral / (end_time - begin_time);
}
} // namespace fuelsim

namespace fuelsim {
bool BackendAccess::uses_thermal(const SteadyProblem& problem) {
    return dynamic_cast<const ThermalBackend*>(problem._impl->_backend.get()) != nullptr;
}

bool BackendAccess::uses_thermal(const TransientProblem& problem) {
    return dynamic_cast<const ThermalBackend*>(problem._impl->_backend.get()) != nullptr;
}

const thermal::SpatialAssembly& BackendAccess::thermal_spatial(const SteadyProblem& problem) {
    return dynamic_cast<const ThermalBackend&>(*problem._impl->_backend).spatial();
}

const thermal::SpatialAssembly& BackendAccess::thermal_spatial(const TransientProblem& problem) {
    return dynamic_cast<const ThermalBackend&>(*problem._impl->_backend).spatial();
}
} // namespace fuelsim
