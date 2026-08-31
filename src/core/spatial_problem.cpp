#include "cartesian3d_assembly.hpp"
#include "core/problem_backend_access.hpp"
#include "fuelsim/core/cartesian3d_hex8.hpp"
#include "fuelsim/core/nonlinear_problem.hpp"
#include "fuelsim/core/spatial_definition.hpp"
#include "fuelsim/core/steady_problem.hpp"
#include "fuelsim/core/transient_problem.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "rz_assembly.hpp"
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
                                    result.surface_heat_input_rate - result.dirichlet_heat_input_rate;
    const double thermal_scale = std::abs(result.stored_heat_rate) + std::abs(result.convection_heat_rate) +
                                 std::abs(result.interface_heat_imbalance) + std::abs(result.generated_heat_rate) +
                                 std::abs(result.surface_heat_input_rate) + std::abs(result.dirichlet_heat_input_rate);
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

void add_trapezoidal_external_work(const NonlinearProblem& problem, const std::vector<double>& current,
    const std::vector<double>& old, const std::vector<double>& current_raw_residual,
    const std::vector<double>& old_raw_residual, const std::vector<double>& current_external_load_residual,
    const std::vector<double>& old_external_load_residual, TransientConservationSummary& result) {
    if (current_raw_residual.size() != problem.dof_count() || old_raw_residual.size() != problem.dof_count() ||
        current_external_load_residual.size() != problem.dof_count() ||
        old_external_load_residual.size() != problem.dof_count())
        throw std::logic_error("Transient trapezoidal work residual layout does not match the problem");
    std::vector<bool> constrained(problem.dof_count(), false);
    for (const DirichletCondition& condition : problem.dirichlet_conditions()) constrained[condition.dof] = true;
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        if (problem.field_layout()[problem.field_index(dof)].category == FieldCategory::thermal) continue;
        const double increment = current[dof] - old[dof];
        result.trapezoidal_pressure_traction_work_increment -=
            0.5 * (old_external_load_residual[dof] + current_external_load_residual[dof]) * increment;
        if (constrained[dof])
            result.trapezoidal_dirichlet_reaction_work_increment +=
                0.5 * (old_raw_residual[dof] + current_raw_residual[dof]) * increment;
    }
}
} // namespace

namespace {
LocalValues rz_local_values(const std::vector<double>& values) {
    if (values.size() != local_dof_count) throw std::invalid_argument("RZ contribution state must contain 12 DOFs");
    LocalValues result{};
    std::copy(values.begin(), values.end(), result.begin());
    return result;
}

LocalValues gather_rz_state(
    const rz::SpatialAssembly& spatial, std::size_t index, const std::vector<double>& global_state) {
    const LocalDofs dofs = spatial.contribution_dofs(index);
    LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local) result[local] = global_state.at(dofs[local]);
    return result;
}
} // namespace

class SpatialProblemStorage {
  public:
    SpatialProblemStorage(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
        : rz(std::make_unique<rz::SpatialAssembly>(std::move(definition), source_mesh)) {
        kernel_data.reserve(rz->region_count());
        for (std::size_t region = 0; region < rz->region_count(); ++region) {
            const RegionDefinition& value = rz->region(region);
            kernel_data.push_back({IsotropicThermoelasticMaterial(value.material), rz->region_heat_source(region), 0.0,
                value.strain_formulation});
        }
    }

    SpatialProblemStorage(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh, bool transient)
        : cartesian(std::make_unique<cartesian::SpatialAssembly>(std::move(definition), source_mesh)) {
        if (!transient)
            for (std::size_t region = 0; region < cartesian->region_count(); ++region) {
                const MaterialFunctionSet& functions = *cartesian->region(region).material.functions;
                if (functions.has_creep() || functions.has_plasticity())
                    throw std::invalid_argument("Steady Cartesian three-dimensional problems support only elasticity");
            }
        if (transient) initialize_cartesian_histories();
    }

    SpatialProblemStorage(SpatialDefinition definition, const UnstructuredHex20Mesh& source_mesh, bool transient)
        : cartesian(std::make_unique<cartesian::SpatialAssembly>(std::move(definition), source_mesh)) {
        if (!transient)
            for (std::size_t region = 0; region < cartesian->region_count(); ++region) {
                const MaterialFunctionSet& functions = *cartesian->region(region).material.functions;
                if (functions.has_creep() || functions.has_plasticity())
                    throw std::invalid_argument("Steady Cartesian three-dimensional problems support only elasticity");
            }
        if (transient) initialize_cartesian_histories();
    }

    void initialize_cartesian_histories() {
        cartesian_material_histories.resize(cartesian->region_count());
        _staged_cartesian_material_histories.resize(cartesian->region_count());
        for (std::size_t region = 0; region < cartesian->region_count(); ++region) {
            const std::size_t points = cartesian->region_material_point_count(region);
            cartesian_material_histories[region].resize(cartesian->region_element_count(region));
            _staged_cartesian_material_histories[region].resize(cartesian->region_element_count(region));
            for (CartesianMaterialHistory& history : cartesian_material_histories[region]) history.resize(points);
            for (CartesianMaterialHistory& history : _staged_cartesian_material_histories[region])
                history.resize(points);
        }
    }

    bool is_cartesian() const noexcept { return cartesian != nullptr; }

    const spatial_detail::SpatialLayout& layout() const noexcept {
        return is_cartesian() ? static_cast<const spatial_detail::SpatialLayout&>(*cartesian) : *rz;
    }

    std::size_t contribution_count() const noexcept {
        return is_cartesian() ? cartesian->contribution_count() : rz->contribution_count();
    }

    std::size_t sparsity_contribution_count() const noexcept {
        return is_cartesian() ? cartesian->sparsity_contribution_count() : rz->sparsity_contribution_count();
    }

    std::pair<std::size_t, std::size_t> contribution_partition(
        std::size_t partition, std::size_t partition_count) const {
        if (is_cartesian()) return cartesian->contribution_partition(partition, partition_count);
        if (partition_count == 0 || partition >= partition_count)
            throw std::out_of_range("Spatial problem contribution partition is out of range");
        return {rz->contribution_count() * partition / partition_count,
            rz->contribution_count() * (partition + 1U) / partition_count};
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
        (void)problem;
        return is_cartesian() ? cartesian->required_state_dofs(first, last) : rz->required_state_dofs(first, last);
    }

    void validate_local_state(
        const NonlinearProblem& problem, std::size_t first, std::size_t last, const std::vector<double>& state) const {
        (void)problem;
        if (is_cartesian())
            cartesian->validate_local_state(first, last, state);
        else
            rz->validate_local_state(first, last, state);
    }

    const std::vector<std::vector<ContactPointHistory>>& committed_contact_histories() const noexcept {
        return is_cartesian() ? cartesian->committed_contact_histories() : rz->committed_contact_histories();
    }

    void commit_contact_state(const std::vector<double>& state) {
        if (is_cartesian())
            cartesian->commit_contact_state(state);
        else
            rz->commit_contact_state(state);
    }

    void restore_contact_state(
        const std::vector<double>& state, std::vector<std::vector<ContactPointHistory>> histories) {
        if (is_cartesian())
            cartesian->restore_contact_state(state, std::move(histories));
        else
            rz->restore_contact_state(state, std::move(histories));
    }

    void contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
        if (is_cartesian()) {
            cartesian->contribution_dofs(index, dofs);
            return;
        }
        const LocalDofs fixed = rz->contribution_dofs(index);
        dofs.assign(fixed.begin(), fixed.end());
    }

    void contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
        if (is_cartesian()) return cartesian->contribution_jacobian_pattern(index, pattern);
        const LocalDofs dofs = rz->contribution_dofs(index);
        pattern.assign(dofs.size() * dofs.size(), 1U);
    }

    void sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
        if (is_cartesian()) {
            cartesian->sparsity_contribution_dofs(index, dofs);
            return;
        }
        const LocalDofs fixed = rz->sparsity_contribution_dofs(index);
        dofs.assign(fixed.begin(), fixed.end());
    }

    void sparsity_contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
        if (is_cartesian()) return cartesian->sparsity_contribution_jacobian_pattern(index, pattern);
        const LocalDofs dofs = rz->sparsity_contribution_dofs(index);
        pattern.assign(dofs.size() * dofs.size(), 1U);
    }

    std::unique_ptr<rz::SpatialAssembly> rz;
    std::unique_ptr<cartesian::SpatialAssembly> cartesian;
    std::vector<Quad4RzData> kernel_data;
    std::vector<std::vector<Quad4MaterialHistory>> material_histories;
    std::vector<std::vector<Quad4MaterialHistory>> _staged_material_histories;
    std::vector<std::vector<CartesianMaterialHistory>> cartesian_material_histories;
    std::vector<std::vector<CartesianMaterialHistory>> _staged_cartesian_material_histories;
    TransientConservationSummary last_conservation_summary;
    std::vector<double> committed_solution, previous_committed_solution, committed_raw_residual,
        committed_external_load_residual;
    double committed_time = 0.0, committed_load_factor = 0.0, active_time_step = 0.0, active_end_time = 0.0,
           active_load_factor = 0.0, previous_committed_time = 0.0;
    std::vector<std::vector<ContactPointHistory>> active_contact_histories;
    bool time_step_active = false, include_thermal_time_term = true, track_previous_committed_solution = false;
};

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh)) {}

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, false)) {}

SteadyProblem::SteadyProblem(SpatialDefinition definition, const UnstructuredHex20Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, false)) {}

SteadyProblem::~SteadyProblem() = default;

const cartesian::SpatialAssembly& BackendAccess::cartesian_spatial(const SteadyProblem& problem) noexcept {
    return *problem._impl->cartesian;
}

rz::SteadyBackendView BackendAccess::steady(const SteadyProblem& problem) noexcept {
    return {*problem._impl->rz, problem._impl->kernel_data};
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
    if (_impl->is_cartesian()) return;
    for (std::size_t region = 0; region < _impl->rz->region_count(); ++region)
        _impl->kernel_data[region].volumetric_heat_source = _impl->rz->region_heat_source(region);
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
    for (Quad4RzData& data : _impl->kernel_data) data.time = value;
    for (std::size_t region = 0; region < _impl->rz->region_count(); ++region)
        _impl->kernel_data[region].volumetric_heat_source = _impl->rz->region_heat_source(region);
}

std::vector<double> SteadyProblem::initial_state() const { return _impl->layout().initial_state(); }

ProblemStateSnapshot SteadyProblem::capture_internal_state() const {
    auto histories =
        std::make_shared<const std::vector<std::vector<ContactPointHistory>>>(_impl->committed_contact_histories());
    return ProblemStateSnapshot(discretization_identity(), histories);
}

void SteadyProblem::restore_internal_state(const ProblemStateSnapshot& snapshot, const std::vector<double>& state) {
    if (snapshot.empty()) throw std::invalid_argument("SteadyProblem cannot restore an empty internal-state snapshot");
    if (snapshot._owner != discretization_identity())
        throw std::invalid_argument("SteadyProblem cannot restore a snapshot from another problem");
    const auto histories =
        std::static_pointer_cast<const std::vector<std::vector<ContactPointHistory>>>(snapshot._state);
    _impl->restore_contact_state(state, *histories);
}

void SteadyProblem::commit_internal_state(const std::vector<double>& state) { _impl->commit_contact_state(state); }

std::size_t SteadyProblem::dof_count() const noexcept { return _impl->layout().dof_count(); }

std::size_t SteadyProblem::contribution_count() const noexcept { return _impl->contribution_count(); }

std::size_t SteadyProblem::sparsity_contribution_count() const noexcept { return _impl->sparsity_contribution_count(); }

std::pair<std::size_t, std::size_t> SteadyProblem::contribution_partition(
    std::size_t partition, std::size_t partition_count) const {
    return _impl->contribution_partition(partition, partition_count);
}

const std::vector<FieldDescriptor>& SteadyProblem::field_layout() const noexcept {
    return _impl->layout().field_layout();
}

const std::vector<DirichletCondition>& SteadyProblem::dirichlet_conditions() const noexcept {
    return _impl->layout().dirichlet_conditions();
}

void SteadyProblem::validate_state(const std::vector<double>& state) const { _impl->validate_state(state); }

std::vector<std::size_t> SteadyProblem::required_state_dofs(std::size_t first, std::size_t last) const {
    return _impl->required_state_dofs(*this, first, last);
}

void SteadyProblem::validate_local_state(std::size_t first, std::size_t last, const std::vector<double>& state) const {
    _impl->validate_local_state(*this, first, last, state);
}

void SteadyProblem::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->contribution_dofs(index, dofs);
}

void SteadyProblem::contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
    _impl->contribution_jacobian_pattern(index, pattern);
}

void SteadyProblem::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->sparsity_contribution_dofs(index, dofs);
}

void SteadyProblem::sparsity_contribution_jacobian_pattern(
    std::size_t index, std::vector<unsigned char>& pattern) const {
    _impl->sparsity_contribution_jacobian_pattern(index, pattern);
}

void SteadyProblem::compute_contribution(std::size_t index, const std::vector<double>& state,
    std::vector<double>& residual, std::vector<double>* jacobian) const {
    if (_impl->is_cartesian()) {
        _impl->cartesian->compute_contribution(index, state, nullptr, nullptr, 0.0, residual, jacobian);
        return;
    }
    const LocalValues local_state = rz_local_values(state);
    LocalJacobian local_jacobian{};
    LocalResidual local_residual{};
    if (index >= _impl->rz->volume_contribution_count())
        local_residual =
            _impl->rz->compute_contribution(index, local_state, jacobian == nullptr ? nullptr : &local_jacobian);
    else {
        const auto location = _impl->rz->element_location(index);
        local_residual = compute_quad4_rz_thermoelastic(_impl->kernel_data[location.first],
            _impl->rz->region_element_geometry(location.first, location.second), local_state,
            jacobian == nullptr ? nullptr : &local_jacobian);
    }
    residual.assign(local_residual.begin(), local_residual.end());
    if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
}

std::vector<double> TransientProblem::accumulate_contribution_conservation(const std::vector<double>& solution,
    TransientConservationSummary& summary, std::vector<double>* external_load_residual) const {
    std::vector<double> raw_residual(dof_count(), 0.0);
    if (external_load_residual != nullptr) external_load_residual->assign(dof_count(), 0.0);
    ContributionWorkspace workspace;
    for (std::size_t entry = 0; entry < contribution_count(); ++entry) {
        evaluate_contribution(entry, solution, workspace, false);
        const SpatialContributionType type =
            _impl->is_cartesian() ? _impl->cartesian->contribution_type(entry) : _impl->rz->contribution_type(entry);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local) {
            const std::size_t dof = workspace.dofs[local];
            const double residual = workspace.residual[local];
            raw_residual[dof] += residual;
            if (field_layout()[field_index(dof)].category == FieldCategory::thermal) {
                if (type == SpatialContributionType::thermal_contact)
                    summary.interface_heat_imbalance += residual;
                else if (type == SpatialContributionType::convection)
                    summary.convection_heat_rate += residual;
                else if (type == SpatialContributionType::heat_flux)
                    summary.surface_heat_input_rate -= residual;
                continue;
            }
            const double work = residual * (solution[dof] - _impl->committed_solution[dof]);
            if (type == SpatialContributionType::volume)
                summary.internal_mechanical_work_increment += work;
            else if (type == SpatialContributionType::mechanical_contact)
                summary.contact_work_increment += work;
            else if (type == SpatialContributionType::pressure || type == SpatialContributionType::traction)
                summary.pressure_traction_work_increment -= work;
            if (external_load_residual != nullptr &&
                (type == SpatialContributionType::pressure || type == SpatialContributionType::traction))
                (*external_load_residual)[dof] += residual;
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

double trapezoidal_stress_strain_inner_product(const AxisymmetricStressValues& old_stress,
    const AxisymmetricStressValues& new_stress, const std::array<double, 4>& strain_increment) noexcept {
    return 0.5 * (stress_strain_inner_product(old_stress, strain_increment) +
                     stress_strain_inner_product(new_stress, strain_increment));
}
} // namespace
} // namespace rz

namespace cartesian {
namespace {
double stress_strain_inner_product(const SymmetricTensor3Values& stress, const std::array<double, 6>& strain) noexcept {
    return stress.xx * strain[0] + stress.yy * strain[1] + stress.zz * strain[2] +
           2.0 * (stress.xy * strain[3] + stress.yz * strain[4] + stress.xz * strain[5]);
}

std::array<double, 6> strain_difference(const std::array<double, 6>& current, const std::array<double, 6>& old) {
    std::array<double, 6> result{};
    for (std::size_t component = 0; component < result.size(); ++component)
        result[component] = current[component] - old[component];
    return result;
}

double trapezoidal_stress_strain_inner_product(const SymmetricTensor3Values& old_stress,
    const SymmetricTensor3Values& new_stress, const std::array<double, 6>& strain_increment) noexcept {
    return 0.5 * (stress_strain_inner_product(old_stress, strain_increment) +
                     stress_strain_inner_product(new_stress, strain_increment));
}

SymmetricTensor3Values rotate_tensor_values(const std::array<double, 6>& tensor, const CartesianRotation& rotation) {
    const SymmetricTensor3 rotated =
        rotate_cartesian_tensor({tensor[0], tensor[1], tensor[2], tensor[3], tensor[4], tensor[5]}, rotation);
    return {rotated.xx.value(), rotated.yy.value(), rotated.zz.value(), rotated.xy.value(), rotated.yz.value(),
        rotated.xz.value()};
}

SymmetricTensor3Values rotate_tensor_values(const SymmetricTensor3Values& tensor, const CartesianRotation& rotation) {
    const std::array<double, 6> values = {tensor.xx, tensor.yy, tensor.zz, tensor.xy, tensor.yz, tensor.xz};
    return rotate_tensor_values(values, rotation);
}

std::array<double, 6> components(const SymmetricTensor3Values& tensor) {
    return {tensor.xx, tensor.yy, tensor.zz, tensor.xy, tensor.yz, tensor.xz};
}
} // namespace
} // namespace cartesian

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

bool valid_material_state(const CartesianMaterialPointState& state) {
    for (std::size_t component = 0; component < state.elastic_strain.size(); ++component)
        if (!std::isfinite(state.elastic_strain[component]) || !std::isfinite(state.plastic_strain[component]) ||
            !std::isfinite(state.creep_strain[component]))
            return false;
    return std::isfinite(state.equivalent_plastic_strain) && state.equivalent_plastic_strain >= 0.0 &&
           std::isfinite(state.equivalent_creep_strain) && state.equivalent_creep_strain >= 0.0 &&
           std::isfinite(state.stress.xx) && std::isfinite(state.stress.yy) && std::isfinite(state.stress.zz) &&
           std::isfinite(state.stress.xy) && std::isfinite(state.stress.yz) && std::isfinite(state.stress.xz);
}
} // namespace

TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredQuad4Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh)) {
    const std::size_t regions = _impl->layout().definition().regions.size();
    _impl->material_histories.resize(regions);
    _impl->_staged_material_histories.resize(regions);
    for (std::size_t region = 0; region < regions; ++region) {
        _impl->material_histories[region].resize(_impl->rz->region_element_count(region));
        _impl->_staged_material_histories[region].resize(_impl->rz->region_element_count(region));
    }
    apply_spatial_controls(0.0, 0.0);
    _impl->committed_solution = _impl->rz->initial_state();
    _impl->committed_raw_residual.assign(dof_count(), 0.0);
    _impl->committed_external_load_residual.assign(dof_count(), 0.0);
    _impl->rz->restore_contact_state(_impl->committed_solution, _impl->rz->committed_contact_histories());
}

TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredHex8Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, true)) {
    apply_spatial_controls(0.0, 0.0);
    _impl->committed_solution = _impl->cartesian->initial_state();
    _impl->committed_raw_residual.assign(dof_count(), 0.0);
    _impl->committed_external_load_residual.assign(dof_count(), 0.0);
    _impl->cartesian->restore_contact_state(_impl->committed_solution, _impl->cartesian->committed_contact_histories());
}

TransientProblem::TransientProblem(SpatialDefinition definition, const UnstructuredHex20Mesh& source_mesh)
    : _impl(std::make_unique<SpatialProblemStorage>(std::move(definition), source_mesh, true)) {
    apply_spatial_controls(0.0, 0.0);
    _impl->committed_solution = _impl->cartesian->initial_state();
    _impl->committed_raw_residual.assign(dof_count(), 0.0);
    _impl->committed_external_load_residual.assign(dof_count(), 0.0);
    _impl->cartesian->restore_contact_state(_impl->committed_solution, _impl->cartesian->committed_contact_histories());
}

TransientProblem::~TransientProblem() = default;

bool TransientProblem::is_cartesian_3d() const noexcept { return _impl->is_cartesian(); }

const cartesian::SpatialAssembly& BackendAccess::cartesian_spatial(const TransientProblem& problem) noexcept {
    return *problem._impl->cartesian;
}

const std::vector<std::vector<CartesianMaterialHistory>>& BackendAccess::cartesian_material_histories(
    const TransientProblem& problem) noexcept {
    return problem._impl->cartesian_material_histories;
}

rz::TransientBackendView BackendAccess::transient(const TransientProblem& problem) noexcept {
    return {*problem._impl->rz, problem._impl->kernel_data, problem._impl->material_histories,
        problem._impl->committed_solution, problem._impl->active_time_step, problem._impl->time_step_active,
        problem._impl->include_thermal_time_term};
}

const SpatialDefinition& TransientProblem::definition() const noexcept { return _impl->layout().definition(); }

std::vector<double> TransientProblem::initial_solution() const { return _impl->layout().initial_state(); }

const std::vector<double>& TransientProblem::committed_solution() const noexcept { return _impl->committed_solution; }

bool TransientProblem::has_previous_committed_solution() const noexcept {
    return !_impl->previous_committed_solution.empty();
}

const std::vector<double>& TransientProblem::previous_committed_solution() const noexcept {
    return _impl->previous_committed_solution;
}

double TransientProblem::previous_committed_time() const noexcept { return _impl->previous_committed_time; }

void TransientProblem::track_previous_committed_solution(bool enabled) {
    if (_impl->time_step_active)
        throw std::logic_error("TransientProblem cannot change predictor tracking during an active time step");
    _impl->track_previous_committed_solution = enabled;
    if (!enabled) {
        _impl->previous_committed_solution.clear();
        _impl->previous_committed_time = 0.0;
    }
}

double TransientProblem::committed_time() const noexcept { return _impl->committed_time; }

double TransientProblem::committed_load_factor() const noexcept { return _impl->committed_load_factor; }

bool TransientProblem::time_step_active() const noexcept { return _impl->time_step_active; }

std::vector<double> TransientProblem::time_events() const {
    std::vector<double> result;
    for (const PiecewiseLinearTimeTable& table : _impl->layout().definition().time_tables)
        result.insert(result.end(), table.times().begin(), table.times().end());
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

RegionStateSummary TransientProblem::summarize_region(std::size_t region) const {
    if (region >= _impl->layout().definition().regions.size())
        throw std::out_of_range("TransientProblem region summary index is out of range");
    const std::size_t node_count =
        _impl->is_cartesian()
            ? (_impl->cartesian->uses_hex20() ? _impl->cartesian->hex20_region_mesh(region).nodes().size()
                                              : _impl->cartesian->region_mesh(region).nodes().size())
            : _impl->rz->region_mesh(region).nodes().size();
    const auto temperature = std::find_if(field_layout().begin(), field_layout().end(),
        [](const FieldDescriptor& field) { return field.category == FieldCategory::thermal; });
    RegionStateSummary result{-std::numeric_limits<double>::infinity(), 0.0, 0.0};
    for (std::size_t node = 0; node < node_count; ++node) {
        if (_impl->is_cartesian() && _impl->cartesian->uses_hex20() &&
            !_impl->cartesian->hex20_region_mesh(region).temperature_nodes().at(node))
            continue;
        result.maximum_temperature = std::max(result.maximum_temperature,
            _impl->committed_solution.at(temperature->begin + _impl->layout().global_temperature_node(region, node)));
    }
    if (_impl->is_cartesian()) {
        for (const CartesianMaterialHistory& element : _impl->cartesian_material_histories[region])
            for (const CartesianMaterialPointState& point : element) {
                result.maximum_equivalent_plastic_strain =
                    std::max(result.maximum_equivalent_plastic_strain, point.equivalent_plastic_strain);
                result.maximum_equivalent_creep_strain =
                    std::max(result.maximum_equivalent_creep_strain, point.equivalent_creep_strain);
            }
        return result;
    }
    for (const Quad4MaterialHistory& element : _impl->material_histories[region])
        for (const MaterialPointState& point : element) {
            result.maximum_equivalent_plastic_strain =
                std::max(result.maximum_equivalent_plastic_strain, point.equivalent_plastic_strain);
            result.maximum_equivalent_creep_strain =
                std::max(result.maximum_equivalent_creep_strain, point.equivalent_creep_strain);
        }
    return result;
}

TransientCommittedState BackendAccess::committed_state(const TransientProblem& problem) {
    const SpatialProblemStorage& storage = *problem._impl;
    return {storage.committed_solution, storage.previous_committed_solution, storage.material_histories,
        storage.cartesian_material_histories, storage.committed_contact_histories(), storage.committed_raw_residual,
        storage.committed_external_load_residual, storage.last_conservation_summary, storage.committed_time,
        storage.committed_load_factor, storage.previous_committed_time};
}

void BackendAccess::restore_committed_state(TransientProblem& problem, TransientCommittedState state) {
    SpatialProblemStorage& storage = *problem._impl;
    if (storage.time_step_active) throw std::logic_error("TransientProblem cannot restore during an active time step");
    const bool previous_valid = state.previous_solution.empty()
                                    ? state.previous_time == 0.0
                                    : state.previous_solution.size() == problem.dof_count() &&
                                          std::isfinite(state.previous_time) && state.previous_time >= 0.0 &&
                                          state.previous_time < state.time &&
                                          std::all_of(
                                              state.previous_solution.begin(), state.previous_solution.end(),
                                              [](double value) { return std::isfinite(value); });
    if (state.solution.size() != problem.dof_count() || !previous_valid || !std::isfinite(state.time) ||
        state.time < 0.0 || !std::isfinite(state.load_factor) || state.load_factor < 0.0 ||
        state.raw_residual.size() != problem.dof_count() ||
        state.external_load_residual.size() != problem.dof_count() ||
        !std::all_of(
            state.raw_residual.begin(), state.raw_residual.end(), [](double value) { return std::isfinite(value); }) ||
        !std::all_of(
            state.external_load_residual.begin(), state.external_load_residual.end(),
            [](double value) { return std::isfinite(value); }))
        throw std::invalid_argument("Transient committed state layout does not match the problem");
    if (storage.is_cartesian()) {
        storage.cartesian->validate_state(state.solution);
        if (!state.material_histories.empty() ||
            state.contact_histories.size() != storage.layout().definition().contacts.size() ||
            state.cartesian_material_histories.size() != storage.cartesian->region_count())
            throw std::invalid_argument("Cartesian transient committed state layout does not match the problem");
        for (std::size_t region = 0; region < storage.cartesian->region_count(); ++region) {
            if (state.cartesian_material_histories[region].size() != storage.cartesian->region_element_count(region))
                throw std::invalid_argument("Cartesian committed element state layout does not match");
            const std::size_t expected_points = storage.cartesian->region_material_point_count(region);
            for (const CartesianMaterialHistory& element : state.cartesian_material_histories[region]) {
                if (element.size() != expected_points)
                    throw std::invalid_argument("Cartesian committed integration-point layout does not match");
                for (const CartesianMaterialPointState& point : element)
                    if (!valid_material_state(point))
                        throw std::invalid_argument("Cartesian committed material state is invalid");
            }
        }
        storage.cartesian->restore_contact_state(state.solution, std::move(state.contact_histories));
        storage.committed_solution = std::move(state.solution);
        storage.previous_committed_solution = std::move(state.previous_solution);
        storage.committed_raw_residual = std::move(state.raw_residual);
        storage.committed_external_load_residual = std::move(state.external_load_residual);
        storage.cartesian_material_histories = std::move(state.cartesian_material_histories);
        storage.last_conservation_summary = state.conservation;
        storage.committed_time = state.time;
        storage.committed_load_factor = state.load_factor;
        storage.previous_committed_time = state.previous_time;
        problem.clear_active_time_step();
        problem.apply_spatial_controls(state.time, state.load_factor);
        return;
    }
    if (state.material_histories.size() != storage.rz->region_count() || !state.cartesian_material_histories.empty() ||
        state.contact_histories.size() != storage.layout().definition().contacts.size())
        throw std::invalid_argument("Transient committed state layout does not match the problem");
    const spatial_detail::SpatialLayout& dofs = *storage.rz;
    for (std::size_t node = 0; node < dofs.node_count(); ++node) {
        const double temperature = state.solution.at(dofs.dof(Field::temperature, node));
        if (!std::isfinite(temperature) || !(temperature > 0.0) ||
            !std::isfinite(state.solution.at(dofs.dof(Field::radial_displacement, node))) ||
            !std::isfinite(state.solution.at(dofs.dof(Field::axial_displacement, node))))
            throw std::invalid_argument("Transient committed nodal state must be finite with positive temperatures");
    }
    for (std::size_t region = 0; region < storage.rz->region_count(); ++region) {
        const std::size_t elements = storage.rz->region_element_count(region);
        if (state.material_histories[region].size() != elements)
            throw std::invalid_argument("Transient committed element state layout does not match");
        for (std::size_t element = 0; element < elements; ++element) {
            for (std::size_t q = 0; q < 4; ++q)
                if (!valid_material_state(state.material_histories[region][element][q]) ||
                    !finite_stress(state.material_histories[region][element][q].stress))
                    throw std::invalid_argument("Transient committed integration-point state is invalid");
        }
    }
    storage.rz->restore_contact_state(state.solution, state.contact_histories);
    storage.committed_solution = std::move(state.solution);
    storage.previous_committed_solution = std::move(state.previous_solution);
    storage.committed_raw_residual = std::move(state.raw_residual);
    storage.committed_external_load_residual = std::move(state.external_load_residual);
    storage.material_histories = std::move(state.material_histories);
    storage.last_conservation_summary = state.conservation;
    storage.committed_time = state.time;
    storage.committed_load_factor = state.load_factor;
    storage.previous_committed_time = state.previous_time;
    problem.clear_active_time_step();
    problem.apply_spatial_controls(storage.committed_time, storage.committed_load_factor);
}

ProblemStateSnapshot TransientProblem::capture_state() const {
    if (_impl->time_step_active) throw std::logic_error("TransientProblem cannot capture an active time step");
    auto state = std::make_shared<const TransientCommittedState>(BackendAccess::committed_state(*this));
    return ProblemStateSnapshot(discretization_identity(), state);
}

void TransientProblem::restore_state(const ProblemStateSnapshot& snapshot) {
    if (snapshot.empty()) throw std::invalid_argument("TransientProblem cannot restore an empty state snapshot");
    if (snapshot._owner != discretization_identity())
        throw std::invalid_argument("TransientProblem cannot restore a snapshot from another problem");
    BackendAccess::restore_committed_state(
        *this, *std::static_pointer_cast<const TransientCommittedState>(snapshot._state));
}

namespace rz {
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

double normalized_time_error(
    const TimeErrorAccumulator& accumulator, double absolute_tolerance, double relative_tolerance) {
    if (accumulator.count == 0) return 0.0;
    const double denominator = absolute_tolerance * std::sqrt(static_cast<double>(accumulator.count)) +
                               relative_tolerance * std::sqrt(accumulator.solution_squared);
    return std::sqrt(accumulator.difference_squared) / denominator;
}

void accumulate_material_time_error(MaterialTimeErrors& errors, const double* full_elastic, const double* half_elastic,
    const double* full_plastic, const double* half_plastic, const double* full_creep, const double* half_creep,
    const double* full_stress, const double* half_stress, std::size_t count, double full_plastic_equivalent,
    double half_plastic_equivalent, double full_creep_equivalent, double half_creep_equivalent) {
    for (std::size_t component = 0; component < count; ++component) {
        accumulate_time_error(errors.elastic, full_elastic[component], half_elastic[component]);
        accumulate_time_error(errors.plastic, full_plastic[component], half_plastic[component]);
        accumulate_time_error(errors.creep, full_creep[component], half_creep[component]);
        accumulate_time_error(errors.stress, full_stress[component], half_stress[component]);
    }
    accumulate_time_error(errors.equivalent_plastic, full_plastic_equivalent, half_plastic_equivalent);
    accumulate_time_error(errors.equivalent_creep, full_creep_equivalent, half_creep_equivalent);
}

void assign_material_time_errors(
    TransientTimeErrorEstimate& result, const MaterialTimeErrors& errors, const TransientTimeOptions& options) {
    const double strain = options.strain_history_time_absolute_tolerance,
                 relative = options.time_error_relative_tolerance;
    result.elastic_strain = normalized_time_error(errors.elastic, strain, relative);
    result.plastic_strain = normalized_time_error(errors.plastic, strain, relative);
    result.creep_strain = normalized_time_error(errors.creep, strain, relative);
    result.equivalent_plastic_strain = normalized_time_error(errors.equivalent_plastic, strain, relative);
    result.equivalent_creep_strain = normalized_time_error(errors.equivalent_creep, strain, relative);
    result.stress = normalized_time_error(errors.stress, options.stress_history_time_absolute_tolerance, relative);
    result.maximum = std::max({result.maximum, result.elastic_strain, result.plastic_strain, result.creep_strain,
        result.equivalent_plastic_strain, result.equivalent_creep_strain, result.stress});
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
    for (std::size_t index = 0; index < 6; ++index) {
        double TransientConservationSummary::* member = transient_conservation_fields[index].member;
        result.*member = 0.5 * (first.*member + second.*member);
    }
    result.global_thermal_balance = result.stored_heat_rate + result.convection_heat_rate +
                                    result.interface_heat_imbalance - result.generated_heat_rate -
                                    result.surface_heat_input_rate - result.dirichlet_heat_input_rate;
    const double thermal_scale = std::abs(result.generated_heat_rate) + std::abs(result.stored_heat_rate) +
                                 std::abs(result.convection_heat_rate) + std::abs(result.interface_heat_imbalance) +
                                 std::abs(result.surface_heat_input_rate) + std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0 ? std::abs(result.global_thermal_balance) / thermal_scale : 0.0;
    result.unconstrained_thermal_residual_l2 =
        std::max(first.unconstrained_thermal_residual_l2, second.unconstrained_thermal_residual_l2);
    for (std::size_t index = 9; index < 13; ++index) {
        double TransientConservationSummary::* member = transient_conservation_fields[index].member;
        result.*member = first.*member + second.*member;
    }
    result.mechanical_work_balance = result.internal_mechanical_work_increment + result.contact_work_increment -
                                     result.pressure_traction_work_increment - result.dirichlet_reaction_work_increment;
    const double mechanical_scale =
        std::abs(result.internal_mechanical_work_increment) + std::abs(result.pressure_traction_work_increment) +
        std::abs(result.dirichlet_reaction_work_increment) + std::abs(result.contact_work_increment);
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

TransientTimeErrorEstimate compare_step_doubling_states(const TransientCommittedState& full_step,
    const TransientCommittedState& two_half_steps, const std::vector<FieldDescriptor>& fields,
    std::size_t expected_dof_count, const TransientTimeOptions& options) {
    if (full_step.solution.size() != two_half_steps.solution.size() || full_step.solution.size() != expected_dof_count)
        throw std::logic_error("step-doubling nodal-state layouts differ");
    if (full_step.material_histories.size() != two_half_steps.material_histories.size())
        throw std::logic_error("step-doubling material-state region layouts differ");
    MaterialTimeErrors material;
    TimeErrorAccumulator contact_friction, contact_normal_multiplier;
    bool contact_state_mismatch = false;
    for (std::size_t region = 0; region < full_step.material_histories.size(); ++region) {
        const auto &full_history = full_step.material_histories[region],
                   &half_history = two_half_steps.material_histories[region];
        if (full_history.size() != half_history.size())
            throw std::logic_error("step-doubling material-state element layouts differ");
        for (std::size_t element = 0; element < full_history.size(); ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const MaterialPointState &full_point = full_history[element][q], &half_point = half_history[element][q];
                const AxisymmetricStressValues &full_value = full_point.stress, &half_value = half_point.stress;
                const double full_stress[] = {full_value.rr, full_value.zz, full_value.hoop, full_value.rz},
                             half_stress[] = {half_value.rr, half_value.zz, half_value.hoop, half_value.rz};
                accumulate_material_time_error(material, full_point.elastic_strain.data(),
                    half_point.elastic_strain.data(), full_point.plastic_strain.data(),
                    half_point.plastic_strain.data(), full_point.creep_strain.data(), half_point.creep_strain.data(),
                    full_stress, half_stress, 4, full_point.equivalent_plastic_strain,
                    half_point.equivalent_plastic_strain, full_point.equivalent_creep_strain,
                    half_point.equivalent_creep_strain);
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
            for (std::size_t component = 0; component < full.cartesian_elastic_tangential_slip.size(); ++component)
                accumulate_time_error(contact_friction, full.cartesian_elastic_tangential_slip[component],
                    half.cartesian_elastic_tangential_slip[component]);
            for (std::size_t component = 0; component < full.cartesian_total_tangential_slip.size(); ++component)
                accumulate_time_error(contact_friction, full.cartesian_total_tangential_slip[component],
                    half.cartesian_total_tangential_slip[component]);
            for (std::size_t component = 0; component < full.cartesian_contact_normal.size(); ++component) {
                accumulate_time_error(contact_friction, full.cartesian_contact_normal[component],
                    half.cartesian_contact_normal[component]);
                accumulate_time_error(contact_friction, full.cartesian_contact_tangent_first[component],
                    half.cartesian_contact_tangent_first[component]);
            }
            accumulate_time_error(contact_normal_multiplier, full.normal_multiplier, half.normal_multiplier);
            contact_state_mismatch =
                contact_state_mismatch || full.sliding != half.sliding ||
                full.cartesian_tangent_basis_initialized != half.cartesian_tangent_basis_initialized;
        }
    }
    TransientTimeErrorEstimate result = nodal_time_error(full_step, two_half_steps, fields, options);
    assign_material_time_errors(result, material, options);
    result.contact_friction =
        contact_state_mismatch ? std::numeric_limits<double>::infinity()
                               : normalized_time_error(contact_friction, options.displacement_time_absolute_tolerance,
                                     options.time_error_relative_tolerance);
    result.contact_normal_multiplier = normalized_time_error(contact_normal_multiplier,
        options.stress_history_time_absolute_tolerance, options.time_error_relative_tolerance);
    result.maximum = std::max({result.maximum, result.contact_friction, result.contact_normal_multiplier});
    for (const TransientFieldTimeError& field : result.nodal_fields)
        result.maximum = std::max(result.maximum, field.value);
    return result;
}
} // namespace
} // namespace rz

TransientTimeErrorEstimate TransientProblem::step_doubling_error(const ProblemStateSnapshot& full_snapshot,
    const ProblemStateSnapshot& half_snapshot, const TransientTimeOptions& options) const {
    if (full_snapshot.empty() || half_snapshot.empty())
        throw std::invalid_argument("step-doubling requires two complete state snapshots");
    if (full_snapshot._owner != discretization_identity() || half_snapshot._owner != discretization_identity())
        throw std::invalid_argument("step-doubling snapshots belong to another problem");
    const TransientCommittedState& full =
        *std::static_pointer_cast<const TransientCommittedState>(full_snapshot._state);
    const TransientCommittedState& half =
        *std::static_pointer_cast<const TransientCommittedState>(half_snapshot._state);
    if (_impl->is_cartesian()) {
        if (full.solution.size() != dof_count() || half.solution.size() != dof_count())
            throw std::logic_error("Cartesian step-doubling snapshot layouts differ");
        TransientTimeErrorEstimate result = rz::nodal_time_error(full, half, field_layout(), options);
        if (full.cartesian_material_histories.size() != half.cartesian_material_histories.size())
            throw std::logic_error("Cartesian step-doubling material region layouts differ");
        rz::MaterialTimeErrors material;
        for (std::size_t region = 0; region < full.cartesian_material_histories.size(); ++region) {
            const auto &full_region = full.cartesian_material_histories[region],
                       &half_region = half.cartesian_material_histories[region];
            if (full_region.size() != half_region.size())
                throw std::logic_error("Cartesian step-doubling material element layouts differ");
            for (std::size_t element = 0; element < full_region.size(); ++element) {
                if (full_region[element].size() != half_region[element].size())
                    throw std::logic_error("Cartesian step-doubling integration-point layouts differ");
                for (std::size_t q = 0; q < full_region[element].size(); ++q) {
                    const CartesianMaterialPointState &first = full_region[element][q],
                                                      &second = half_region[element][q];
                    const std::array<double, 6> first_stress = {first.stress.xx, first.stress.yy, first.stress.zz,
                        first.stress.xy, first.stress.yz, first.stress.xz};
                    const std::array<double, 6> second_stress = {second.stress.xx, second.stress.yy, second.stress.zz,
                        second.stress.xy, second.stress.yz, second.stress.xz};
                    rz::accumulate_material_time_error(material, first.elastic_strain.data(),
                        second.elastic_strain.data(), first.plastic_strain.data(), second.plastic_strain.data(),
                        first.creep_strain.data(), second.creep_strain.data(), first_stress.data(),
                        second_stress.data(), 6, first.equivalent_plastic_strain, second.equivalent_plastic_strain,
                        first.equivalent_creep_strain, second.equivalent_creep_strain);
                }
            }
        }
        rz::assign_material_time_errors(result, material, options);
        return result;
    }
    return rz::compare_step_doubling_states(full, half, field_layout(), dof_count(), options);
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
    _impl->include_thermal_time_term = input.include_thermal_time_term;
    if (uses_augmented_contact()) _impl->active_contact_histories = _impl->committed_contact_histories();
    try {
        apply_spatial_controls(input.end_time, input.load_factor);
    } catch (...) {
        if (!_impl->active_contact_histories.empty())
            _impl->restore_contact_state(_impl->committed_solution, std::move(_impl->active_contact_histories));
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
        std::vector<double> external_load_residual;
        const std::vector<double> raw_residual =
            accumulate_contribution_conservation(converged_solution, conservation, &external_load_residual);
        auto& staged = _impl->_staged_cartesian_material_histories;
        for (std::size_t region = 0; region < _impl->cartesian->region_count(); ++region) {
            const std::size_t offset = _impl->cartesian->region_element_offset(region);
            for (std::size_t element = 0; element < _impl->cartesian->region_element_count(region); ++element) {
                if (_impl->cartesian->uses_hex20()) {
                    const Hex20LocalValues current =
                        _impl->cartesian->hex20_volume_state(offset + element, converged_solution);
                    const Hex20LocalValues old =
                        _impl->cartesian->hex20_volume_state(offset + element, _impl->committed_solution);
                    const Hex20Geometry& geometry = _impl->cartesian->hex20_region_element_geometry(region, element);
                    CartesianMaterialHistory update = _impl->cartesian->transient_update(region, element, current, old,
                        _impl->cartesian_material_histories[region][element], _impl->active_time_step);
                    for (const Hex20ThermalQuadraturePoint& point : geometry.thermal_points) {
                        double current_temperature = 0.0, old_temperature = 0.0;
                        for (std::size_t node = 0; node < hex20_temperature_node_count; ++node) {
                            current_temperature += point.temperature_shape[node] * current[node];
                            old_temperature += point.temperature_shape[node] * old[node];
                        }
                        if (_impl->include_thermal_time_term)
                            conservation.stored_heat_rate +=
                                point.weighted_measure *
                                _impl->cartesian->heat_capacity(region, current_temperature, point.position) *
                                (current_temperature - old_temperature) / _impl->active_time_step;
                        conservation.generated_heat_rate +=
                            point.weighted_measure * _impl->cartesian->region_heat_source(region);
                    }
                    for (std::size_t q = 0; q < geometry.mechanical_points.size(); ++q) {
                        const Hex20MechanicalQuadraturePoint& point = geometry.mechanical_points[q];
                        const CartesianMaterialPointState& old_history =
                            _impl->cartesian_material_histories[region][element][q];
                        const CartesianMaterialPointState& new_history = update[q];
                        conservation.elastic_energy_change +=
                            0.5 * point.weighted_measure *
                            (cartesian::stress_strain_inner_product(new_history.stress, new_history.elastic_strain) -
                                cartesian::stress_strain_inner_product(old_history.stress, old_history.elastic_strain));
                        conservation.plastic_dissipation_increment +=
                            point.weighted_measure *
                            cartesian::trapezoidal_stress_strain_inner_product(old_history.stress, new_history.stress,
                                cartesian::strain_difference(new_history.plastic_strain, old_history.plastic_strain));
                        conservation.creep_dissipation_increment +=
                            point.weighted_measure *
                            cartesian::trapezoidal_stress_strain_inner_product(old_history.stress, new_history.stress,
                                cartesian::strain_difference(new_history.creep_strain, old_history.creep_strain));
                    }
                    staged[region][element] = std::move(update);
                    continue;
                }
                const Hex8LocalValues current = _impl->cartesian->volume_state(offset + element, converged_solution);
                const Hex8LocalValues old = _impl->cartesian->volume_state(offset + element, _impl->committed_solution);
                const Hex8Geometry& geometry = _impl->cartesian->region_element_geometry(region, element);
                const bool reduced =
                    _impl->cartesian->region(region).hex8_element_formulation == Hex8ElementFormulation::c3d8rt;
                const auto& capacity_points = reduced ? geometry.reduced_capacity_points : geometry.capacity_points;
                CartesianMaterialHistory update = _impl->cartesian->transient_update(region, element, current, old,
                    _impl->cartesian_material_histories[region][element], _impl->active_time_step);
                if (_impl->include_thermal_time_term)
                    for (std::size_t node = 0; node < hex8_node_count; ++node) {
                        const Hex8CapacityPoint& point = capacity_points[node];
                        conservation.stored_heat_rate +=
                            point.weighted_measure *
                            _impl->cartesian->heat_capacity(region, current[node], point.position) *
                            (current[node] - old[node]) / _impl->active_time_step;
                    }
                if (reduced)
                    conservation.generated_heat_rate +=
                        geometry.reduced_body_source_measure * _impl->cartesian->region_heat_source(region);
                else
                    for (const Hex8CapacityPoint& point : capacity_points)
                        conservation.generated_heat_rate +=
                            point.weighted_measure * _impl->cartesian->region_heat_source(region);
                if (reduced) {
                    const double current_hourglass =
                        _impl->cartesian->mechanical_hourglass_energy(region, element, current);
                    const double old_hourglass = _impl->cartesian->mechanical_hourglass_energy(region, element, old);
                    conservation.mechanical_hourglass_energy += current_hourglass;
                    conservation.mechanical_hourglass_energy_change += current_hourglass - old_hourglass;
                }
                double finite_current_volume = 0.0, finite_old_volume = 0.0;
                if (_impl->cartesian->region(region).strain_formulation == StrainFormulation::finite) {
                    Hex8LocalAdValues active_current{}, active_old{};
                    for (std::size_t local = 0; local < current.size(); ++local) {
                        active_current[local] = current[local];
                        active_old[local] = old[local];
                    }
                    for (const Hex8QuadraturePoint& point : geometry.points) {
                        finite_current_volume += evaluate_cartesian_incremental_kinematics(
                            point, active_current, old, StrainFormulation::finite)
                                                     .current_weighted_measure.value();
                        finite_old_volume +=
                            evaluate_cartesian_incremental_kinematics(point, active_old, old, StrainFormulation::finite)
                                .current_weighted_measure.value();
                    }
                }
                for (std::size_t q = 0; q < update.size(); ++q) {
                    const Hex8QuadraturePoint& point = reduced ? geometry.reduced_point : geometry.points[q];
                    const CartesianMaterialPointState &old_history =
                                                          _impl->cartesian_material_histories[region][element][q],
                                                      &new_history = update[q];
                    double current_measure = point.weighted_measure, old_measure = point.weighted_measure;
                    SymmetricTensor3Values diagnostic_new_stress = new_history.stress;
                    std::array<double, 6> diagnostic_new_plastic = new_history.plastic_strain;
                    std::array<double, 6> diagnostic_new_creep = new_history.creep_strain;
                    if (_impl->cartesian->region(region).strain_formulation == StrainFormulation::finite) {
                        Hex8LocalAdValues active_current{};
                        for (std::size_t local = 0; local < current.size(); ++local)
                            active_current[local] = current[local];
                        const CartesianKinematics kinematics = evaluate_cartesian_incremental_kinematics(
                            point, active_current, old, StrainFormulation::finite);
                        current_measure = point.weighted_measure / geometry.reference_volume * finite_current_volume;
                        old_measure = point.weighted_measure / geometry.reference_volume * finite_old_volume;
                        const CartesianRotation inverse_rotation = {kinematics.rotation.xx, kinematics.rotation.yx,
                            kinematics.rotation.zx, kinematics.rotation.xy, kinematics.rotation.yy,
                            kinematics.rotation.zy, kinematics.rotation.xz, kinematics.rotation.yz,
                            kinematics.rotation.zz};
                        diagnostic_new_stress = cartesian::rotate_tensor_values(new_history.stress, inverse_rotation);
                        diagnostic_new_plastic = cartesian::components(
                            cartesian::rotate_tensor_values(new_history.plastic_strain, inverse_rotation));
                        diagnostic_new_creep = cartesian::components(
                            cartesian::rotate_tensor_values(new_history.creep_strain, inverse_rotation));
                    }
                    conservation.elastic_energy_change +=
                        0.5 * (current_measure * cartesian::stress_strain_inner_product(
                                                     new_history.stress, new_history.elastic_strain) -
                                  old_measure * cartesian::stress_strain_inner_product(
                                                    old_history.stress, old_history.elastic_strain));
                    conservation.plastic_dissipation_increment +=
                        current_measure *
                        cartesian::trapezoidal_stress_strain_inner_product(old_history.stress, diagnostic_new_stress,
                            cartesian::strain_difference(diagnostic_new_plastic, old_history.plastic_strain));
                    conservation.creep_dissipation_increment +=
                        current_measure *
                        cartesian::trapezoidal_stress_strain_inner_product(old_history.stress, diagnostic_new_stress,
                            cartesian::strain_difference(diagnostic_new_creep, old_history.creep_strain));
                }
                staged[region][element] = std::move(update);
            }
        }
        add_trapezoidal_external_work(*this, converged_solution, _impl->committed_solution, raw_residual,
            _impl->committed_raw_residual, external_load_residual, _impl->committed_external_load_residual,
            conservation);
        finalize_conservation(*this, converged_solution, _impl->committed_solution, raw_residual, conservation);
        conservation.friction_dissipation_increment = _impl->cartesian->commit_contact_state(converged_solution);
        _impl->last_conservation_summary = conservation;
        _impl->cartesian_material_histories.swap(staged);
        _impl->committed_raw_residual = raw_residual;
        _impl->committed_external_load_residual = std::move(external_load_residual);
    } else {
        TransientConservationSummary conservation;
        std::vector<double> external_load_residual;
        const std::vector<double> raw_residual =
            accumulate_contribution_conservation(converged_solution, conservation, &external_load_residual);
        const std::size_t regions = _impl->rz->region_count();
        auto& staged = _impl->_staged_material_histories;
        for (std::size_t region = 0; region < regions; ++region) {
            const std::size_t offset = _impl->rz->region_element_offset(region);
            for (std::size_t element = 0; element < staged[region].size(); ++element) {
                const LocalValues state = gather_rz_state(*_impl->rz, offset + element, converged_solution);
                const LocalValues committed_state =
                    gather_rz_state(*_impl->rz, offset + element, _impl->committed_solution);
                const Quad4RzGeometry& geometry = _impl->rz->region_element_geometry(region, element);
                Quad4MaterialHistory update = compute_quad4_rz_transient_update(_impl->kernel_data[region], geometry,
                    state, committed_state, _impl->material_histories[region][element], _impl->active_time_step);
                for (std::size_t q = 0; q < geometry.points.size(); ++q) {
                    const RzQuadraturePoint& point = geometry.points[q];
                    double current_temperature = 0.0, old_temperature = 0.0;
                    for (std::size_t node = 0; node < quad4_node_count; ++node) {
                        current_temperature += point.shape[node] * state[node];
                        old_temperature += point.shape[node] * committed_state[node];
                    }
                    const Quad4RzData& kernel_data = _impl->kernel_data[region];
                    const MaterialFunctionContext context = {
                        kernel_data.time, point.radius, 0.0, point.axial_coordinate};
                    const double heat_capacity =
                        kernel_data.material.heat_capacity(current_temperature, context).value();
                    if (_impl->include_thermal_time_term)
                        conservation.stored_heat_rate += point.weighted_measure * heat_capacity *
                                                         (current_temperature - old_temperature) /
                                                         _impl->active_time_step;
                    conservation.generated_heat_rate += point.weighted_measure * kernel_data.volumetric_heat_source;
                    const MaterialPointState &old_history = _impl->material_histories[region][element][q],
                                             &new_history = update[q];
                    const AxisymmetricStressValues &old_stress = old_history.stress, &new_stress = new_history.stress;
                    conservation.elastic_energy_change +=
                        0.5 * point.weighted_measure *
                        (rz::stress_strain_inner_product(new_stress, new_history.elastic_strain) -
                            rz::stress_strain_inner_product(old_stress, old_history.elastic_strain));
                    conservation.plastic_dissipation_increment +=
                        point.weighted_measure *
                        rz::trapezoidal_stress_strain_inner_product(old_stress, new_stress,
                            rz::strain_difference(new_history.plastic_strain, old_history.plastic_strain));
                    conservation.creep_dissipation_increment +=
                        point.weighted_measure *
                        rz::trapezoidal_stress_strain_inner_product(old_stress, new_stress,
                            rz::strain_difference(new_history.creep_strain, old_history.creep_strain));
                }
                staged[region][element] = std::move(update);
            }
        }
        add_trapezoidal_external_work(*this, converged_solution, _impl->committed_solution, raw_residual,
            _impl->committed_raw_residual, external_load_residual, _impl->committed_external_load_residual,
            conservation);
        finalize_conservation(*this, converged_solution, _impl->committed_solution, raw_residual, conservation);
        _impl->last_conservation_summary = conservation;
        _impl->rz->commit_contact_state(converged_solution);
        _impl->material_histories.swap(staged);
        _impl->committed_raw_residual = raw_residual;
        _impl->committed_external_load_residual = std::move(external_load_residual);
    }
    if (_impl->track_previous_committed_solution) {
        _impl->previous_committed_solution = _impl->committed_solution;
        _impl->previous_committed_time = _impl->committed_time;
    }
    _impl->committed_solution = converged_solution;
    _impl->committed_time = _impl->active_end_time;
    _impl->committed_load_factor = _impl->active_load_factor;
    clear_active_time_step();
}

void TransientProblem::rollback_time_step() noexcept {
    if (!_impl->time_step_active) return;
    apply_spatial_controls(_impl->committed_time, _impl->committed_load_factor);
    if (!_impl->active_contact_histories.empty())
        _impl->restore_contact_state(_impl->committed_solution, std::move(_impl->active_contact_histories));
    clear_active_time_step();
}

void TransientProblem::apply_spatial_controls(double time, double load_factor) {
    _impl->set_time(time);
    _impl->set_load_factor(load_factor);
    if (_impl->is_cartesian()) return;
    for (Quad4RzData& kernel_data : _impl->kernel_data) kernel_data.time = time;
    for (std::size_t region = 0; region < _impl->rz->region_count(); ++region)
        _impl->kernel_data[region].volumetric_heat_source = _impl->rz->region_heat_source(region);
}

void TransientProblem::clear_active_time_step() noexcept {
    _impl->active_time_step = 0.0;
    _impl->active_end_time = _impl->committed_time;
    _impl->active_load_factor = _impl->committed_load_factor;
    _impl->active_contact_histories.clear();
    _impl->include_thermal_time_term = true;
    _impl->time_step_active = false;
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

std::size_t TransientProblem::sparsity_contribution_count() const noexcept {
    return _impl->sparsity_contribution_count();
}

std::pair<std::size_t, std::size_t> TransientProblem::contribution_partition(
    std::size_t partition, std::size_t partition_count) const {
    return _impl->contribution_partition(partition, partition_count);
}

const std::vector<FieldDescriptor>& TransientProblem::field_layout() const noexcept {
    return _impl->layout().field_layout();
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

void TransientProblem::validate_local_state(
    std::size_t first, std::size_t last, const std::vector<double>& state) const {
    require_active_time_step();
    _impl->validate_local_state(*this, first, last, state);
}

void TransientProblem::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->contribution_dofs(index, dofs);
}

void TransientProblem::contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
    _impl->contribution_jacobian_pattern(index, pattern);
}

void TransientProblem::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _impl->sparsity_contribution_dofs(index, dofs);
}

void TransientProblem::sparsity_contribution_jacobian_pattern(
    std::size_t index, std::vector<unsigned char>& pattern) const {
    _impl->sparsity_contribution_jacobian_pattern(index, pattern);
}

void TransientProblem::compute_contribution(std::size_t index, const std::vector<double>& state,
    std::vector<double>& residual, std::vector<double>* jacobian) const {
    require_active_time_step();
    if (_impl->is_cartesian()) {
        const CartesianMaterialHistory* history = nullptr;
        if (index < _impl->cartesian->volume_contribution_count()) {
            const auto location = _impl->cartesian->element_location(index);
            history = &_impl->cartesian_material_histories[location.first][location.second];
        }
        _impl->cartesian->compute_contribution(index, state, &_impl->committed_solution, history,
            _impl->active_time_step, residual, jacobian, _impl->include_thermal_time_term);
        return;
    }
    const rz::TransientBackendView backend = BackendAccess::transient(*this);
    const LocalValues local_state = rz_local_values(state);
    LocalJacobian local_jacobian{};
    LocalResidual local_residual{};
    if (index >= backend.spatial.volume_contribution_count())
        local_residual =
            backend.spatial.compute_contribution(index, local_state, jacobian == nullptr ? nullptr : &local_jacobian);
    else {
        const auto location = backend.spatial.element_location(index);
        local_residual = compute_quad4_rz_transient(backend.kernel_data[location.first],
            backend.spatial.region_element_geometry(location.first, location.second), local_state,
            gather_rz_state(backend.spatial, index, backend.committed_solution),
            backend.histories[location.first][location.second], backend.active_time_step,
            jacobian == nullptr ? nullptr : &local_jacobian, backend.include_thermal_time_term);
    }
    residual.assign(local_residual.begin(), local_residual.end());
    if (jacobian != nullptr) jacobian->assign(local_jacobian.begin(), local_jacobian.end());
}

void TransientProblem::require_active_time_step() const {
    if (!_impl->time_step_active)
        throw std::logic_error("TransientProblem residual evaluation requires an active time step");
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
