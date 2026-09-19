#include "c3d20rt.hpp"
#include "c3d20t.hpp"
#include "c3d8_types.hpp"
#include "c3d8rt.hpp"
#include "c3d8t.hpp"
#include "cartesian3d_assembly.hpp"
#include "cax2t_gps.hpp"
#include "cax4rt.hpp"
#include "cax4t.hpp"
#include "cax8rt.hpp"
#include "cax8t.hpp"
#include "contact_types.hpp"
#include "core/cax_evaluation.hpp"
#include "core/element_evaluation.hpp"
#include "core/element_region_data.hpp"
#include "core/nonlinear_problem.hpp"
#include "core/problem_backend_access.hpp"
#include "core/spatial_definition.hpp"
#include "core/steady_problem.hpp"
#include "core/transient_problem.hpp"
#include "cpeg8t.hpp"
#include "plane_assembly.hpp"
#include "rz_assembly.hpp"
#include "solver/solve_workflows.hpp"
#include "thermal_assembly.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "spatial_backend_common.hpp"

namespace fuelsim {
std::string set_commit_failure(std::vector<double>& values, const std::exception_ptr& failure) {
    std::string failure_message;
    if (failure) {
        try {
            std::rethrow_exception(failure);
        } catch (const std::domain_error& error) {
            values[0] = 1.0;
            failure_message = error.what();
        } catch (const std::overflow_error& error) {
            values[0] = 1.0;
            failure_message = error.what();
        } catch (const std::exception& error) {
            values[1] = 1.0;
            failure_message = error.what();
        } catch (...) {
            values[1] = 1.0;
            failure_message = "unknown partitioned commit failure";
        }
    }
    return failure_message;
}

void sum_commit_partitions(std::vector<double>& values,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    std::string failure_message) {
    const std::size_t size = values.size();
    sum_partitions(values);
    if (values.size() != size)
        throw std::logic_error("Partition sum changed the transient commit buffer size");
    if (values[0] != 0.0 || values[1] != 0.0) {
        if (failure_message.empty())
            failure_message = "transient commit failed on another partition";
        if (values[1] != 0.0)
            throw std::runtime_error(failure_message);
        throw std::domain_error(failure_message);
    }
}

void append_commit_diagnostics(std::vector<double>& values,
    const std::vector<double>& residual,
    const std::vector<double>& external,
    const TransientConservationSummary& conservation) {
    values.insert(values.end(), residual.begin(), residual.end());
    values.insert(values.end(), external.begin(), external.end());
    for (const auto& field : transient_conservation_fields)
        values.push_back(conservation.*field.member);
}

std::size_t read_commit_diagnostics(const std::vector<double>& values,
    std::vector<double>& residual,
    std::vector<double>& external,
    TransientConservationSummary& conservation) {
    std::size_t entry = 2;
    for (double& value : residual)
        value = values[entry++];
    for (double& value : external)
        value = values[entry++];
    for (const auto& field : transient_conservation_fields)
        conservation.*field.member = values[entry++];
    return entry;
}

void synchronize_commit_diagnostics(std::vector<double>& residual,
    std::vector<double>& external,
    TransientConservationSummary& conservation,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    const std::exception_ptr& failure) {
    std::vector<double> values;
    values.reserve(2 + residual.size() + external.size() + transient_conservation_fields.size());
    values.resize(2, 0.0);
    const auto message = set_commit_failure(values, failure);
    append_commit_diagnostics(values, residual, external, conservation);
    sum_commit_partitions(values, sum_partitions, message);
    (void)read_commit_diagnostics(values, residual, external, conservation);
}

void synchronize_cartesian_commit(std::vector<std::vector<CartesianMaterialHistory>>& histories,
    std::vector<double>& residual,
    std::vector<double>& external,
    TransientConservationSummary& conservation,
    std::size_t first,
    std::size_t last,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    const std::exception_ptr& failure) {
    // Only staged data is exchanged. No committed field or history is changed
    // until every partition has completed successfully.
    std::size_t point_count = 0;
    for (const auto& region : histories)
        for (const auto& history : region)
            point_count += history.size();
    std::vector<double> values;
    values.reserve(2 + residual.size() + external.size() + transient_conservation_fields.size() + 26 * point_count);
    values.resize(2, 0.0);
    const std::string failure_message = set_commit_failure(values, failure);
    append_commit_diagnostics(values, residual, external, conservation);
    std::size_t element = 0;
    for (const auto& region : histories)
        for (const auto& history : region) {
            const bool owned = element >= first && element < last && !failure;
            for (const auto& point : history) {
                for (const auto* tensor : {&point.elastic_strain, &point.plastic_strain, &point.creep_strain})
                    for (const double component : *tensor)
                        values.push_back(owned ? component : 0.0);
                values.push_back(owned ? point.equivalent_plastic_strain : 0.0);
                values.push_back(owned ? point.equivalent_creep_strain : 0.0);
                for (const double component : {point.stress.xx,
                         point.stress.yy,
                         point.stress.zz,
                         point.stress.xy,
                         point.stress.yz,
                         point.stress.xz})
                    values.push_back(owned ? component : 0.0);
            }
            ++element;
        }
    sum_commit_partitions(values, sum_partitions, failure_message);
    std::size_t entry = read_commit_diagnostics(values, residual, external, conservation);
    for (auto& region : histories)
        for (auto& history : region)
            for (auto& point : history) {
                for (auto* tensor : {&point.elastic_strain, &point.plastic_strain, &point.creep_strain})
                    for (double& component : *tensor)
                        component = values[entry++];
                point.equivalent_plastic_strain = values[entry++];
                point.equivalent_creep_strain = values[entry++];
                for (double* component : {&point.stress.xx,
                         &point.stress.yy,
                         &point.stress.zz,
                         &point.stress.xy,
                         &point.stress.yz,
                         &point.stress.xz})
                    *component = values[entry++];
            }
}

// Both axisymmetric topologies use the same concrete four-component material
// state. Include inactive slots so their required zero values survive exchange.
void synchronize_rz_commit(const std::vector<MaterialPointState*>& points,
    std::size_t points_per_element,
    std::vector<double>& residual,
    std::vector<double>& external,
    TransientConservationSummary& conservation,
    std::size_t first,
    std::size_t last,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    const std::exception_ptr& failure) {
    std::vector<double> values;
    values.reserve(2 + residual.size() + external.size() + transient_conservation_fields.size() + 18 * points.size());
    values.resize(2, 0.0);
    const std::string failure_message = set_commit_failure(values, failure);
    append_commit_diagnostics(values, residual, external, conservation);
    for (std::size_t index = 0; index < points.size(); ++index) {
        const auto& point = *points[index];
        const std::size_t element = index / points_per_element;
        const bool owned = element >= first && element < last && !failure;
        for (const auto* tensor : {&point.elastic_strain, &point.plastic_strain, &point.creep_strain})
            for (const double component : *tensor)
                values.push_back(owned ? component : 0.0);
        values.push_back(owned ? point.equivalent_plastic_strain : 0.0);
        values.push_back(owned ? point.equivalent_creep_strain : 0.0);
        for (const double component : {point.stress.rr, point.stress.zz, point.stress.hoop, point.stress.rz})
            values.push_back(owned ? component : 0.0);
    }
    sum_commit_partitions(values, sum_partitions, failure_message);
    std::size_t entry = read_commit_diagnostics(values, residual, external, conservation);
    for (auto* point : points) {
        for (auto* tensor : {&point->elastic_strain, &point->plastic_strain, &point->creep_strain})
            for (double& component : *tensor)
                component = values[entry++];
        point->equivalent_plastic_strain = values[entry++];
        point->equivalent_creep_strain = values[entry++];
        for (double* component : {&point->stress.rr, &point->stress.zz, &point->stress.hoop, &point->stress.rz})
            *component = values[entry++];
    }
}

double body_point_work(const RegionDefinition& region,
    const MaterialFunctionContext& context,
    double reference_measure,
    const std::array<double, 3>& displacement_increment) {
    if (region.body_acceleration == std::array<double, 3>{})
        return 0.0;
    const IsotropicThermoelasticMaterial material(region.material);
    double specific_work = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        specific_work += region.body_acceleration[component] * displacement_increment[component];
    return reference_measure * material.initial_density(region.initial_temperature, context) * specific_work;
}

void finalize_conservation(const NonlinearProblem& problem,
    const std::vector<double>& current,
    const std::vector<double>& old,
    const std::vector<double>& current_raw_residual,
    const std::vector<double>& old_raw_residual,
    const std::vector<double>& current_external_load_residual,
    const std::vector<double>& old_external_load_residual,
    TransientConservationSummary& result) {
    if (current_raw_residual.size() != problem.dof_count() || old_raw_residual.size() != problem.dof_count()
        || current_external_load_residual.size() != problem.dof_count()
        || old_external_load_residual.size() != problem.dof_count())
        throw std::logic_error("Transient trapezoidal work residual layout does not match the problem");
    std::vector<bool> constrained(problem.dof_count(), false);
    for (const DirichletCondition& condition : problem.dirichlet_conditions())
        constrained[condition.dof] = true;
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        if (problem.field_layout()[problem.field_index(dof)].category == FieldCategory::thermal)
            continue;
        const double increment = current[dof] - old[dof];
        result.trapezoidal_pressure_traction_work_increment -=
            0.5 * (old_external_load_residual[dof] + current_external_load_residual[dof]) * increment;
        if (constrained[dof])
            result.trapezoidal_dirichlet_reaction_work_increment +=
                0.5 * (old_raw_residual[dof] + current_raw_residual[dof]) * increment;
    }

    for (const DirichletCondition& condition : problem.dirichlet_conditions()) {
        if (problem.field_layout()[problem.field_index(condition.dof)].category == FieldCategory::thermal)
            result.dirichlet_heat_input_rate += current_raw_residual[condition.dof];
        else
            result.dirichlet_reaction_work_increment +=
                current_raw_residual[condition.dof] * (current[condition.dof] - old[condition.dof]);
    }
    double thermal_residual_squared = 0.0, mechanical_residual_squared = 0.0;
    for (std::size_t dof = 0; dof < problem.dof_count(); ++dof) {
        if (constrained[dof])
            continue;
        double& norm_squared = problem.field_layout()[problem.field_index(dof)].category == FieldCategory::thermal
                                   ? thermal_residual_squared
                                   : mechanical_residual_squared;
        norm_squared += current_raw_residual[dof] * current_raw_residual[dof];
    }
    result.unconstrained_thermal_residual_l2 = std::sqrt(thermal_residual_squared);
    result.unconstrained_mechanical_residual_l2 = std::sqrt(mechanical_residual_squared);
    result.global_thermal_balance = result.stored_heat_rate + result.convection_heat_rate
                                    + result.interface_heat_imbalance - result.generated_heat_rate
                                    - result.surface_heat_input_rate - result.dirichlet_heat_input_rate;
    const double thermal_scale = std::abs(result.stored_heat_rate) + std::abs(result.convection_heat_rate)
                                 + std::abs(result.interface_heat_imbalance) + std::abs(result.generated_heat_rate)
                                 + std::abs(result.surface_heat_input_rate)
                                 + std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0 ? std::abs(result.global_thermal_balance) / thermal_scale : 0.0;
    result.internal_mechanical_work_increment += result.body_force_work_increment;
    result.mechanical_work_balance = result.internal_mechanical_work_increment + result.contact_work_increment
                                     - result.pressure_traction_work_increment
                                     - result.dirichlet_reaction_work_increment - result.body_force_work_increment;
    const double mechanical_scale =
        std::abs(result.body_force_work_increment) + std::abs(result.internal_mechanical_work_increment)
        + std::abs(result.contact_work_increment) + std::abs(result.pressure_traction_work_increment)
        + std::abs(result.dirichlet_reaction_work_increment);
    result.relative_mechanical_work_balance =
        mechanical_scale > 0.0 ? std::abs(result.mechanical_work_balance) / mechanical_scale : 0.0;
}

Cax4LocalValues rz_local_values(const std::vector<double>& values) {
    if (values.size() != cax4_local_dof_count)
        throw std::invalid_argument("RZ contribution state must contain 12 DOFs");
    Cax4LocalValues result{};
    std::copy(values.begin(), values.end(), result.begin());
    return result;
}

Cax4LocalValues
gather_rz_state(const rz::SpatialAssembly& spatial, std::size_t index, const std::vector<double>& global_state) {
    const Cax4LocalDofs dofs = spatial.contribution_dofs(index);
    Cax4LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local)
        result[local] = global_state.at(dofs[local]);
    return result;
}

namespace rz {
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

double trapezoidal_stress_strain_inner_product(const AxisymmetricStressValues& old_stress,
    const AxisymmetricStressValues& new_stress,
    const std::array<double, 4>& strain_increment) noexcept {
    return 0.5
           * (stress_strain_inner_product(old_stress, strain_increment)
               + stress_strain_inner_product(new_stress, strain_increment));
}

void accumulate_material_conservation(TransientConservationSummary& result,
    const MaterialPointState& old,
    const MaterialPointState& current,
    double measure) {
    result.elastic_energy_change += 0.5 * measure
                                    * (stress_strain_inner_product(current.stress, current.elastic_strain)
                                        - stress_strain_inner_product(old.stress, old.elastic_strain));
    result.plastic_dissipation_increment += measure
                                            * trapezoidal_stress_strain_inner_product(old.stress,
                                                current.stress,
                                                strain_difference(current.plastic_strain, old.plastic_strain));
    result.creep_dissipation_increment += measure
                                          * trapezoidal_stress_strain_inner_product(old.stress,
                                              current.stress,
                                              strain_difference(current.creep_strain, old.creep_strain));
}
} // namespace rz

namespace cartesian {
double stress_strain_inner_product(const SymmetricTensor3Values& stress, const std::array<double, 6>& strain) noexcept {
    return stress.xx * strain[0] + stress.yy * strain[1] + stress.zz * strain[2]
           + 2.0 * (stress.xy * strain[3] + stress.yz * strain[4] + stress.xz * strain[5]);
}

std::array<double, 6> strain_difference(const std::array<double, 6>& current, const std::array<double, 6>& old) {
    std::array<double, 6> result{};
    for (std::size_t component = 0; component < result.size(); ++component)
        result[component] = current[component] - old[component];
    return result;
}

double trapezoidal_stress_strain_inner_product(const SymmetricTensor3Values& old_stress,
    const SymmetricTensor3Values& new_stress,
    const std::array<double, 6>& strain_increment) noexcept {
    return 0.5
           * (stress_strain_inner_product(old_stress, strain_increment)
               + stress_strain_inner_product(new_stress, strain_increment));
}

SymmetricTensor3Values rotate_tensor_values(const std::array<double, 6>& tensor, const CartesianRotation& rotation) {
    const SymmetricTensor3 rotated =
        rotate_cartesian_tensor({tensor[0], tensor[1], tensor[2], tensor[3], tensor[4], tensor[5]}, rotation);
    return {rotated.xx.value(),
        rotated.yy.value(),
        rotated.zz.value(),
        rotated.xy.value(),
        rotated.yz.value(),
        rotated.xz.value()};
}

SymmetricTensor3Values rotate_tensor_values(const SymmetricTensor3Values& tensor, const CartesianRotation& rotation) {
    const std::array<double, 6> values = {tensor.xx, tensor.yy, tensor.zz, tensor.xy, tensor.yz, tensor.xz};
    return rotate_tensor_values(values, rotation);
}

std::array<double, 6> components(const SymmetricTensor3Values& tensor) {
    return {tensor.xx, tensor.yy, tensor.zz, tensor.xy, tensor.yz, tensor.xz};
}
} // namespace cartesian

bool finite_stress(const AxisymmetricStressValues& stress) {
    return std::isfinite(stress.rr) && std::isfinite(stress.zz) && std::isfinite(stress.hoop)
           && std::isfinite(stress.rz);
}

bool valid_material_state(const MaterialPointState& state) {
    for (std::size_t component = 0; component < 4; ++component)
        if (!std::isfinite(state.elastic_strain[component]) || !std::isfinite(state.plastic_strain[component])
            || !std::isfinite(state.creep_strain[component]))
            return false;
    return std::isfinite(state.equivalent_plastic_strain) && state.equivalent_plastic_strain >= 0.0
           && std::isfinite(state.equivalent_creep_strain) && state.equivalent_creep_strain >= 0.0;
}

bool valid_material_state(const CartesianMaterialPointState& state) {
    for (std::size_t component = 0; component < state.elastic_strain.size(); ++component)
        if (!std::isfinite(state.elastic_strain[component]) || !std::isfinite(state.plastic_strain[component])
            || !std::isfinite(state.creep_strain[component]))
            return false;
    return std::isfinite(state.equivalent_plastic_strain) && state.equivalent_plastic_strain >= 0.0
           && std::isfinite(state.equivalent_creep_strain) && state.equivalent_creep_strain >= 0.0
           && std::isfinite(state.stress.xx) && std::isfinite(state.stress.yy) && std::isfinite(state.stress.zz)
           && std::isfinite(state.stress.xy) && std::isfinite(state.stress.yz) && std::isfinite(state.stress.xz);
}

std::pair<std::size_t, std::size_t> SpatialBackend::contribution_partition(std::size_t partition,
    std::size_t partition_count) const {
    if (partition_count == 0 || partition >= partition_count)
        throw std::out_of_range("Spatial problem contribution partition is out of range");
    return {contribution_count() * partition / partition_count,
        contribution_count() * (partition + 1U) / partition_count};
}

void SpatialBackend::contribution_jacobian_pattern(std::size_t index, std::vector<unsigned char>& pattern) const {
    std::vector<std::size_t> dofs;
    contribution_dofs(index, dofs);
    pattern.assign(dofs.size() * dofs.size(), 1U);
}

void SpatialBackend::sparsity_contribution_jacobian_pattern(std::size_t index,
    std::vector<unsigned char>& pattern) const {
    std::vector<std::size_t> dofs;
    sparsity_contribution_dofs(index, dofs);
    pattern.assign(dofs.size() * dofs.size(), 1U);
}

const std::vector<std::vector<ContactPointHistory>>& SpatialBackend::committed_contact_histories() const noexcept {
    static const std::vector<std::vector<ContactPointHistory>> empty;
    return empty;
}

void SpatialBackend::commit_contact_state(const std::vector<double>&) {
}

void SpatialBackend::restore_contact_state(const std::vector<double>& state,
    std::vector<std::vector<ContactPointHistory>> histories) {
    if (!histories.empty())
        throw std::invalid_argument("This backend has no mechanical contact history");
    validate_state(state);
}

AugmentedContactUpdate SpatialBackend::update_augmented_contact_multipliers(const std::vector<double>&, std::size_t) {
    throw std::logic_error("This spatial backend does not support augmented contact");
}

void SpatialBackend::capture_active_contact_state() {
    if (uses_augmented_contact())
        _active_contact_histories = committed_contact_histories();
}

void SpatialBackend::restore_active_contact_state() {
    if (!_active_contact_histories.empty())
        restore_contact_state(_state.committed_solution, std::move(_active_contact_histories));
}

void SpatialBackend::prepare_contacts(const std::vector<double>& state,
    std::exception_ptr& failure,
    const std::function<void(std::vector<double>&)>& sum_partitions) {
    if (failure)
        return;
    try {
        check_commit_test_hook(CommitPreparationStage::material);
        stage_contact_history(state);
        check_commit_test_hook(CommitPreparationStage::contact);
    } catch (...) {
        if (!sum_partitions)
            throw;
        failure = std::current_exception();
    }
}

} // namespace fuelsim
