#include "spatial_backend_common.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
void ThermalBackend::accumulate_time_error(const TransientCommittedState&,
    const TransientCommittedState&,
    const TransientTimeOptions&,
    TransientTimeErrorEstimate&) const {
}

ThermalBackend::ThermalBackend(SpatialTimeState& state,
    SpatialDefinition definition,
    const UnstructuredQuad4Mesh& mesh,
    bool transient)
    : SpatialBackend(state), _spatial(std::move(definition), mesh) {
    if (transient && _spatial.has_pellet_response())
        throw std::invalid_argument("Pellet responses support steady conduction only");
}

ThermalBackend::ThermalBackend(SpatialTimeState& state,
    SpatialDefinition definition,
    const UnstructuredQuad8Mesh& mesh,
    bool transient)
    : SpatialBackend(state), _spatial(std::move(definition), mesh) {
    if (transient && _spatial.has_pellet_response())
        throw std::invalid_argument("Pellet responses support steady conduction only");
}

ThermalBackend::ThermalBackend(SpatialTimeState& state,
    SpatialDefinition definition,
    const UnstructuredHex8Mesh& mesh,
    bool transient)
    : SpatialBackend(state), _spatial(std::move(definition), mesh) {
    if (transient && _spatial.has_pellet_response())
        throw std::invalid_argument("Pellet responses support steady conduction only");
}

ThermalBackend::ThermalBackend(SpatialTimeState& state,
    SpatialDefinition definition,
    const UnstructuredHex20Mesh& mesh,
    bool transient)
    : SpatialBackend(state), _spatial(std::move(definition), mesh) {
    if (transient && _spatial.has_pellet_response())
        throw std::invalid_argument("Pellet responses support steady conduction only");
}

void ThermalBackend::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _spatial.contribution_dofs(index, dofs);
}

void ThermalBackend::set_time(double value) {
    _spatial.set_time(value);
}

void ThermalBackend::set_load_factor(double value) {
    _spatial.set_load_factor(value);
}

void ThermalBackend::set_heat_source_interval(double, double) {
}

void ThermalBackend::set_region_strain_formulation(std::size_t region, StrainFormulation formulation) {
    _spatial.set_region_strain_formulation(region, formulation);
}

void ThermalBackend::compute_contribution(std::size_t index,
    const std::vector<double>& state,
    std::vector<double>& residual,
    std::vector<double>* jacobian,
    bool transient) const {
    auto result = _spatial.evaluate(index,
        state,
        _state.committed_solution,
        transient && _state.include_thermal_time_term ? _state.active_time_step : 0.0,
        _state.committed_time,
        jacobian != nullptr);
    residual = std::move(result.residual);
    if (jacobian)
        *jacobian = std::move(result.jacobian);
    return;
}

void ThermalBackend::prepare_time_step(const NonlinearProblem& problem,
    const std::vector<double>& converged_solution,
    std::size_t first_contribution,
    std::size_t last_contribution,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    TransientConservationSummary& conservation,
    std::vector<double>& raw_residual,
    std::vector<double>& external_load_residual,
    std::exception_ptr partition_failure) {
    try {
        if (!partition_failure)
            for (std::size_t index = first_contribution;
                index < std::min(last_contribution, _spatial.volume_contribution_count());
                ++index) {
                std::vector<std::size_t> dofs;
                _spatial.contribution_dofs(index, dofs);
                std::vector<double> local;
                for (auto n : dofs)
                    local.push_back(converged_solution[n]);
                const auto result = _spatial.evaluate(index,
                    local,
                    _state.committed_solution,
                    _state.include_thermal_time_term ? _state.active_time_step : 0.0,
                    _state.committed_time,
                    false);
                conservation.generated_heat_rate += result.generated_heat_rate;
                conservation.stored_heat_rate += result.stored_heat_rate;
            }
    } catch (...) {
        if (!sum_partitions)
            throw;
        partition_failure = std::current_exception();
    }
    prepare_contacts(converged_solution, partition_failure, sum_partitions);
    if (sum_partitions)
        synchronize_commit_diagnostics(raw_residual,
            external_load_residual,
            conservation,
            sum_partitions,
            partition_failure);
    finalize_conservation(problem,
        converged_solution,
        _state.committed_solution,
        raw_residual,
        _state.committed_raw_residual,
        external_load_residual,
        _state.committed_external_load_residual,
        conservation);
}

void ThermalBackend::export_histories(TransientCommittedState& state) const {
    state.contact_histories = committed_contact_histories();
}

void ThermalBackend::restore_histories(TransientCommittedState& state) {
    if (!state.material_histories.empty() || !state.cartesian_material_histories.empty()
        || !state.quad8_material_histories.empty() || !state.radial_material_histories.empty())
        throw std::invalid_argument("Unexpected material history for this spatial backend");
    validate_state(state.solution);
    restore_contact_state(state.solution, std::move(state.contact_histories));
}

RegionStateSummary ThermalBackend::summarize_region(std::size_t region) const {
    if (region >= _spatial.region_count())
        throw std::out_of_range("TransientProblem region summary index is out of range");
    const auto& fields = _spatial.field_layout();
    const auto temperature = std::find_if(fields.begin(), fields.end(), [](const FieldDescriptor& field) {
        return field.category == FieldCategory::thermal;
    });
    RegionStateSummary result{-std::numeric_limits<double>::infinity(), 0.0, 0.0};
    for (std::size_t node = 0; node < _spatial.source_nodes(region).size(); ++node) {
        result.maximum_temperature = std::max(result.maximum_temperature,
            _state.committed_solution.at(temperature->begin + _spatial.global_temperature_node(region, node)));
    }
    return result;
}

std::vector<double> ThermalBackend::committed_creep_rates() const {
    return {};
}

void ThermalBackend::publish_material_history() noexcept {
}
} // namespace fuelsim
