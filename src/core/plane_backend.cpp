#include "spatial_backend_common.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
void PlaneBackend::accumulate_time_error(const TransientCommittedState& full,
    const TransientCommittedState& half,
    const TransientTimeOptions& options,
    TransientTimeErrorEstimate& result) const {
    accumulate_cartesian_material_time_error(full, half, options, result);
    accumulate_cartesian_contact_time_error(full, half, options, result);
}

PlaneBackend::PlaneBackend(SpatialTimeState& state,
    SpatialDefinition definition,
    const UnstructuredPlaneQuad8Mesh& mesh,
    bool)
    : SpatialBackend(state), _spatial(std::move(definition), mesh) {
    {
        _histories.resize(_spatial.region_count());
        for (std::size_t r = 0; r < _spatial.region_count(); ++r)
            _histories[r].resize(_spatial.region_element_count(r), CartesianMaterialHistory(9));
    }
}

void PlaneBackend::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _spatial.contribution_dofs(index, dofs);
}

void PlaneBackend::set_time(double value) {
    _spatial.set_time(value);
    _time = value;
}

void PlaneBackend::set_load_factor(double value) {
    _spatial.set_load_factor(value);
}

void PlaneBackend::set_heat_source_interval(double, double) {
}

void PlaneBackend::set_region_strain_formulation(std::size_t region, StrainFormulation formulation) {
    _spatial.set_region_strain_formulation(region, formulation);
}

elements::Cpeg8Result PlaneBackend::evaluate_plane(std::size_t index,
    const std::vector<double>& local,
    bool transient,
    bool jacobian,
    bool history) const {
    const auto [r, e] = _spatial.element_location(index);
    std::vector<std::size_t> dofs;
    _spatial.contribution_dofs(index, dofs);
    if (local.size() != 23)
        throw std::invalid_argument("CPEG8T volume state must contain 23 values");
    elements::Cpeg8Values current{}, old{};
    std::copy(local.begin(), local.end(), current.begin());
    const auto& region = _spatial.region(r);
    elements::Cpeg8History committed{};
    if (transient) {
        for (std::size_t i = 0; i < 23; ++i)
            old[i] = _state.committed_solution[dofs[i]];
        std::copy(_histories[r][e].begin(), _histories[r][e].end(), committed.begin());
    } else
        for (std::size_t i = 0; i < 4; ++i)
            old[i] = region.initial_temperature;
    const IsotropicThermoelasticMaterial material(region.material);
    elements::Cpeg8Input input{material, _spatial.geometry(r, e), current, old, transient ? &committed : nullptr};
    input.time = transient ? _state.active_end_time : _time;
    input.time_step = transient ? _state.active_time_step : 0.0;
    input.include_thermal_time_term = transient && _state.include_thermal_time_term;
    input.strain_formulation = region.strain_formulation;
    input.volumetric_heat_source =
        transient && region.heat_source_time_evaluation == HeatSourceTimeEvaluation::interval_average
            ? _spatial.region_heat_source_average(r, _state.committed_time, _state.active_end_time)
            : _spatial.region_heat_source(r);
    input.initial_temperature = region.initial_temperature;
    input.body_acceleration = {region.body_acceleration[0], region.body_acceleration[1]};
    return elements::evaluate_cpeg8t(input, {true, jacobian, history, false});
}

void PlaneBackend::compute_contribution(std::size_t index,
    const std::vector<double>& state,
    std::vector<double>& residual,
    std::vector<double>* jacobian,
    bool transient) const {
    if (index >= _spatial.contact_offset()) {
        const auto evaluated = _spatial.compute_contact(index, state, jacobian != nullptr);
        residual.assign(evaluated.residual.begin(), evaluated.residual.end());
        if (jacobian)
            jacobian->assign(evaluated.jacobian.begin(), evaluated.jacobian.end());
        return;
    }
    const auto evaluated = index < _spatial.volume_contribution_count()
                               ? evaluate_plane(index, state, transient, jacobian != nullptr, false)
                               : _spatial.compute_boundary(index, state, jacobian != nullptr);
    residual.assign(evaluated.residual.begin(), evaluated.residual.end());
    if (jacobian)
        jacobian->assign(evaluated.jacobian.begin(), evaluated.jacobian.end());
    return;
}

void PlaneBackend::prepare_time_step(const NonlinearProblem& problem,
    const std::vector<double>& converged_solution,
    std::size_t first_contribution,
    std::size_t last_contribution,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    TransientConservationSummary& conservation,
    std::vector<double>& raw_residual,
    std::vector<double>& external_load_residual,
    std::exception_ptr partition_failure) {
    _staged = _histories;
    auto& staged = _staged;
    try {
        if (!partition_failure)
            for (std::size_t index = first_contribution;
                index < std::min(last_contribution, _spatial.volume_contribution_count());
                ++index) {
                std::vector<std::size_t> dofs;
                _spatial.contribution_dofs(index, dofs);
                std::vector<double> local(dofs.size());
                for (std::size_t i = 0; i < dofs.size(); ++i)
                    local[i] = converged_solution[dofs[i]];
                const auto result = evaluate_plane(index, local, true, false, true);
                const auto [r, e] = _spatial.element_location(index);
                staged[r][e].assign(result.history.begin(), result.history.end());
                conservation.generated_heat_rate += result.generated_heat_rate;
                conservation.stored_heat_rate += result.stored_heat_rate;
                conservation.elastic_energy_change += result.elastic_energy_change;
                conservation.plastic_dissipation_increment += result.plastic_dissipation_increment;
                conservation.creep_dissipation_increment += result.creep_dissipation_increment;
                const auto& geometry = _spatial.geometry(r, e);
                for (std::size_t q = 0; q < result.history.size(); ++q) {
                    const auto& point = geometry.points[q];
                    std::array<double, 3> increment{};
                    for (std::size_t component = 0; component < 2; ++component)
                        for (std::size_t node = 0; node < 8; ++node) {
                            const auto dof = dofs[4 + 8 * component + node];
                            increment[component] +=
                                point.shape[node] * (converged_solution[dof] - _state.committed_solution[dof]);
                        }
                    conservation.body_force_work_increment +=
                        body_point_work(_spatial.region(r), {0.0, point.x, point.y, 0.0}, point.measure, increment);
                }
            }
    } catch (...) {
        if (!sum_partitions)
            throw;
        partition_failure = std::current_exception();
    }
    prepare_contacts(converged_solution, partition_failure, sum_partitions);
    if (sum_partitions)
        synchronize_cartesian_commit(staged,
            raw_residual,
            external_load_residual,
            conservation,
            first_contribution,
            last_contribution,
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

void PlaneBackend::export_histories(TransientCommittedState& state) const {
    state.cartesian_material_histories = _histories;
    state.contact_histories = committed_contact_histories();
}

void PlaneBackend::restore_histories(TransientCommittedState& state) {
    if (!state.material_histories.empty() || !state.quad8_material_histories.empty()
        || !state.radial_material_histories.empty())
        throw std::invalid_argument("Unexpected material history for this spatial backend");
    if (state.cartesian_material_histories.size() != _histories.size())
        throw std::invalid_argument("Committed material region layout mismatch");
    for (std::size_t r = 0; r < _histories.size(); ++r) {
        if (state.cartesian_material_histories[r].size() != _histories[r].size())
            throw std::invalid_argument("Committed material element layout mismatch");
        for (std::size_t e = 0; e < _histories[r].size(); ++e) {
            if (state.cartesian_material_histories[r][e].size() != _histories[r][e].size())
                throw std::invalid_argument("Committed material point layout mismatch");
            for (const auto& point : state.cartesian_material_histories[r][e])
                if (!valid_material_state(point))
                    throw std::invalid_argument("Invalid committed material state");
        }
    }
    validate_state(state.solution);
    restore_contact_state(state.solution, std::move(state.contact_histories));
    _histories = std::move(state.cartesian_material_histories);
}

RegionStateSummary PlaneBackend::summarize_region(std::size_t region) const {
    if (region >= _spatial.region_count())
        throw std::out_of_range("TransientProblem region summary index is out of range");
    const auto& fields = _spatial.field_layout();
    const auto temperature = std::find_if(fields.begin(), fields.end(), [](const FieldDescriptor& field) {
        return field.category == FieldCategory::thermal;
    });
    RegionStateSummary result{-std::numeric_limits<double>::infinity(), 0.0, 0.0};
    for (std::size_t node = 0; node < _spatial.source_nodes(region).size(); ++node) {
        if (!_spatial.thermal_nodes(region)[node])
            continue;
        result.maximum_temperature = std::max(result.maximum_temperature,
            _state.committed_solution.at(temperature->begin + _spatial.global_temperature_node(region, node)));
    }
    for (const auto& element : _histories[region])
        for (const auto& point : element) {
            result.maximum_equivalent_plastic_strain =
                std::max(result.maximum_equivalent_plastic_strain, point.equivalent_plastic_strain);
            result.maximum_equivalent_creep_strain =
                std::max(result.maximum_equivalent_creep_strain, point.equivalent_creep_strain);
        }
    return result;
}

std::vector<double> PlaneBackend::committed_creep_rates() const {
    std::vector<double> rates;
    bool has_creep = false;
    for (std::size_t r = 0; r < _spatial.definition().regions.size(); ++r) {
        const auto& region = _spatial.definition().regions[r];
        const IsotropicThermoelasticMaterial material(region.material);
        has_creep = has_creep || material.functions().has_creep();
        for (std::size_t e = 0; e < _spatial.region_element_count(r); ++e) {
            const auto index = _spatial.region_element_offset(r) + e;
            std::vector<std::size_t> dofs;
            contribution_dofs(index, dofs);
            std::vector<double> local(dofs.size());
            for (std::size_t i = 0; i < dofs.size(); ++i)
                local[i] = _state.committed_solution[dofs[i]];
            std::vector<double> values;
            const double time = _state.committed_time;

            elements::Cpeg8Values state{};
            std::copy(local.begin(), local.end(), state.begin());
            values = elements::cpeg8t_creep_rates(material, _spatial.geometry(r, e), state, _histories[r][e], time);
            rates.insert(rates.end(), values.begin(), values.end());
        }
    }
    if (!has_creep)
        throw std::invalid_argument("Creep-rate time control requires a creep material");
    return rates;
}

void PlaneBackend::publish_material_history() noexcept {
    _histories.swap(_staged);
}
} // namespace fuelsim
