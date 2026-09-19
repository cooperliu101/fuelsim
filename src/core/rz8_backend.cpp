#include "cax8rt.hpp"
#include "cax8t.hpp"
#include "core/cax_evaluation.hpp"
#include "spatial_backend_common.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
void Rz8Backend::accumulate_time_error(const TransientCommittedState& full,
    const TransientCommittedState& half,
    const TransientTimeOptions& options,
    TransientTimeErrorEstimate& result) const {
    if (full.quad8_material_histories.size() != half.quad8_material_histories.size())
        throw std::logic_error("step-doubling material region layouts differ");
    MaterialTimeErrors material;
    for (std::size_t r = 0; r < full.quad8_material_histories.size(); ++r) {
        const auto& first = full.quad8_material_histories[r];
        const auto& second = half.quad8_material_histories[r];
        if (first.size() != second.size())
            throw std::logic_error("step-doubling material element layouts differ");
        for (std::size_t e = 0; e < first.size(); ++e) {
            const std::size_t count = _spatial.region_element_geometry(r, e).point_count;
            for (std::size_t q = 0; q < count; ++q)
                accumulate_material_time_error(material, first[e][q], second[e][q]);
        }
    }
    assign_material_time_errors(result, material, options);
    accumulate_axisymmetric_contact_time_error(full, half, options, result, false);
}

Rz8Backend::Rz8Backend(SpatialTimeState& state, SpatialDefinition definition, const UnstructuredQuad8Mesh& mesh, bool)
    : SpatialBackend(state), _spatial(std::move(definition), mesh) {

    for (std::size_t r = 0; r < _spatial.region_count(); ++r) {
        const auto& value = _spatial.region(r);
        _kernel_data.push_back({IsotropicThermoelasticMaterial(value.material),
            _spatial.region_heat_source(r),
            0.0,
            value.strain_formulation,
            value.rz_element_formulation,
            value.initial_temperature,
            value.body_acceleration});
    }
    {
        _histories.resize(_spatial.region_count());
        for (std::size_t r = 0; r < _spatial.region_count(); ++r)
            _histories[r].resize(_spatial.region_element_count(r));
    }
}

void Rz8Backend::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _spatial.contribution_dofs(index, dofs);
}

void Rz8Backend::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _spatial.sparsity_contribution_dofs(index, dofs);
}

void Rz8Backend::set_time(double value) {
    _spatial.set_time(value);
    for (auto& data : _kernel_data)
        data.time = value;
    for (std::size_t r = 0; r < _spatial.region_count(); ++r)
        _kernel_data[r].volumetric_heat_source = _spatial.region_heat_source(r);
}

void Rz8Backend::set_load_factor(double value) {
    _spatial.set_load_factor(value);
    for (std::size_t r = 0; r < _spatial.region_count(); ++r)
        _kernel_data[r].volumetric_heat_source = _spatial.region_heat_source(r);
}

void Rz8Backend::set_heat_source_interval(double begin, double end) {
    for (std::size_t r = 0; r < _spatial.region_count(); ++r)
        _kernel_data[r].volumetric_heat_source = _spatial.region_heat_source_average(r, begin, end);
}

void Rz8Backend::set_region_strain_formulation(std::size_t region, StrainFormulation formulation) {
    _spatial.set_region_strain_formulation(region, formulation);
    _kernel_data[region].strain_formulation = formulation;
}

void Rz8Backend::commit_contact_state(const std::vector<double>& state) {
    (void)_spatial.commit_contact_state(state);
}

void Rz8Backend::compute_contribution(std::size_t index,
    const std::vector<double>& state,
    std::vector<double>& residual,
    std::vector<double>* jacobian,
    bool transient) const {
    if (index >= _spatial.volume_contribution_count())
        return _spatial.compute_boundary(index, state, residual, jacobian);
    if (state.size() != 20)
        throw std::invalid_argument("CAX8T volume state must contain 20 degrees of freedom");
    Quad8RzValues local{}, old{};
    std::copy(state.begin(), state.end(), local.begin());
    if (transient) {
        std::vector<std::size_t> dofs;
        _spatial.contribution_dofs(index, dofs);
        for (std::size_t i = 0; i < 20; ++i)
            old[i] = _state.committed_solution[dofs[i]];
    }
    const auto [r, e] = _spatial.element_location(index);
    const auto result = compute_cax8(_kernel_data[r],
        _spatial.region_element_geometry(r, e),
        local,
        old,
        transient ? &_histories[r][e] : nullptr,
        transient ? _state.active_time_step : 0.0,
        transient && _state.include_thermal_time_term,
        {true, jacobian != nullptr, false, false});
    residual.assign(result.residual.begin(), result.residual.end());
    if (jacobian)
        jacobian->assign(result.jacobian.begin(), result.jacobian.end());
    return;
}

void Rz8Backend::prepare_time_step(const NonlinearProblem& problem,
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
    std::vector<std::size_t> dofs;
    try {
        if (!partition_failure) {
            for (std::size_t index = 0; index < _spatial.volume_contribution_count(); ++index) {
                if (index < first_contribution || index >= last_contribution)
                    continue;
                const auto [r, e] = _spatial.element_location(index);
                _spatial.contribution_dofs(index, dofs);
                Quad8RzValues current{}, old{};
                for (std::size_t i = 0; i < 20; ++i) {
                    current[i] = converged_solution[dofs[i]];
                    old[i] = _state.committed_solution[dofs[i]];
                }
                const auto& geometry = _spatial.region_element_geometry(r, e);
                if (_spatial.region(r).body_acceleration != std::array<double, 3>{})
                    for (std::size_t q = 0; q < geometry.point_count; ++q) {
                        const auto& point = geometry.points[q];
                        std::array<double, 3> increment{};
                        for (std::size_t node = 0; node < 8; ++node) {
                            increment[0] += point.shape[node] * (current[4 + node] - old[4 + node]);
                            increment[2] += point.shape[node] * (current[12 + node] - old[12 + node]);
                        }
                        conservation.body_force_work_increment += body_point_work(_spatial.region(r),
                            {0.0, point.radius, 0.0, point.axial_coordinate},
                            point.weighted_measure,
                            increment);
                    }

                const auto update = compute_cax8(_kernel_data[r],
                    geometry,
                    current,
                    old,
                    &_histories[r][e],
                    _state.active_time_step,
                    _state.include_thermal_time_term,
                    {true, false, true, false});
                conservation.stored_heat_rate += update.stored_heat_rate;
                conservation.generated_heat_rate += update.generated_heat_rate;
                for (std::size_t q = 0; q < geometry.point_count; ++q) {
                    rz::accumulate_material_conservation(conservation,
                        _histories[r][e][q],
                        update.history[q],
                        geometry.points[q].weighted_measure);
                }
                staged[r][e] = update.history;
            }
        }
    } catch (...) {
        if (!sum_partitions)
            throw;
        partition_failure = std::current_exception();
    }
    prepare_contacts(converged_solution, partition_failure, sum_partitions);
    if (sum_partitions) {
        std::vector<MaterialPointState*> points;
        for (auto& region : staged)
            for (auto& history : region)
                for (auto& point : history)
                    points.push_back(&point);
        synchronize_rz_commit(points,
            9,
            raw_residual,
            external_load_residual,
            conservation,
            first_contribution,
            last_contribution,
            sum_partitions,
            partition_failure);
    }
    finalize_conservation(problem,
        converged_solution,
        _state.committed_solution,
        raw_residual,
        _state.committed_raw_residual,
        external_load_residual,
        _state.committed_external_load_residual,
        conservation);
}

void Rz8Backend::export_histories(TransientCommittedState& state) const {
    state.quad8_material_histories = _histories;
    state.contact_histories = committed_contact_histories();
}

void Rz8Backend::restore_histories(TransientCommittedState& state) {
    if (!state.material_histories.empty() || !state.cartesian_material_histories.empty()
        || !state.radial_material_histories.empty())
        throw std::invalid_argument("Unexpected material history for this spatial backend");
    if (state.quad8_material_histories.size() != _histories.size())
        throw std::invalid_argument("Committed material region layout mismatch");
    for (std::size_t r = 0; r < _histories.size(); ++r) {
        if (state.quad8_material_histories[r].size() != _histories[r].size())
            throw std::invalid_argument("Committed material element layout mismatch");
        for (std::size_t e = 0; e < _histories[r].size(); ++e) {
            if (state.quad8_material_histories[r][e].size() != _histories[r][e].size())
                throw std::invalid_argument("Committed material point layout mismatch");
            for (const auto& point : state.quad8_material_histories[r][e])
                if (!valid_material_state(point))
                    throw std::invalid_argument("Invalid committed material state");
        }
    }
    validate_state(state.solution);
    restore_contact_state(state.solution, std::move(state.contact_histories));
    _histories = std::move(state.quad8_material_histories);
}

RegionStateSummary Rz8Backend::summarize_region(std::size_t region) const {
    if (region >= _spatial.region_count())
        throw std::out_of_range("TransientProblem region summary index is out of range");
    const auto& fields = _spatial.field_layout();
    const auto temperature = std::find_if(fields.begin(), fields.end(), [](const FieldDescriptor& field) {
        return field.category == FieldCategory::thermal;
    });
    RegionStateSummary result{-std::numeric_limits<double>::infinity(), 0.0, 0.0};
    for (std::size_t node = 0; node < _spatial.region_mesh(region).nodes().size(); ++node) {
        if (!_spatial.region_mesh(region).temperature_nodes()[node])
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

std::vector<double> Rz8Backend::committed_creep_rates() const {
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

            Quad8RzValues state{};
            std::copy(local.begin(), local.end(), state.begin());
            const auto& geometry = _spatial.region_element_geometry(r, e);
            const auto& history = _histories[r][e];
            values = region.rz_element_formulation == RzElementFormulation::cax8t
                         ? elements::cax8t_creep_rates(material, geometry, state, history, time)
                         : elements::cax8rt_creep_rates(material, geometry, state, history, time);
            rates.insert(rates.end(), values.begin(), values.end());
        }
    }
    if (!has_creep)
        throw std::invalid_argument("Creep-rate time control requires a creep material");
    return rates;
}

void Rz8Backend::publish_material_history() noexcept {
    _histories.swap(_staged);
}
} // namespace fuelsim
