#include "cax4rt.hpp"
#include "cax4t.hpp"
#include "core/cax_evaluation.hpp"
#include "spatial_backend_common.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
Rz4Backend::Rz4Backend(SpatialTimeState& state,
    SpatialDefinition definition,
    const UnstructuredQuad4Mesh& mesh,
    bool transient)
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
    if (transient) {
        _histories.resize(_spatial.region_count());
        for (std::size_t r = 0; r < _spatial.region_count(); ++r)
            _histories[r].resize(_spatial.region_element_count(r));
        _staged = _histories;
    }
}

void Rz4Backend::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    const auto fixed = _spatial.contribution_dofs(index);
    dofs.assign(fixed.begin(), fixed.end());
}

void Rz4Backend::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    const auto fixed = _spatial.sparsity_contribution_dofs(index);
    dofs.assign(fixed.begin(), fixed.end());
}

void Rz4Backend::set_time(double value) {
    _spatial.set_time(value);
    for (auto& data : _kernel_data)
        data.time = value;
    for (std::size_t r = 0; r < _spatial.region_count(); ++r)
        _kernel_data[r].volumetric_heat_source = _spatial.region_heat_source(r);
}

void Rz4Backend::set_load_factor(double value) {
    _spatial.set_load_factor(value);
    for (std::size_t r = 0; r < _spatial.region_count(); ++r)
        _kernel_data[r].volumetric_heat_source = _spatial.region_heat_source(r);
}

void Rz4Backend::set_heat_source_interval(double begin, double end) {
    for (std::size_t r = 0; r < _spatial.region_count(); ++r)
        _kernel_data[r].volumetric_heat_source = _spatial.region_heat_source_average(r, begin, end);
}

void Rz4Backend::set_region_strain_formulation(std::size_t region, StrainFormulation formulation) {
    _spatial.set_region_strain_formulation(region, formulation);
    _kernel_data[region].strain_formulation = formulation;
}

void Rz4Backend::commit_contact_state(const std::vector<double>& state) {
    (void)_spatial.commit_contact_state(state);
}

void Rz4Backend::compute_contribution(std::size_t index,
    const std::vector<double>& state,
    std::vector<double>& residual,
    std::vector<double>* jacobian,
    bool transient) const {
    const Cax4LocalValues local_state = rz_local_values(state);
    Cax4LocalJacobian local_jacobian{};
    Cax4LocalResidual local_residual{};
    if (index >= _spatial.volume_contribution_count())
        local_residual =
            _spatial.compute_contribution(index, local_state, jacobian == nullptr ? nullptr : &local_jacobian);
    else {
        const auto location = _spatial.element_location(index);
        const auto volume = evaluate_cax4(_kernel_data[location.first],
            _spatial.region_element_geometry(location.first, location.second),
            local_state,
            transient ? gather_rz_state(_spatial, index, _state.committed_solution) : Cax4LocalValues{},
            transient ? &_histories[location.first][location.second] : nullptr,
            transient ? _state.active_time_step : 0.0,
            transient && _state.include_thermal_time_term,
            {true, jacobian != nullptr, false, false});
        local_residual = volume.residual;
        local_jacobian = volume.jacobian;
    }
    residual.assign(local_residual.begin(), local_residual.end());
    if (jacobian != nullptr)
        jacobian->assign(local_jacobian.begin(), local_jacobian.end());
}

void Rz4Backend::prepare_time_step(const NonlinearProblem& problem,
    const std::vector<double>& converged_solution,
    std::size_t first_contribution,
    std::size_t last_contribution,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    TransientConservationSummary& conservation,
    std::vector<double>& raw_residual,
    std::vector<double>& external_load_residual,
    std::exception_ptr partition_failure) {
    const std::size_t regions = _spatial.region_count();
    auto& staged = _staged;
    try {
        if (!partition_failure) {
            for (std::size_t region = 0; region < regions; ++region) {
                const std::size_t offset = _spatial.region_element_offset(region);
                for (std::size_t element = 0; element < staged[region].size(); ++element) {
                    if (offset + element < first_contribution || offset + element >= last_contribution)
                        continue;
                    const Cax4LocalValues state = gather_rz_state(_spatial, offset + element, converged_solution);
                    const Cax4LocalValues committed_state =
                        gather_rz_state(_spatial, offset + element, _state.committed_solution);
                    const Quad4RzGeometry& geometry = _spatial.region_element_geometry(region, element);
                    const auto& element_data = _kernel_data[region];
                    const bool reduced = element_data.element_formulation == RzElementFormulation::cax4rt;
                    auto result = evaluate_cax4(element_data,
                        geometry,
                        state,
                        committed_state,
                        &_histories[region][element],
                        _state.active_time_step,
                        _state.include_thermal_time_term,
                        {true, false, true, false});
                    auto update = std::move(result.history);
                    if (element_data.body_acceleration != std::array<double, 3>{})
                        for (const auto& point : geometry.points) {
                            std::array<double, 3> increment{};
                            for (std::size_t node = 0; node < 4; ++node) {
                                increment[0] += point.shape[node] * (state[4 + node] - committed_state[4 + node]);
                                increment[2] += point.shape[node] * (state[8 + node] - committed_state[8 + node]);
                            }
                            conservation.body_force_work_increment += body_point_work(_spatial.region(region),
                                {0.0, point.radius, 0.0, point.axial_coordinate},
                                point.weighted_measure,
                                increment);
                        }
                    conservation.stored_heat_rate += result.stored_heat_rate;
                    conservation.generated_heat_rate += result.generated_heat_rate;
                    if (reduced) {
                        const double current_hourglass =
                            fuelsim::elements::cax4rt_hourglass_energy({_kernel_data[region].material,
                                geometry,
                                state,
                                state,
                                nullptr,
                                0.0,
                                _kernel_data[region].time,
                                _kernel_data[region].volumetric_heat_source,
                                _kernel_data[region].strain_formulation,
                                false,
                                _kernel_data[region].initial_temperature,
                                _kernel_data[region].body_acceleration});
                        const double old_hourglass =
                            fuelsim::elements::cax4rt_hourglass_energy({_kernel_data[region].material,
                                geometry,
                                committed_state,
                                committed_state,
                                nullptr,
                                0.0,
                                _kernel_data[region].time,
                                _kernel_data[region].volumetric_heat_source,
                                _kernel_data[region].strain_formulation,
                                false,
                                _kernel_data[region].initial_temperature,
                                _kernel_data[region].body_acceleration});
                        conservation.mechanical_hourglass_energy += current_hourglass;
                        conservation.mechanical_hourglass_energy_change += current_hourglass - old_hourglass;
                    }
                    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
                        const RzQuadraturePoint& point = geometry.points[q];
                        const MaterialPointState &old_history = _histories[region][element][reduced ? 0 : q],
                                                 &new_history = update[reduced ? 0 : q];
                        rz::accumulate_material_conservation(conservation,
                            old_history,
                            new_history,
                            point.weighted_measure);
                    }
                    staged[region][element] = std::move(update);
                }
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
            4,
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

void Rz4Backend::export_histories(TransientCommittedState& state) const {
    state.material_histories = _histories;
    state.contact_histories = committed_contact_histories();
}

void Rz4Backend::restore_histories(TransientCommittedState& state) {
    if (!state.cartesian_material_histories.empty() || !state.quad8_material_histories.empty()
        || !state.radial_material_histories.empty())
        throw std::invalid_argument("Unexpected material history for this spatial backend");
    if (state.material_histories.size() != _histories.size())
        throw std::invalid_argument("Committed material region layout mismatch");
    for (std::size_t r = 0; r < _histories.size(); ++r) {
        if (state.material_histories[r].size() != _histories[r].size())
            throw std::invalid_argument("Committed material element layout mismatch");
        for (std::size_t e = 0; e < _histories[r].size(); ++e) {
            if (state.material_histories[r][e].size() != _histories[r][e].size())
                throw std::invalid_argument("Committed material point layout mismatch");
            for (const auto& point : state.material_histories[r][e])
                if (!valid_material_state(point) || !finite_stress(point.stress))
                    throw std::invalid_argument("Invalid committed material state");
        }
    }
    for (std::size_t node = 0; node < _spatial.node_count(); ++node) {
        const double temperature = state.solution.at(_spatial.dof(Field::temperature, node));
        if (!std::isfinite(temperature) || !(temperature > 0.0)
            || !std::isfinite(state.solution.at(_spatial.dof(Field::radial_displacement, node)))
            || !std::isfinite(state.solution.at(_spatial.dof(Field::axial_displacement, node))))
            throw std::invalid_argument("Transient committed nodal state must be finite with positive temperatures");
    }
    restore_contact_state(state.solution, std::move(state.contact_histories));
    _histories = std::move(state.material_histories);
}

RegionStateSummary Rz4Backend::summarize_region(std::size_t region) const {
    if (region >= _spatial.region_count())
        throw std::out_of_range("TransientProblem region summary index is out of range");
    const auto& fields = _spatial.field_layout();
    const auto temperature = std::find_if(fields.begin(), fields.end(), [](const FieldDescriptor& field) {
        return field.category == FieldCategory::thermal;
    });
    RegionStateSummary result{-std::numeric_limits<double>::infinity(), 0.0, 0.0};
    for (std::size_t node = 0; node < _spatial.region_mesh(region).nodes().size(); ++node) {
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

std::vector<double> Rz4Backend::committed_creep_rates() const {
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

            const auto state = rz_local_values(local);
            const auto& geometry = _spatial.region_element_geometry(r, e);
            const auto& history = _histories[r][e];
            values =
                region.rz_element_formulation == RzElementFormulation::cax4t
                    ? elements::cax4t_creep_rates(material, geometry, state, history, time)
                    : elements::cax4rt_creep_rates(material, geometry, state, history, time, region.strain_formulation);
            rates.insert(rates.end(), values.begin(), values.end());
        }
    }
    if (!has_creep)
        throw std::invalid_argument("Creep-rate time control requires a creep material");
    return rates;
}

void Rz4Backend::publish_material_history() noexcept {
    _histories.swap(_staged);
}
} // namespace fuelsim
