#include "c3d20rt.hpp"
#include "c3d20t.hpp"
#include "c3d8rt.hpp"
#include "c3d8t.hpp"
#include "spatial_backend_common.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
void CartesianBackend::accumulate_time_error(const TransientCommittedState& full,
    const TransientCommittedState& half,
    const TransientTimeOptions& options,
    TransientTimeErrorEstimate& result) const {
    accumulate_cartesian_material_time_error(full, half, options, result);
    accumulate_cartesian_contact_time_error(full, half, options, result);
}

CartesianBackend::CartesianBackend(SpatialTimeState& state,
    SpatialDefinition definition,
    const UnstructuredHex8Mesh& mesh,
    bool transient)
    : SpatialBackend(state), _spatial(std::move(definition), mesh) {

    if (!transient) {
        for (std::size_t region = 0; region < _spatial.region_count(); ++region) {
            const MaterialFunctionSet& functions = *_spatial.region(region).material.functions;
            if (functions.has_creep() || functions.has_plasticity())
                throw std::invalid_argument("Steady Cartesian three-dimensional problems support only elasticity");
        }
        return;
    }
    _histories.resize(_spatial.region_count());
    _staged.resize(_spatial.region_count());
    for (std::size_t region = 0; region < _spatial.region_count(); ++region) {
        const std::size_t points = _spatial.region_material_point_count(region);
        _histories[region].resize(_spatial.region_element_count(region));
        _staged[region].resize(_spatial.region_element_count(region));
        for (CartesianMaterialHistory& history : _histories[region])
            history.resize(points);
        for (CartesianMaterialHistory& history : _staged[region])
            history.resize(points);
    }
}

CartesianBackend::CartesianBackend(SpatialTimeState& state,
    SpatialDefinition definition,
    const UnstructuredHex20Mesh& mesh,
    bool transient)
    : SpatialBackend(state), _spatial(std::move(definition), mesh) {

    if (!transient) {
        for (std::size_t region = 0; region < _spatial.region_count(); ++region) {
            const MaterialFunctionSet& functions = *_spatial.region(region).material.functions;
            if (functions.has_creep() || functions.has_plasticity())
                throw std::invalid_argument("Steady Cartesian three-dimensional problems support only elasticity");
        }
        return;
    }
    _histories.resize(_spatial.region_count());
    _staged.resize(_spatial.region_count());
    for (std::size_t region = 0; region < _spatial.region_count(); ++region) {
        const std::size_t points = _spatial.region_material_point_count(region);
        _histories[region].resize(_spatial.region_element_count(region));
        _staged[region].resize(_spatial.region_element_count(region));
        for (CartesianMaterialHistory& history : _histories[region])
            history.resize(points);
        for (CartesianMaterialHistory& history : _staged[region])
            history.resize(points);
    }
}

void CartesianBackend::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _spatial.contribution_dofs(index, dofs);
}

void CartesianBackend::sparsity_contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _spatial.sparsity_contribution_dofs(index, dofs);
}

void CartesianBackend::set_time(double value) {
    _spatial.set_time(value);
}

void CartesianBackend::set_load_factor(double value) {
    _spatial.set_load_factor(value);
}

void CartesianBackend::set_heat_source_interval(double begin, double end) {
    _spatial.set_heat_source_interval(begin, end);
}

void CartesianBackend::set_region_strain_formulation(std::size_t region, StrainFormulation formulation) {
    _spatial.set_region_strain_formulation(region, formulation);
}

void CartesianBackend::commit_contact_state(const std::vector<double>& state) {
    (void)_spatial.commit_contact_state(state);
}

void CartesianBackend::compute_contribution(std::size_t index,
    const std::vector<double>& state,
    std::vector<double>& residual,
    std::vector<double>* jacobian,
    bool transient) const {
    const CartesianMaterialHistory* history = nullptr;
    if (transient && index < _spatial.volume_contribution_count()) {
        const auto location = _spatial.element_location(index);
        history = &_histories[location.first][location.second];
    }
    _spatial.compute_contribution(index,
        state,
        transient ? &_state.committed_solution : nullptr,
        history,
        transient ? _state.active_time_step : 0.0,
        residual,
        jacobian,
        transient ? _state.include_thermal_time_term : true);
    return;
}

void CartesianBackend::prepare_time_step(const NonlinearProblem& problem,
    const std::vector<double>& converged_solution,
    std::size_t first_contribution,
    std::size_t last_contribution,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    TransientConservationSummary& conservation,
    std::vector<double>& raw_residual,
    std::vector<double>& external_load_residual,
    std::exception_ptr partition_failure) {
    auto& staged = _staged;
    try {
        if (!partition_failure)
            for (std::size_t region = 0; region < _spatial.region_count(); ++region) {
                const std::size_t offset = _spatial.region_element_offset(region);
                for (std::size_t element = 0; element < _spatial.region_element_count(region); ++element) {
                    if (offset + element < first_contribution || offset + element >= last_contribution)
                        continue;
                    if (_spatial.uses_hex20()) {
                        const Hex20LocalValues current =
                            _spatial.hex20_volume_state(offset + element, converged_solution);
                        const Hex20LocalValues old =
                            _spatial.hex20_volume_state(offset + element, _state.committed_solution);
                        const Hex20Geometry& geometry = _spatial.hex20_region_element_geometry(region, element);
                        if (_spatial.region(region).body_acceleration != std::array<double, 3>{})
                            for (const auto& point : geometry.mechanical_points) {
                                std::array<double, 3> increment{};
                                for (std::size_t component = 0; component < 3; ++component)
                                    for (std::size_t node = 0; node < 20; ++node)
                                        increment[component] +=
                                            point.displacement_shape[node]
                                            * (current[8 + component * 20 + node] - old[8 + component * 20 + node]);
                                conservation.body_force_work_increment += body_point_work(_spatial.region(region),
                                    {0.0, point.position.x, point.position.y, point.position.z},
                                    point.weighted_measure,
                                    increment);
                            }
                        CartesianMaterialHistory update = _spatial.transient_update(region,
                            element,
                            current,
                            old,
                            _histories[region][element],
                            _state.active_time_step);
                        const bool finite = _spatial.region(region).strain_formulation == StrainFormulation::finite;
                        const double heat_source =
                            _spatial.region_heat_source_average(region, _state.committed_time, _state.active_end_time);
                        const auto accumulate_stored_heat = [&](const std::array<double, 8>& shape,
                                                                const CartesianPoint3& position,
                                                                double reference_measure) {
                            if (!_state.include_thermal_time_term)
                                return;
                            double current_temperature = 0.0, old_temperature = 0.0;
                            for (std::size_t node = 0; node < hex20_temperature_node_count; ++node) {
                                current_temperature += shape[node] * current[node];
                                old_temperature += shape[node] * old[node];
                            }
                            conservation.stored_heat_rate +=
                                reference_measure
                                * _spatial.reference_heat_capacity(region, current_temperature, position)
                                * (current_temperature - old_temperature) / _state.active_time_step;
                        };
                        if (!finite) {
                            for (const Hex20ThermalQuadraturePoint& point : geometry.thermal_points) {
                                accumulate_stored_heat(point.temperature_shape, point.position, point.weighted_measure);
                                conservation.generated_heat_rate += point.weighted_measure * heat_source;
                            }
                        } else if (heat_source != 0.0) {
                            // The finite source uses the linear corner geometry, independently of the midside
                            // motion.
                            const Hex20RegionMesh& mesh = _spatial.hex20_region_mesh(region);
                            Hex8Coordinates current_corners{};
                            for (std::size_t node = 0; node < hex20_temperature_node_count; ++node) {
                                const CartesianPoint3& reference = mesh.nodes()[mesh.elements()[element].nodes[node]];
                                current_corners[node] = {reference.x + current[8 + node],
                                    reference.y + current[28 + node],
                                    reference.z + current[48 + node]};
                            }
                            // Two-point Gauss integration is exact for this trilinear geometric volume, as is
                            // the local source's three-point rule in C3D20T (two-point rule in C3D20RT).
                            conservation.generated_heat_rate +=
                                elements::make_c3d8t_geometry(current_corners).reference_volume * heat_source;
                        }
                        for (std::size_t q = 0; q < geometry.mechanical_points.size(); ++q) {
                            const Hex20MechanicalQuadraturePoint& point = geometry.mechanical_points[q];
                            if (finite)
                                accumulate_stored_heat(point.temperature_shape, point.position, point.weighted_measure);
                            const CartesianMaterialPointState& old_history = _histories[region][element][q];
                            const CartesianMaterialPointState& new_history = update[q];
                            conservation.elastic_energy_change +=
                                0.5 * point.weighted_measure
                                * (cartesian::stress_strain_inner_product(new_history.stress,
                                       new_history.elastic_strain)
                                    - cartesian::stress_strain_inner_product(old_history.stress,
                                        old_history.elastic_strain));
                            conservation.plastic_dissipation_increment +=
                                point.weighted_measure
                                * cartesian::trapezoidal_stress_strain_inner_product(old_history.stress,
                                    new_history.stress,
                                    cartesian::strain_difference(new_history.plastic_strain,
                                        old_history.plastic_strain));
                            conservation.creep_dissipation_increment +=
                                point.weighted_measure
                                * cartesian::trapezoidal_stress_strain_inner_product(old_history.stress,
                                    new_history.stress,
                                    cartesian::strain_difference(new_history.creep_strain, old_history.creep_strain));
                        }
                        staged[region][element] = std::move(update);
                        continue;
                    }
                    const Hex8LocalValues current = _spatial.volume_state(offset + element, converged_solution);
                    const Hex8LocalValues old = _spatial.volume_state(offset + element, _state.committed_solution);
                    const Hex8Geometry& geometry = _spatial.region_element_geometry(region, element);
                    if (_spatial.region(region).body_acceleration != std::array<double, 3>{})
                        for (const auto& point : geometry.points) {
                            std::array<double, 3> increment{};
                            for (std::size_t component = 0; component < 3; ++component)
                                for (std::size_t node = 0; node < 8; ++node)
                                    increment[component] +=
                                        point.shape[node]
                                        * (current[8 + component * 8 + node] - old[8 + component * 8 + node]);
                            conservation.body_force_work_increment += body_point_work(_spatial.region(region),
                                {0.0, point.position.x, point.position.y, point.position.z},
                                point.weighted_measure,
                                increment);
                        }
                    const bool reduced =
                        _spatial.region(region).hex8_element_formulation == Hex8ElementFormulation::c3d8rt;
                    const auto& capacity_points = reduced ? geometry.reduced_capacity_points : geometry.capacity_points;
                    auto element_result = _spatial.transient_update(region,
                        element,
                        current,
                        old,
                        _histories[region][element],
                        _state.active_time_step);
                    auto& update = element_result.history;
                    if (_state.include_thermal_time_term)
                        for (std::size_t node = 0; node < hex8_node_count; ++node) {
                            const Hex8CapacityPoint& point = capacity_points[node];
                            conservation.stored_heat_rate +=
                                point.weighted_measure
                                * _spatial.reference_heat_capacity(region, current[node], point.position)
                                * (current[node] - old[node]) / _state.active_time_step;
                        }
                    const double heat_source =
                        _spatial.region_heat_source_average(region, _state.committed_time, _state.active_end_time);
                    if (heat_source != 0.0) {
                        double source_measure =
                            reduced ? geometry.reduced_body_source_measure : geometry.reference_volume;
                        if (_spatial.region(region).strain_formulation == StrainFormulation::finite) {
                            if (reduced) {
                                Hex8Coordinates current_coordinates{};
                                for (std::size_t node = 0; node < hex8_node_count; ++node) {
                                    const CartesianPoint3& reference = geometry.capacity_points[node].position;
                                    current_coordinates[node] = {reference.x + current[8 + node],
                                        reference.y + current[16 + node],
                                        reference.z + current[24 + node]};
                                }
                                source_measure =
                                    elements::make_c3d8rt_geometry(current_coordinates).reduced_body_source_measure;
                            } else {
                                source_measure = element_result.current_volume;
                            }
                        }
                        conservation.generated_heat_rate += source_measure * heat_source;
                    }
                    if (reduced) {
                        const double current_hourglass = _spatial.mechanical_hourglass_energy(region, element, current);
                        const double old_hourglass = _spatial.mechanical_hourglass_energy(region, element, old);
                        conservation.mechanical_hourglass_energy += current_hourglass;
                        conservation.mechanical_hourglass_energy_change += current_hourglass - old_hourglass;
                    }
                    for (std::size_t q = 0; q < update.size(); ++q) {
                        const Hex8QuadraturePoint& point = reduced ? geometry.reduced_point : geometry.points[q];
                        const CartesianMaterialPointState &old_history = _histories[region][element][q],
                                                          &new_history = update[q];
                        double current_measure = point.weighted_measure, old_measure = point.weighted_measure;
                        SymmetricTensor3Values diagnostic_new_stress = new_history.stress;
                        std::array<double, 6> diagnostic_new_plastic = new_history.plastic_strain;
                        std::array<double, 6> diagnostic_new_creep = new_history.creep_strain;
                        if (_spatial.region(region).strain_formulation == StrainFormulation::finite) {
                            current_measure =
                                point.weighted_measure / geometry.reference_volume * element_result.current_volume;
                            old_measure =
                                point.weighted_measure / geometry.reference_volume * element_result.committed_volume;
                            const auto& rotation = element_result.incremental_rotations[q];
                            const CartesianRotation inverse_rotation = {rotation[0],
                                rotation[3],
                                rotation[6],
                                rotation[1],
                                rotation[4],
                                rotation[7],
                                rotation[2],
                                rotation[5],
                                rotation[8]};
                            diagnostic_new_stress =
                                cartesian::rotate_tensor_values(new_history.stress, inverse_rotation);
                            diagnostic_new_plastic = cartesian::components(
                                cartesian::rotate_tensor_values(new_history.plastic_strain, inverse_rotation));
                            diagnostic_new_creep = cartesian::components(
                                cartesian::rotate_tensor_values(new_history.creep_strain, inverse_rotation));
                        }
                        conservation.elastic_energy_change +=
                            0.5
                            * (current_measure
                                    * cartesian::stress_strain_inner_product(new_history.stress,
                                        new_history.elastic_strain)
                                - old_measure
                                      * cartesian::stress_strain_inner_product(old_history.stress,
                                          old_history.elastic_strain));
                        conservation.plastic_dissipation_increment +=
                            current_measure
                            * cartesian::trapezoidal_stress_strain_inner_product(old_history.stress,
                                diagnostic_new_stress,
                                cartesian::strain_difference(diagnostic_new_plastic, old_history.plastic_strain));
                        conservation.creep_dissipation_increment +=
                            current_measure
                            * cartesian::trapezoidal_stress_strain_inner_product(old_history.stress,
                                diagnostic_new_stress,
                                cartesian::strain_difference(diagnostic_new_creep, old_history.creep_strain));
                    }
                    staged[region][element] = std::move(update);
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

void CartesianBackend::export_histories(TransientCommittedState& state) const {
    state.cartesian_material_histories = _histories;
    state.contact_histories = committed_contact_histories();
}

void CartesianBackend::restore_histories(TransientCommittedState& state) {
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

RegionStateSummary CartesianBackend::summarize_region(std::size_t region) const {
    if (region >= _spatial.region_count())
        throw std::out_of_range("TransientProblem region summary index is out of range");
    const auto& fields = _spatial.field_layout();
    const auto temperature = std::find_if(fields.begin(), fields.end(), [](const FieldDescriptor& field) {
        return field.category == FieldCategory::thermal;
    });
    RegionStateSummary result{-std::numeric_limits<double>::infinity(), 0.0, 0.0};
    for (std::size_t node = 0; node < (_spatial.uses_hex20() ? _spatial.hex20_region_mesh(region).nodes().size()
                                                             : _spatial.region_mesh(region).nodes().size());
        ++node) {
        if (_spatial.uses_hex20() && !_spatial.hex20_region_mesh(region).temperature_nodes().at(node))
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

std::vector<double> CartesianBackend::committed_creep_rates() const {
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

            const auto& history = _histories[r][e];
            if (_spatial.uses_hex20()) {
                Hex20LocalValues state{};
                std::copy(local.begin(), local.end(), state.begin());
                const auto& geometry = _spatial.hex20_region_element_geometry(r, e);
                values = region.hex20_element_formulation == Hex20ElementFormulation::c3d20t
                             ? elements::c3d20t_creep_rates(material, geometry, state, history, time)
                             : elements::c3d20rt_creep_rates(material, geometry, state, history, time);
            } else {
                Hex8LocalValues state{};
                std::copy(local.begin(), local.end(), state.begin());
                const auto& geometry = _spatial.region_element_geometry(r, e);
                if (region.hex8_element_formulation == Hex8ElementFormulation::c3d8t) {
                    const auto full = elements::c3d8t_creep_rates(material, geometry, state, history, time);
                    values.assign(full.begin(), full.end());
                } else {
                    values = elements::c3d8rt_creep_rates(material,
                        geometry,
                        state,
                        history,
                        time,
                        region.strain_formulation);
                }
            }
            rates.insert(rates.end(), values.begin(), values.end());
        }
    }
    if (!has_creep)
        throw std::invalid_argument("Creep-rate time control requires a creep material");
    return rates;
}

void CartesianBackend::publish_material_history() noexcept {
    _histories.swap(_staged);
}
} // namespace fuelsim
