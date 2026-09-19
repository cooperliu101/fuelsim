#include "spatial_backend_common.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
void RadialBackend::accumulate_time_error(const TransientCommittedState& full,
    const TransientCommittedState& half,
    const TransientTimeOptions& options,
    TransientTimeErrorEstimate& result) const {
    if (full.radial_material_histories.size() != half.radial_material_histories.size())
        throw std::logic_error("step-doubling material region layouts differ");
    MaterialTimeErrors material;
    for (std::size_t r = 0; r < full.radial_material_histories.size(); ++r) {
        const auto& first = full.radial_material_histories[r];
        const auto& second = half.radial_material_histories[r];
        if (first.size() != second.size())
            throw std::logic_error("step-doubling material element layouts differ");
        for (std::size_t e = 0; e < first.size(); ++e) {
            const std::size_t count = cax2t_gps_material_point_count;
            for (std::size_t q = 0; q < count; ++q)
                accumulate_material_time_error(material, first[e][q], second[e][q]);
        }
    }
    assign_material_time_errors(result, material, options);
    accumulate_axisymmetric_contact_time_error(full, half, options, result, true);
}

RadialBackend::RadialBackend(SpatialTimeState& state,
    SpatialDefinition definition,
    const UnstructuredBar2Mesh& mesh,
    bool)
    : SpatialBackend(state), _spatial(std::move(definition), mesh) {
    {
        _histories.resize(_spatial.region_count());
        for (std::size_t r = 0; r < _spatial.region_count(); ++r)
            _histories[r].resize(_spatial.region_element_count(r));
    }
}

void RadialBackend::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    _spatial.contribution_dofs(index, dofs);
}

void RadialBackend::set_time(double value) {
    _spatial.set_time(value);
    _time = value;
}

void RadialBackend::set_load_factor(double value) {
    _spatial.set_load_factor(value);
}

void RadialBackend::set_heat_source_interval(double begin, double end) {
    _spatial.set_heat_source_interval(begin, end);
}

void RadialBackend::set_region_strain_formulation(std::size_t region, StrainFormulation formulation) {
    _spatial.set_region_strain_formulation(region, formulation);
}

void RadialBackend::commit_contact_state(const std::vector<double>& state) {
    (void)_spatial.commit_contact_state(state);
}

void RadialBackend::compute_contribution(std::size_t index,
    const std::vector<double>& state,
    std::vector<double>& residual,
    std::vector<double>* jacobian,
    bool transient) const {
    const Cax2tGpsMaterialHistory* history = nullptr;
    if (transient && index < _spatial.volume_contribution_count()) {
        const auto [region, element] = _spatial.element_location(index);
        history = &_histories[region][element];
    }
    radial::LocalContribution contribution;
    _spatial.compute_contribution(index,
        state,
        _state.committed_solution,
        history,
        transient ? _state.active_time_step : 0.0,
        transient ? _state.active_end_time : _time,
        transient && _state.include_thermal_time_term,
        jacobian != nullptr,
        contribution);
    residual = std::move(contribution.residual);
    if (jacobian)
        *jacobian = std::move(contribution.jacobian);
    return;
}

void RadialBackend::prepare_time_step(const NonlinearProblem& problem,
    const std::vector<double>& converged_solution,
    std::size_t first_contribution,
    std::size_t last_contribution,
    const std::function<void(std::vector<double>&)>& sum_partitions,
    TransientConservationSummary& conservation,
    std::vector<double>& raw_residual,
    std::vector<double>& external_load_residual,
    std::exception_ptr partition_failure) {
    (void)first_contribution;
    (void)last_contribution;
    (void)sum_partitions;
    (void)partition_failure;

    _staged = _histories;
    auto& staged = _staged;
    constexpr double gauss = 0.577350269189625764509148780501957456;
    constexpr double pi = 3.141592653589793238462643383279502884;
    const std::array<double, 2> stations = {-gauss, gauss};
    for (std::size_t region = 0; region < staged.size(); ++region)
        for (std::size_t element = 0; element < staged[region].size(); ++element) {
            const auto update = _spatial.evaluate_volume(region,
                element,
                converged_solution,
                _state.committed_solution,
                &_histories[region][element],
                _state.active_time_step,
                _state.active_end_time,
                _state.include_thermal_time_term,
                false);
            conservation.stored_heat_rate += update.stored_heat_rate;
            conservation.generated_heat_rate += update.generated_heat_rate;
            const auto& geometry = _spatial.region_element_geometry(region, element);
            const bool finite = _spatial.region(region).strain_formulation == StrainFormulation::finite;
            Cax2tGpsLocalValues current{}, old{};
            if (finite || _spatial.region(region).body_acceleration != std::array<double, 3>{}) {
                const std::size_t index = _spatial.region_element_offset(region) + element;
                current = _spatial.volume_state(index, converged_solution);
                old = _spatial.volume_state(index, _state.committed_solution);
            }
            for (std::size_t q = 0; q < stations.size(); ++q) {
                const std::array<double, 2> shape = {0.5 * (1.0 - stations[q]), 0.5 * (1.0 + stations[q])};
                const double radius = shape[0] * geometry.radii[0] + shape[1] * geometry.radii[1];
                const double thickness = geometry.radii[1] - geometry.radii[0];
                const double height = geometry.z_upper - geometry.z_lower;
                const std::array<double, 3> increment{shape[0] * (current[2] - old[2])
                                                          + shape[1] * (current[3] - old[3]),
                    0.0,
                    0.5 * (current[4] - old[4] + current[5] - old[5])};
                conservation.body_force_work_increment += body_point_work(_spatial.region(region),
                    {0.0, radius, 0.0, 0.5 * (geometry.z_lower + geometry.z_upper)},
                    pi * radius * thickness * height,
                    increment);
                double current_measure = pi * radius * thickness * height;
                double committed_measure = current_measure;
                if (finite) {
                    current_measure *= (1.0 + (current[3] - current[2]) / thickness)
                                       * (1.0 + (current[5] - current[4]) / height)
                                       * (1.0 + (current[2] + current[3]) / (geometry.radii[0] + geometry.radii[1]));
                    committed_measure *= (1.0 + (old[3] - old[2]) / thickness) * (1.0 + (old[5] - old[4]) / height)
                                         * (1.0 + (old[2] + old[3]) / (geometry.radii[0] + geometry.radii[1]));
                }
                const auto& previous = _histories[region][element][q];
                rz::accumulate_material_conservation(conservation, previous, update.history[q], current_measure);
                // The diagonal kinematics have no objective rotation. The
                // elastic energies use their respective assumed-strain
                // mechanical measures; dissipation uses the current
                // mechanical measure and trapezoidal stress.
                if (finite)
                    conservation.elastic_energy_change +=
                        0.5 * (current_measure - committed_measure)
                        * rz::stress_strain_inner_product(previous.stress, previous.elastic_strain);
            }
            staged[region][element] = update.history;
        }
    prepare_contacts(converged_solution, partition_failure, sum_partitions);
    finalize_conservation(problem,
        converged_solution,
        _state.committed_solution,
        raw_residual,
        _state.committed_raw_residual,
        external_load_residual,
        _state.committed_external_load_residual,
        conservation);
}

void RadialBackend::export_histories(TransientCommittedState& state) const {
    state.radial_material_histories = _histories;
    state.contact_histories = committed_contact_histories();
}

void RadialBackend::restore_histories(TransientCommittedState& state) {
    if (!state.material_histories.empty() || !state.cartesian_material_histories.empty()
        || !state.quad8_material_histories.empty())
        throw std::invalid_argument("Unexpected material history for this spatial backend");
    if (state.radial_material_histories.size() != _histories.size())
        throw std::invalid_argument("Committed material region layout mismatch");
    for (std::size_t r = 0; r < _histories.size(); ++r) {
        if (state.radial_material_histories[r].size() != _histories[r].size())
            throw std::invalid_argument("Committed material element layout mismatch");
        for (std::size_t e = 0; e < _histories[r].size(); ++e) {
            if (state.radial_material_histories[r][e].size() != _histories[r][e].size())
                throw std::invalid_argument("Committed material point layout mismatch");
            for (const auto& point : state.radial_material_histories[r][e])
                if (!valid_material_state(point) || !finite_stress(point.stress))
                    throw std::invalid_argument("Invalid committed material state");
        }
    }
    validate_state(state.solution);
    restore_contact_state(state.solution, std::move(state.contact_histories));
    _histories = std::move(state.radial_material_histories);
}

RegionStateSummary RadialBackend::summarize_region(std::size_t region) const {
    if (region >= _spatial.region_count())
        throw std::out_of_range("TransientProblem region summary index is out of range");
    const auto& fields = _spatial.field_layout();
    const auto temperature = std::find_if(fields.begin(), fields.end(), [](const FieldDescriptor& field) {
        return field.category == FieldCategory::thermal;
    });
    RegionStateSummary result{-std::numeric_limits<double>::infinity(), 0.0, 0.0};
    for (std::size_t node = 0; node < _spatial.region_source_node_ids(region).size(); ++node) {
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

std::vector<double> RadialBackend::committed_creep_rates() const {
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

            Cax2tGpsLocalValues state{};
            std::copy(local.begin(), local.end(), state.begin());
            values = elements::cax2t_gps_creep_rates(material,
                _spatial.region_element_geometry(r, e),
                state,
                _histories[r][e],
                time);
            rates.insert(rates.end(), values.begin(), values.end());
        }
    }
    if (!has_creep)
        throw std::invalid_argument("Creep-rate time control requires a creep material");
    return rates;
}

void RadialBackend::publish_material_history() noexcept {
    _histories.swap(_staged);
}
} // namespace fuelsim
