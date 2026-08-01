#include "fuelsim/transient_problem.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

RegionInelasticSummary
summarize_history(const std::vector<Quad4MaterialHistory>& history) noexcept {
    RegionInelasticSummary summary{0.0, 0.0};
    for (const Quad4MaterialHistory& element : history) {
        for (const MaterialPointState& point : element) {
            summary.maximum_equivalent_plastic_strain =
                std::max(summary.maximum_equivalent_plastic_strain,
                         point.equivalent_plastic_strain);
            summary.maximum_equivalent_creep_strain =
                std::max(summary.maximum_equivalent_creep_strain,
                         point.equivalent_creep_strain);
        }
    }
    return summary;
}

void validate_definition(const TransientProblemDefinition& definition) {
    if (definition.regions.size() != definition.spatial.regions.size())
        throw std::invalid_argument(
            "TransientProblem requires one transient material per region");
    for (std::size_t region = 0; region < definition.regions.size(); ++region) {
        if (definition.regions[region].region !=
            definition.spatial.regions[region].name)
            throw std::invalid_argument(
                "TransientProblem region material order does not match the "
                "spatial regions");
    }
}

} // namespace

TransientProblem::TransientProblem(TransientProblemDefinition definition,
                                   const UnstructuredQuad4Mesh& source_mesh)
    : _definition(std::move(definition)),
      _spatial_model(_definition.spatial, source_mesh),
      _committed_solution(_spatial_model.initial_state()), _committed_time(0.0),
      _committed_load_factor(0.0), _active_time_step(0.0),
      _active_end_time(0.0), _active_load_factor(0.0),
      _time_step_active(false) {
    validate_definition(_definition);
    _region_kernels.reserve(region_count());
    _material_histories.resize(region_count());
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        _region_kernels.emplace_back(
            IsotropicInelasticMaterial(
                _definition.spatial.regions[region_value].material,
                _definition.regions[region_value].material),
            0.0);
        _material_histories[region_value].resize(
            _spatial_model.region_element_count(region_value));
    }
    _spatial_model.set_load_factor(0.0);
}

const TransientProblemDefinition&
TransientProblem::definition() const noexcept {
    return _definition;
}

const DofMap& TransientProblem::dof_map() const noexcept {
    return _spatial_model.dof_map();
}

std::size_t TransientProblem::region_count() const noexcept {
    return _definition.spatial.regions.size();
}

const RegionDefinition& TransientProblem::region(std::size_t index) const {
    return _spatial_model.region(index);
}

const StructuredRzMesh& TransientProblem::region_mesh(std::size_t index) const {
    return _spatial_model.region_mesh(index);
}

const Quad4RzTransientKernel&
TransientProblem::region_kernel(std::size_t index) const {
    return _region_kernels.at(index);
}

const Quad4RzGeometry&
TransientProblem::region_element_geometry(std::size_t region_value,
                                          std::size_t element_index) const {
    return _spatial_model.region_element_geometry(region_value, element_index);
}

const std::vector<double>&
TransientProblem::committed_solution() const noexcept {
    return _committed_solution;
}

double TransientProblem::committed_time() const noexcept {
    return _committed_time;
}

double TransientProblem::committed_load_factor() const noexcept {
    return _committed_load_factor;
}

bool TransientProblem::time_step_active() const noexcept {
    return _time_step_active;
}

double TransientProblem::active_time_step() const {
    require_active_time_step();
    return _active_time_step;
}

double TransientProblem::active_end_time() const {
    require_active_time_step();
    return _active_end_time;
}

void TransientProblem::begin_time_step(const TransientStepInput& input) {
    if (_time_step_active)
        throw std::logic_error(
            "TransientProblem already has an active time step");
    if (!std::isfinite(input.end_time) || input.end_time <= _committed_time)
        throw std::invalid_argument(
            "TransientProblem end time must exceed committed time");
    if (!std::isfinite(input.load_factor) || input.load_factor < 0.0)
        throw std::invalid_argument(
            "TransientProblem load factor must be finite and nonnegative");
    _active_time_step = input.end_time - _committed_time;
    _active_end_time = input.end_time;
    _active_load_factor = input.load_factor;
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        _region_kernels[region_value].set_volumetric_heat_source(
            input.load_factor *
            _definition.spatial.regions[region_value].volumetric_heat_source);
    }
    _time_step_active = true;
}

void TransientProblem::commit_time_step(
    const std::vector<double>& converged_solution) {
    require_active_time_step();
    if (converged_solution.size() != dof_count())
        throw std::invalid_argument(
            "TransientProblem committed solution size mismatch");
    if (!std::all_of(converged_solution.begin(), converged_solution.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::domain_error(
            "TransientProblem committed solution must be finite");
    for (std::size_t node = 0; node < dof_map().node_count(); ++node) {
        if (!(converged_solution[dof_map().temperature(node)] > 0.0))
            throw std::domain_error(
                "TransientProblem committed temperatures must be positive");
    }

    std::vector<std::vector<Quad4MaterialHistory>> staged(region_count());
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        staged[region_value].resize(
            _spatial_model.region_element_count(region_value));
        const std::size_t offset =
            _spatial_model.region_element_offset(region_value);
        for (std::size_t element = 0; element < staged[region_value].size();
             ++element) {
            const LocalValues state =
                contribution_state(offset + element, converged_solution);
            staged[region_value][element] =
                _region_kernels[region_value].trial_state_values(
                    region_element_geometry(region_value, element), state,
                    _material_histories[region_value][element],
                    _active_time_step);
        }
    }

    _material_histories.swap(staged);
    _committed_solution = converged_solution;
    _committed_time = _active_end_time;
    _committed_load_factor = _active_load_factor;
    _active_time_step = 0.0;
    _active_end_time = _committed_time;
    _time_step_active = false;
}

void TransientProblem::rollback_time_step() noexcept {
    if (!_time_step_active)
        return;
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        _region_kernels[region_value].set_volumetric_heat_source(
            _committed_load_factor *
            _definition.spatial.regions[region_value].volumetric_heat_source);
    }
    _active_time_step = 0.0;
    _active_end_time = _committed_time;
    _active_load_factor = _committed_load_factor;
    _time_step_active = false;
}

const Quad4MaterialHistory&
TransientProblem::material_history(std::size_t region_value,
                                   std::size_t element_index) const {
    return _material_histories.at(region_value).at(element_index);
}

RegionInelasticSummary
TransientProblem::summarize_region_history(std::size_t region_value) const {
    return summarize_history(_material_histories.at(region_value));
}

InterfaceSummary
TransientProblem::summarize_interface(std::size_t contact_index,
                                      const std::vector<double>& state) const {
    return _spatial_model.summarize_interface(contact_index, state);
}

std::size_t TransientProblem::dof_count() const noexcept {
    return _spatial_model.dof_count();
}

std::size_t TransientProblem::contribution_count() const noexcept {
    return _spatial_model.contribution_count();
}

const std::vector<DirichletCondition>&
TransientProblem::dirichlet_conditions() const noexcept {
    return _spatial_model.dirichlet_conditions();
}

LocalDofs
TransientProblem::contribution_dofs(std::size_t contribution_index) const {
    return _spatial_model.contribution_dofs(contribution_index);
}

LocalResidual
TransientProblem::contribution_residual(std::size_t contribution_index,
                                        const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _spatial_model.volume_contribution_count()) {
        const auto location =
            _spatial_model.element_location(contribution_index);
        return _region_kernels[location.first].residual(
            region_element_geometry(location.first, location.second), state,
            committed_element_temperature(contribution_index),
            _material_histories[location.first][location.second],
            _active_time_step);
    }
    return _spatial_model.contribution_residual(contribution_index, state);
}

LocalSystem
TransientProblem::linearize_contribution(std::size_t contribution_index,
                                         const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _spatial_model.volume_contribution_count()) {
        const auto location =
            _spatial_model.element_location(contribution_index);
        return _region_kernels[location.first].linearize(
            region_element_geometry(location.first, location.second), state,
            committed_element_temperature(contribution_index),
            _material_histories[location.first][location.second],
            _active_time_step);
    }
    return _spatial_model.linearize_contribution(contribution_index, state);
}

void TransientProblem::add_state_independent_residual(
    std::vector<double>& residual) const {
    _spatial_model.add_external_residual(residual);
}

Quad4TemperatureHistory TransientProblem::committed_element_temperature(
    std::size_t contribution_index) const {
    const LocalValues committed =
        contribution_state(contribution_index, _committed_solution);
    Quad4TemperatureHistory temperature{};
    std::copy_n(committed.begin(), temperature.size(), temperature.begin());
    return temperature;
}

void TransientProblem::require_active_time_step() const {
    if (!_time_step_active)
        throw std::logic_error("TransientProblem residual evaluation requires "
                               "an active time step");
}

} // namespace fuelsim
