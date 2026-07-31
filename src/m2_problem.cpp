#include "fuelsim/m2_problem.hpp"

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

} // namespace

M2Problem::M2Problem(M2Parameters parameters)
    : M2Problem(
          parameters,
          StructuredRzMesh::make_annulus(0.0, parameters.base.fuel_radius,
                                         parameters.base.fuel_length,
                                         parameters.base.fuel_radial_elements,
                                         parameters.base.axial_elements),
          StructuredRzMesh::make_annulus(
              parameters.base.cladding_inner_radius,
              parameters.base.cladding_outer_radius,
              parameters.base.cladding_length,
              parameters.base.cladding_radial_elements,
              parameters.base.axial_elements)) {}

M2Problem::M2Problem(M2Parameters parameters, StructuredRzMesh fuel_mesh,
                     StructuredRzMesh cladding_mesh)
    : _parameters(parameters),
      _base_problem(parameters.base, std::move(fuel_mesh),
                    std::move(cladding_mesh)),
      _fuel_kernel(
          IsotropicInelasticMaterial(parameters.base.fuel, parameters.fuel),
          0.0),
      _cladding_kernel(IsotropicInelasticMaterial(parameters.base.cladding,
                                                  parameters.cladding),
                       0.0),
      _fuel_material_history(_base_problem.fuel_element_count()),
      _cladding_material_history(_base_problem.cladding_element_count()),
      _committed_solution(_base_problem.initial_state()), _committed_time(0.0),
      _committed_heat_source(0.0), _active_time_step(0.0),
      _active_end_time(0.0), _active_heat_source(0.0),
      _time_step_active(false) {}

const M2Parameters& M2Problem::parameters() const noexcept {
    return _parameters;
}

const M1Problem& M2Problem::base_problem() const noexcept {
    return _base_problem;
}

const Quad4RzTransientKernel& M2Problem::fuel_kernel() const noexcept {
    return _fuel_kernel;
}

const Quad4RzTransientKernel& M2Problem::cladding_kernel() const noexcept {
    return _cladding_kernel;
}

const std::vector<double>& M2Problem::committed_solution() const noexcept {
    return _committed_solution;
}

double M2Problem::committed_time() const noexcept {
    return _committed_time;
}

double M2Problem::committed_heat_source() const noexcept {
    return _committed_heat_source;
}

bool M2Problem::time_step_active() const noexcept {
    return _time_step_active;
}

double M2Problem::active_time_step() const {
    require_active_time_step();
    return _active_time_step;
}

double M2Problem::active_end_time() const {
    require_active_time_step();
    return _active_end_time;
}

void M2Problem::begin_time_step(const M2TimeStepInput& input) {
    if (_time_step_active)
        throw std::logic_error("M2Problem already has an active time step");
    if (!std::isfinite(input.end_time) || !(input.end_time > _committed_time))
        throw std::invalid_argument(
            "M2Problem time-step end_time must be finite and greater than "
            "committed_time");
    if (!std::isfinite(input.volumetric_heat_source) ||
        !(input.volumetric_heat_source >= 0.0))
        throw std::invalid_argument(
            "M2Problem volumetric heat source must be finite and "
            "nonnegative");

    _active_time_step = input.end_time - _committed_time;
    _active_end_time = input.end_time;
    _active_heat_source = input.volumetric_heat_source;
    _fuel_kernel.set_volumetric_heat_source(_active_heat_source);
    _time_step_active = true;
}

void M2Problem::commit_time_step(
    const std::vector<double>& converged_solution) {
    require_active_time_step();
    if (converged_solution.size() != dof_count())
        throw std::invalid_argument(
            "M2Problem committed solution size mismatch");
    if (!std::all_of(converged_solution.begin(), converged_solution.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::domain_error(
            "M2Problem committed solution must contain finite values");
    for (std::size_t node = 0; node < _base_problem.dof_map().node_count();
         ++node) {
        const double temperature =
            converged_solution[_base_problem.dof_map().temperature(node)];
        if (!(temperature > 0.0))
            throw std::domain_error(
                "M2Problem committed nodal temperatures must be positive");
    }

    std::vector<Quad4MaterialHistory> staged_fuel(
        _fuel_material_history.size());
    for (std::size_t element = 0; element < staged_fuel.size(); ++element) {
        const LocalValues state =
            contribution_state(element, converged_solution);
        staged_fuel[element] = _fuel_kernel.trial_state_values(
            _base_problem.fuel_element_geometry(element), state,
            _fuel_material_history[element], _active_time_step);
    }

    std::vector<Quad4MaterialHistory> staged_cladding(
        _cladding_material_history.size());
    const std::size_t cladding_offset = _base_problem.fuel_element_count();
    for (std::size_t element = 0; element < staged_cladding.size(); ++element) {
        const LocalValues state =
            contribution_state(cladding_offset + element, converged_solution);
        staged_cladding[element] = _cladding_kernel.trial_state_values(
            _base_problem.cladding_element_geometry(element), state,
            _cladding_material_history[element], _active_time_step);
    }

    std::vector<double> staged_solution = converged_solution;
    _fuel_material_history.swap(staged_fuel);
    _cladding_material_history.swap(staged_cladding);
    _committed_solution.swap(staged_solution);
    _committed_time = _active_end_time;
    _committed_heat_source = _active_heat_source;
    _active_time_step = 0.0;
    _active_end_time = _committed_time;
    _time_step_active = false;
}

void M2Problem::rollback_time_step() noexcept {
    if (!_time_step_active)
        return;
    _fuel_kernel.set_volumetric_heat_source(_committed_heat_source);
    _active_time_step = 0.0;
    _active_end_time = _committed_time;
    _active_heat_source = _committed_heat_source;
    _time_step_active = false;
}

const Quad4MaterialHistory&
M2Problem::fuel_material_history(std::size_t element_index) const {
    return _fuel_material_history.at(element_index);
}

const Quad4MaterialHistory&
M2Problem::cladding_material_history(std::size_t element_index) const {
    return _cladding_material_history.at(element_index);
}

RegionInelasticSummary M2Problem::summarize_fuel_history() const noexcept {
    return summarize_history(_fuel_material_history);
}

RegionInelasticSummary M2Problem::summarize_cladding_history() const noexcept {
    return summarize_history(_cladding_material_history);
}

InterfaceSummary
M2Problem::summarize_interface(const std::vector<double>& state) const {
    return _base_problem.summarize_interface(state);
}

std::size_t M2Problem::dof_count() const noexcept {
    return _base_problem.dof_count();
}

std::size_t M2Problem::contribution_count() const noexcept {
    return _base_problem.contribution_count();
}

const std::vector<DirichletCondition>&
M2Problem::dirichlet_conditions() const noexcept {
    return _base_problem.dirichlet_conditions();
}

LocalDofs M2Problem::contribution_dofs(std::size_t contribution_index) const {
    return _base_problem.contribution_dofs(contribution_index);
}

LocalResidual M2Problem::contribution_residual(std::size_t contribution_index,
                                               const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _base_problem.fuel_element_count()) {
        return _fuel_kernel.residual(
            _base_problem.fuel_element_geometry(contribution_index), state,
            committed_element_temperature(contribution_index),
            _fuel_material_history.at(contribution_index), _active_time_step);
    }

    const std::size_t cladding_index =
        contribution_index - _base_problem.fuel_element_count();
    if (cladding_index < _base_problem.cladding_element_count()) {
        return _cladding_kernel.residual(
            _base_problem.cladding_element_geometry(cladding_index), state,
            committed_element_temperature(contribution_index),
            _cladding_material_history.at(cladding_index), _active_time_step);
    }

    return _base_problem.contribution_residual(contribution_index, state);
}

LocalSystem M2Problem::linearize_contribution(std::size_t contribution_index,
                                              const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _base_problem.fuel_element_count()) {
        return _fuel_kernel.linearize(
            _base_problem.fuel_element_geometry(contribution_index), state,
            committed_element_temperature(contribution_index),
            _fuel_material_history.at(contribution_index), _active_time_step);
    }

    const std::size_t cladding_index =
        contribution_index - _base_problem.fuel_element_count();
    if (cladding_index < _base_problem.cladding_element_count()) {
        return _cladding_kernel.linearize(
            _base_problem.cladding_element_geometry(cladding_index), state,
            committed_element_temperature(contribution_index),
            _cladding_material_history.at(cladding_index), _active_time_step);
    }

    return _base_problem.linearize_contribution(contribution_index, state);
}

void M2Problem::add_state_independent_residual(
    std::vector<double>& residual) const {
    (void)residual;
}

Quad4TemperatureHistory
M2Problem::committed_element_temperature(std::size_t contribution_index) const {
    const LocalValues committed =
        contribution_state(contribution_index, _committed_solution);
    Quad4TemperatureHistory temperature{};
    std::copy_n(committed.begin(), temperature.size(), temperature.begin());
    return temperature;
}

void M2Problem::require_active_time_step() const {
    if (!_time_step_active)
        throw std::logic_error(
            "M2Problem residual evaluation requires an active time step");
}

} // namespace fuelsim
