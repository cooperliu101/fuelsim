#include "fuelsim/transient_fuel_cladding_problem.hpp"

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

TransientFuelCladdingProblem::TransientFuelCladdingProblem(
    TransientFuelCladdingParameters parameters)
    : TransientFuelCladdingProblem(
          parameters,
          StructuredRzMesh::make_annulus(0.0, parameters.steady.fuel_radius,
                                         parameters.steady.fuel_length,
                                         parameters.steady.fuel_radial_elements,
                                         parameters.steady.axial_elements),
          StructuredRzMesh::make_annulus(
              parameters.steady.cladding_inner_radius,
              parameters.steady.cladding_outer_radius,
              parameters.steady.cladding_length,
              parameters.steady.cladding_radial_elements,
              parameters.steady.axial_elements)) {}

TransientFuelCladdingProblem::TransientFuelCladdingProblem(
    TransientFuelCladdingParameters parameters, StructuredRzMesh fuel_mesh,
    StructuredRzMesh cladding_mesh)
    : _parameters(parameters),
      _steady_problem(parameters.steady, std::move(fuel_mesh),
                      std::move(cladding_mesh)),
      _fuel_kernel(
          IsotropicInelasticMaterial(parameters.steady.fuel, parameters.fuel),
          0.0),
      _cladding_kernel(IsotropicInelasticMaterial(parameters.steady.cladding,
                                                  parameters.cladding),
                       0.0),
      _fuel_material_history(_steady_problem.fuel_element_count()),
      _cladding_material_history(_steady_problem.cladding_element_count()),
      _committed_solution(_steady_problem.initial_state()),
      _committed_time(0.0), _committed_heat_source(0.0), _active_time_step(0.0),
      _active_end_time(0.0), _active_heat_source(0.0),
      _time_step_active(false) {}

const TransientFuelCladdingParameters&
TransientFuelCladdingProblem::parameters() const noexcept {
    return _parameters;
}

const SteadyFuelCladdingProblem&
TransientFuelCladdingProblem::steady_problem() const noexcept {
    return _steady_problem;
}

const Quad4RzTransientKernel&
TransientFuelCladdingProblem::fuel_kernel() const noexcept {
    return _fuel_kernel;
}

const Quad4RzTransientKernel&
TransientFuelCladdingProblem::cladding_kernel() const noexcept {
    return _cladding_kernel;
}

const std::vector<double>&
TransientFuelCladdingProblem::committed_solution() const noexcept {
    return _committed_solution;
}

double TransientFuelCladdingProblem::committed_time() const noexcept {
    return _committed_time;
}

double TransientFuelCladdingProblem::committed_heat_source() const noexcept {
    return _committed_heat_source;
}

bool TransientFuelCladdingProblem::time_step_active() const noexcept {
    return _time_step_active;
}

double TransientFuelCladdingProblem::active_time_step() const {
    require_active_time_step();
    return _active_time_step;
}

double TransientFuelCladdingProblem::active_end_time() const {
    require_active_time_step();
    return _active_end_time;
}

void TransientFuelCladdingProblem::begin_time_step(
    const TransientStepInput& input) {
    if (_time_step_active)
        throw std::logic_error(
            "TransientFuelCladdingProblem already has an active time step");
    if (!std::isfinite(input.end_time) || !(input.end_time > _committed_time))
        throw std::invalid_argument("TransientFuelCladdingProblem time-step "
                                    "end_time must be finite and greater than "
                                    "committed_time");
    if (!std::isfinite(input.volumetric_heat_source) ||
        !(input.volumetric_heat_source >= 0.0))
        throw std::invalid_argument("TransientFuelCladdingProblem volumetric "
                                    "heat source must be finite and "
                                    "nonnegative");

    _active_time_step = input.end_time - _committed_time;
    _active_end_time = input.end_time;
    _active_heat_source = input.volumetric_heat_source;
    _fuel_kernel.set_volumetric_heat_source(_active_heat_source);
    _time_step_active = true;
}

void TransientFuelCladdingProblem::commit_time_step(
    const std::vector<double>& converged_solution) {
    require_active_time_step();
    if (converged_solution.size() != dof_count())
        throw std::invalid_argument(
            "TransientFuelCladdingProblem committed solution size mismatch");
    if (!std::all_of(converged_solution.begin(), converged_solution.end(),
                     [](double value) { return std::isfinite(value); }))
        throw std::domain_error("TransientFuelCladdingProblem committed "
                                "solution must contain finite values");
    for (std::size_t node = 0; node < _steady_problem.dof_map().node_count();
         ++node) {
        const double temperature =
            converged_solution[_steady_problem.dof_map().temperature(node)];
        if (!(temperature > 0.0))
            throw std::domain_error("TransientFuelCladdingProblem committed "
                                    "nodal temperatures must be positive");
    }

    std::vector<Quad4MaterialHistory> staged_fuel(
        _fuel_material_history.size());
    for (std::size_t element = 0; element < staged_fuel.size(); ++element) {
        const LocalValues state =
            contribution_state(element, converged_solution);
        staged_fuel[element] = _fuel_kernel.trial_state_values(
            _steady_problem.fuel_element_geometry(element), state,
            _fuel_material_history[element], _active_time_step);
    }

    std::vector<Quad4MaterialHistory> staged_cladding(
        _cladding_material_history.size());
    const std::size_t cladding_offset = _steady_problem.fuel_element_count();
    for (std::size_t element = 0; element < staged_cladding.size(); ++element) {
        const LocalValues state =
            contribution_state(cladding_offset + element, converged_solution);
        staged_cladding[element] = _cladding_kernel.trial_state_values(
            _steady_problem.cladding_element_geometry(element), state,
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

void TransientFuelCladdingProblem::rollback_time_step() noexcept {
    if (!_time_step_active)
        return;
    _fuel_kernel.set_volumetric_heat_source(_committed_heat_source);
    _active_time_step = 0.0;
    _active_end_time = _committed_time;
    _active_heat_source = _committed_heat_source;
    _time_step_active = false;
}

const Quad4MaterialHistory& TransientFuelCladdingProblem::fuel_material_history(
    std::size_t element_index) const {
    return _fuel_material_history.at(element_index);
}

const Quad4MaterialHistory&
TransientFuelCladdingProblem::cladding_material_history(
    std::size_t element_index) const {
    return _cladding_material_history.at(element_index);
}

RegionInelasticSummary
TransientFuelCladdingProblem::summarize_fuel_history() const noexcept {
    return summarize_history(_fuel_material_history);
}

RegionInelasticSummary
TransientFuelCladdingProblem::summarize_cladding_history() const noexcept {
    return summarize_history(_cladding_material_history);
}

InterfaceSummary TransientFuelCladdingProblem::summarize_interface(
    const std::vector<double>& state) const {
    return _steady_problem.summarize_interface(state);
}

std::size_t TransientFuelCladdingProblem::dof_count() const noexcept {
    return _steady_problem.dof_count();
}

std::size_t TransientFuelCladdingProblem::contribution_count() const noexcept {
    return _steady_problem.contribution_count();
}

const std::vector<DirichletCondition>&
TransientFuelCladdingProblem::dirichlet_conditions() const noexcept {
    return _steady_problem.dirichlet_conditions();
}

LocalDofs TransientFuelCladdingProblem::contribution_dofs(
    std::size_t contribution_index) const {
    return _steady_problem.contribution_dofs(contribution_index);
}

LocalResidual TransientFuelCladdingProblem::contribution_residual(
    std::size_t contribution_index, const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _steady_problem.fuel_element_count()) {
        return _fuel_kernel.residual(
            _steady_problem.fuel_element_geometry(contribution_index), state,
            committed_element_temperature(contribution_index),
            _fuel_material_history.at(contribution_index), _active_time_step);
    }

    const std::size_t cladding_index =
        contribution_index - _steady_problem.fuel_element_count();
    if (cladding_index < _steady_problem.cladding_element_count()) {
        return _cladding_kernel.residual(
            _steady_problem.cladding_element_geometry(cladding_index), state,
            committed_element_temperature(contribution_index),
            _cladding_material_history.at(cladding_index), _active_time_step);
    }

    return _steady_problem.contribution_residual(contribution_index, state);
}

LocalSystem TransientFuelCladdingProblem::linearize_contribution(
    std::size_t contribution_index, const LocalValues& state) const {
    require_active_time_step();
    if (contribution_index < _steady_problem.fuel_element_count()) {
        return _fuel_kernel.linearize(
            _steady_problem.fuel_element_geometry(contribution_index), state,
            committed_element_temperature(contribution_index),
            _fuel_material_history.at(contribution_index), _active_time_step);
    }

    const std::size_t cladding_index =
        contribution_index - _steady_problem.fuel_element_count();
    if (cladding_index < _steady_problem.cladding_element_count()) {
        return _cladding_kernel.linearize(
            _steady_problem.cladding_element_geometry(cladding_index), state,
            committed_element_temperature(contribution_index),
            _cladding_material_history.at(cladding_index), _active_time_step);
    }

    return _steady_problem.linearize_contribution(contribution_index, state);
}

void TransientFuelCladdingProblem::add_state_independent_residual(
    std::vector<double>& residual) const {
    (void)residual;
}

Quad4TemperatureHistory
TransientFuelCladdingProblem::committed_element_temperature(
    std::size_t contribution_index) const {
    const LocalValues committed =
        contribution_state(contribution_index, _committed_solution);
    Quad4TemperatureHistory temperature{};
    std::copy_n(committed.begin(), temperature.size(), temperature.begin());
    return temperature;
}

void TransientFuelCladdingProblem::require_active_time_step() const {
    if (!_time_step_active)
        throw std::logic_error("TransientFuelCladdingProblem residual "
                               "evaluation requires an active time step");
}

} // namespace fuelsim
