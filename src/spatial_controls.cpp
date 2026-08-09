#include "spatial_assembly.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim {

void SpatialAssembly::set_load_factor(double value) {
    if (!std::isfinite(value) || value < 0.0)
        throw std::invalid_argument(
            "SpatialAssembly load factor must be finite and nonnegative");
    _load_factor = value;
    refresh_controlled_values();
}

double SpatialAssembly::load_factor() const noexcept {
    return _load_factor;
}

void SpatialAssembly::set_time(double value) {
    if (!std::isfinite(value) || value < 0.0)
        throw std::invalid_argument(
            "SpatialAssembly time must be finite and nonnegative");
    _time = value;
    refresh_controlled_values();
}

double SpatialAssembly::function_value(const std::string& name) const {
    const auto found = std::find_if(
        _definition.time_tables.begin(), _definition.time_tables.end(),
        [&name](const PiecewiseLinearTimeTable& table) {
            return table.name() == name;
        });
    if (found == _definition.time_tables.end())
        throw std::invalid_argument("Unknown time-table function: " + name);
    return found->value(_time);
}

double SpatialAssembly::load_multiplier(bool scale_with_load,
                                      const std::string& function) const {
    if (!function.empty())
        return function_value(function);
    return scale_with_load ? _load_factor : 1.0;
}

void SpatialAssembly::refresh_controlled_values() {
    for (const ControlledDirichlet& controlled :
         _controlled_dirichlet_conditions) {
        const auto condition = std::lower_bound(
            _dirichlet_conditions.begin(), _dirichlet_conditions.end(),
            controlled.dof,
            [](const DirichletCondition& candidate, std::size_t dof) {
                return candidate.dof < dof;
            });
        if (condition == _dirichlet_conditions.end() ||
            condition->dof != controlled.dof)
            throw std::logic_error(
                "SpatialAssembly controlled Dirichlet mapping is invalid");
        condition->value =
            load_multiplier(controlled.scale_with_load, controlled.function) *
            controlled.value;
    }
    for (std::size_t load = 0; load < _convection_loads.size(); ++load) {
        const ConvectionLoad& convection = _convection_loads[load];
        const double coefficient_multiplier =
            convection.coefficient_function.empty()
                ? 1.0
                : function_value(convection.coefficient_function);
        const double ambient_multiplier =
            convection.ambient_temperature_function.empty()
                ? 1.0
                : function_value(convection.ambient_temperature_function);
        _convection_kernels[load].set_properties(
            {coefficient_multiplier * convection.heat_transfer_coefficient,
             ambient_multiplier * convection.ambient_temperature});
    }
    for (std::size_t load = 0; load < _pressure_loads.size(); ++load) {
        const PressureLoad& pressure = _pressure_loads[load];
        const double value =
            load_multiplier(pressure.scale_with_load, pressure.function) *
            pressure.pressure;
        if (value < 0.0)
            throw std::domain_error(
                "Pressure time function produced a negative load");
        PressureProperties properties = _pressure_kernels[load].properties();
        properties.pressure = value;
        _pressure_kernels[load].set_properties(properties);
    }
    for (std::size_t load = 0; load < _traction_loads.size(); ++load) {
        const TractionLoad& traction = _traction_loads[load];
        TractionProperties properties = _traction_kernels[load].properties();
        properties.traction =
            load_multiplier(traction.scale_with_load, traction.function) *
            traction.traction;
        _traction_kernels[load].set_properties(properties);
    }
}


} // namespace fuelsim
