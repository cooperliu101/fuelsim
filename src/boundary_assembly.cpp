#include "spatial_assembly.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

void validate_dirichlet_conditions(
    std::vector<DirichletCondition>& conditions) {
    std::sort(conditions.begin(), conditions.end(),
              [](const DirichletCondition& lhs,
                 const DirichletCondition& rhs) {
                  return lhs.dof < rhs.dof;
              });
    for (std::size_t index = 1; index < conditions.size(); ++index) {
        const DirichletCondition& previous = conditions[index - 1];
        const DirichletCondition& current = conditions[index];
        if (previous.dof != current.dof)
            continue;
        if (previous.value != current.value)
            throw std::invalid_argument(
                "BoundaryAssembly has conflicting Dirichlet conditions");
        throw std::invalid_argument(
            "BoundaryAssembly has duplicate Dirichlet conditions");
    }
}

} // namespace

void BoundaryAssembly::build(const UnstructuredQuad4Mesh& source_mesh,
                             const SpatialLayout& layout) {
    for (const BoundaryConditionDefinition& definition :
         layout.definition().boundary_conditions) {
        if (definition.name.empty() || definition.boundary.empty() ||
            !std::isfinite(definition.value))
            throw std::invalid_argument("Boundary-condition names, boundaries, "
                                        "and values must be valid");
        const SpatialLayout::ResolvedBoundary resolved =
            layout.resolve_boundary(source_mesh, definition.boundary);
        if (definition.scale_with_load && !definition.function.empty())
            throw std::invalid_argument(
                "Boundary condition cannot combine scale_with_load and a "
                "time function: " +
                definition.name);
        if (definition.type == BoundaryConditionType::dirichlet) {
            for (const std::size_t local_node : resolved.boundary.nodes) {
                const std::size_t dof = layout.dof_map().dof(
                    definition.field,
                    layout.global_node(resolved.region, local_node));
                _dirichlet_conditions.push_back(
                    {dof, load_multiplier(definition.scale_with_load,
                                          definition.function, layout) *
                              definition.value});
                if (definition.scale_with_load || !definition.function.empty())
                    _controlled_dirichlet_conditions.push_back(
                        {dof, definition.value, definition.scale_with_load,
                         definition.function});
            }
        } else if (definition.type == BoundaryConditionType::pressure) {
            if (definition.value < 0.0)
                throw std::invalid_argument(
                    "Pressure boundary conditions must be nonnegative");
            const std::size_t load = _pressure_loads.size();
            _pressure_loads.push_back({definition.value,
                                       definition.scale_with_load,
                                       definition.function});
            const bool displaced =
                layout.region(resolved.region).strain_formulation ==
                StrainFormulation::finite;
            _pressure_kernels.emplace_back(PressureProperties{
                load_multiplier(definition.scale_with_load,
                                definition.function, layout) *
                    definition.value,
                displaced});
            const RegionMesh& mesh = layout.region_mesh(resolved.region);
            for (const Line2BoundaryElement& edge :
                 resolved.boundary.elements) {
                const auto parent = layout.edge_parent(resolved.region, edge);
                const Quad4Element& element =
                    mesh.elements().at(parent.first);
                std::array<std::size_t, 4> nodes{};
                for (std::size_t node = 0; node < nodes.size(); ++node)
                    nodes[node] = layout.global_node(
                        resolved.region, element.nodes[node]);
                _pressure_contributions.push_back(
                    {load, nodes,
                     make_line2_rz_pressure_geometry(
                         {{mesh.nodes().at(element.nodes[parent.second[0]]),
                           mesh.nodes().at(element.nodes[parent.second[1]])}},
                         parent.second)});
            }
        } else if (definition.type == BoundaryConditionType::traction) {
            if (definition.field == Field::temperature)
                throw std::invalid_argument(
                    "Traction requires a displacement field: " +
                    definition.boundary);
            if (definition.use_displaced_geometry &&
                layout.region(resolved.region).strain_formulation !=
                    StrainFormulation::finite)
                throw std::invalid_argument(
                    "Current-configuration traction requires finite strain: " +
                    definition.name);
            const std::size_t load = _traction_loads.size();
            _traction_loads.push_back({definition.value,
                                       definition.scale_with_load,
                                       definition.function});
            _traction_kernels.emplace_back(TractionProperties{
                definition.field == Field::radial_displacement
                    ? TractionComponent::radial
                    : TractionComponent::axial,
                load_multiplier(definition.scale_with_load,
                                definition.function, layout) *
                    definition.value,
                definition.use_displaced_geometry});
            const RegionMesh& mesh = layout.region_mesh(resolved.region);
            for (const Line2BoundaryElement& edge :
                 resolved.boundary.elements) {
                const auto parent = layout.edge_parent(resolved.region, edge);
                const Quad4Element& element =
                    mesh.elements().at(parent.first);
                std::array<std::size_t, 4> nodes{};
                for (std::size_t node = 0; node < nodes.size(); ++node)
                    nodes[node] = layout.global_node(
                        resolved.region, element.nodes[node]);
                _traction_contributions.push_back(
                    {load, nodes,
                     make_line2_rz_traction_geometry(
                         {{mesh.nodes().at(element.nodes[parent.second[0]]),
                           mesh.nodes().at(element.nodes[parent.second[1]])}},
                         parent.second)});
            }
        } else {
            if (!(definition.heat_transfer_coefficient > 0.0) ||
                !(definition.ambient_temperature > 0.0))
                throw std::invalid_argument(
                    "Convection coefficient and ambient temperature must be "
                    "positive: " +
                    definition.name);
            const std::size_t load = _convection_loads.size();
            _convection_loads.push_back(
                {definition.heat_transfer_coefficient,
                 definition.ambient_temperature,
                 definition.coefficient_function,
                 definition.ambient_temperature_function});
            _convection_kernels.emplace_back(
                ConvectionProperties{definition.heat_transfer_coefficient,
                                     definition.ambient_temperature});
            const RegionMesh& mesh = layout.region_mesh(resolved.region);
            for (const Line2BoundaryElement& edge :
                 resolved.boundary.elements) {
                const auto parent = layout.edge_parent(resolved.region, edge);
                const Quad4Element& element =
                    mesh.elements().at(parent.first);
                std::array<std::size_t, 4> nodes{};
                for (std::size_t node = 0; node < nodes.size(); ++node)
                    nodes[node] = layout.global_node(
                        resolved.region, element.nodes[node]);
                _convection_contributions.push_back(
                    {load, nodes,
                     make_line2_rz_convection_geometry(
                         {{mesh.nodes().at(element.nodes[parent.second[0]]),
                           mesh.nodes().at(element.nodes[parent.second[1]])}},
                         parent.second)});
            }
        }
    }
    validate_dirichlet_conditions(_dirichlet_conditions);
    refresh_controlled_values(layout);
}

void BoundaryAssembly::set_load_factor(double value,
                                       const SpatialLayout& layout) {
    if (!std::isfinite(value) || value < 0.0)
        throw std::invalid_argument(
            "SpatialAssembly load factor must be finite and nonnegative");
    _load_factor = value;
    refresh_controlled_values(layout);
}

double BoundaryAssembly::load_factor() const noexcept {
    return _load_factor;
}

void BoundaryAssembly::set_time(double value, const SpatialLayout& layout) {
    if (!std::isfinite(value) || value < 0.0)
        throw std::invalid_argument(
            "SpatialAssembly time must be finite and nonnegative");
    _time = value;
    refresh_controlled_values(layout);
}

double BoundaryAssembly::function_value(const std::string& name,
                                        const SpatialLayout& layout) const {
    const auto found = std::find_if(
        layout.definition().time_tables.begin(),
        layout.definition().time_tables.end(),
        [&name](const PiecewiseLinearTimeTable& table) {
            return table.name() == name;
        });
    if (found == layout.definition().time_tables.end())
        throw std::invalid_argument("Unknown time-table function: " + name);
    return found->value(_time);
}

double BoundaryAssembly::load_multiplier(bool scale_with_load,
                                         const std::string& function,
                                         const SpatialLayout& layout) const {
    if (!function.empty())
        return function_value(function, layout);
    return scale_with_load ? _load_factor : 1.0;
}

void BoundaryAssembly::refresh_controlled_values(const SpatialLayout& layout) {
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
            load_multiplier(controlled.scale_with_load, controlled.function,
                            layout) *
            controlled.value;
    }
    for (std::size_t load = 0; load < _convection_loads.size(); ++load) {
        const ConvectionLoad& convection = _convection_loads[load];
        const double coefficient_multiplier =
            convection.coefficient_function.empty()
                ? 1.0
                : function_value(convection.coefficient_function, layout);
        const double ambient_multiplier =
            convection.ambient_temperature_function.empty()
                ? 1.0
                : function_value(convection.ambient_temperature_function,
                                 layout);
        _convection_kernels[load].set_properties(
            {coefficient_multiplier * convection.heat_transfer_coefficient,
             ambient_multiplier * convection.ambient_temperature});
    }
    for (std::size_t load = 0; load < _pressure_loads.size(); ++load) {
        const PressureLoad& pressure = _pressure_loads[load];
        const double value =
            load_multiplier(pressure.scale_with_load, pressure.function,
                            layout) *
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
            load_multiplier(traction.scale_with_load, traction.function,
                            layout) *
            traction.traction;
        _traction_kernels[load].set_properties(properties);
    }
}

double BoundaryAssembly::region_heat_source(
    std::size_t region_index, const SpatialLayout& layout) const {
    const RegionDefinition& region =
        layout.definition().regions.at(region_index);
    const double multiplier =
        region.heat_source_function.empty()
            ? _load_factor
            : function_value(region.heat_source_function, layout);
    return multiplier * region.volumetric_heat_source;
}

std::size_t BoundaryAssembly::pressure_contribution_count() const noexcept {
    return _pressure_contributions.size();
}

std::size_t BoundaryAssembly::traction_contribution_count() const noexcept {
    return _traction_contributions.size();
}

std::size_t BoundaryAssembly::convection_contribution_count() const noexcept {
    return _convection_contributions.size();
}

LocalDofs BoundaryAssembly::contribution_dofs(
    SpatialContributionType type, std::size_t index,
    const SpatialLayout& layout) const {
    switch (type) {
    case SpatialContributionType::pressure:
        return layout.dof_map().local_dofs(
            _pressure_contributions.at(index).nodes);
    case SpatialContributionType::traction:
        return layout.dof_map().local_dofs(
            _traction_contributions.at(index).nodes);
    case SpatialContributionType::convection:
        return layout.dof_map().local_dofs(
            _convection_contributions.at(index).nodes);
    default:
        throw std::invalid_argument(
            "BoundaryAssembly requires a boundary contribution type");
    }
}

LocalResidual BoundaryAssembly::contribution_residual(
    SpatialContributionType type, std::size_t index,
    const LocalValues& state) const {
    switch (type) {
    case SpatialContributionType::pressure: {
        const PressureContribution& contribution =
            _pressure_contributions.at(index);
        return _pressure_kernels[contribution.load].residual(
            contribution.geometry, state);
    }
    case SpatialContributionType::traction: {
        const TractionContribution& contribution =
            _traction_contributions.at(index);
        return _traction_kernels[contribution.load].residual(
            contribution.geometry, state);
    }
    case SpatialContributionType::convection: {
        const ConvectionContribution& contribution =
            _convection_contributions.at(index);
        return _convection_kernels[contribution.load].residual(
            contribution.geometry, state);
    }
    default:
        throw std::invalid_argument(
            "BoundaryAssembly requires a boundary contribution type");
    }
}

LocalSystem BoundaryAssembly::linearize_contribution(
    SpatialContributionType type, std::size_t index,
    const LocalValues& state) const {
    switch (type) {
    case SpatialContributionType::pressure: {
        const PressureContribution& contribution =
            _pressure_contributions.at(index);
        return _pressure_kernels[contribution.load].linearize(
            contribution.geometry, state);
    }
    case SpatialContributionType::traction: {
        const TractionContribution& contribution =
            _traction_contributions.at(index);
        return _traction_kernels[contribution.load].linearize(
            contribution.geometry, state);
    }
    case SpatialContributionType::convection: {
        const ConvectionContribution& contribution =
            _convection_contributions.at(index);
        return _convection_kernels[contribution.load].linearize(
            contribution.geometry, state);
    }
    default:
        throw std::invalid_argument(
            "BoundaryAssembly requires a boundary contribution type");
    }
}

const std::vector<DirichletCondition>&
BoundaryAssembly::dirichlet_conditions() const noexcept {
    return _dirichlet_conditions;
}

void SpatialAssembly::set_load_factor(double value) {
    _boundary.set_load_factor(value, _layout);
}

double SpatialAssembly::load_factor() const noexcept {
    return _boundary.load_factor();
}

void SpatialAssembly::set_time(double value) {
    _boundary.set_time(value, _layout);
}


} // namespace fuelsim
