#include "plane_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

namespace fuelsim::plane {
SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredPlaneQuad8Mesh& mesh)
    : SpatialLayout(definition,
          spatial_detail::resolve_block_ids(definition, mesh),
          spatial_detail::DofLayout::cartesian_3d) {
    std::vector<std::size_t> counts, elements;
    std::map<std::size_t, std::size_t> connected_controls;
    const auto& sections = _definition.generalized_plane_strain;
    if (sections.empty())
        throw std::invalid_argument("CPEG8T requires generalized plane strain sections");
    _region_sections.assign(region_count(), sections.size());
    _section_origins.resize(sections.size());
    std::vector<double> areas(sections.size(), 0.0);
    std::set<std::string> section_names;
    for (std::size_t c = 0; c < sections.size(); ++c) {
        const auto& section = sections[c];
        if (section.name.empty() || section.name.size() > 38 || !section_names.insert(section.name).second
            || section.blocks.empty() || !std::isfinite(section.initial_thickness) || !(section.initial_thickness > 0))
            throw std::invalid_argument(
                "CPEG8T section requires a unique name of 1-38 characters, blocks and positive thickness");
        for (std::size_t i = 0; i < 3; ++i) {
            if ((!section.functions[i].empty() && !section.prescribed[i])
                || (section.prescribed[i] && !std::isfinite(*section.prescribed[i])))
                throw std::invalid_argument("Invalid CPEG8T section constraint");
        }
        for (const auto& block : section.blocks) {
            const auto id = mesh.element_block(block).id;
            const auto found = std::find(_block_ids.begin(), _block_ids.end(), id);
            if (found == _block_ids.end())
                throw std::invalid_argument("CPEG8T section block is outside selected regions: " + block);
            const auto r = static_cast<std::size_t>(found - _block_ids.begin());
            if (_region_sections[r] != sections.size())
                throw std::invalid_argument("CPEG8T block belongs to multiple section entries: " + block);
            _region_sections[r] = c;
        }
    }
    if (std::find(_region_sections.begin(), _region_sections.end(), sections.size()) != _region_sections.end())
        throw std::invalid_argument("Every CPEG8T region must belong to a section");
    for (std::size_t r = 0; r < region_count(); ++r) {
        _source_elements.emplace_back();
        _geometry.emplace_back();
        std::set<std::size_t> used, thermal;
        for (std::size_t e = 0; e < mesh.elements().size(); ++e) {
            if (mesh.element_block_ids()[e] != _block_ids[r])
                continue;
            const auto& element = mesh.elements()[e];
            _source_elements.back().push_back(e);
            const auto control = _region_sections[r];
            elements::Cpeg8Coordinates coordinates;
            for (std::size_t n = 0; n < 8; ++n) {
                const auto node = element.nodes[n];
                used.insert(node);
                if (n < 4)
                    thermal.insert(node);
                coordinates[n] = mesh.nodes()[node];
                const auto inserted = connected_controls.emplace(node, control);
                if (!inserted.second && inserted.first->second != control)
                    throw std::invalid_argument("Connected CPEG8T elements must share their section");
            }
            const double thickness = sections[control].initial_thickness;
            _geometry.back().push_back(elements::make_cpeg8t_geometry(coordinates, thickness));
            for (const auto& point : _geometry.back().back().points) {
                const double area = point.measure / thickness;
                areas[control] += area;
                _section_origins[control][0] += area * point.x;
                _section_origins[control][1] += area * point.y;
            }
        }
        if (used.empty())
            throw std::invalid_argument("CPEG8T region contains no elements");
        _source_nodes.emplace_back(used.begin(), used.end());
        _thermal_nodes.emplace_back();
        for (auto node : _source_nodes.back())
            _thermal_nodes.back().push_back(thermal.count(node) != 0);
        counts.push_back(used.size());
        elements.push_back(_geometry.back().size());
    }
    for (std::size_t c = 0; c < sections.size(); ++c)
        for (auto& coordinate : _section_origins[c])
            coordinate /= areas[c];
    for (std::size_t r = 0; r < region_count(); ++r)
        for (auto& geometry : _geometry[r])
            geometry.reference_point = _section_origins[_region_sections[r]];
    initialize_counts(counts, elements);
    initialize_mixed_shared_nodes(_source_nodes, _thermal_nodes);
    _field_layout.resize(3);
    const std::array<const char*, 3> names{"u3", "rotation_x", "rotation_y"};
    for (std::size_t i = 0; i < 3; ++i) {
        const auto start = _field_layout.back().end;
        _field_layout.push_back({names[i], start, start + sections.size(), FieldCategory::mechanical});
    }
    for (std::size_t r = 0; r < region_count(); ++r)
        for (auto source_element : _source_elements[r]) {
            const auto& element = mesh.elements()[source_element];
            std::array<std::size_t, 23> dofs{};
            for (std::size_t n = 0; n < 8; ++n) {
                const auto it = std::lower_bound(_source_nodes[r].begin(), _source_nodes[r].end(), element.nodes[n]);
                const auto local = static_cast<std::size_t>(it - _source_nodes[r].begin());
                if (n < 4)
                    dofs[n] = dof(Field::temperature, global_temperature_node(r, local));
                dofs[4 + n] = dof(Field::displacement_x, global_node(r, local));
                dofs[12 + n] = dof(Field::displacement_y, global_node(r, local));
            }
            for (std::size_t i = 0; i < 3; ++i)
                dofs[20 + i] = section_dof(i, _region_sections[r]);
            _dofs.push_back(dofs);
        }
    for (std::size_t b = 0; b < _definition.boundary_conditions.size(); ++b) {
        const auto& boundary = _definition.boundary_conditions[b];
        if (boundary.type != BoundaryConditionType::dirichlet) {
            if (boundary.type == BoundaryConditionType::axial_force)
                throw std::invalid_argument("CPEG8T reference forces require a separate control load definition");
            const auto& set = mesh.side_set(boundary.boundary);
            for (const auto& side : set.sides) {
                bool found = false;
                for (std::size_t r = 0; r < region_count(); ++r) {
                    const auto it = std::find(_source_elements[r].begin(), _source_elements[r].end(), side.element);
                    if (it == _source_elements[r].end())
                        continue;
                    _boundaries.push_back(
                        {b, r, static_cast<std::size_t>(it - _source_elements[r].begin()), side.local_side});
                    record_configuration_warning(boundary, region(r));
                    found = true;
                }
                if (!found)
                    throw std::invalid_argument("CPEG8T boundary lies outside selected regions");
            }
            continue;
        }
        const auto set = std::find_if(mesh.node_sets().begin(), mesh.node_sets().end(), [&](const NodeSet& value) {
            return value.name == boundary.boundary;
        });
        if (set == mesh.node_sets().end())
            throw std::invalid_argument("Unknown CPEG8T node set: " + boundary.boundary);
        std::set<std::size_t> selected;
        for (auto source_node : set->nodes) {
            if (boundary.field != Field::temperature && boundary.field != Field::displacement_x
                && boundary.field != Field::displacement_y)
                throw std::invalid_argument("CPEG8T field must be temperature, displacement_x or displacement_y");
            for (std::size_t r = 0; r < region_count(); ++r) {
                const auto it = std::lower_bound(_source_nodes[r].begin(), _source_nodes[r].end(), source_node);
                if (it == _source_nodes[r].end() || *it != source_node)
                    continue;
                const auto local = static_cast<std::size_t>(it - _source_nodes[r].begin());
                if (boundary.field == Field::temperature && !_thermal_nodes[r][local])
                    continue;
                selected.insert(dof(boundary.field,
                    boundary.field == Field::temperature ? global_temperature_node(r, local) : global_node(r, local)));
            }
        }
        if (selected.empty())
            throw std::invalid_argument("CPEG8T boundary contains no active field nodes");
        for (auto index : selected)
            add_dirichlet(index, b);
    }
    for (std::size_t c = 0; c < sections.size(); ++c)
        for (std::size_t i = 0; i < 3; ++i)
            if (sections[c].prescribed[i])
                _dirichlet_conditions.push_back({section_dof(i, c), 0.0});
    spatial_detail::validate_dirichlet_conditions(_dirichlet_conditions,
        "Conflicting CPEG8T conditions",
        "Duplicate CPEG8T conditions");
    refresh_constraints();
    build_contacts(mesh);
    validate_state(initial_state());
}

std::size_t SpatialAssembly::section_dof(std::size_t component, std::size_t control) const {
    if (component >= 3 || control >= section_count())
        throw std::out_of_range("CPEG8T section control index");
    return _field_layout.at(3 + component).begin + control;
}

void SpatialAssembly::refresh_constraints() {
    refresh_dirichlet_values();
    for (std::size_t c = 0; c < section_count(); ++c) {
        const auto& section = _definition.generalized_plane_strain[c];
        for (std::size_t i = 0; i < 3; ++i) {
            if (!section.prescribed[i])
                continue;
            const auto index = section_dof(i, c);
            const auto condition = std::lower_bound(_dirichlet_conditions.begin(),
                _dirichlet_conditions.end(),
                index,
                [](const DirichletCondition& candidate, std::size_t dof) { return candidate.dof < dof; });
            condition->value = spatial_detail::controlled_value(_definition,
                _time,
                _load_factor,
                *section.prescribed[i],
                false,
                section.functions[i]);
        }
    }
}

void SpatialAssembly::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    if (index >= contact_offset()) {
        dofs = _contact_constraints.at(index - contact_offset()).dofs;
        return;
    }
    if (index >= volume_contribution_count()) {
        const auto& boundary = _boundaries.at(index - volume_contribution_count());
        index = region_element_offset(boundary.region) + boundary.element;
    }
    const auto& values = _dofs.at(index);
    dofs.assign(values.begin(), values.end());
}

std::vector<std::size_t> SpatialAssembly::required_state_dofs(std::size_t first, std::size_t last) const {
    if (first > last || last > contribution_count())
        throw std::out_of_range("CPEG8T contribution range");
    std::set<std::size_t> selected;
    for (auto index = first; index < last; ++index) {
        std::vector<std::size_t> dofs;
        contribution_dofs(index, dofs);
        selected.insert(dofs.begin(), dofs.end());
    }
    return {selected.begin(), selected.end()};
}

void SpatialAssembly::validate_local_state(std::size_t first,
    std::size_t last,
    const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("CPEG8T state size mismatch");
    for (auto index : required_state_dofs(first, last))
        if (!std::isfinite(state[index]))
            throw std::domain_error("CPEG8T trial state must be finite");
    for (auto index = first; index < std::min(last, volume_contribution_count()); ++index) {
        const auto [r, e] = element_location(index);
        elements::Cpeg8Values local{};
        for (std::size_t i = 0; i < 23; ++i)
            local[i] = state[_dofs[index][i]];
        const IsotropicThermoelasticMaterial material(region(r).material);
        elements::Cpeg8Input input{material, geometry(r, e), local, local};
        input.strain_formulation = region(r).strain_formulation;
        (void)elements::evaluate_cpeg8t(input, {false, false, false, false});
    }
    validate_contacts(first, last, state);
}

SpatialContributionType SpatialAssembly::contribution_type(std::size_t index) const {
    if (index >= contact_offset())
        return _contact_constraints.at(index - contact_offset()).geometry.temperature_nodes.empty()
                   ? SpatialContributionType::mechanical_contact
                   : SpatialContributionType::thermal_contact;
    if (index < volume_contribution_count())
        return SpatialContributionType::volume;
    const auto& boundary =
        _definition.boundary_conditions.at(_boundaries.at(index - volume_contribution_count()).definition);
    switch (boundary.type) {
    case BoundaryConditionType::pressure:
        return SpatialContributionType::pressure;
    case BoundaryConditionType::traction:
        return SpatialContributionType::traction;
    case BoundaryConditionType::heat_flux:
        return SpatialContributionType::heat_flux;
    case BoundaryConditionType::convection:
        return SpatialContributionType::convection;
    default:
        throw std::logic_error("Invalid CPEG8T boundary kind");
    }
}

elements::Cpeg8Result
SpatialAssembly::compute_boundary(std::size_t index, const std::vector<double>& local, bool jacobian) const {
    const auto& side = _boundaries.at(index - volume_contribution_count());
    const auto& boundary = _definition.boundary_conditions.at(side.definition);
    if (local.size() != 23)
        throw std::invalid_argument("CPEG8T boundary requires a parent-element state");
    elements::Cpeg8Values values{};
    std::copy(local.begin(), local.end(), values.begin());
    const double load = spatial_detail::controlled_value(_definition,
        _time,
        _load_factor,
        boundary.value,
        boundary.scale_with_load,
        boundary.function);
    elements::Cpeg8BoundaryInput input{geometry(side.region, side.element),
        values,
        side.side,
        elements::Cpeg8BoundaryKind::pressure,
        load};
    input.current = boundary_uses_displaced_geometry(boundary, region(side.region));
    if (boundary.type == BoundaryConditionType::traction) {
        if (boundary.field != Field::displacement_x && boundary.field != Field::displacement_y)
            throw std::invalid_argument("CPEG8T traction must use an in-plane component");
        input.kind = boundary.field == Field::displacement_x ? elements::Cpeg8BoundaryKind::traction_x
                                                             : elements::Cpeg8BoundaryKind::traction_y;
    } else if (boundary.type == BoundaryConditionType::heat_flux)
        input.kind = elements::Cpeg8BoundaryKind::heat_flux;
    else if (boundary.type == BoundaryConditionType::convection) {
        input.kind = elements::Cpeg8BoundaryKind::convection;
        const auto data = convection_values(boundary);
        input.value = data.coefficient;
        input.ambient = data.ambient;
    }
    return elements::evaluate_cpeg8t_boundary(input, jacobian);
}

void SpatialAssembly::validate_state(const std::vector<double>& state) const {
    validate_local_state(0, contribution_count(), state);
}
} // namespace fuelsim::plane
