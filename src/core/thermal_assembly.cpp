#include "thermal_assembly.hpp"
#include "dc3d20.hpp"
#include "dc3d8.hpp"
#include "dcax4.hpp"
#include "dcax8.hpp"
#include "thermal_boundary.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim::thermal {
SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad4Mesh& mesh)
    : SpatialLayout(definition,
          spatial_detail::resolve_block_ids(definition, mesh),
          spatial_detail::DofLayout::axisymmetric_rz),
      _axisymmetric(true) {
    for (const auto& node : mesh.nodes())
        _coordinates.push_back({node.r, node.z, 0.0});
    for (const auto& element : mesh.elements())
        _connectivity.emplace_back(element.nodes.begin(), element.nodes.end());
    initialize(mesh, ThermalElement::dcax4);
}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredQuad8Mesh& mesh)
    : SpatialLayout(definition,
          spatial_detail::resolve_block_ids(definition, mesh),
          spatial_detail::DofLayout::axisymmetric_rz),
      _axisymmetric(true) {
    for (const auto& node : mesh.nodes())
        _coordinates.push_back({node.r, node.z, 0.0});
    for (const auto& element : mesh.elements())
        _connectivity.emplace_back(element.nodes.begin(), element.nodes.end());
    initialize(mesh, ThermalElement::dcax8);
}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredHex8Mesh& mesh)
    : SpatialLayout(definition,
          spatial_detail::resolve_block_ids(definition, mesh),
          spatial_detail::DofLayout::cartesian_3d),
      _axisymmetric(false) {
    for (const auto& node : mesh.nodes())
        _coordinates.push_back(node);
    for (const auto& element : mesh.elements())
        _connectivity.emplace_back(element.nodes.begin(), element.nodes.end());
    initialize(mesh, ThermalElement::dc3d8);
}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition, const UnstructuredHex20Mesh& mesh)
    : SpatialLayout(definition,
          spatial_detail::resolve_block_ids(definition, mesh),
          spatial_detail::DofLayout::cartesian_3d),
      _axisymmetric(false) {
    for (const auto& node : mesh.nodes())
        _coordinates.push_back(node);
    for (const auto& element : mesh.elements())
        _connectivity.emplace_back(element.nodes.begin(), element.nodes.end());
    initialize(mesh, ThermalElement::dc3d20);
}

void SpatialAssembly::initialize(const UnstructuredMeshMetadata& mesh, ThermalElement topology) {
    if (_definition.physics != Physics::thermal)
        throw std::invalid_argument("Thermal assembly requires thermal physics");
    const auto missing = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> node_counts, element_counts, source_to_global(_coordinates.size(), missing);
    std::vector<std::size_t> element_regions(_connectivity.size(), missing);
    for (std::size_t r = 0; r < region_count(); ++r) {
        const auto& region = _definition.regions[r];
        if (region.thermal_element != topology || region.strain_formulation != StrainFormulation::small
            || !region.material.functions || !region.material.functions->thermal.function)
            throw std::invalid_argument("Thermal region topology or material mismatch");
        std::vector<std::size_t> nodes;
        std::size_t count = 0;
        for (std::size_t e = 0; e < _connectivity.size(); ++e)
            if (mesh.element_block_ids()[e] == _block_ids[r]) {
                element_regions[e] = r;
                ++count;
                nodes.insert(nodes.end(), _connectivity[e].begin(), _connectivity[e].end());
            }
        std::sort(nodes.begin(), nodes.end());
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        node_counts.push_back(nodes.size());
        element_counts.push_back(count);
        _source_nodes.push_back(std::move(nodes));
    }
    initialize_counts(node_counts, element_counts);
    initialize_shared_nodes(_source_nodes);
    for (std::size_t r = 0; r < region_count(); ++r)
        for (std::size_t n = 0; n < _source_nodes[r].size(); ++n)
            source_to_global[_source_nodes[r][n]] = global_temperature_node(r, n);
    for (std::size_t r = 0; r < region_count(); ++r)
        for (std::size_t e = 0; e < _connectivity.size(); ++e) {
            if (element_regions[e] != r)
                continue;
            Contribution c;
            c.region = r;
            for (auto n : _connectivity[e])
                c.dofs.push_back(source_to_global[n]);
            if (topology == ThermalElement::dcax4) {
                std::array<CartesianPoint3, 4> coordinates{};
                for (std::size_t n = 0; n < 4; ++n)
                    coordinates[n] = _coordinates[_connectivity[e][n]];
                c.geometry = elements::make_dcax4_geometry(coordinates);
            } else if (topology == ThermalElement::dc3d20) {
                std::array<CartesianPoint3, 20> coordinates{};
                for (std::size_t n = 0; n < 20; ++n)
                    coordinates[n] = _coordinates[_connectivity[e][n]];
                c.geometry = elements::make_dc3d20_geometry(coordinates);
            } else {
                std::array<CartesianPoint3, 8> coordinates{};
                for (std::size_t n = 0; n < 8; ++n)
                    coordinates[n] = _coordinates[_connectivity[e][n]];
                c.geometry = _axisymmetric ? elements::make_dcax8_geometry(coordinates)
                                           : elements::make_dc3d8_geometry(coordinates);
            }
            _source_elements.push_back(e);
            _contributions.push_back(std::move(c));
        }
    constexpr std::array<std::array<std::size_t, 8>, 6> faces{{{0, 1, 5, 4, 8, 13, 16, 12},
        {1, 2, 6, 5, 9, 14, 17, 13},
        {2, 3, 7, 6, 10, 15, 18, 14},
        {0, 4, 7, 3, 12, 19, 15, 11},
        {0, 3, 2, 1, 11, 10, 9, 8},
        {4, 5, 6, 7, 16, 17, 18, 19}}};
    const auto side_nodes = [&](const ElementSide& side) {
        if (element_regions.at(side.element) == missing)
            throw std::invalid_argument("Thermal surface belongs to an inactive region");
        const auto& connectivity = _connectivity[side.element];
        std::vector<std::size_t> nodes;
        if (_axisymmetric) {
            nodes = {connectivity.at(side.local_side), connectivity.at((side.local_side + 1) % 4)};
            if (connectivity.size() == 8)
                nodes.push_back(connectivity.at(4 + side.local_side));
        } else
            for (std::size_t n = 0; n < (connectivity.size() == 20 ? 8U : 4U); ++n)
                nodes.push_back(connectivity.at(faces.at(side.local_side)[n]));
        return nodes;
    };
    for (std::size_t b = 0; b < _definition.boundary_conditions.size(); ++b) {
        const auto& boundary = _definition.boundary_conditions[b];
        if (boundary.type != BoundaryConditionType::dirichlet && boundary.type != BoundaryConditionType::heat_flux
            && boundary.type != BoundaryConditionType::convection)
            throw std::invalid_argument("Thermal physics does not accept mechanical boundary conditions");
        if (boundary.type == BoundaryConditionType::dirichlet && boundary.field != Field::temperature)
            throw std::invalid_argument("Thermal Dirichlet boundary requires temperature");
        std::vector<std::size_t> fixed;
        if (boundary.type == BoundaryConditionType::dirichlet)
            for (const auto& nodeset : mesh.node_sets())
                if (nodeset.name == boundary.boundary)
                    fixed = nodeset.nodes;
        if (fixed.empty())
            for (const auto& side : mesh.side_set(boundary.boundary).sides) {
                if (element_regions[side.element] == missing)
                    throw std::invalid_argument("Thermal boundary belongs to an inactive region");
                Contribution c;
                c.region = element_regions[side.element];
                c.boundary = b;
                c.type = boundary.type == BoundaryConditionType::convection ? SpatialContributionType::convection
                                                                            : SpatialContributionType::heat_flux;
                const auto nodes = side_nodes(side);
                if (boundary.type == BoundaryConditionType::dirichlet) {
                    fixed.insert(fixed.end(), nodes.begin(), nodes.end());
                    continue;
                }
                std::vector<CartesianPoint3> coordinates;
                for (auto n : nodes) {
                    c.dofs.push_back(source_to_global[n]);
                    coordinates.push_back(_coordinates[n]);
                }
                c.geometry = elements::make_thermal_boundary_geometry(coordinates, _axisymmetric);
                _contributions.push_back(std::move(c));
            }
        std::sort(fixed.begin(), fixed.end());
        fixed.erase(std::unique(fixed.begin(), fixed.end()), fixed.end());
        for (auto n : fixed) {
            if (source_to_global.at(n) == missing)
                throw std::invalid_argument("Thermal temperature boundary contains an inactive node");
            add_dirichlet(source_to_global[n], b);
        }
    }
    for (std::size_t c = 0; c < _definition.contacts.size(); ++c) {
        const auto& contact = _definition.contacts[c];
        if (!contact.thermal || contact.mechanical || contact.gap_conductance_pressure_derivative != 0
            || contact.thermal_discretization != ThermalContactDiscretization::surface_to_surface)
            throw std::invalid_argument(
                "Thermal interfaces require pressure-independent surface-to-surface heat transfer");
        std::vector<std::vector<std::size_t>> primary_nodes;
        std::vector<std::vector<CartesianPoint3>> primary;
        for (const auto& side : mesh.side_set(contact.primary).sides) {
            primary_nodes.push_back(side_nodes(side));
            std::vector<CartesianPoint3> coordinates;
            for (auto n : primary_nodes.back())
                coordinates.push_back(_coordinates[n]);
            primary.push_back(std::move(coordinates));
        }
        for (const auto& side : mesh.side_set(contact.secondary).sides) {
            const auto nodes = side_nodes(side);
            std::vector<CartesianPoint3> coordinates;
            for (auto n : nodes)
                coordinates.push_back(_coordinates[n]);
            for (auto point : elements::make_thermal_interface_points(coordinates, primary, _axisymmetric)) {
                Contribution entry;
                entry.type = SpatialContributionType::thermal_contact;
                entry.boundary = c;
                for (auto n : nodes)
                    entry.dofs.push_back(source_to_global[n]);
                for (auto n : primary_nodes[point.primary])
                    entry.dofs.push_back(source_to_global[n]);
                auto unique = entry.dofs;
                std::sort(unique.begin(), unique.end());
                if (std::adjacent_find(unique.begin(), unique.end()) != unique.end())
                    throw std::invalid_argument("Thermal interface sides must use independent source nodes");
                entry.interface = std::move(point);
                _contributions.push_back(std::move(entry));
            }
        }
    }
    refresh_dirichlet_values();
}

SpatialContributionType SpatialAssembly::contribution_type(std::size_t index) const {
    return _contributions.at(index).type;
}

void SpatialAssembly::contribution_dofs(std::size_t index, std::vector<std::size_t>& dofs) const {
    dofs = _contributions.at(index).dofs;
}

void SpatialAssembly::validate_state(const std::vector<double>& state) const {
    validate_local_state(0, contribution_count(), state);
}

void SpatialAssembly::validate_local_state(std::size_t first,
    std::size_t last,
    const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument("Thermal state size mismatch");
    for (auto n : required_state_dofs(first, last))
        if (!std::isfinite(state[n]))
            throw std::domain_error("Thermal temperature must be finite");
}

std::vector<std::size_t> SpatialAssembly::required_state_dofs(std::size_t first, std::size_t last) const {
    if (first > last || last > contribution_count())
        throw std::invalid_argument("Thermal contribution range mismatch");
    std::vector<std::size_t> result;
    for (std::size_t i = first; i < last; ++i)
        result.insert(result.end(), _contributions[i].dofs.begin(), _contributions[i].dofs.end());
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

elements::ThermalResult SpatialAssembly::evaluate(std::size_t index,
    const std::vector<double>& state,
    const std::vector<double>& previous,
    double step,
    double begin,
    bool jacobian) const {
    const auto& c = _contributions.at(index);
    if (c.type == SpatialContributionType::thermal_contact) {
        const auto& contact = _definition.contacts[c.boundary];
        const GapHeatProperties material{contact.gap_conductivity,
            contact.minimum_gap,
            contact.gap_heat_conductance_law,
            contact.gap_conductance,
            contact.gap_conductance_clearance_derivative,
            0.0,
            contact.gap_conductance_temperature_derivative,
            contact.gap_conductance_reference_temperature,
            0.0};
        return elements::evaluate_thermal_interface(c.interface, material, state, jacobian);
    }
    if (c.type != SpatialContributionType::volume) {
        const auto& b = _definition.boundary_conditions[c.boundary];
        const auto convection = convection_values(b);
        const double flux =
            spatial_detail::controlled_value(_definition, _time, _load_factor, b.value, b.scale_with_load, b.function);
        return elements::evaluate_thermal_boundary(c.geometry,
            state,
            c.type == SpatialContributionType::heat_flux ? flux : 0.0,
            c.type == SpatialContributionType::convection ? convection.coefficient : 0.0,
            convection.ambient,
            jacobian);
    }
    const auto& r = region(c.region);
    std::vector<double> old;
    if (step > 0)
        for (auto n : c.dofs)
            old.push_back(previous.at(n));
    elements::ThermalInput input{r.material.functions->thermal,
        c.geometry,
        state,
        old,
        r.initial_temperature,
        _time,
        step,
        !previous.empty() && _time > begin ? region_heat_source_average(c.region, begin, _time)
                                           : region_heat_source(c.region)};
    switch (r.thermal_element) {
    case ThermalElement::dcax4:
        return elements::evaluate_dcax4(input, jacobian);
    case ThermalElement::dcax8:
        return elements::evaluate_dcax8(input, jacobian);
    case ThermalElement::dc3d8:
        return elements::evaluate_dc3d8(input, jacobian);
    case ThermalElement::dc3d20:
        return elements::evaluate_dc3d20(input, jacobian);
    }
    throw std::logic_error("Unknown thermal element");
}
} // namespace fuelsim::thermal

namespace fuelsim::thermal {
InterfaceSummary SpatialAssembly::summarize_interface(std::size_t contact, const std::vector<double>& state) const {
    InterfaceSummary summary;
    for (std::size_t index = 0; index < _contributions.size(); ++index) {
        const auto& c = _contributions[index];
        if (c.type != SpatialContributionType::thermal_contact || c.boundary != contact)
            continue;
        std::vector<double> local;
        for (auto n : c.dofs)
            local.push_back(state.at(n));
        const auto result = evaluate(index, local, {}, 0.0, 0.0, false);
        summary.minimum_gap = std::min(summary.minimum_gap, c.interface.gap);
        for (std::size_t i = 0; i < c.interface.secondary_shape.size(); ++i)
            summary.total_heat_rate += result.residual[i];
        ++summary.projected_contact_nodes;
    }
    return summary;
}
} // namespace fuelsim::thermal
