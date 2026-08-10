#include "assembly.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {

namespace {

Quad4Coordinates element_coordinates(const RegionMesh& mesh,
                                     const Quad4Element& element) {
    Quad4Coordinates coordinates{};
    for (std::size_t node = 0; node < element.nodes.size(); ++node)
        coordinates[node] = mesh.nodes().at(element.nodes[node]);
    return coordinates;
}

std::size_t checked_layout_node_count(const std::vector<RegionMesh>& meshes) {
    std::size_t result = 0;
    for (const RegionMesh& mesh : meshes) {
        if (mesh.nodes().size() >
            std::numeric_limits<std::size_t>::max() - result)
            throw std::length_error("Spatial layout node count overflows");
        result += mesh.nodes().size();
    }
    return result;
}

} // namespace

SpatialLayout::SpatialLayout(SpatialDefinition definition,
                             std::vector<std::int64_t> block_ids,
                             std::vector<RegionMesh> meshes)
    : _definition(std::move(definition)), _block_ids(std::move(block_ids)),
      _meshes(std::move(meshes)),
      _dof_map(checked_layout_node_count(_meshes)) {
    _node_offsets.reserve(_meshes.size() + 1);
    _element_offsets.reserve(_meshes.size() + 1);
    _node_offsets.push_back(0);
    _element_offsets.push_back(0);
    for (const RegionMesh& mesh : _meshes) {
        _node_offsets.push_back(_node_offsets.back() + mesh.nodes().size());
        _element_offsets.push_back(_element_offsets.back() +
                                   mesh.elements().size());
    }
}

const SpatialDefinition& SpatialLayout::definition() const noexcept {
    return _definition;
}

const DofMap& SpatialLayout::dof_map() const noexcept {
    return _dof_map;
}

std::size_t SpatialLayout::region_count() const noexcept {
    return _definition.regions.size();
}

const RegionDefinition& SpatialLayout::region(std::size_t index) const {
    return _definition.regions.at(index);
}

const RegionMesh& SpatialLayout::region_mesh(std::size_t index) const {
    return _meshes.at(index);
}

std::size_t SpatialLayout::region_node_offset(std::size_t index) const {
    if (index >= region_count())
        throw std::out_of_range("Spatial layout region index is out of range");
    return _node_offsets[index];
}

std::size_t SpatialLayout::region_element_offset(std::size_t index) const {
    if (index >= region_count())
        throw std::out_of_range("Spatial layout region index is out of range");
    return _element_offsets[index];
}

std::size_t SpatialLayout::volume_contribution_count() const noexcept {
    return _element_offsets.back();
}

const Quad4RzGeometry& SpatialLayout::element_geometry(
    std::size_t region_value, std::size_t element) const {
    return _region_geometries.at(region_value).at(element);
}

SpatialLayout::ResolvedBoundary SpatialLayout::resolve_boundary(
    const UnstructuredQuad4Mesh& source_mesh,
    const std::string& name) const {
    const std::int64_t block_id = source_mesh.side_set_block_id(name);
    const auto found =
        std::find(_block_ids.begin(), _block_ids.end(), block_id);
    if (found == _block_ids.end())
        throw std::invalid_argument(
            "Boundary belongs to an undeclared block: " + name);
    const std::size_t region =
        static_cast<std::size_t>(found - _block_ids.begin());
    return {region, _meshes[region].map_side_set(source_mesh, name)};
}

std::size_t SpatialLayout::global_node(std::size_t region,
                                       std::size_t local_node) const {
    if (local_node >= region_mesh(region).nodes().size())
        throw std::out_of_range("Spatial layout local node is out of range");
    return _node_offsets.at(region) + local_node;
}

std::pair<std::size_t, std::array<std::size_t, 2>>
SpatialLayout::edge_parent(std::size_t region,
                           const Line2BoundaryElement& edge) const {
    const RegionMesh& mesh = region_mesh(region);
    std::size_t parent = mesh.elements().size();
    std::array<std::size_t, 2> local_nodes{};
    for (std::size_t element_index = 0;
         element_index < mesh.elements().size(); ++element_index) {
        const Quad4Element& element = mesh.elements()[element_index];
        std::array<std::size_t, 2> candidate{};
        bool contains = false;
        for (std::size_t side = 0; side < element.nodes.size(); ++side) {
            const std::size_t next = (side + 1U) % element.nodes.size();
            const bool forward = element.nodes[side] == edge.nodes[0] &&
                                 element.nodes[next] == edge.nodes[1];
            const bool reverse = element.nodes[side] == edge.nodes[1] &&
                                 element.nodes[next] == edge.nodes[0];
            if (forward || reverse) {
                candidate = {{side, next}};
                contains = true;
                break;
            }
        }
        if (!contains)
            continue;
        if (parent != mesh.elements().size())
            throw std::invalid_argument(
                "Boundary edge has more than one adjacent region element");
        parent = element_index;
        local_nodes = candidate;
    }
    if (parent == mesh.elements().size())
        throw std::invalid_argument(
            "Boundary edge has no adjacent region element");
    return {parent, local_nodes};
}

void SpatialLayout::build_volume_geometries() {
    _region_geometries.resize(region_count());
    for (std::size_t region_index = 0; region_index < region_count();
         ++region_index) {
        const RegionMesh& mesh = _meshes[region_index];
        std::vector<Quad4RzGeometry>& geometries =
            _region_geometries[region_index];
        geometries.reserve(mesh.elements().size());
        for (const Quad4Element& element : mesh.elements())
            geometries.push_back(
                make_quad4_rz_geometry(element_coordinates(mesh, element)));
    }
}

const SpatialDefinition& SpatialAssembly::definition() const noexcept {
    return _layout.definition();
}

const DofMap& SpatialAssembly::dof_map() const noexcept {
    return _layout.dof_map();
}

std::size_t SpatialAssembly::region_count() const noexcept {
    return _layout.region_count();
}

std::size_t SpatialAssembly::region_index(const std::string& name) const {
    for (std::size_t region = 0; region < region_count(); ++region) {
        if (_layout._definition.regions[region].name == name)
            return region;
    }
    throw std::invalid_argument("Unknown region: " + name);
}

const RegionDefinition& SpatialAssembly::region(std::size_t index) const {
    return _layout.region(index);
}

const RegionMesh& SpatialAssembly::region_mesh(std::size_t index) const {
    return _layout.region_mesh(index);
}

std::size_t SpatialAssembly::region_node_offset(std::size_t index) const {
    return _layout.region_node_offset(index);
}

std::size_t SpatialAssembly::region_element_count(std::size_t index) const {
    return _layout.region_mesh(index).elements().size();
}

std::size_t SpatialAssembly::region_element_offset(std::size_t index) const {
    return _layout.region_element_offset(index);
}

std::size_t SpatialAssembly::volume_contribution_count() const noexcept {
    return _layout.volume_contribution_count();
}

SpatialAssembly::ContributionRanges
SpatialAssembly::contribution_ranges() const noexcept {
    const std::size_t thermal_begin = volume_contribution_count();
    const std::size_t mechanical_begin =
        thermal_begin + _contact._thermal_contributions.size();
    const std::size_t pressure_begin =
        mechanical_begin + _contact._mechanical_contributions.size();
    const std::size_t traction_begin =
        pressure_begin + _boundary.pressure_contribution_count();
    const std::size_t convection_begin =
        traction_begin + _boundary.traction_contribution_count();
    return {thermal_begin,
            mechanical_begin,
            pressure_begin,
            traction_begin,
            convection_begin,
            convection_begin + _boundary.convection_contribution_count()};
}

SpatialAssembly::ContributionLocation SpatialAssembly::locate_contribution(
    std::size_t contribution_index) const {
    const ContributionRanges ranges = contribution_ranges();
    if (contribution_index >= ranges.end)
        throw std::out_of_range(
            "SpatialAssembly contribution index is out of range");
    if (contribution_index < ranges.thermal_begin)
        return {SpatialContributionType::volume, contribution_index};
    if (contribution_index < ranges.mechanical_begin)
        return {SpatialContributionType::thermal_contact,
                contribution_index - ranges.thermal_begin};
    if (contribution_index < ranges.pressure_begin)
        return {SpatialContributionType::mechanical_contact,
                contribution_index - ranges.mechanical_begin};
    if (contribution_index < ranges.traction_begin)
        return {SpatialContributionType::pressure,
                contribution_index - ranges.pressure_begin};
    if (contribution_index < ranges.convection_begin)
        return {SpatialContributionType::traction,
                contribution_index - ranges.traction_begin};
    return {SpatialContributionType::convection,
            contribution_index - ranges.convection_begin};
}

SpatialContributionType
SpatialAssembly::contribution_type(std::size_t contribution_index) const {
    return locate_contribution(contribution_index).type;
}

const Quad4RzGeometry&
SpatialAssembly::region_element_geometry(std::size_t region_value,
                                         std::size_t element_index) const {
    return _layout.element_geometry(region_value, element_index);
}

double SpatialAssembly::region_heat_source(std::size_t region_value) const {
    return _boundary.region_heat_source(region_value, _layout);
}

// Boundary-condition assembly.
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

struct BoundaryEdgeData final {
    std::array<std::size_t, 4> nodes;
    std::array<RzPoint, 2> coordinates;
    std::array<std::size_t, 2> local_nodes;
};

BoundaryEdgeData boundary_edge_data(
    const SpatialLayout& layout,
    const SpatialLayout::ResolvedBoundary& boundary,
    const Line2BoundaryElement& edge) {
    const RegionMesh& mesh = layout.region_mesh(boundary.region);
    const auto parent = layout.edge_parent(boundary.region, edge);
    const Quad4Element& element = mesh.elements().at(parent.first);
    BoundaryEdgeData result{{}, {}, parent.second};
    for (std::size_t node = 0; node < result.nodes.size(); ++node)
        result.nodes[node] =
            layout.global_node(boundary.region, element.nodes[node]);
    for (std::size_t node = 0; node < result.coordinates.size(); ++node)
        result.coordinates[node] =
            mesh.nodes().at(element.nodes[parent.second[node]]);
    return result;
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
            for (const Line2BoundaryElement& edge :
                 resolved.boundary.elements) {
                const BoundaryEdgeData data =
                    boundary_edge_data(layout, resolved, edge);
                _pressure_contributions.push_back(
                    {load, data.nodes, make_line2_rz_pressure_geometry(
                                           data.coordinates, data.local_nodes)});
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
            for (const Line2BoundaryElement& edge :
                 resolved.boundary.elements) {
                const BoundaryEdgeData data =
                    boundary_edge_data(layout, resolved, edge);
                _traction_contributions.push_back(
                    {load, data.nodes, make_line2_rz_traction_geometry(
                                           data.coordinates, data.local_nodes)});
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
            for (const Line2BoundaryElement& edge :
                 resolved.boundary.elements) {
                const BoundaryEdgeData data =
                    boundary_edge_data(layout, resolved, edge);
                _convection_contributions.push_back(
                    {load, data.nodes, make_line2_rz_convection_geometry(
                                           data.coordinates, data.local_nodes)});
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

// Contact assembly and state transactions.

std::size_t SpatialAssembly::contact_count() const noexcept {
    return _contact.contact_count(_layout);
}

const ContactDefinition& SpatialAssembly::contact(std::size_t index) const {
    return _contact.contact(index, _layout);
}

const std::vector<std::vector<ContactPointHistory>>&
SpatialAssembly::committed_contact_histories() const noexcept {
    return _contact.committed_histories();
}

bool SpatialAssembly::uses_augmented_contact() const noexcept {
    return _contact.uses_augmented_contact(_layout);
}

AugmentedContactUpdate SpatialAssembly::update_augmented_contact_multipliers(
    const std::vector<double>& state, std::size_t completed_updates) {
    return _contact.update_augmented_multipliers(*this, state,
                                                 completed_updates);
}

void SpatialAssembly::commit_contact_state(const std::vector<double>& state) {
    _contact.commit_state(*this, state);
}

void SpatialAssembly::restore_contact_state(
    const std::vector<double>& state,
    std::vector<std::vector<ContactPointHistory>> histories) {
    _contact.restore_state(*this, state, std::move(histories));
}

std::size_t
ContactAssembly::contact_count(const SpatialLayout& layout) const noexcept {
    return layout.definition().contacts.size();
}

const ContactDefinition&
ContactAssembly::contact(std::size_t index, const SpatialLayout& layout) const {
    return layout.definition().contacts.at(index);
}

const std::vector<std::vector<ContactPointHistory>>&
ContactAssembly::committed_histories() const noexcept {
    return _contact_histories;
}

bool ContactAssembly::uses_augmented_contact(
    const SpatialLayout& layout) const noexcept {
    return std::any_of(
        layout.definition().contacts.begin(), layout.definition().contacts.end(),
        [](const ContactDefinition& contact) {
            return contact.mechanical &&
                   contact.mechanical_formulation ==
                       MechanicalContactFormulation::augmented_lagrangian;
        });
}

AugmentedContactUpdate ContactAssembly::update_augmented_multipliers(
    SpatialAssembly& assembly, const std::vector<double>& state,
    std::size_t completed_updates) {
    if (state.size() != assembly.dof_count())
        throw std::invalid_argument(
            "SpatialAssembly augmented-contact state size mismatch");
    AugmentedContactUpdate result;
    result.penetration_tolerance = std::numeric_limits<double>::infinity();
    std::vector<std::vector<ContactPointHistory>> staged = _contact_histories;
    for (std::size_t contact_value = 0;
         contact_value < contact_count(assembly._layout);
         ++contact_value) {
        const ContactDefinition& definition =
            assembly._layout.definition().contacts[contact_value];
        if (!definition.mechanical ||
            definition.mechanical_formulation !=
                MechanicalContactFormulation::augmented_lagrangian)
            continue;
        result.penetration_tolerance = std::min(
            result.penetration_tolerance, definition.penetration_tolerance);
        const std::vector<ContactNodeSummary> nodes =
            assembly.summarize_contact_nodes(contact_value, state);
        double contact_penetration = 0.0;
        double contact_constraint_violation = 0.0;
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            contact_penetration =
                std::max(contact_penetration, std::max(-nodes[node].gap, 0.0));
            const bool captured =
                staged[contact_value][node].normal_multiplier > 0.0 ||
                nodes[node].gap <= 0.0;
            if (captured)
                contact_constraint_violation = std::max(
                    contact_constraint_violation, std::abs(nodes[node].gap));
        }
        result.maximum_penetration =
            std::max(result.maximum_penetration, contact_penetration);
        result.maximum_constraint_violation = std::max(
            result.maximum_constraint_violation, contact_constraint_violation);
        if (contact_constraint_violation <= definition.penetration_tolerance)
            continue;
        result.converged = false;
        if (completed_updates >= definition.maximum_augmented_iterations) {
            result.update_allowed = false;
            continue;
        }
        for (std::size_t node = 0; node < nodes.size(); ++node) {
            ContactPointHistory& history = staged[contact_value][node];
            history.normal_multiplier =
                std::max(0.0, history.normal_multiplier -
                                  definition.penalty * nodes[node].gap);
        }
    }
    if (!std::isfinite(result.penetration_tolerance))
        result.penetration_tolerance = 0.0;
    if (!result.converged && result.update_allowed)
        _contact_histories.swap(staged);
    return result;
}

void ContactAssembly::commit_state(SpatialAssembly& assembly,
                                   const std::vector<double>& state) {
    if (state.size() != assembly.dof_count())
        throw std::invalid_argument(
            "SpatialAssembly committed contact state size mismatch");
    assembly.update_mechanical_candidates(state);
    std::vector<std::vector<ContactPointHistory>> staged = _contact_histories;
    std::vector<std::vector<bool>> updated(contact_count(assembly._layout));
    for (std::size_t contact_value = 0;
         contact_value < contact_count(assembly._layout);
         ++contact_value)
        updated[contact_value].resize(_contact_histories[contact_value].size(),
                                      false);

    const std::size_t first_mechanical =
        assembly.contribution_ranges().mechanical_begin;
    for (std::size_t contribution = 0;
         contribution < _mechanical_contributions.size(); ++contribution) {
        const MechanicalContribution& candidate =
            _mechanical_contributions[contribution];
        if (!candidate.active)
            continue;
        const LocalValues local_state =
            assembly.contribution_state(first_mechanical + contribution, state);
        const LocalValues committed_state = assembly.contribution_state(
            first_mechanical + contribution, _committed_contact_solution);
        const ContactPointValue value =
            _mechanical_kernels[candidate.contact].value(
                candidate.geometry, local_state, committed_state,
                _contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected)
            continue;
        const ContactPointHistory trial =
            _mechanical_kernels[candidate.contact].trial_history(
                candidate.geometry, local_state, committed_state,
                _contact_histories[candidate.contact][candidate.secondary]);
        if (updated[candidate.contact][candidate.secondary]) {
            const ContactPointHistory& prior =
                staged[candidate.contact][candidate.secondary];
            const double scale =
                std::max({1.0, std::abs(prior.elastic_tangential_slip),
                          std::abs(trial.elastic_tangential_slip)});
            if (std::abs(prior.elastic_tangential_slip -
                         trial.elastic_tangential_slip) >
                    32.0 * std::numeric_limits<double>::epsilon() * scale ||
                prior.sliding != trial.sliding)
                throw std::logic_error(
                    "Mechanical half-edge contributions disagree on the "
                    "contact-node friction history");
            continue;
        }
        staged[candidate.contact][candidate.secondary] = trial;
        updated[candidate.contact][candidate.secondary] = true;
    }
    for (std::size_t contact_value = 0;
         contact_value < contact_count(assembly._layout);
         ++contact_value) {
        if (!assembly._layout.definition().contacts[contact_value].mechanical)
            continue;
        if (std::find(updated[contact_value].begin(),
                      updated[contact_value].end(),
                      false) != updated[contact_value].end())
            throw std::domain_error(
                "Cannot commit friction history for an unprojected contact "
                "node");
    }
    _contact_histories.swap(staged);
    _committed_contact_solution = state;
}

void ContactAssembly::restore_state(
    const SpatialAssembly& assembly,
    const std::vector<double>& state,
    std::vector<std::vector<ContactPointHistory>> histories) {
    if (state.size() != assembly.dof_count() ||
        histories.size() != contact_count(assembly._layout))
        throw std::invalid_argument(
            "SpatialAssembly restored contact state layout mismatch");
    for (std::size_t contact_value = 0;
         contact_value < contact_count(assembly._layout);
         ++contact_value) {
        if (histories[contact_value].size() !=
            _contact_histories[contact_value].size())
            throw std::invalid_argument(
                "SpatialAssembly restored contact history layout mismatch");
        for (const ContactPointHistory& history : histories[contact_value]) {
            if (!std::isfinite(history.elastic_tangential_slip) ||
                !std::isfinite(history.normal_multiplier) ||
                history.normal_multiplier < 0.0)
                throw std::invalid_argument(
                    "SpatialAssembly restored contact history is invalid");
        }
    }
    _committed_contact_solution = state;
    _contact_histories = std::move(histories);
}

void SpatialAssembly::update_mechanical_candidates(
    const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly contact-search state size mismatch");
    const ContributionRanges ranges = contribution_ranges();
    const GlobalStateView state_view(state);
    update_mechanical_candidates(ranges.mechanical_begin,
                                 ranges.pressure_begin, state_view);
}

void SpatialAssembly::update_mechanical_candidates(
    std::size_t contribution_begin, std::size_t contribution_end,
    const GlobalStateView& state) const {
    if (state.global_size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly contact-search shadow state size mismatch");
    const ContributionRanges ranges = contribution_ranges();
    for (const MechanicalContribution& contribution :
         _contact._mechanical_contributions)
        contribution.active = false;
    for (std::vector<bool>& nodes : _contact._projected_mechanical_nodes)
        std::fill(nodes.begin(), nodes.end(), false);
    const std::vector<std::vector<bool>> touched =
        touched_mechanical_nodes(contribution_begin, contribution_end);
    if (touched.empty())
        return;

    std::vector<std::vector<double>> minimum_distance(contact_count());
    std::vector<std::vector<std::size_t>> selected_primary(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        const std::size_t node_count =
            _contact._contact_histories[contact_value].size();
        minimum_distance[contact_value].assign(
            node_count, std::numeric_limits<double>::infinity());
        selected_primary[contact_value].assign(
            node_count, std::numeric_limits<std::size_t>::max());
    }
    std::vector<bool> projected(_contact._mechanical_contributions.size(),
                                false);
    for (std::size_t contribution = 0;
         contribution < _contact._mechanical_contributions.size(); ++contribution) {
        const MechanicalContribution& candidate =
            _contact._mechanical_contributions[contribution];
        if (!touched[candidate.contact][candidate.secondary])
            continue;
        const LocalValues local_state =
            contribution_state(ranges.mechanical_begin + contribution, state);
        const LocalValues committed_state = contribution_state(
            ranges.mechanical_begin + contribution,
            _contact._committed_contact_solution);
        const ContactPointValue value =
            _contact._mechanical_kernels[candidate.contact].value(
                candidate.geometry, local_state, committed_state,
                _contact._contact_histories[candidate.contact][candidate.secondary]);
        if (!value.projected)
            continue;
        projected[contribution] = true;
        const double distance = std::abs(value.gap);
        if (distance <
                minimum_distance[candidate.contact][candidate.secondary] ||
            (distance ==
                 minimum_distance[candidate.contact][candidate.secondary] &&
             candidate.primary <
                 selected_primary[candidate.contact][candidate.secondary])) {
            minimum_distance[candidate.contact][candidate.secondary] = distance;
            selected_primary[candidate.contact][candidate.secondary] =
                candidate.primary;
        }
    }

    for (std::size_t contribution = 0;
         contribution < _contact._mechanical_contributions.size(); ++contribution) {
        if (!projected[contribution])
            continue;
        const MechanicalContribution& candidate =
            _contact._mechanical_contributions[contribution];
        if (candidate.primary !=
            selected_primary[candidate.contact][candidate.secondary])
            continue;
        candidate.active = true;
        _contact._projected_mechanical_nodes[candidate.contact][candidate.secondary] =
            true;
    }
}

std::vector<std::vector<bool>> SpatialAssembly::touched_mechanical_nodes(
    std::size_t contribution_begin, std::size_t contribution_end) const {
    const ContributionRanges ranges = contribution_ranges();
    const std::size_t begin =
        std::max(contribution_begin, ranges.mechanical_begin);
    const std::size_t end = std::min(contribution_end, ranges.pressure_begin);
    if (begin >= end)
        return {};

    std::vector<std::vector<bool>> result(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value)
        result[contact_value].resize(
            _contact._contact_histories[contact_value].size(), false);
    for (std::size_t full = begin; full < end; ++full) {
        const MechanicalContribution& candidate =
            _contact._mechanical_contributions[full - ranges.mechanical_begin];
        result[candidate.contact][candidate.secondary] = true;
    }
    return result;
}

std::vector<ContactNodeSummary>
SpatialAssembly::summarize_contact_nodes(std::size_t contact_value,
                                       const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly contact summary state size mismatch");
    update_mechanical_candidates(state);
    const ResolvedBoundary& secondary = _contact._secondary_boundaries.at(contact_value);
    const RegionMesh& mesh = _layout._meshes[secondary.region];
    std::vector<ContactNodeSummary> result;
    result.reserve(secondary.boundary.nodes.size());
    for (std::size_t node : secondary.boundary.nodes) {
        result.push_back({mesh.nodes().at(node).r, mesh.nodes().at(node).z,
                          false, std::numeric_limits<std::size_t>::max(),
                          std::numeric_limits<double>::infinity(), 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0, 0.0, false});
    }
    const std::size_t first_mechanical =
        contribution_ranges().mechanical_begin;
    for (std::size_t contribution = 0;
         contribution < _contact._mechanical_contributions.size(); ++contribution) {
        const MechanicalContribution& candidate =
            _contact._mechanical_contributions[contribution];
        if (candidate.contact != contact_value)
            continue;
        if (!candidate.active)
            continue;
        const LocalValues local_state =
            contribution_state(first_mechanical + contribution, state);
        const LocalValues committed_state = contribution_state(
            first_mechanical + contribution, _contact._committed_contact_solution);
        const std::size_t secondary_index = candidate.secondary;
        const ContactPointValue value =
            _contact._mechanical_kernels[contact_value].value(
                candidate.geometry, local_state, committed_state,
                _contact._contact_histories[contact_value][secondary_index]);
        if (!value.projected)
            continue;
        ContactNodeSummary& node = result.at(secondary_index);
        const std::size_t primary_segment = candidate.primary;
        if (node.projected && node.primary_segment != primary_segment)
            throw std::logic_error(
                "Mechanical contact node has more than one active primary "
                "segment");
        node.projected = true;
        node.primary_segment = primary_segment;
        node.gap = std::min(node.gap, value.gap);
        node.tributary_area += value.tributary_area;
        node.tributary_length += value.tributary_length;
        node.contact_force += value.contact_force;
        node.tangential_force += value.tangential_force;
        node.elastic_tangential_slip = value.elastic_tangential_slip;
        node.sliding = value.sliding;
    }
    for (ContactNodeSummary& node : result) {
        if (node.tributary_area > 0.0) {
            node.pressure = node.contact_force / node.tributary_area;
            node.tangential_traction =
                node.tangential_force / node.tributary_area;
        }
    }
    return result;
}

void SpatialAssembly::validate_state(const std::vector<double>& state) const {
    if (state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly global state has the wrong size");
    update_mechanical_candidates(state);
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        if (!_layout._definition.contacts[contact_value].mechanical)
            continue;
        const std::size_t unprojected = static_cast<std::size_t>(std::count(
            _contact._projected_mechanical_nodes[contact_value].begin(),
            _contact._projected_mechanical_nodes[contact_value].end(), false));
        if (unprojected != 0)
            throw std::domain_error(
                "Mechanical contact '" +
                _layout._definition.contacts[contact_value].name +
                "' lost projection for " + std::to_string(unprojected) +
                " secondary nodes after searching the complete primary chain");
    }
}

void SpatialAssembly::validate_local_state(std::size_t contribution_begin,
                                           std::size_t contribution_end,
                                           const GlobalStateView& state) const {
    if (contribution_begin > contribution_end ||
        contribution_end > contribution_count())
        throw std::out_of_range(
            "SpatialAssembly contribution range is out of bounds");
    if (state.global_size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly local state has the wrong global size");
    update_mechanical_candidates(contribution_begin, contribution_end, state);
    const std::vector<std::vector<bool>> touched =
        touched_mechanical_nodes(contribution_begin, contribution_end);
    if (touched.empty())
        return;
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        if (!_layout._definition.contacts[contact_value].mechanical)
            continue;
        std::size_t unprojected = 0;
        for (std::size_t node = 0; node < touched[contact_value].size();
             ++node) {
            if (touched[contact_value][node] &&
                !_contact._projected_mechanical_nodes[contact_value][node])
                ++unprojected;
        }
        if (unprojected != 0)
            throw std::domain_error(
                "Mechanical contact '" +
                _layout._definition.contacts[contact_value].name +
                "' lost projection for " + std::to_string(unprojected) +
                " locally owned secondary nodes after searching the complete "
                "primary chain");
    }
}

std::vector<std::size_t>
SpatialAssembly::contact_secondary_source_nodes(std::size_t contact_value) const {
    const ResolvedBoundary& secondary = _contact._secondary_boundaries.at(contact_value);
    const RegionMesh& mesh = _layout._meshes.at(secondary.region);
    std::vector<std::size_t> result;
    result.reserve(secondary.boundary.nodes.size());
    for (const std::size_t node : secondary.boundary.nodes)
        result.push_back(mesh.source_node_ids().at(node));
    return result;
}

InterfaceSummary
SpatialAssembly::summarize_interface(std::size_t contact_value,
                                   const std::vector<double>& state) const {
    if (contact_value >= contact_count())
        throw std::out_of_range("SpatialAssembly contact index is out of range");
    InterfaceSummary summary = {
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        0.0,
        0.0,
        0.0,
        0.0,
        0,
        0,
        0,
        0.0,
    };
    const std::size_t first_thermal = contribution_ranges().thermal_begin;
    bool has_thermal = false;
    for (std::size_t contribution = 0;
         contribution < _contact._thermal_contributions.size(); ++contribution) {
        const ThermalContribution& candidate =
            _contact._thermal_contributions[contribution];
        if (candidate.contact != contact_value)
            continue;
        has_thermal = true;
        const LocalValues local_state =
            contribution_state(first_thermal + contribution, state);
        const HeatQuadratureValues values =
            _contact._thermal_kernels[contact_value].quadrature_values(
                candidate.geometry, local_state);
        for (const HeatQuadratureValue& value : values) {
            summary.minimum_gap = std::min(summary.minimum_gap, value.gap);
            summary.maximum_gap = std::max(summary.maximum_gap, value.gap);
            summary.total_heat_rate += value.weighted_measure * value.heat_flux;
        }
    }
    if (!has_thermal) {
        summary.minimum_gap = 0.0;
        summary.maximum_gap = 0.0;
    }

    bool has_mechanical = false;
    for (const MechanicalContribution& contribution :
         _contact._mechanical_contributions) {
        if (contribution.contact == contact_value) {
            has_mechanical = true;
            break;
        }
    }
    if (has_mechanical) {
        for (const ContactNodeSummary& node :
             summarize_contact_nodes(contact_value, state)) {
            if (!node.projected) {
                ++summary.unprojected_contact_nodes;
                continue;
            }
            ++summary.projected_contact_nodes;
            summary.minimum_contact_gap =
                std::min(summary.minimum_contact_gap, node.gap);
            summary.maximum_contact_pressure =
                std::max(summary.maximum_contact_pressure, node.pressure);
            summary.total_contact_force += node.contact_force;
            summary.total_tangential_force += node.tangential_force;
            if (node.pressure > 0.0) {
                ++summary.active_contact_nodes;
                summary.active_contact_length += node.tributary_length;
            }
        }
    } else {
        summary.minimum_contact_gap = 0.0;
    }
    return summary;
}

// Spatial assembly construction and contribution dispatch.
namespace {

Line2InterfaceSideCoordinates
edge_coordinates(const RegionMesh& mesh, const Line2BoundaryElement& edge) {
    return {{mesh.nodes().at(edge.nodes[0]), mesh.nodes().at(edge.nodes[1])}};
}

double planar_quad_area(const RegionMesh& mesh, const Quad4Element& element) {
    double twice_area = 0.0;
    for (std::size_t node = 0; node < element.nodes.size(); ++node) {
        const RzPoint& current = mesh.nodes().at(element.nodes[node]);
        const RzPoint& next =
            mesh.nodes().at(element.nodes[(node + 1U) % element.nodes.size()]);
        twice_area += current.r * next.z - next.r * current.z;
    }
    return 0.5 * std::abs(twice_area);
}

double minimum_boundary_normal_length(const SpatialLayout& layout,
                                      std::size_t region,
                                      const RegionBoundary& boundary) {
    const RegionMesh& mesh = layout.region_mesh(region);
    double result = std::numeric_limits<double>::infinity();
    for (const Line2BoundaryElement& edge : boundary.elements) {
        const auto parent = layout.edge_parent(region, edge);
        const Line2InterfaceSideCoordinates coordinates =
            edge_coordinates(mesh, edge);
        const double edge_length =
            std::hypot(coordinates[1].r - coordinates[0].r,
                       coordinates[1].z - coordinates[0].z);
        const double area =
            planar_quad_area(mesh, mesh.elements().at(parent.first));
        const double normal_length = area / edge_length;
        if (!std::isfinite(normal_length) || !(normal_length > 0.0))
            throw std::invalid_argument(
                "Contact boundary has a nonpositive characteristic element "
                "length");
        result = std::min(result, normal_length);
    }
    return result;
}

RzPoint element_centroid(const RegionMesh& mesh, const Quad4Element& element) {
    RzPoint centroid{0.0, 0.0};
    for (const std::size_t node : element.nodes) {
        centroid.r += mesh.nodes().at(node).r;
        centroid.z += mesh.nodes().at(node).z;
    }
    centroid.r /= static_cast<double>(element.nodes.size());
    centroid.z /= static_cast<double>(element.nodes.size());
    return centroid;
}

std::vector<RzPoint> boundary_parent_centroids(
    const SpatialLayout& layout, std::size_t region,
    const RegionBoundary& boundary) {
    const RegionMesh& mesh = layout.region_mesh(region);
    std::vector<RzPoint> result;
    result.reserve(boundary.elements.size());
    for (const Line2BoundaryElement& edge : boundary.elements)
        result.push_back(element_centroid(
            mesh, mesh.elements().at(layout.edge_parent(region, edge).first)));
    return result;
}

std::array<std::size_t, 4> interface_nodes(
    const SpatialLayout& layout, std::size_t secondary_region,
    const Line2BoundaryElement& secondary_edge, std::size_t primary_region,
    const Line2BoundaryElement& primary_edge) {
    return {layout.global_node(secondary_region, secondary_edge.nodes[0]),
            layout.global_node(secondary_region, secondary_edge.nodes[1]),
            layout.global_node(primary_region, primary_edge.nodes[0]),
            layout.global_node(primary_region, primary_edge.nodes[1])};
}

// Signed distance of `point` from the primary line along the base normal
// (tangent_z, -tangent_r)/length. This is the same normal convention as the
// raw gap in reference_normal_orientation in interface.cpp, so the sign of
// this value for a point that rides on the line is exactly the raw-gap sign
// the point would have if it were moved to that side of the line.
double
primary_line_side(const RzPoint& point,
                  const Line2InterfaceSideCoordinates& primary_coordinates) {
    const double tangent_r =
        primary_coordinates[1].r - primary_coordinates[0].r;
    const double tangent_z =
        primary_coordinates[1].z - primary_coordinates[0].z;
    const double length = std::hypot(tangent_r, tangent_z);
    return ((point.r - primary_coordinates[0].r) * tangent_z -
            (point.z - primary_coordinates[0].z) * tangent_r) /
           length;
}

// Zero-gap orientation hint from the material-side topology: the signed side
// of the secondary parent-element centroid relative to the primary line. If
// both parent centroids fall on the same strict side of the line, the two
// material bodies overlap next to the interface and no meaningful hint
// exists; returning zero keeps the explicit construction error for any
// secondary point that actually rides on the segment. If either centroid
// lies on the line, the corresponding element is degenerate and the same
// explicit error remains in force.
double zero_gap_orientation_hint(
    const RzPoint& secondary_centroid, const RzPoint& primary_centroid,
    const Line2InterfaceSideCoordinates& primary_coordinates) {
    const double secondary_side =
        primary_line_side(secondary_centroid, primary_coordinates);
    const double primary_side =
        primary_line_side(primary_centroid, primary_coordinates);
    if (secondary_side == 0.0 || primary_side == 0.0)
        return 0.0;
    const bool same_side = (secondary_side > 0.0 && primary_side > 0.0) ||
                           (secondary_side < 0.0 && primary_side < 0.0);
    return same_side ? 0.0 : secondary_side;
}

void validate_definitions(const SpatialDefinition& definition) {
    if (definition.regions.empty())
        throw std::invalid_argument(
            "SpatialAssembly requires at least one region");
    for (std::size_t region = 0; region < definition.regions.size(); ++region) {
        const RegionDefinition& value = definition.regions[region];
        if (value.name.empty() || (value.block.empty() && value.block_id < 0) ||
            (!value.block.empty() && value.block_id >= 0))
            throw std::invalid_argument(
                "SpatialAssembly regions require a name and exactly one block "
                "selector");
        if (!std::isfinite(value.volumetric_heat_source) ||
            value.volumetric_heat_source < 0.0)
            throw std::invalid_argument("SpatialAssembly region heat sources "
                                        "must be finite and nonnegative");
        if (!std::isfinite(value.initial_temperature) ||
            !(value.initial_temperature > 0.0))
            throw std::invalid_argument("SpatialAssembly region temperatures "
                                        "must be finite and positive");
        for (std::size_t previous = 0; previous < region; ++previous) {
            if (definition.regions[previous].name == value.name)
                throw std::invalid_argument("Duplicate region name: " +
                                            value.name);
            if ((!value.block.empty() &&
                 definition.regions[previous].block == value.block) ||
                (value.block_id >= 0 &&
                 definition.regions[previous].block_id == value.block_id))
                throw std::invalid_argument(
                    "Duplicate region block: " +
                    (value.block.empty() ? std::to_string(value.block_id)
                                         : value.block));
        }
    }
    for (std::size_t table = 0; table < definition.time_tables.size();
         ++table) {
        for (std::size_t previous = 0; previous < table; ++previous) {
            if (definition.time_tables[previous].name() ==
                definition.time_tables[table].name())
                throw std::invalid_argument(
                    "Duplicate time-table name: " +
                    definition.time_tables[table].name());
        }
    }
    for (std::size_t contact = 0; contact < definition.contacts.size();
         ++contact) {
        const ContactDefinition& value = definition.contacts[contact];
        if (value.name.empty() || value.primary.empty() ||
            value.secondary.empty())
            throw std::invalid_argument(
                "SpatialAssembly contact names and boundaries must not be empty");
        if (!value.thermal && !value.mechanical)
            throw std::invalid_argument(
                "Contact must enable thermal or mechanical coupling: " +
                value.name);
        if (value.primary == value.secondary)
            throw std::invalid_argument(
                "Contact primary and secondary must differ: " + value.name);
        if (value.thermal &&
            (!(value.gap_conductivity > 0.0) || !(value.minimum_gap > 0.0)))
            throw std::invalid_argument(
                "Thermal contact parameters must be positive: " + value.name);
        if (value.mechanical && !value.automatic_penalty &&
            (!std::isfinite(value.penalty) || !(value.penalty > 0.0)))
            throw std::invalid_argument(
                "Mechanical contact penalty must be positive: " + value.name);
        if (value.mechanical && (!std::isfinite(value.penalty_factor) ||
                                 !(value.penalty_factor > 0.0)))
            throw std::invalid_argument(
                "Mechanical contact penalty factor must be finite and "
                "positive: " +
                value.name);
        if (value.mechanical &&
            value.mechanical_formulation ==
                MechanicalContactFormulation::augmented_lagrangian &&
            (!std::isfinite(value.penetration_tolerance) ||
             !(value.penetration_tolerance > 0.0) ||
             value.maximum_augmented_iterations == 0))
            throw std::invalid_argument(
                "Augmented contact requires a positive penetration tolerance "
                "and iteration limit: " +
                value.name);
        if (!std::isfinite(value.friction_coefficient) ||
            value.friction_coefficient < 0.0)
            throw std::invalid_argument(
                "Mechanical contact friction coefficient must be finite and "
                "nonnegative: " +
                value.name);
        for (std::size_t previous = 0; previous < contact; ++previous) {
            if (definition.contacts[previous].name == value.name)
                throw std::invalid_argument("Duplicate contact name: " +
                                            value.name);
            if (value.mechanical && definition.contacts[previous].mechanical &&
                definition.contacts[previous].secondary == value.secondary)
                throw std::invalid_argument(
                    "Mechanical secondary boundary is reused: " +
                    value.secondary);
        }
    }
}

RegionBoundary ordered_connected_boundary(const RegionMesh& mesh,
                                          RegionBoundary boundary,
                                          const std::string& name) {
    if (boundary.elements.empty())
        throw std::invalid_argument("Contact side set is empty: " + name);

    std::vector<std::size_t> degree(mesh.nodes().size(), 0);
    for (const Line2BoundaryElement& edge : boundary.elements) {
        if (edge.nodes[0] == edge.nodes[1])
            throw std::invalid_argument("Contact side set has a zero edge: " +
                                        name);
        for (std::size_t node : edge.nodes) {
            if (++degree.at(node) > 2)
                throw std::invalid_argument(
                    "Contact side set must not branch: " + name);
        }
    }
    std::vector<std::size_t> endpoints;
    for (std::size_t node = 0; node < degree.size(); ++node) {
        if (degree[node] == 1)
            endpoints.push_back(node);
    }
    if (endpoints.size() != 2)
        throw std::invalid_argument(
            "Contact side set must be one open connected chain: " + name);

    const auto coordinate_less = [&](std::size_t lhs, std::size_t rhs) {
        const RzPoint& left = mesh.nodes().at(lhs);
        const RzPoint& right = mesh.nodes().at(rhs);
        if (left.r != right.r)
            return left.r < right.r;
        if (left.z != right.z)
            return left.z < right.z;
        return lhs < rhs;
    };
    std::size_t current = coordinate_less(endpoints[0], endpoints[1])
                              ? endpoints[0]
                              : endpoints[1];
    std::vector<bool> used(boundary.elements.size(), false);
    std::vector<Line2BoundaryElement> ordered;
    std::vector<std::size_t> nodes = {current};
    ordered.reserve(boundary.elements.size());
    nodes.reserve(boundary.elements.size() + 1);
    for (std::size_t position = 0; position < boundary.elements.size();
         ++position) {
        std::size_t selected = boundary.elements.size();
        Line2BoundaryElement oriented{};
        for (std::size_t edge = 0; edge < boundary.elements.size(); ++edge) {
            if (used[edge])
                continue;
            const Line2BoundaryElement& candidate = boundary.elements[edge];
            if (candidate.nodes[0] == current) {
                selected = edge;
                oriented = candidate;
                break;
            }
            if (candidate.nodes[1] == current) {
                selected = edge;
                oriented = {{{candidate.nodes[1], candidate.nodes[0]}}};
                break;
            }
        }
        if (selected == boundary.elements.size())
            throw std::invalid_argument(
                "Contact side set must be one connected chain: " + name);
        used[selected] = true;
        ordered.push_back(oriented);
        current = oriented.nodes[1];
        nodes.push_back(current);
    }
    boundary.nodes = std::move(nodes);
    boundary.elements = std::move(ordered);
    return boundary;
}

double projection_fraction(const RzPoint& point,
                           const Line2InterfaceSideCoordinates& segment) {
    const double dr = segment[1].r - segment[0].r;
    const double dz = segment[1].z - segment[0].z;
    return ((point.r - segment[0].r) * dr + (point.z - segment[0].z) * dz) /
           (dr * dr + dz * dz);
}

bool projection_interval(const Line2InterfaceSideCoordinates& secondary,
                         const Line2InterfaceSideCoordinates& primary,
                         double& lower, double& upper) {
    const RzPoint midpoint = {
        0.5 * (secondary[0].r + secondary[1].r),
        0.5 * (secondary[0].z + secondary[1].z),
    };
    const RzPoint half = {
        0.5 * (secondary[1].r - secondary[0].r),
        0.5 * (secondary[1].z - secondary[0].z),
    };
    const double center = projection_fraction(midpoint, primary);
    const double slope =
        projection_fraction({midpoint.r + half.r, midpoint.z + half.z},
                            primary) -
        center;
    lower = -1.0;
    upper = 1.0;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
    if (std::abs(slope) <= tolerance) {
        return center >= -tolerance && center <= 1.0 + tolerance;
    }
    double first = (0.0 - center) / slope;
    double second = (1.0 - center) / slope;
    if (first > second)
        std::swap(first, second);
    lower = std::max(lower, first);
    upper = std::min(upper, second);
    lower = std::max(-1.0, std::min(1.0, lower));
    upper = std::max(-1.0, std::min(1.0, upper));
    return upper - lower > tolerance;
}

} // namespace

std::vector<std::int64_t>
SpatialAssembly::resolve_block_ids(const SpatialDefinition& definition,
                                 const UnstructuredQuad4Mesh& source_mesh) {
    std::vector<std::int64_t> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions) {
        const std::int64_t id =
            region.block_id >= 0 ? region.block_id
                                 : source_mesh.element_block(region.block).id;
        const bool known = std::any_of(
            source_mesh.element_blocks().begin(),
            source_mesh.element_blocks().end(),
            [id](const ElementBlockInfo& block) { return block.id == id; });
        if (!known)
            throw std::invalid_argument("Unknown element block ID: " +
                                        std::to_string(id));
        if (std::find(result.begin(), result.end(), id) != result.end())
            throw std::invalid_argument("Duplicate resolved region block ID: " +
                                        std::to_string(id));
        result.push_back(id);
    }
    return result;
}

std::vector<RegionMesh>
SpatialAssembly::build_meshes(const SpatialDefinition& definition,
                            const UnstructuredQuad4Mesh& source_mesh) {
    std::vector<RegionMesh> result;
    result.reserve(definition.regions.size());
    for (const RegionDefinition& region : definition.regions) {
        if (region.block_id >= 0)
            result.push_back(RegionMesh::from_unstructured_block(
                source_mesh, region.block_id));
        else
            result.push_back(
                RegionMesh::from_unstructured_block(source_mesh, region.block));
    }
    return result;
}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition,
                             const UnstructuredQuad4Mesh& source_mesh)
    : SpatialAssembly(definition, source_mesh,
                    resolve_block_ids(definition, source_mesh),
                    build_meshes(definition, source_mesh)) {}

SpatialAssembly::SpatialAssembly(SpatialDefinition definition,
                             const UnstructuredQuad4Mesh& source_mesh,
                             std::vector<std::int64_t> block_ids,
                             std::vector<RegionMesh> meshes)
    : _layout(std::move(definition), std::move(block_ids), std::move(meshes)) {
    validate_definitions(_layout._definition);

    _layout.build_volume_geometries();
    build_contacts(source_mesh);
    _boundary.build(source_mesh, _layout);
    _contact._contact_histories.resize(contact_count());
    _contact._projected_mechanical_nodes.resize(contact_count());
    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        _contact._contact_histories[contact_value].resize(
            _contact._secondary_boundaries[contact_value].boundary.nodes.size());
        _contact._projected_mechanical_nodes[contact_value].resize(
            _contact._secondary_boundaries[contact_value].boundary.nodes.size(), false);
    }
    _contact._committed_contact_solution = initial_state();
    update_mechanical_candidates(_contact._committed_contact_solution);
}

std::vector<double> SpatialAssembly::initial_state() const {
    std::vector<double> result(dof_count(), 0.0);
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        const std::size_t offset = _layout._node_offsets[region_value];
        for (std::size_t node = 0; node < _layout._meshes[region_value].nodes().size();
             ++node) {
            result[_layout._dof_map.temperature(offset + node)] =
                _layout._definition.regions[region_value].initial_temperature;
        }
    }
    for (const DirichletCondition& condition :
         _boundary.dirichlet_conditions())
        result[condition.dof] = condition.value;
    return result;
}

std::size_t SpatialAssembly::dof_count() const noexcept {
    return _layout._dof_map.dof_count();
}

std::size_t SpatialAssembly::contribution_count() const noexcept {
    return contribution_ranges().end;
}

const std::vector<DirichletCondition>&
SpatialAssembly::dirichlet_conditions() const noexcept {
    return _boundary.dirichlet_conditions();
}

std::vector<std::size_t>
SpatialAssembly::required_state_dofs(std::size_t contribution_begin,
                                     std::size_t contribution_end) const {
    if (contribution_begin > contribution_end ||
        contribution_end > contribution_count())
        throw std::out_of_range(
            "SpatialAssembly contribution range is out of bounds");
    std::vector<std::size_t> result;
    result.reserve((contribution_end - contribution_begin) * local_dof_count);
    for (std::size_t contribution = contribution_begin;
         contribution < contribution_end; ++contribution) {
        const LocalDofs dofs = contribution_dofs(contribution);
        for (const std::size_t dof : dofs) {
            if (dof >= dof_count())
                throw std::out_of_range(
                    "SpatialAssembly contribution has an invalid DOF");
            result.push_back(dof);
        }
    }
    const std::vector<std::vector<bool>> touched =
        touched_mechanical_nodes(contribution_begin, contribution_end);
    if (!touched.empty()) {
        const ContributionRanges ranges = contribution_ranges();
        for (std::size_t contribution = 0;
             contribution < _contact._mechanical_contributions.size();
             ++contribution) {
            const MechanicalContribution& candidate =
                _contact._mechanical_contributions[contribution];
            if (!touched[candidate.contact][candidate.secondary])
                continue;
            const LocalDofs dofs =
                contribution_dofs(ranges.mechanical_begin + contribution);
            result.insert(result.end(), dofs.begin(), dofs.end());
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

std::pair<std::size_t, std::size_t>
SpatialAssembly::element_location(std::size_t contribution_index) const {
    if (contribution_index >= _layout._element_offsets.back())
        throw std::out_of_range(
            "SpatialAssembly volume contribution is out of range");
    const auto upper = std::upper_bound(
        _layout._element_offsets.begin(), _layout._element_offsets.end(), contribution_index);
    const std::size_t region_value =
        static_cast<std::size_t>(upper - _layout._element_offsets.begin() - 1);
    return {region_value, contribution_index - _layout._element_offsets[region_value]};
}

LocalDofs
SpatialAssembly::contribution_dofs(std::size_t contribution_index) const {
    const ContributionLocation contribution =
        locate_contribution(contribution_index);
    switch (contribution.type) {
    case SpatialContributionType::volume: {
        const auto location = element_location(contribution.local_index);
        const Quad4Element& element =
            _layout._meshes[location.first].elements().at(location.second);
        std::array<std::size_t, 4> nodes{};
        for (std::size_t node = 0; node < nodes.size(); ++node)
            nodes[node] =
                _layout.global_node(location.first, element.nodes[node]);
        return _layout._dof_map.local_dofs(nodes);
    }
    case SpatialContributionType::thermal_contact:
        return _layout._dof_map.local_dofs(
            _contact._thermal_contributions.at(contribution.local_index).nodes);
    case SpatialContributionType::mechanical_contact:
        return _layout._dof_map.local_dofs(
            _contact._mechanical_contributions.at(contribution.local_index).nodes);
    case SpatialContributionType::pressure:
    case SpatialContributionType::traction:
    case SpatialContributionType::convection:
        return _boundary.contribution_dofs(
            contribution.type, contribution.local_index, _layout);
    }
    throw std::logic_error("SpatialAssembly contribution type is invalid");
}

LocalValues SpatialAssembly::contribution_state(
    std::size_t contribution_index,
    const std::vector<double>& global_state) const {
    if (global_state.size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly global state has the wrong size");
    return contribution_state(contribution_index,
                              GlobalStateView(global_state));
}

LocalValues SpatialAssembly::contribution_state(
    std::size_t contribution_index,
    const GlobalStateView& global_state) const {
    if (global_state.global_size() != dof_count())
        throw std::invalid_argument(
            "SpatialAssembly shadow state has the wrong global size");
    const LocalDofs dofs = contribution_dofs(contribution_index);
    LocalValues result{};
    for (std::size_t local = 0; local < dofs.size(); ++local)
        result[local] = global_state.value(dofs[local]);
    return result;
}

LocalResidual
SpatialAssembly::contribution_residual(std::size_t contribution_index,
                                       const LocalValues& state) const {
    const ContributionLocation location =
        locate_contribution(contribution_index);
    switch (location.type) {
    case SpatialContributionType::volume:
        throw std::logic_error(
            "SpatialAssembly does not own volume residual physics");
    case SpatialContributionType::thermal_contact: {
        const ThermalContribution& contribution =
            _contact._thermal_contributions[location.local_index];
        return _contact._thermal_kernels[contribution.contact].residual(
            contribution.geometry, state);
    }
    case SpatialContributionType::mechanical_contact: {
        const MechanicalContribution& contribution =
            _contact._mechanical_contributions.at(location.local_index);
        if (!contribution.active)
            return {};
        return _contact._mechanical_kernels[contribution.contact].residual(
            contribution.geometry, state,
            contribution_state(contribution_index,
                               _contact._committed_contact_solution),
            _contact._contact_histories[contribution.contact][contribution.secondary]);
    }
    case SpatialContributionType::pressure:
    case SpatialContributionType::traction:
    case SpatialContributionType::convection:
        return _boundary.contribution_residual(
            location.type, location.local_index, state);
    }
    throw std::logic_error("SpatialAssembly contribution type is invalid");
}

LocalSystem
SpatialAssembly::linearize_contribution(std::size_t contribution_index,
                                        const LocalValues& state) const {
    const ContributionLocation location =
        locate_contribution(contribution_index);
    switch (location.type) {
    case SpatialContributionType::volume:
        throw std::logic_error(
            "SpatialAssembly does not own volume Jacobian physics");
    case SpatialContributionType::thermal_contact: {
        const ThermalContribution& contribution =
            _contact._thermal_contributions[location.local_index];
        return _contact._thermal_kernels[contribution.contact].linearize(
            contribution.geometry, state);
    }
    case SpatialContributionType::mechanical_contact: {
        const MechanicalContribution& contribution =
            _contact._mechanical_contributions.at(location.local_index);
        if (!contribution.active)
            return {};
        return _contact._mechanical_kernels[contribution.contact].linearize(
            contribution.geometry, state,
            contribution_state(contribution_index,
                               _contact._committed_contact_solution),
            _contact._contact_histories[contribution.contact][contribution.secondary]);
    }
    case SpatialContributionType::pressure:
    case SpatialContributionType::traction:
    case SpatialContributionType::convection:
        return _boundary.linearize_contribution(
            location.type, location.local_index, state);
    }
    throw std::logic_error("SpatialAssembly contribution type is invalid");
}

void SpatialAssembly::build_contacts(const UnstructuredQuad4Mesh& source_mesh) {
    _contact._primary_boundaries.reserve(contact_count());
    _contact._secondary_boundaries.reserve(contact_count());
    _contact._thermal_kernels.reserve(contact_count());
    _contact._mechanical_kernels.reserve(contact_count());

    for (std::size_t contact_value = 0; contact_value < contact_count();
         ++contact_value) {
        ContactDefinition& contact_definition =
            _layout._definition.contacts[contact_value];
        ResolvedBoundary primary =
            _layout.resolve_boundary(source_mesh, contact_definition.primary);
        ResolvedBoundary secondary =
            _layout.resolve_boundary(source_mesh,
                                     contact_definition.secondary);
        if (primary.region == secondary.region)
            throw std::invalid_argument("Self-contact is not supported: " +
                                        contact_definition.name);
        const RegionMesh& primary_mesh = _layout._meshes[primary.region];
        const RegionMesh& secondary_mesh = _layout._meshes[secondary.region];
        primary.boundary = ordered_connected_boundary(
            primary_mesh, std::move(primary.boundary),
            contact_definition.primary);
        secondary.boundary = ordered_connected_boundary(
            secondary_mesh, std::move(secondary.boundary),
            contact_definition.secondary);

        // Parent-element centroids give the material side of each boundary
        // edge; they form the zero-gap orientation hint passed to every
        // thermal and mechanical contact geometry below.
        const std::vector<RzPoint> primary_parent_centroids =
            boundary_parent_centroids(_layout, primary.region,
                                      primary.boundary);
        const std::vector<RzPoint> secondary_parent_centroids =
            boundary_parent_centroids(_layout, secondary.region,
                                      secondary.boundary);

        if (contact_definition.mechanical &&
            contact_definition.automatic_penalty) {
            const double primary_length =
                minimum_boundary_normal_length(_layout, primary.region,
                                               primary.boundary);
            const double secondary_length = minimum_boundary_normal_length(
                _layout, secondary.region, secondary.boundary);
            const double primary_modulus =
                _layout._definition.regions[primary.region].material.young_modulus;
            const double secondary_modulus =
                _layout._definition.regions[secondary.region].material.young_modulus;
            contact_definition.penalty = contact_definition.penalty_factor /
                                         (primary_length / primary_modulus +
                                          secondary_length / secondary_modulus);
            if (!std::isfinite(contact_definition.penalty) ||
                !(contact_definition.penalty > 0.0))
                throw std::overflow_error(
                    "Automatic contact penalty is not finite and positive: " +
                    contact_definition.name);
        }

        _contact._thermal_kernels.emplace_back(GapHeatProperties{
            contact_definition.thermal ? contact_definition.gap_conductivity
                                       : 1.0,
            contact_definition.thermal ? contact_definition.minimum_gap : 1.0});
        _contact._mechanical_kernels.emplace_back(NormalContactProperties{
            contact_definition.mechanical ? contact_definition.penalty : 1.0,
            contact_definition.mechanical
                ? contact_definition.friction_coefficient
                : 0.0,
            contact_definition.mechanical &&
                contact_definition.mechanical_formulation ==
                    MechanicalContactFormulation::augmented_lagrangian});

        if (contact_definition.thermal) {
            for (std::size_t edge_index = 0;
                 edge_index < secondary.boundary.elements.size();
                 ++edge_index) {
                const Line2BoundaryElement& secondary_edge =
                    secondary.boundary.elements[edge_index];
                const Line2InterfaceSideCoordinates secondary_coordinates =
                    edge_coordinates(secondary_mesh, secondary_edge);
                struct ThermalProjection final {
                    std::size_t primary_edge;
                    double lower;
                    double upper;
                };
                std::vector<ThermalProjection> projections;
                for (std::size_t primary_edge = 0;
                     primary_edge < primary.boundary.elements.size();
                     ++primary_edge) {
                    double lower = 0.0;
                    double upper = 0.0;
                    if (projection_interval(
                            secondary_coordinates,
                            edge_coordinates(
                                primary_mesh,
                                primary.boundary.elements[primary_edge]),
                            lower, upper))
                        projections.push_back({primary_edge, lower, upper});
                }
                std::sort(projections.begin(), projections.end(),
                          [](const ThermalProjection& lhs,
                             const ThermalProjection& rhs) {
                              return lhs.lower < rhs.lower;
                          });
                double covered = -1.0;
                constexpr double coverage_tolerance = 1.0e-10;
                for (const ThermalProjection& projection : projections) {
                    if (projection.lower > covered + coverage_tolerance ||
                        projection.lower < covered - coverage_tolerance)
                        throw std::invalid_argument(
                            "Secondary thermal edge projection has a gap or "
                            "overlap on the primary surface: " +
                            contact_definition.name);
                    covered = projection.upper;
                }
                if (covered < 1.0 - coverage_tolerance)
                    throw std::invalid_argument(
                        "Secondary thermal edge is outside the primary "
                        "surface projection: " +
                        contact_definition.name);
                for (const ThermalProjection& projection : projections) {
                    const Line2BoundaryElement& primary_edge =
                        primary.boundary.elements[projection.primary_edge];
                    const Line2InterfaceSideCoordinates primary_coordinates =
                        edge_coordinates(primary_mesh, primary_edge);
                    _contact._thermal_contributions.push_back(
                        {contact_value,
                         interface_nodes(_layout, secondary.region,
                                         secondary_edge, primary.region,
                                         primary_edge),
                         make_line2_rz_heat_geometry(
                             secondary_coordinates, primary_coordinates,
                             projection.lower, projection.upper,
                             zero_gap_orientation_hint(
                                 secondary_parent_centroids[edge_index],
                                 primary_parent_centroids[projection
                                                              .primary_edge],
                                 primary_coordinates))});
                }
            }
        }

        if (contact_definition.mechanical) {
            for (std::size_t edge_index = 0;
                 edge_index < secondary.boundary.elements.size();
                 ++edge_index) {
                const Line2BoundaryElement& secondary_edge =
                    secondary.boundary.elements[edge_index];
                const Line2InterfaceSideCoordinates secondary_coordinates =
                    edge_coordinates(secondary_mesh, secondary_edge);
                for (std::size_t secondary_node = 0;
                     secondary_node < line2_interface_side_node_count;
                     ++secondary_node) {
                    for (std::size_t candidate = 0;
                         candidate < primary.boundary.elements.size();
                         ++candidate) {
                        const Line2BoundaryElement& primary_edge =
                            primary.boundary.elements[candidate];
                        const Line2InterfaceSideCoordinates
                            primary_coordinates =
                                edge_coordinates(primary_mesh, primary_edge);
                        _contact._mechanical_contributions.push_back(
                            {contact_value,
                             interface_nodes(_layout, secondary.region,
                                             secondary_edge, primary.region,
                                             primary_edge),
                             make_node_to_line_rz_contact_geometry(
                                 secondary_coordinates, primary_coordinates,
                                 secondary_node, candidate == 0,
                                 candidate + 1 ==
                                     primary.boundary.elements.size(),
                                 zero_gap_orientation_hint(
                                     secondary_parent_centroids[edge_index],
                                     primary_parent_centroids[candidate],
                                     primary_coordinates)),
                             edge_index + secondary_node,
                             candidate});
                    }
                }
            }
        }
        _contact._primary_boundaries.push_back(std::move(primary));
        _contact._secondary_boundaries.push_back(std::move(secondary));
    }
}

} // namespace fuelsim
