#include "spatial_assembly.hpp"

#include <algorithm>
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
    return _layout._definition;
}

const DofMap& SpatialAssembly::dof_map() const noexcept {
    return _layout._dof_map;
}

std::size_t SpatialAssembly::region_count() const noexcept {
    return _layout._definition.regions.size();
}

std::size_t SpatialAssembly::region_index(const std::string& name) const {
    for (std::size_t region = 0; region < region_count(); ++region) {
        if (_layout._definition.regions[region].name == name)
            return region;
    }
    throw std::invalid_argument("Unknown region: " + name);
}

const RegionDefinition& SpatialAssembly::region(std::size_t index) const {
    return _layout._definition.regions.at(index);
}

const RegionMesh& SpatialAssembly::region_mesh(std::size_t index) const {
    return _layout._meshes.at(index);
}

std::size_t SpatialAssembly::region_node_offset(std::size_t index) const {
    if (index >= region_count())
        throw std::out_of_range("SpatialAssembly region index is out of range");
    return _layout._node_offsets[index];
}

std::size_t SpatialAssembly::region_element_count(std::size_t index) const {
    return _layout._meshes.at(index).elements().size();
}

std::size_t SpatialAssembly::region_element_offset(std::size_t index) const {
    if (index >= region_count())
        throw std::out_of_range("SpatialAssembly region index is out of range");
    return _layout._element_offsets[index];
}

std::size_t SpatialAssembly::volume_contribution_count() const noexcept {
    return _layout._element_offsets.back();
}

SpatialAssembly::ContributionRanges
SpatialAssembly::contribution_ranges() const noexcept {
    const std::size_t thermal_begin = volume_contribution_count();
    const std::size_t mechanical_begin =
        thermal_begin + _contact._thermal_contributions.size();
    const std::size_t pressure_begin =
        mechanical_begin + _contact._mechanical_contributions.size();
    const std::size_t traction_begin =
        pressure_begin + _boundary._pressure_contributions.size();
    const std::size_t convection_begin =
        traction_begin + _boundary._traction_contributions.size();
    return {thermal_begin,
            mechanical_begin,
            pressure_begin,
            traction_begin,
            convection_begin,
            convection_begin + _boundary._convection_contributions.size()};
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
    return _layout._region_geometries.at(region_value).at(element_index);
}

double SpatialAssembly::region_heat_source(std::size_t region_value) const {
    return _boundary.region_heat_source(region_value, _layout);
}


} // namespace fuelsim
