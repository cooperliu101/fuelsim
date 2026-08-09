#include "spatial_assembly.hpp"

#include <algorithm>
#include <stdexcept>

namespace fuelsim {

const SpatialDefinition& SpatialAssembly::definition() const noexcept {
    return _definition;
}

const DofMap& SpatialAssembly::dof_map() const noexcept {
    return _dof_map;
}

std::size_t SpatialAssembly::region_count() const noexcept {
    return _definition.regions.size();
}

std::size_t SpatialAssembly::region_index(const std::string& name) const {
    for (std::size_t region = 0; region < region_count(); ++region) {
        if (_definition.regions[region].name == name)
            return region;
    }
    throw std::invalid_argument("Unknown region: " + name);
}

const RegionDefinition& SpatialAssembly::region(std::size_t index) const {
    return _definition.regions.at(index);
}

const RegionMesh& SpatialAssembly::region_mesh(std::size_t index) const {
    return _meshes.at(index);
}

std::size_t SpatialAssembly::region_node_offset(std::size_t index) const {
    if (index >= region_count())
        throw std::out_of_range("SpatialAssembly region index is out of range");
    return _node_offsets[index];
}

std::size_t SpatialAssembly::region_element_count(std::size_t index) const {
    return _meshes.at(index).elements().size();
}

std::size_t SpatialAssembly::region_element_offset(std::size_t index) const {
    if (index >= region_count())
        throw std::out_of_range("SpatialAssembly region index is out of range");
    return _element_offsets[index];
}

std::size_t SpatialAssembly::volume_contribution_count() const noexcept {
    return _element_offsets.back();
}

SpatialAssembly::ContributionRanges
SpatialAssembly::contribution_ranges() const noexcept {
    const std::size_t thermal_begin = volume_contribution_count();
    const std::size_t mechanical_begin =
        thermal_begin + _thermal_contributions.size();
    const std::size_t pressure_begin =
        mechanical_begin + _mechanical_contributions.size();
    const std::size_t traction_begin =
        pressure_begin + _pressure_contributions.size();
    const std::size_t convection_begin =
        traction_begin + _traction_contributions.size();
    return {thermal_begin,
            mechanical_begin,
            pressure_begin,
            traction_begin,
            convection_begin,
            convection_begin + _convection_contributions.size()};
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
    return _region_geometries.at(region_value).at(element_index);
}

double SpatialAssembly::region_heat_source(std::size_t region_value) const {
    const RegionDefinition& value = _definition.regions.at(region_value);
    const double multiplier = value.heat_source_function.empty()
                                  ? _load_factor
                                  : function_value(value.heat_source_function);
    return multiplier * value.volumetric_heat_source;
}


} // namespace fuelsim
