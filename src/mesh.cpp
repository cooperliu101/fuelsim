#include "fuelsim/mesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

bool same_coordinate(double lhs, double rhs) {
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= 1.0e-12 * scale;
}

std::size_t coordinate_index(const std::vector<double>& coordinates,
                             double value) {
    for (std::size_t index = 0; index < coordinates.size(); ++index) {
        if (same_coordinate(coordinates[index], value))
            return index;
    }
    throw std::invalid_argument(
        "StructuredRzMesh node does not lie on the reconstructed grid");
}

std::vector<double>
sorted_unique_coordinates(const UnstructuredQuad4Mesh& source,
                          const std::vector<bool>& used_nodes, bool radial) {
    std::vector<double> coordinates;
    for (std::size_t node = 0; node < source.nodes().size(); ++node) {
        if (used_nodes[node])
            coordinates.push_back(radial ? source.nodes()[node].r
                                         : source.nodes()[node].z);
    }
    std::sort(coordinates.begin(), coordinates.end());
    coordinates.erase(
        std::unique(coordinates.begin(), coordinates.end(), same_coordinate),
        coordinates.end());
    return coordinates;
}

} // namespace

UnstructuredQuad4Mesh::UnstructuredQuad4Mesh(
    std::vector<RzPoint> nodes, std::vector<Quad4Element> elements,
    std::vector<std::int64_t> element_block_ids,
    std::vector<ElementBlockInfo> element_blocks,
    std::vector<NodeSet> node_sets, std::vector<SideSet> side_sets)
    : _nodes(std::move(nodes)), _elements(std::move(elements)),
      _element_block_ids(std::move(element_block_ids)),
      _element_blocks(std::move(element_blocks)),
      _node_sets(std::move(node_sets)), _side_sets(std::move(side_sets)) {
    if (_nodes.empty())
        throw std::invalid_argument(
            "UnstructuredQuad4Mesh requires at least one node");
    if (_elements.empty())
        throw std::invalid_argument(
            "UnstructuredQuad4Mesh requires at least one element");
    if (_elements.size() != _element_block_ids.size())
        throw std::invalid_argument(
            "UnstructuredQuad4Mesh block ID count must match element count");

    for (const RzPoint& node : _nodes) {
        if (!std::isfinite(node.r) || !std::isfinite(node.z) || node.r < 0.0)
            throw std::invalid_argument(
                "UnstructuredQuad4Mesh requires finite RZ coordinates and "
                "nonnegative radius");
    }
    if (_element_blocks.empty())
        throw std::invalid_argument(
            "UnstructuredQuad4Mesh requires element block metadata");
    for (std::size_t block = 0; block < _element_blocks.size(); ++block) {
        if (_element_blocks[block].id < 0)
            throw std::invalid_argument(
                "UnstructuredQuad4Mesh block IDs must be nonnegative");
        for (std::size_t previous = 0; previous < block; ++previous) {
            if (_element_blocks[previous].id == _element_blocks[block].id)
                throw std::invalid_argument(
                    "UnstructuredQuad4Mesh block IDs must be unique");
            if (!_element_blocks[block].name.empty() &&
                _element_blocks[previous].name == _element_blocks[block].name)
                throw std::invalid_argument(
                    "UnstructuredQuad4Mesh block names must be unique");
        }
    }
    for (std::size_t index = 0; index < _elements.size(); ++index) {
        const bool known_block = std::any_of(
            _element_blocks.begin(), _element_blocks.end(),
            [id = _element_block_ids[index]](const ElementBlockInfo& block) {
                return block.id == id;
            });
        if (!known_block)
            throw std::invalid_argument(
                "UnstructuredQuad4Mesh element references an unknown block");
        for (const std::size_t node : _elements[index].nodes) {
            if (node >= _nodes.size())
                throw std::out_of_range(
                    "UnstructuredQuad4Mesh connectivity is out of range");
        }
    }

    for (std::size_t set = 0; set < _node_sets.size(); ++set) {
        const NodeSet& node_set = _node_sets[set];
        if (node_set.id < 0)
            throw std::invalid_argument(
                "UnstructuredQuad4Mesh node set IDs must be nonnegative");
        for (std::size_t previous_index = 0; previous_index < set;
             ++previous_index) {
            const NodeSet& previous = _node_sets[previous_index];
            if (previous.id == node_set.id ||
                (!node_set.name.empty() && previous.name == node_set.name))
                throw std::invalid_argument(
                    "UnstructuredQuad4Mesh node set IDs and names must be "
                    "unique");
        }
        for (const std::size_t node : node_set.nodes) {
            if (node >= _nodes.size())
                throw std::out_of_range(
                    "UnstructuredQuad4Mesh node set is out of range");
        }
    }
    for (std::size_t set = 0; set < _side_sets.size(); ++set) {
        const SideSet& side_set = _side_sets[set];
        if (side_set.id < 0)
            throw std::invalid_argument(
                "UnstructuredQuad4Mesh side set IDs must be nonnegative");
        for (std::size_t previous_index = 0; previous_index < set;
             ++previous_index) {
            const SideSet& previous = _side_sets[previous_index];
            if (previous.id == side_set.id ||
                (!side_set.name.empty() && previous.name == side_set.name))
                throw std::invalid_argument(
                    "UnstructuredQuad4Mesh side set IDs and names must be "
                    "unique");
        }
        for (const ElementSide& side : side_set.sides) {
            if (side.element >= _elements.size() || side.local_side >= 4U)
                throw std::out_of_range(
                    "UnstructuredQuad4Mesh side set is out of range");
        }
    }
}

const std::vector<RzPoint>& UnstructuredQuad4Mesh::nodes() const noexcept {
    return _nodes;
}

const std::vector<Quad4Element>&
UnstructuredQuad4Mesh::elements() const noexcept {
    return _elements;
}

const std::vector<std::int64_t>&
UnstructuredQuad4Mesh::element_block_ids() const noexcept {
    return _element_block_ids;
}

const std::vector<ElementBlockInfo>&
UnstructuredQuad4Mesh::element_blocks() const noexcept {
    return _element_blocks;
}

const std::vector<NodeSet>& UnstructuredQuad4Mesh::node_sets() const noexcept {
    return _node_sets;
}

const std::vector<SideSet>& UnstructuredQuad4Mesh::side_sets() const noexcept {
    return _side_sets;
}

const ElementBlockInfo&
UnstructuredQuad4Mesh::element_block(const std::string& name) const {
    const auto block =
        std::find_if(_element_blocks.begin(), _element_blocks.end(),
                     [&name](const ElementBlockInfo& candidate) {
                         return candidate.name == name;
                     });
    if (block == _element_blocks.end())
        throw std::invalid_argument("Unknown element block: " + name);
    return *block;
}

const NodeSet& UnstructuredQuad4Mesh::node_set(const std::string& name) const {
    const auto set = std::find_if(
        _node_sets.begin(), _node_sets.end(),
        [&name](const NodeSet& candidate) { return candidate.name == name; });
    if (set == _node_sets.end())
        throw std::invalid_argument("Unknown node set: " + name);
    return *set;
}

const SideSet& UnstructuredQuad4Mesh::side_set(const std::string& name) const {
    const auto set = std::find_if(
        _side_sets.begin(), _side_sets.end(),
        [&name](const SideSet& candidate) { return candidate.name == name; });
    if (set == _side_sets.end())
        throw std::invalid_argument("Unknown side set: " + name);
    return *set;
}

std::int64_t
UnstructuredQuad4Mesh::side_set_block_id(const std::string& name) const {
    const SideSet& set = side_set(name);
    if (set.sides.empty())
        throw std::invalid_argument("Side set is empty: " + name);
    const std::int64_t block_id =
        _element_block_ids.at(set.sides.front().element);
    for (const ElementSide& side : set.sides) {
        if (_element_block_ids.at(side.element) != block_id)
            throw std::invalid_argument("Side set crosses element blocks: " +
                                        name);
    }
    return block_id;
}

RegionMesh
RegionMesh::from_unstructured_block(const UnstructuredQuad4Mesh& source,
                                    const std::string& block_name) {
    return from_unstructured_block(source, source.element_block(block_name).id);
}

RegionMesh
RegionMesh::from_unstructured_block(const UnstructuredQuad4Mesh& source,
                                    std::int64_t block_id) {
    const bool known_block = std::any_of(
        source.element_blocks().begin(), source.element_blocks().end(),
        [block_id](const ElementBlockInfo& block) {
            return block.id == block_id;
        });
    if (!known_block)
        throw std::invalid_argument("Unknown element block ID: " +
                                    std::to_string(block_id));

    const std::size_t invalid = std::numeric_limits<std::size_t>::max();
    RegionMesh mesh;
    mesh._block_id = block_id;
    mesh._source_node_to_local.assign(source.nodes().size(), invalid);
    mesh._source_element_to_local.assign(source.elements().size(), invalid);

    std::vector<bool> used_nodes(source.nodes().size(), false);
    for (std::size_t source_element = 0;
         source_element < source.elements().size(); ++source_element) {
        if (source.element_block_ids()[source_element] != block_id)
            continue;
        mesh._source_element_to_local[source_element] =
            mesh._source_element_ids.size();
        mesh._source_element_ids.push_back(source_element);
        for (std::size_t node : source.elements()[source_element].nodes)
            used_nodes[node] = true;
    }
    if (mesh._source_element_ids.empty())
        throw std::invalid_argument("Element block is empty: " +
                                    std::to_string(block_id));

    for (std::size_t source_node = 0; source_node < source.nodes().size();
         ++source_node) {
        if (!used_nodes[source_node])
            continue;
        mesh._source_node_to_local[source_node] = mesh._nodes.size();
        mesh._source_node_ids.push_back(source_node);
        mesh._nodes.push_back(source.nodes()[source_node]);
    }

    mesh._elements.reserve(mesh._source_element_ids.size());
    for (std::size_t source_element : mesh._source_element_ids) {
        Quad4Element element{};
        for (std::size_t node = 0; node < element.nodes.size(); ++node) {
            const std::size_t local = mesh._source_node_to_local.at(
                source.elements()[source_element].nodes[node]);
            if (local == invalid)
                throw std::logic_error(
                    "RegionMesh connectivity crosses element blocks");
            element.nodes[node] = local;
        }
        mesh._elements.push_back(element);
    }
    return mesh;
}

const std::vector<RzPoint>& RegionMesh::nodes() const noexcept {
    return _nodes;
}

const std::vector<Quad4Element>& RegionMesh::elements() const noexcept {
    return _elements;
}

const std::vector<std::size_t>& RegionMesh::source_node_ids() const noexcept {
    return _source_node_ids;
}

const std::vector<std::size_t>&
RegionMesh::source_element_ids() const noexcept {
    return _source_element_ids;
}

std::int64_t RegionMesh::block_id() const noexcept {
    return _block_id;
}

RegionBoundary
RegionMesh::map_side_set(const UnstructuredQuad4Mesh& source,
                         const std::string& side_set_name) const {
    if (source.side_set_block_id(side_set_name) != _block_id)
        throw std::invalid_argument(
            "Side set belongs to an unexpected block: " + side_set_name);

    const std::size_t invalid = std::numeric_limits<std::size_t>::max();
    std::vector<Line2BoundaryElement> elements;
    std::vector<RzPoint> adjacent_centroids;
    const SideSet& side_set = source.side_set(side_set_name);
    elements.reserve(side_set.sides.size());
    adjacent_centroids.reserve(side_set.sides.size());
    for (const ElementSide& side : side_set.sides) {
        if (_source_element_to_local.at(side.element) == invalid)
            throw std::invalid_argument("Side set is outside its region: " +
                                        side_set_name);
        const Quad4Element& source_element = source.elements().at(side.element);
        const std::size_t first_source =
            source_element.nodes.at(side.local_side);
        const std::size_t second_source =
            source_element.nodes.at((side.local_side + 1U) % 4U);
        const std::size_t first = _source_node_to_local.at(first_source);
        const std::size_t second = _source_node_to_local.at(second_source);
        if (first == invalid || second == invalid)
            throw std::logic_error("RegionMesh side-set node mapping failed");
        elements.push_back({{{first, second}}});
        RzPoint centroid{0.0, 0.0};
        for (std::size_t node : source_element.nodes) {
            centroid.r += source.nodes().at(node).r / 4.0;
            centroid.z += source.nodes().at(node).z / 4.0;
        }
        adjacent_centroids.push_back(centroid);
    }

    const RzPoint& reference_point = _nodes.at(elements.front().nodes[0]);
    const auto all_on = [&](bool radial) {
        const double coordinate =
            radial ? reference_point.r : reference_point.z;
        return std::all_of(elements.begin(), elements.end(),
                           [&](const Line2BoundaryElement& edge) {
                               return std::all_of(
                                   edge.nodes.begin(), edge.nodes.end(),
                                   [&](std::size_t node) {
                                       const RzPoint& point = _nodes.at(node);
                                       return same_coordinate(radial ? point.r
                                                                     : point.z,
                                                              coordinate);
                                   });
                           });
    };
    const auto material_side = [&](bool radial) {
        const double coordinate =
            radial ? reference_point.r : reference_point.z;
        int side = 0;
        for (const RzPoint& centroid : adjacent_centroids) {
            const double value = radial ? centroid.r : centroid.z;
            if (same_coordinate(value, coordinate))
                return 0;
            const int current = value > coordinate ? 1 : -1;
            if (side != 0 && current != side)
                return 0;
            side = current;
        }
        return side;
    };

    RegionBoundaryKind kind = RegionBoundaryKind::general;
    bool sort_by_axial = false;
    const int radial_material_side = all_on(true) ? material_side(true) : 0;
    const int axial_material_side = all_on(false) ? material_side(false) : 0;
    if (radial_material_side > 0) {
        kind = RegionBoundaryKind::radial_inner;
        sort_by_axial = true;
    } else if (radial_material_side < 0) {
        kind = RegionBoundaryKind::radial_outer;
        sort_by_axial = true;
    } else if (axial_material_side > 0) {
        kind = RegionBoundaryKind::bottom;
    } else if (axial_material_side < 0) {
        kind = RegionBoundaryKind::top;
    }

    if (kind != RegionBoundaryKind::general) {
        for (Line2BoundaryElement& edge : elements) {
            const RzPoint& first = _nodes.at(edge.nodes[0]);
            const RzPoint& second = _nodes.at(edge.nodes[1]);
            if ((sort_by_axial && first.z > second.z) ||
                (!sort_by_axial && first.r > second.r))
                std::swap(edge.nodes[0], edge.nodes[1]);
        }
        std::sort(elements.begin(), elements.end(),
                  [&](const Line2BoundaryElement& lhs,
                      const Line2BoundaryElement& rhs) {
                      const RzPoint& lhs_point = _nodes.at(lhs.nodes[0]);
                      const RzPoint& rhs_point = _nodes.at(rhs.nodes[0]);
                      return sort_by_axial ? lhs_point.z < rhs_point.z
                                           : lhs_point.r < rhs_point.r;
                  });
    }

    std::vector<std::size_t> nodes;
    for (const Line2BoundaryElement& edge : elements)
        nodes.insert(nodes.end(), edge.nodes.begin(), edge.nodes.end());
    std::sort(nodes.begin(), nodes.end(),
              [&](std::size_t lhs, std::size_t rhs) {
                  if (kind == RegionBoundaryKind::radial_inner ||
                      kind == RegionBoundaryKind::radial_outer)
                      return _nodes[lhs].z < _nodes[rhs].z;
                  if (kind == RegionBoundaryKind::bottom ||
                      kind == RegionBoundaryKind::top)
                      return _nodes[lhs].r < _nodes[rhs].r;
                  return _source_node_ids[lhs] < _source_node_ids[rhs];
              });
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    return {kind, std::move(nodes), std::move(elements)};
}

StructuredRzMesh
StructuredRzMesh::from_unstructured_block(const UnstructuredQuad4Mesh& source,
                                          const std::string& block_name) {
    return from_unstructured_block(source, source.element_block(block_name).id,
                                   {});
}

StructuredRzMesh StructuredRzMesh::from_unstructured_block(
    const UnstructuredQuad4Mesh& source, const std::string& block_name,
    const RzBoundaryNames& boundary_names) {
    return from_unstructured_block(source, source.element_block(block_name).id,
                                   boundary_names);
}

StructuredRzMesh StructuredRzMesh::from_unstructured_block(
    const UnstructuredQuad4Mesh& source, std::int64_t block_id,
    const RzBoundaryNames& boundary_names) {
    const bool known_block = std::any_of(
        source.element_blocks().begin(), source.element_blocks().end(),
        [block_id](const ElementBlockInfo& block) {
            return block.id == block_id;
        });
    if (!known_block)
        throw std::invalid_argument("Unknown element block ID: " +
                                    std::to_string(block_id));
    std::vector<bool> used_nodes(source.nodes().size(), false);
    std::size_t block_element_count = 0;
    for (std::size_t element = 0; element < source.elements().size();
         ++element) {
        if (source.element_block_ids()[element] != block_id)
            continue;
        ++block_element_count;
        for (const std::size_t node : source.elements()[element].nodes)
            used_nodes[node] = true;
    }
    if (block_element_count == 0)
        throw std::invalid_argument("Element block is empty: " +
                                    std::to_string(block_id));

    const std::vector<double> radial_coordinates =
        sorted_unique_coordinates(source, used_nodes, true);
    const std::vector<double> axial_coordinates =
        sorted_unique_coordinates(source, used_nodes, false);
    if (radial_coordinates.size() < 2 || axial_coordinates.size() < 2)
        throw std::invalid_argument(
            "StructuredRzMesh block requires at least one element per axis");

    const std::size_t radial_nodes = radial_coordinates.size();
    const std::size_t axial_nodes = axial_coordinates.size();
    if (radial_nodes > std::numeric_limits<std::size_t>::max() / axial_nodes)
        throw std::length_error("StructuredRzMesh node count overflows");
    const std::size_t expected_node_count = radial_nodes * axial_nodes;
    const std::size_t radial_elements = radial_nodes - 1;
    const std::size_t axial_elements = axial_nodes - 1;
    if (radial_elements >
        std::numeric_limits<std::size_t>::max() / axial_elements)
        throw std::length_error("StructuredRzMesh element count overflows");
    if (block_element_count != radial_elements * axial_elements)
        throw std::invalid_argument(
            "Element block is not a complete structured Quad4 grid");

    StructuredRzMesh mesh;
    mesh._inner_radius = radial_coordinates.front();
    mesh._outer_radius = radial_coordinates.back();
    mesh._length = axial_coordinates.back() - axial_coordinates.front();
    mesh._radial_elements = radial_elements;
    mesh._axial_elements = axial_elements;
    mesh._nodes.resize(expected_node_count);

    const std::size_t invalid_node = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> source_to_local(source.nodes().size(),
                                             invalid_node);
    std::vector<bool> local_node_found(expected_node_count, false);
    for (std::size_t source_node = 0; source_node < source.nodes().size();
         ++source_node) {
        if (!used_nodes[source_node])
            continue;
        const RzPoint& point = source.nodes()[source_node];
        const std::size_t radial =
            coordinate_index(radial_coordinates, point.r);
        const std::size_t axial = coordinate_index(axial_coordinates, point.z);
        const std::size_t local = axial * radial_nodes + radial;
        if (local_node_found[local])
            throw std::invalid_argument(
                "Element block contains duplicate grid coordinates");
        local_node_found[local] = true;
        source_to_local[source_node] = local;
        mesh._nodes[local] = point;
    }
    if (std::find(local_node_found.begin(), local_node_found.end(), false) !=
        local_node_found.end())
        throw std::invalid_argument(
            "Element block is missing a structured grid node");

    std::vector<bool> cell_found(radial_elements * axial_elements, false);
    mesh._elements.reserve(block_element_count);
    for (std::size_t source_element = 0;
         source_element < source.elements().size(); ++source_element) {
        if (source.element_block_ids()[source_element] != block_id)
            continue;
        Quad4Element element{};
        std::size_t minimum_radial = radial_elements;
        std::size_t maximum_radial = 0;
        std::size_t minimum_axial = axial_elements;
        std::size_t maximum_axial = 0;
        for (std::size_t local_node = 0; local_node < 4U; ++local_node) {
            const std::size_t mapped = source_to_local.at(
                source.elements()[source_element].nodes[local_node]);
            if (mapped == invalid_node)
                throw std::invalid_argument(
                    "Element block connectivity crosses block boundaries");
            element.nodes[local_node] = mapped;
            const std::size_t radial = mapped % radial_nodes;
            const std::size_t axial = mapped / radial_nodes;
            minimum_radial = std::min(minimum_radial, radial);
            maximum_radial = std::max(maximum_radial, radial);
            minimum_axial = std::min(minimum_axial, axial);
            maximum_axial = std::max(maximum_axial, axial);
        }
        if (maximum_radial != minimum_radial + 1 ||
            maximum_axial != minimum_axial + 1)
            throw std::invalid_argument(
                "Element block contains a non-grid-aligned Quad4");
        const std::size_t cell =
            minimum_axial * radial_elements + minimum_radial;
        if (cell_found[cell])
            throw std::invalid_argument(
                "Element block contains duplicate structured cells");
        cell_found[cell] = true;
        mesh._elements.push_back(element);
    }

    const auto map_boundary = [&](const std::string& name, bool sort_by_axial,
                                  double fixed_coordinate, bool check_radial) {
        std::vector<std::size_t> local_nodes;
        std::vector<std::size_t> source_nodes;
        if (name.empty()) {
            for (std::size_t source_node = 0;
                 source_node < source.nodes().size(); ++source_node) {
                if (!used_nodes[source_node])
                    continue;
                const RzPoint& point = source.nodes()[source_node];
                const double coordinate = check_radial ? point.r : point.z;
                if (same_coordinate(coordinate, fixed_coordinate))
                    source_nodes.push_back(source_node);
            }
        } else {
            source_nodes = source.node_set(name).nodes;
        }
        local_nodes.reserve(source_nodes.size());
        for (const std::size_t source_node : source_nodes) {
            const std::size_t local = source_to_local.at(source_node);
            if (local == invalid_node)
                throw std::invalid_argument(
                    "Node set crosses element blocks: " + name);
            const RzPoint& point = source.nodes()[source_node];
            const double coordinate = check_radial ? point.r : point.z;
            if (!same_coordinate(coordinate, fixed_coordinate))
                throw std::invalid_argument(
                    "Node set is not on the expected RZ boundary: " + name);
            local_nodes.push_back(local);
        }
        std::sort(local_nodes.begin(), local_nodes.end(),
                  [&](std::size_t lhs, std::size_t rhs) {
                      const RzPoint& lhs_point = mesh._nodes[lhs];
                      const RzPoint& rhs_point = mesh._nodes[rhs];
                      return sort_by_axial ? lhs_point.z < rhs_point.z
                                           : lhs_point.r < rhs_point.r;
                  });
        if (std::adjacent_find(local_nodes.begin(), local_nodes.end()) !=
            local_nodes.end())
            throw std::invalid_argument("Node set contains duplicates: " +
                                        name);
        return local_nodes;
    };

    mesh._radial_inner_nodes = map_boundary(boundary_names.radial_inner, true,
                                            mesh._inner_radius, true);
    mesh._radial_outer_nodes = map_boundary(boundary_names.radial_outer, true,
                                            mesh._outer_radius, true);
    mesh._bottom_nodes = map_boundary(boundary_names.bottom, false,
                                      axial_coordinates.front(), false);
    mesh._top_nodes = map_boundary(boundary_names.top, false,
                                   axial_coordinates.back(), false);
    if (mesh._radial_inner_nodes.size() != axial_nodes ||
        mesh._radial_outer_nodes.size() != axial_nodes ||
        mesh._bottom_nodes.size() != radial_nodes ||
        mesh._top_nodes.size() != radial_nodes)
        throw std::invalid_argument(
            "StructuredRzMesh boundary node counts do not match the grid");

    const auto validate_side_set = [&](const std::string& name,
                                       std::size_t expected_sides,
                                       double fixed_coordinate,
                                       bool check_radial,
                                       const std::vector<std::size_t>&
                                           expected_boundary_nodes) {
        const SideSet& source_set = source.side_set(name);
        if (source_set.sides.size() != expected_sides)
            throw std::invalid_argument(
                "StructuredRzMesh side set count does not match the grid: " +
                name);
        for (const ElementSide& side : source_set.sides) {
            if (source.element_block_ids().at(side.element) != block_id)
                throw std::invalid_argument(
                    "Side set crosses element blocks: " + name);
            const Quad4Element& element = source.elements().at(side.element);
            const std::array<std::size_t, 2> side_nodes = {
                element.nodes[side.local_side],
                element.nodes[(side.local_side + 1U) % 4U]};
            for (const std::size_t source_node : side_nodes) {
                const RzPoint& point = source.nodes().at(source_node);
                const double coordinate = check_radial ? point.r : point.z;
                const std::size_t local_node = source_to_local.at(source_node);
                if (!same_coordinate(coordinate, fixed_coordinate) ||
                    std::find(expected_boundary_nodes.begin(),
                              expected_boundary_nodes.end(),
                              local_node) == expected_boundary_nodes.end())
                    throw std::invalid_argument(
                        "Side set is not on the expected RZ boundary: " + name);
            }
        }
    };
    if (!boundary_names.radial_inner.empty())
        validate_side_set(boundary_names.radial_inner, axial_elements,
                          mesh._inner_radius, true, mesh._radial_inner_nodes);
    if (!boundary_names.radial_outer.empty())
        validate_side_set(boundary_names.radial_outer, axial_elements,
                          mesh._outer_radius, true, mesh._radial_outer_nodes);
    if (!boundary_names.bottom.empty())
        validate_side_set(boundary_names.bottom, radial_elements,
                          axial_coordinates.front(), false, mesh._bottom_nodes);
    if (!boundary_names.top.empty())
        validate_side_set(boundary_names.top, radial_elements,
                          axial_coordinates.back(), false, mesh._top_nodes);

    const auto make_edges = [](const std::vector<std::size_t>& nodes) {
        std::vector<Line2BoundaryElement> edges;
        edges.reserve(nodes.size() - 1);
        for (std::size_t edge = 0; edge + 1 < nodes.size(); ++edge)
            edges.push_back({{{nodes[edge], nodes[edge + 1]}}});
        return edges;
    };
    mesh._radial_inner_elements = make_edges(mesh._radial_inner_nodes);
    mesh._radial_outer_elements = make_edges(mesh._radial_outer_nodes);
    mesh._bottom_elements = make_edges(mesh._bottom_nodes);
    mesh._top_elements = make_edges(mesh._top_nodes);
    return mesh;
}

StructuredRzBoundary
StructuredRzMesh::map_side_set(const UnstructuredQuad4Mesh& source,
                               const std::string& side_set_name,
                               std::int64_t expected_block_id) const {
    if (source.side_set_block_id(side_set_name) != expected_block_id)
        throw std::invalid_argument(
            "Side set belongs to an unexpected block: " + side_set_name);

    const SideSet& source_set = source.side_set(side_set_name);
    const auto local_node = [&](std::size_t source_node) {
        const RzPoint& source_point = source.nodes().at(source_node);
        const auto found = std::find_if(
            _nodes.begin(), _nodes.end(), [&](const RzPoint& point) {
                return same_coordinate(point.r, source_point.r) &&
                       same_coordinate(point.z, source_point.z);
            });
        if (found == _nodes.end())
            throw std::invalid_argument("Side set node is outside its block: " +
                                        side_set_name);
        return static_cast<std::size_t>(found - _nodes.begin());
    };

    std::vector<Line2BoundaryElement> elements;
    elements.reserve(source_set.sides.size());
    for (const ElementSide& side : source_set.sides) {
        const Quad4Element& element = source.elements().at(side.element);
        std::array<std::size_t, 2> nodes = {
            local_node(element.nodes.at(side.local_side)),
            local_node(element.nodes.at((side.local_side + 1U) % 4U)),
        };
        const RzPoint& first = _nodes.at(nodes[0]);
        const RzPoint& second = _nodes.at(nodes[1]);
        if (first.z > second.z ||
            (same_coordinate(first.z, second.z) && first.r > second.r))
            std::swap(nodes[0], nodes[1]);
        elements.push_back({nodes});
    }

    const auto all_on = [&](double coordinate, bool radial) {
        return std::all_of(elements.begin(), elements.end(),
                           [&](const Line2BoundaryElement& edge) {
                               return std::all_of(
                                   edge.nodes.begin(), edge.nodes.end(),
                                   [&](std::size_t node) {
                                       const RzPoint& point = _nodes.at(node);
                                       return same_coordinate(radial ? point.r
                                                                     : point.z,
                                                              coordinate);
                                   });
                           });
    };

    const auto axial_bounds = std::minmax_element(
        _nodes.begin(), _nodes.end(),
        [](const RzPoint& lhs, const RzPoint& rhs) { return lhs.z < rhs.z; });
    const double bottom = axial_bounds.first->z;
    const double top = axial_bounds.second->z;

    BoundaryId id = BoundaryId::radial_inner;
    bool sort_by_axial = true;
    if (all_on(_inner_radius, true)) {
        id = BoundaryId::radial_inner;
    } else if (all_on(_outer_radius, true)) {
        id = BoundaryId::radial_outer;
    } else if (all_on(bottom, false)) {
        id = BoundaryId::bottom;
        sort_by_axial = false;
    } else if (all_on(top, false)) {
        id = BoundaryId::top;
        sort_by_axial = false;
    } else {
        throw std::invalid_argument(
            "Side set is not on one structured RZ boundary: " + side_set_name);
    }

    std::sort(
        elements.begin(), elements.end(),
        [&](const Line2BoundaryElement& lhs, const Line2BoundaryElement& rhs) {
            const RzPoint& lhs_point = _nodes.at(lhs.nodes[0]);
            const RzPoint& rhs_point = _nodes.at(rhs.nodes[0]);
            return sort_by_axial ? lhs_point.z < rhs_point.z
                                 : lhs_point.r < rhs_point.r;
        });

    std::vector<std::size_t> nodes;
    for (const Line2BoundaryElement& edge : elements) {
        nodes.insert(nodes.end(), edge.nodes.begin(), edge.nodes.end());
    }
    std::sort(nodes.begin(), nodes.end(),
              [&](std::size_t lhs, std::size_t rhs) {
                  const RzPoint& lhs_point = _nodes.at(lhs);
                  const RzPoint& rhs_point = _nodes.at(rhs);
                  return sort_by_axial ? lhs_point.z < rhs_point.z
                                       : lhs_point.r < rhs_point.r;
              });
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    return {id, std::move(nodes), std::move(elements)};
}

StructuredRzMesh StructuredRzMesh::make_annulus(double inner_radius,
                                                double outer_radius,
                                                double length,
                                                std::size_t radial_elements,
                                                std::size_t axial_elements) {
    if (!(inner_radius >= 0.0))
        throw std::invalid_argument(
            "StructuredRzMesh inner_radius must be nonnegative");
    if (!(outer_radius > inner_radius))
        throw std::invalid_argument(
            "StructuredRzMesh outer_radius must exceed inner_radius");
    if (!(length > 0.0))
        throw std::invalid_argument("StructuredRzMesh length must be positive");
    if (radial_elements == 0 || axial_elements == 0)
        throw std::invalid_argument(
            "StructuredRzMesh element counts must be positive");

    const std::size_t radial_nodes = radial_elements + 1;
    const std::size_t axial_nodes = axial_elements + 1;
    if (radial_nodes > std::numeric_limits<std::size_t>::max() / axial_nodes)
        throw std::length_error("StructuredRzMesh node count overflows");

    StructuredRzMesh mesh;
    mesh._inner_radius = inner_radius;
    mesh._outer_radius = outer_radius;
    mesh._length = length;
    mesh._radial_elements = radial_elements;
    mesh._axial_elements = axial_elements;

    mesh._nodes.reserve(radial_nodes * axial_nodes);
    for (std::size_t iz = 0; iz < axial_nodes; ++iz) {
        const double z = length * static_cast<double>(iz) /
                         static_cast<double>(axial_elements);
        for (std::size_t ir = 0; ir < radial_nodes; ++ir) {
            const double fraction =
                static_cast<double>(ir) / static_cast<double>(radial_elements);
            const double r =
                inner_radius + (outer_radius - inner_radius) * fraction;
            mesh._nodes.push_back({r, z});
        }
    }

    mesh._elements.reserve(radial_elements * axial_elements);
    for (std::size_t iz = 0; iz < axial_elements; ++iz) {
        for (std::size_t ir = 0; ir < radial_elements; ++ir) {
            mesh._elements.push_back(
                {{mesh.node_id(ir, iz), mesh.node_id(ir + 1, iz),
                  mesh.node_id(ir + 1, iz + 1), mesh.node_id(ir, iz + 1)}});
        }
    }

    mesh._radial_inner_nodes.reserve(axial_nodes);
    mesh._radial_outer_nodes.reserve(axial_nodes);
    for (std::size_t iz = 0; iz < axial_nodes; ++iz) {
        mesh._radial_inner_nodes.push_back(mesh.node_id(0, iz));
        mesh._radial_outer_nodes.push_back(mesh.node_id(radial_elements, iz));
    }

    mesh._bottom_nodes.reserve(radial_nodes);
    mesh._top_nodes.reserve(radial_nodes);
    for (std::size_t ir = 0; ir < radial_nodes; ++ir) {
        mesh._bottom_nodes.push_back(mesh.node_id(ir, 0));
        mesh._top_nodes.push_back(mesh.node_id(ir, axial_elements));
    }

    mesh._radial_inner_elements.reserve(axial_elements);
    mesh._radial_outer_elements.reserve(axial_elements);
    for (std::size_t iz = 0; iz < axial_elements; ++iz) {
        mesh._radial_inner_elements.push_back(
            {{{mesh.node_id(0, iz), mesh.node_id(0, iz + 1)}}});
        mesh._radial_outer_elements.push_back(
            {{{mesh.node_id(radial_elements, iz),
               mesh.node_id(radial_elements, iz + 1)}}});
    }

    mesh._bottom_elements.reserve(radial_elements);
    mesh._top_elements.reserve(radial_elements);
    for (std::size_t ir = 0; ir < radial_elements; ++ir) {
        mesh._bottom_elements.push_back(
            {{{mesh.node_id(ir, 0), mesh.node_id(ir + 1, 0)}}});
        mesh._top_elements.push_back(
            {{{mesh.node_id(ir, axial_elements),
               mesh.node_id(ir + 1, axial_elements)}}});
    }

    return mesh;
}

const std::vector<RzPoint>& StructuredRzMesh::nodes() const noexcept {
    return _nodes;
}

const std::vector<Quad4Element>& StructuredRzMesh::elements() const noexcept {
    return _elements;
}

const std::vector<std::size_t>&
StructuredRzMesh::boundary_nodes(BoundaryId boundary) const {
    switch (boundary) {
    case BoundaryId::radial_inner:
        return _radial_inner_nodes;
    case BoundaryId::radial_outer:
        return _radial_outer_nodes;
    case BoundaryId::bottom:
        return _bottom_nodes;
    case BoundaryId::top:
        return _top_nodes;
    }
    throw std::invalid_argument("Unknown StructuredRzMesh boundary");
}

const std::vector<Line2BoundaryElement>&
StructuredRzMesh::boundary_elements(BoundaryId boundary) const {
    switch (boundary) {
    case BoundaryId::radial_inner:
        return _radial_inner_elements;
    case BoundaryId::radial_outer:
        return _radial_outer_elements;
    case BoundaryId::bottom:
        return _bottom_elements;
    case BoundaryId::top:
        return _top_elements;
    }
    throw std::invalid_argument("Unknown StructuredRzMesh boundary");
}

std::size_t StructuredRzMesh::radial_elements() const noexcept {
    return _radial_elements;
}

std::size_t StructuredRzMesh::axial_elements() const noexcept {
    return _axial_elements;
}

std::size_t StructuredRzMesh::node_id(std::size_t radial_index,
                                      std::size_t axial_index) const {
    if (radial_index > _radial_elements || axial_index > _axial_elements)
        throw std::out_of_range("StructuredRzMesh node index is out of range");
    return axial_index * (_radial_elements + 1) + radial_index;
}

double StructuredRzMesh::inner_radius() const noexcept {
    return _inner_radius;
}

double StructuredRzMesh::outer_radius() const noexcept {
    return _outer_radius;
}

double StructuredRzMesh::length() const noexcept {
    return _length;
}

} // namespace fuelsim
