#include "fuelsim/dof_map.hpp"
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

// Field-major degree-of-freedom mapping.

DofMap::DofMap(std::size_t node_count) : _node_count(node_count) {
    if (node_count == 0)
        throw std::invalid_argument("DofMap node_count must be positive");
    if (node_count > std::numeric_limits<std::size_t>::max() / 3)
        throw std::length_error("DofMap DOF count overflows");
}

std::size_t DofMap::node_count() const noexcept {
    return _node_count;
}

std::size_t DofMap::dof_count() const noexcept {
    return 3 * _node_count;
}

std::size_t DofMap::dof(Field field, std::size_t node) const {
    if (node >= _node_count)
        throw std::out_of_range("DofMap node index is out of range");

    switch (field) {
    case Field::temperature:
        return node;
    case Field::radial_displacement:
        return _node_count + node;
    case Field::axial_displacement:
        return 2 * _node_count + node;
    }
    throw std::invalid_argument("Unknown DofMap field");
}

std::size_t DofMap::temperature(std::size_t node) const {
    return dof(Field::temperature, node);
}

std::size_t DofMap::radial_displacement(std::size_t node) const {
    return dof(Field::radial_displacement, node);
}

std::size_t DofMap::axial_displacement(std::size_t node) const {
    return dof(Field::axial_displacement, node);
}

LocalDofs
DofMap::local_dofs(const std::array<std::size_t, 4>& global_nodes) const {
    LocalDofs result{};
    for (std::size_t node = 0; node < global_nodes.size(); ++node) {
        result[node] = temperature(global_nodes[node]);
        result[global_nodes.size() + node] =
            radial_displacement(global_nodes[node]);
        result[2 * global_nodes.size() + node] =
            axial_displacement(global_nodes[node]);
    }
    return result;
}

LocalDofs DofMap::element_dofs(const Quad4Element& element) const {
    return local_dofs(element.nodes);
}

} // namespace fuelsim
