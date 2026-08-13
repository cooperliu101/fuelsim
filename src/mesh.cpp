#include "fuelsim/mesh.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
namespace fuelsim {
namespace {
bool same_coordinate(double lhs, double rhs) {
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= 1.0e-12 * scale;
}
} // namespace
UnstructuredMeshMetadata::UnstructuredMeshMetadata(std::size_t node_count, std::size_t element_count,
    std::size_t sides_per_element, std::vector<std::int64_t> element_block_ids,
    std::vector<ElementBlockInfo> element_blocks, std::vector<NodeSet> node_sets, std::vector<SideSet> side_sets,
    const std::string& geometry_name)
    : _element_block_ids(std::move(element_block_ids)), _element_blocks(std::move(element_blocks)),
      _node_sets(std::move(node_sets)), _side_sets(std::move(side_sets)) {
    if (node_count == 0 || element_count == 0)
        throw std::invalid_argument(geometry_name + " requires at least one node and one element");
    if (element_count != _element_block_ids.size())
        throw std::invalid_argument(geometry_name + " block ID count must match element count");
    if (_element_blocks.empty()) throw std::invalid_argument(geometry_name + " requires element block metadata");
    for (std::size_t block = 0; block < _element_blocks.size(); ++block) {
        if (_element_blocks[block].id < 0)
            throw std::invalid_argument(geometry_name + " block IDs must be nonnegative");
        for (std::size_t previous = 0; previous < block; ++previous)
            if (_element_blocks[previous].id == _element_blocks[block].id ||
                (!_element_blocks[block].name.empty() && _element_blocks[previous].name == _element_blocks[block].name))
                throw std::invalid_argument(geometry_name + " block IDs and names must be unique");
    }
    for (const std::int64_t id : _element_block_ids)
        if (std::none_of(_element_blocks.begin(), _element_blocks.end(),
                [id](const ElementBlockInfo& block) { return block.id == id; }))
            throw std::invalid_argument(geometry_name + " element references an unknown block");
    for (std::size_t set = 0; set < _node_sets.size(); ++set) {
        const NodeSet& current = _node_sets[set];
        if (current.id < 0) throw std::invalid_argument(geometry_name + " node set IDs must be nonnegative");
        for (std::size_t previous = 0; previous < set; ++previous)
            if (_node_sets[previous].id == current.id ||
                (!current.name.empty() && _node_sets[previous].name == current.name))
                throw std::invalid_argument(geometry_name + " node set IDs and names must be unique");
        for (const std::size_t node : current.nodes)
            if (node >= node_count) throw std::out_of_range(geometry_name + " node set is out of range");
    }
    for (std::size_t set = 0; set < _side_sets.size(); ++set) {
        const SideSet& current = _side_sets[set];
        if (current.id < 0) throw std::invalid_argument(geometry_name + " side set IDs must be nonnegative");
        for (std::size_t previous = 0; previous < set; ++previous)
            if (_side_sets[previous].id == current.id ||
                (!current.name.empty() && _side_sets[previous].name == current.name))
                throw std::invalid_argument(geometry_name + " side set IDs and names must be unique");
        for (const ElementSide& side : current.sides)
            if (side.element >= element_count || side.local_side >= sides_per_element)
                throw std::out_of_range(geometry_name + " side set is out of range");
    }
}
const ElementBlockInfo& UnstructuredMeshMetadata::element_block(const std::string& name) const {
    const auto block = std::find_if(_element_blocks.begin(), _element_blocks.end(),
        [&name](const ElementBlockInfo& candidate) { return candidate.name == name; });
    if (block == _element_blocks.end()) throw std::invalid_argument("Unknown element block: " + name);
    return *block;
}
const SideSet& UnstructuredMeshMetadata::side_set(const std::string& name) const {
    const auto set = std::find_if(
        _side_sets.begin(), _side_sets.end(), [&name](const SideSet& candidate) { return candidate.name == name; });
    if (set == _side_sets.end()) throw std::invalid_argument("Unknown side set: " + name);
    return *set;
}
std::int64_t UnstructuredMeshMetadata::side_set_block_id(const std::string& name) const {
    const SideSet& set = side_set(name);
    if (set.sides.empty()) throw std::invalid_argument("Side set is empty: " + name);
    const std::int64_t block_id = _element_block_ids.at(set.sides.front().element);
    for (const ElementSide& side : set.sides)
        if (_element_block_ids.at(side.element) != block_id)
            throw std::invalid_argument("Side set crosses element blocks: " + name);
    return block_id;
}
UnstructuredQuad4Mesh::UnstructuredQuad4Mesh(std::vector<RzPoint> nodes, std::vector<Quad4Element> elements,
    std::vector<std::int64_t> element_block_ids, std::vector<ElementBlockInfo> element_blocks,
    std::vector<NodeSet> node_sets, std::vector<SideSet> side_sets)
    : UnstructuredMeshMetadata(nodes.size(), elements.size(), 4, std::move(element_block_ids),
          std::move(element_blocks), std::move(node_sets), std::move(side_sets), "UnstructuredQuad4Mesh"),
      _nodes(std::move(nodes)), _elements(std::move(elements)) {
    for (const RzPoint& node : _nodes)
        if (!std::isfinite(node.r) || !std::isfinite(node.z) || node.r < 0.0)
            throw std::invalid_argument("UnstructuredQuad4Mesh requires finite RZ coordinates and nonnegative radius");
    for (const Quad4Element& element : _elements) {
        for (const std::size_t node : element.nodes)
            if (node >= _nodes.size()) throw std::out_of_range("UnstructuredQuad4Mesh connectivity is out of range");
    }
}
UnstructuredHex8Mesh::UnstructuredHex8Mesh(std::vector<CartesianPoint3> nodes, std::vector<Hex8Element> elements,
    std::vector<std::int64_t> element_block_ids, std::vector<ElementBlockInfo> element_blocks,
    std::vector<NodeSet> node_sets, std::vector<SideSet> side_sets)
    : UnstructuredMeshMetadata(nodes.size(), elements.size(), 6, std::move(element_block_ids),
          std::move(element_blocks), std::move(node_sets), std::move(side_sets), "UnstructuredHex8Mesh"),
      _nodes(std::move(nodes)), _elements(std::move(elements)) {
    for (const CartesianPoint3& node : _nodes)
        if (!std::isfinite(node.x) || !std::isfinite(node.y) || !std::isfinite(node.z))
            throw std::invalid_argument("UnstructuredHex8Mesh requires finite Cartesian coordinates");
    for (const Hex8Element& element : _elements) {
        for (const std::size_t node : element.nodes)
            if (node >= _nodes.size()) throw std::out_of_range("UnstructuredHex8Mesh connectivity is out of range");
    }
}
RegionMeshMapping::RegionMeshMapping(
    const UnstructuredMeshMetadata& source, std::size_t node_count, std::int64_t block_id)
    : _block_id(block_id), _source_node_to_local(node_count, invalid_index),
      _source_element_to_local(source.element_block_ids().size(), invalid_index) {
    const bool known_block = std::any_of(source.element_blocks().begin(), source.element_blocks().end(),
        [block_id](const ElementBlockInfo& block) { return block.id == block_id; });
    if (!known_block) throw std::invalid_argument("Unknown element block ID: " + std::to_string(block_id));
    for (std::size_t source_element = 0; source_element < source.element_block_ids().size(); ++source_element) {
        if (source.element_block_ids()[source_element] != block_id) continue;
        _source_element_to_local[source_element] = _source_element_ids.size();
        _source_element_ids.push_back(source_element);
    }
    if (_source_element_ids.empty()) throw std::invalid_argument("Element block is empty: " + std::to_string(block_id));
}
void RegionMeshMapping::select_nodes(const std::vector<bool>& used_nodes) {
    if (used_nodes.size() != _source_node_to_local.size())
        throw std::logic_error("Region node selection has the wrong size");
    for (std::size_t source_node = 0; source_node < used_nodes.size(); ++source_node) {
        if (!used_nodes[source_node]) continue;
        _source_node_to_local[source_node] = _source_node_ids.size();
        _source_node_ids.push_back(source_node);
    }
}
Hex8RegionMesh::Hex8RegionMesh(const UnstructuredHex8Mesh& source, std::int64_t block_id)
    : RegionMeshMapping(source, source.nodes().size(), block_id) {}
Hex8RegionMesh Hex8RegionMesh::from_unstructured_block(const UnstructuredHex8Mesh& source, std::int64_t block_id) {
    Hex8RegionMesh mesh(source, block_id);
    std::vector<bool> used_nodes(source.nodes().size(), false);
    for (const std::size_t source_element : mesh._source_element_ids)
        for (std::size_t node : source.elements()[source_element].nodes) used_nodes[node] = true;
    mesh.select_nodes(used_nodes);
    for (const std::size_t source_node : mesh._source_node_ids) mesh._nodes.push_back(source.nodes()[source_node]);
    mesh._elements.reserve(mesh._source_element_ids.size());
    for (std::size_t source_element : mesh._source_element_ids) {
        Hex8Element element{};
        for (std::size_t node = 0; node < element.nodes.size(); ++node) {
            const std::size_t local = mesh._source_node_to_local.at(source.elements()[source_element].nodes[node]);
            if (local == invalid_index) throw std::logic_error("Hex8RegionMesh connectivity crosses element blocks");
            element.nodes[node] = local;
        }
        mesh._elements.push_back(element);
    }
    return mesh;
}
Hex8RegionBoundary Hex8RegionMesh::map_side_set(
    const UnstructuredHex8Mesh& source, const std::string& side_set_name) const {
    if (source.side_set_block_id(side_set_name) != _block_id)
        throw std::invalid_argument("Side set belongs to an unexpected block: " + side_set_name);
    static constexpr std::array<std::array<std::size_t, 4>, 6> face_nodes = {{
        {{0, 1, 5, 4}},
        {{1, 2, 6, 5}},
        {{2, 3, 7, 6}},
        {{0, 4, 7, 3}},
        {{0, 3, 2, 1}},
        {{4, 5, 6, 7}},
    }};
    Hex8RegionBoundary result;
    for (const ElementSide& side : source.side_set(side_set_name).sides) {
        const std::size_t local_element = _source_element_to_local.at(side.element);
        if (local_element == invalid_index)
            throw std::invalid_argument("Side set is outside its region: " + side_set_name);
        Quad4FaceElement face{{}, local_element, side.local_side};
        for (std::size_t node = 0; node < face.nodes.size(); ++node) {
            const std::size_t source_node = source.elements()[side.element].nodes[face_nodes[side.local_side][node]],
                              local_node = _source_node_to_local.at(source_node);
            if (local_node == invalid_index) throw std::logic_error("Hex8RegionMesh side-set node mapping failed");
            face.nodes[node] = local_node;
            result.nodes.push_back(local_node);
        }
        result.faces.push_back(face);
    }
    std::sort(result.nodes.begin(), result.nodes.end(),
        [&](std::size_t lhs, std::size_t rhs) { return _source_node_ids[lhs] < _source_node_ids[rhs]; });
    result.nodes.erase(std::unique(result.nodes.begin(), result.nodes.end()), result.nodes.end());
    return result;
}
RegionMesh::RegionMesh(const UnstructuredQuad4Mesh& source, std::int64_t block_id)
    : RegionMeshMapping(source, source.nodes().size(), block_id) {}
RegionMesh RegionMesh::from_unstructured_block(const UnstructuredQuad4Mesh& source, std::int64_t block_id) {
    RegionMesh mesh(source, block_id);
    std::vector<bool> used_nodes(source.nodes().size(), false);
    for (const std::size_t source_element : mesh._source_element_ids)
        for (std::size_t node : source.elements()[source_element].nodes) used_nodes[node] = true;
    mesh.select_nodes(used_nodes);
    for (const std::size_t source_node : mesh._source_node_ids) mesh._nodes.push_back(source.nodes()[source_node]);
    mesh._elements.reserve(mesh._source_element_ids.size());
    for (std::size_t source_element : mesh._source_element_ids) {
        Quad4Element element{};
        for (std::size_t node = 0; node < element.nodes.size(); ++node) {
            const std::size_t local = mesh._source_node_to_local.at(source.elements()[source_element].nodes[node]);
            if (local == invalid_index) throw std::logic_error("RegionMesh connectivity crosses element blocks");
            element.nodes[node] = local;
        }
        mesh._elements.push_back(element);
    }
    return mesh;
}
RegionBoundary RegionMesh::map_side_set(const UnstructuredQuad4Mesh& source, const std::string& side_set_name) const {
    if (source.side_set_block_id(side_set_name) != _block_id)
        throw std::invalid_argument("Side set belongs to an unexpected block: " + side_set_name);
    std::vector<Line2BoundaryElement> elements;
    std::vector<RzPoint> adjacent_centroids;
    const SideSet& side_set = source.side_set(side_set_name);
    elements.reserve(side_set.sides.size());
    adjacent_centroids.reserve(side_set.sides.size());
    for (const ElementSide& side : side_set.sides) {
        if (_source_element_to_local.at(side.element) == invalid_index)
            throw std::invalid_argument("Side set is outside its region: " + side_set_name);
        const Quad4Element& source_element = source.elements().at(side.element);
        const std::size_t first_source = source_element.nodes.at(side.local_side),
                          second_source = source_element.nodes.at((side.local_side + 1U) % 4U),
                          first = _source_node_to_local.at(first_source),
                          second = _source_node_to_local.at(second_source);
        if (first == invalid_index || second == invalid_index)
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
        const double coordinate = radial ? reference_point.r : reference_point.z;
        return std::all_of(elements.begin(), elements.end(), [&](const Line2BoundaryElement& edge) {
            return std::all_of(edge.nodes.begin(), edge.nodes.end(), [&](std::size_t node) {
                const RzPoint& point = _nodes.at(node);
                return same_coordinate(radial ? point.r : point.z, coordinate);
            });
        });
    };
    const auto material_side = [&](bool radial) {
        const double coordinate = radial ? reference_point.r : reference_point.z;
        int side = 0;
        for (const RzPoint& centroid : adjacent_centroids) {
            const double value = radial ? centroid.r : centroid.z;
            if (same_coordinate(value, coordinate)) return 0;
            const int current = value > coordinate ? 1 : -1;
            if (side != 0 && current != side) return 0;
            side = current;
        }
        return side;
    };
    RegionBoundaryKind kind = RegionBoundaryKind::general;
    bool sort_by_axial = false;
    const int radial_material_side = all_on(true) ? material_side(true) : 0,
              axial_material_side = all_on(false) ? material_side(false) : 0;
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
            const RzPoint &first = _nodes.at(edge.nodes[0]), &second = _nodes.at(edge.nodes[1]);
            if ((sort_by_axial && first.z > second.z) || (!sort_by_axial && first.r > second.r))
                std::swap(edge.nodes[0], edge.nodes[1]);
        }
        std::sort(
            elements.begin(), elements.end(), [&](const Line2BoundaryElement& lhs, const Line2BoundaryElement& rhs) {
                const RzPoint &lhs_point = _nodes.at(lhs.nodes[0]), &rhs_point = _nodes.at(rhs.nodes[0]);
                return sort_by_axial ? lhs_point.z < rhs_point.z : lhs_point.r < rhs_point.r;
            });
    }
    std::vector<std::size_t> nodes;
    for (const Line2BoundaryElement& edge : elements) nodes.insert(nodes.end(), edge.nodes.begin(), edge.nodes.end());
    std::sort(nodes.begin(), nodes.end(), [&](std::size_t lhs, std::size_t rhs) {
        if (kind == RegionBoundaryKind::radial_inner || kind == RegionBoundaryKind::radial_outer)
            return _nodes[lhs].z < _nodes[rhs].z;
        if (kind == RegionBoundaryKind::bottom || kind == RegionBoundaryKind::top) return _nodes[lhs].r < _nodes[rhs].r;
        return _source_node_ids[lhs] < _source_node_ids[rhs];
    });
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    return {kind, std::move(nodes), std::move(elements)};
}
} // namespace fuelsim
