#include "fuelsim/core/mesh.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fuelsim {
UnstructuredQuad8Mesh::UnstructuredQuad8Mesh(std::vector<RzPoint> nodes,
    std::vector<Quad8Element> elements,
    std::vector<std::int64_t> ids,
    std::vector<ElementBlockInfo> blocks,
    std::vector<NodeSet> node_sets,
    std::vector<SideSet> side_sets)
    : UnstructuredMeshMetadata(nodes.size(),
          elements.size(),
          4,
          std::move(ids),
          std::move(blocks),
          std::move(node_sets),
          std::move(side_sets),
          "UnstructuredQuad8Mesh"),
      _nodes(std::move(nodes)), _elements(std::move(elements)) {
    for (const auto& node : _nodes)
        if (!std::isfinite(node.r) || !std::isfinite(node.z) || node.r < 0.0)
            throw std::invalid_argument("QUAD8 requires finite coordinates and nonnegative radius");
    std::vector<unsigned char> roles(_nodes.size(), 0);
    for (const auto& element : _elements) {
        auto sorted = element.nodes;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            throw std::invalid_argument("QUAD8 connectivity contains repeated nodes");
        for (std::size_t n = 0; n < 8; ++n) {
            const auto node = element.nodes[n];
            if (node >= _nodes.size())
                throw std::out_of_range("QUAD8 connectivity is out of range");
            roles[node] |= n < 4 ? 1U : 2U;
            if (roles[node] == 3U)
                throw std::invalid_argument("QUAD8 node cannot be both a temperature corner and an edge midpoint");
        }
    }
}

Quad8RegionMesh::Quad8RegionMesh(const UnstructuredQuad8Mesh& source, std::int64_t block)
    : RegionMeshMapping(source, source.nodes().size(), block) {
}

Quad8RegionMesh Quad8RegionMesh::from_unstructured_block(const UnstructuredQuad8Mesh& source, std::int64_t block) {
    Quad8RegionMesh mesh(source, block);
    std::vector<bool> used(source.nodes().size(), false), thermal(source.nodes().size(), false);
    for (auto e : mesh._source_element_ids) {
        const auto& element = source.elements()[e];
        for (auto n : element.nodes)
            used[n] = true;
        for (std::size_t n = 0; n < 4; ++n)
            thermal[element.nodes[n]] = true;
    }
    mesh.select_nodes(used);
    for (auto n : mesh._source_node_ids) {
        mesh._nodes.push_back(source.nodes()[n]);
        mesh._temperature_nodes.push_back(thermal[n]);
    }
    for (auto e : mesh._source_element_ids) {
        Quad8Element element{};
        for (std::size_t n = 0; n < 8; ++n)
            element.nodes[n] = mesh._source_node_to_local[source.elements()[e].nodes[n]];
        mesh._elements.push_back(element);
    }
    return mesh;
}

Quad8RegionBoundary Quad8RegionMesh::map_side_set(const UnstructuredQuad8Mesh& source, const std::string& name) const {
    if (source.side_set_block_id(name) != _block_id)
        throw std::invalid_argument("QUAD8 side set belongs to an unexpected region: " + name);
    Quad8RegionBoundary result;
    for (const auto& side : source.side_set(name).sides) {
        const auto e = _source_element_to_local.at(side.element), s = side.local_side;
        if (e == invalid_index)
            throw std::invalid_argument("QUAD8 side is outside its region");
        const auto& n = _elements[e].nodes;
        const Line3BoundaryElement edge{{n[s], n[(s + 1) % 4], n[4 + s]}, e, s};
        result.elements.push_back(edge);
        result.temperature_nodes.insert(result.temperature_nodes.end(), edge.nodes.begin(), edge.nodes.begin() + 2);
        result.displacement_nodes.insert(result.displacement_nodes.end(), edge.nodes.begin(), edge.nodes.end());
    }
    std::sort(result.temperature_nodes.begin(), result.temperature_nodes.end());
    result.temperature_nodes.erase(std::unique(result.temperature_nodes.begin(), result.temperature_nodes.end()),
        result.temperature_nodes.end());
    std::sort(result.displacement_nodes.begin(), result.displacement_nodes.end());
    result.displacement_nodes.erase(std::unique(result.displacement_nodes.begin(), result.displacement_nodes.end()),
        result.displacement_nodes.end());
    return result;
}
} // namespace fuelsim
