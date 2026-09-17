#include "core/plane_mesh.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fuelsim {
UnstructuredPlaneQuad8Mesh::UnstructuredPlaneQuad8Mesh(std::vector<std::array<double, 2>> nodes,
    std::vector<PlaneQuad8Element> elements,
    std::vector<std::int64_t> block_ids,
    std::vector<ElementBlockInfo> blocks,
    std::vector<NodeSet> node_sets,
    std::vector<SideSet> side_sets)
    : UnstructuredMeshMetadata(nodes.size(),
          elements.size(),
          4,
          std::move(block_ids),
          std::move(blocks),
          std::move(node_sets),
          std::move(side_sets),
          "Generalized plane strain QUAD8"),
      _nodes(std::move(nodes)), _elements(std::move(elements)) {
    for (const auto& node : _nodes)
        for (double value : node)
            if (!std::isfinite(value))
                throw std::invalid_argument("Plane QUAD8 coordinates must be finite");
    std::vector<unsigned char> roles(_nodes.size(), 0);
    for (const auto& element : _elements) {
        auto sorted = element.nodes;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            throw std::invalid_argument("Plane QUAD8 element has repeated nodes");
        elements::Cpeg8Coordinates coordinates;
        for (std::size_t n = 0; n < 8; ++n) {
            if (element.nodes[n] >= _nodes.size())
                throw std::invalid_argument("Plane QUAD8 connectivity is out of range");
            roles[element.nodes[n]] |= n < 4 ? 1U : 2U;
            coordinates[n] = _nodes[element.nodes[n]];
        }
        (void)elements::make_cpeg8t_geometry(coordinates, 1.0);
    }
    for (unsigned char role : roles)
        if (role != 0 && role != 1 && role != 2)
            throw std::invalid_argument("Plane QUAD8 corner and midpoint roles must be distinct");
}
} // namespace fuelsim
