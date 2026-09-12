#include "core/mesh.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace fuelsim {
UnstructuredBar2Mesh::UnstructuredBar2Mesh(std::vector<RzPoint> nodes,
    std::vector<Bar2Element> elements,
    std::vector<std::int64_t> element_block_ids,
    std::vector<ElementBlockInfo> element_blocks,
    std::vector<NodeSet> node_sets,
    std::vector<SideSet> side_sets)
    : UnstructuredMeshMetadata(nodes.size(),
          elements.size(),
          2,
          std::move(element_block_ids),
          std::move(element_blocks),
          std::move(node_sets),
          std::move(side_sets),
          "UnstructuredBar2Mesh"),
      _nodes(std::move(nodes)), _elements(std::move(elements)) {
    for (const auto& node : _nodes)
        if (!std::isfinite(node.r) || !std::isfinite(node.z) || node.r < 0.0)
            throw std::invalid_argument("BAR2 requires finite RZ coordinates and nonnegative radius");

    constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();
    std::vector<unsigned char> roles(_nodes.size(), 0);
    std::vector<std::size_t> predecessor(_nodes.size(), invalid), successor(_nodes.size(), invalid);
    std::vector<std::array<std::size_t, 2>> radial_layers(_nodes.size(), {{invalid, invalid}});
    std::map<std::array<std::size_t, 2>, std::vector<std::array<double, 2>>> layer_intervals;
    for (const auto& element : _elements) {
        for (const auto node : element.nodes) {
            if (node >= _nodes.size())
                throw std::out_of_range("BAR2 radial connectivity is out of range");
            roles[node] |= 1U;
        }
        for (const auto node : element.axial_nodes) {
            if (node >= _nodes.size())
                throw std::out_of_range("BAR2 axial connectivity is out of range");
            roles[node] |= 2U;
        }
        for (const auto node : element.nodes)
            if (roles[node] == 3U)
                throw std::invalid_argument("BAR2 source node cannot carry both radial and axial fields");
        for (const auto node : element.axial_nodes)
            if (roles[node] == 3U)
                throw std::invalid_argument("BAR2 source node cannot carry both radial and axial fields");

        const auto lower = element.axial_nodes[0], upper = element.axial_nodes[1];
        const double lower_z = _nodes[lower].z, upper_z = _nodes[upper].z;
        const double height = upper_z - lower_z;
        if (!(height > 0.0) || !std::isfinite(height))
            throw std::invalid_argument("BAR2 axial controls must have finite positive ordered height");
        const auto& inner = _nodes[element.nodes[0]];
        const auto& outer = _nodes[element.nodes[1]];
        if (!(outer.r > inner.r))
            throw std::invalid_argument("BAR2 radial nodes must have strictly increasing radii");
        const double midpoint = lower_z + 0.5 * height;
        const double scale = std::max({height, std::abs(lower_z), std::abs(upper_z)});
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
        if (std::abs(inner.z - midpoint) > tolerance || std::abs(outer.z - midpoint) > tolerance)
            throw std::invalid_argument("BAR2 radial nodes must lie at the axial control midpoint");

        if ((successor[lower] != invalid && successor[lower] != upper)
            || (predecessor[upper] != invalid && predecessor[upper] != lower))
            throw std::invalid_argument("BAR2 axial controls must form independent unbranched chains");
        successor[lower] = upper;
        predecessor[upper] = lower;
        for (const auto node : element.nodes) {
            if (radial_layers[node][0] != invalid && radial_layers[node] != element.axial_nodes)
                throw std::invalid_argument("BAR2 radial source node cannot be shared across axial layers");
            radial_layers[node] = element.axial_nodes;
        }
        layer_intervals[element.axial_nodes].push_back({{inner.r, outer.r}});
    }
    for (auto& layer : layer_intervals) {
        auto& intervals = layer.second;
        std::sort(intervals.begin(), intervals.end());
        for (std::size_t i = 1; i < intervals.size(); ++i)
            if (intervals[i][0] < intervals[i - 1][1])
                throw std::invalid_argument("BAR2 radial elements in one axial layer must not overlap");
    }
    for (std::size_t node = 0; node < roles.size(); ++node) {
        if (roles[node] == 1U)
            _radial_source_node_ids.push_back(node);
        else if (roles[node] == 2U)
            _axial_source_node_ids.push_back(node);
    }
}
} // namespace fuelsim
