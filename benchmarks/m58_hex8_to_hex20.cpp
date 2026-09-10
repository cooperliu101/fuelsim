#include "io/results_io.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
fuelsim::UnstructuredHex20Mesh upgrade(const fuelsim::UnstructuredHex8Mesh& source) {
    static constexpr std::array<std::array<std::size_t, 2>, 12> edge_corners = {{{{0, 1}},
        {{1, 2}},
        {{2, 3}},
        {{3, 0}},
        {{0, 4}},
        {{1, 5}},
        {{2, 6}},
        {{3, 7}},
        {{4, 5}},
        {{5, 6}},
        {{6, 7}},
        {{7, 4}}}};

    std::vector<fuelsim::CartesianPoint3> nodes = source.nodes();
    std::vector<fuelsim::Hex20Element> elements;
    elements.reserve(source.elements().size());
    std::map<std::array<std::size_t, 2>, std::size_t> edge_nodes;
    for (const fuelsim::Hex8Element& source_element : source.elements()) {
        fuelsim::Hex20Element element{};
        std::copy(source_element.nodes.begin(), source_element.nodes.end(), element.nodes.begin());
        for (std::size_t edge = 0; edge < edge_corners.size(); ++edge) {
            const std::size_t first = source_element.nodes[edge_corners[edge][0]],
                              second = source_element.nodes[edge_corners[edge][1]];
            const std::array<std::size_t, 2> key = {std::min(first, second), std::max(first, second)};
            const auto inserted = edge_nodes.emplace(key, nodes.size());
            if (inserted.second) {
                const fuelsim::CartesianPoint3& a = source.nodes().at(first);
                const fuelsim::CartesianPoint3& b = source.nodes().at(second);
                nodes.push_back({0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z)});
            }
            element.nodes[8 + edge] = inserted.first->second;
        }
        elements.push_back(element);
    }
    return fuelsim::UnstructuredHex20Mesh(std::move(nodes),
        std::move(elements),
        source.element_block_ids(),
        source.element_blocks(),
        source.node_sets(),
        source.side_sets());
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_m58_hex8_to_hex20 <input-hex8.e> <output-hex20.e>\n";
        return 2;
    }
    try {
        const fuelsim::UnstructuredHex8Mesh source = fuelsim::read_exodus_hex8(argv[1]);
        const fuelsim::UnstructuredHex20Mesh converted = upgrade(source);
        fuelsim::write_exodus_hex20(argv[2], converted);
        const std::size_t corner_nodes = source.nodes().size();
        const std::size_t coupled_dofs = corner_nodes + 3U * converted.nodes().size();
        std::cout << "m58_hex8_elements=" << source.elements().size() << '\n'
                  << "m58_hex20_nodes=" << converted.nodes().size() << '\n'
                  << "m58_hex20_corner_temperature_nodes=" << corner_nodes << '\n'
                  << "m58_hex20_coupled_dofs=" << coupled_dofs << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "M5.8 HEX8-to-HEX20 conversion failed: " << error.what() << '\n';
        return 1;
    }
}
