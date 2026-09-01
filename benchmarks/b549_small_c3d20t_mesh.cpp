#include "fuelsim/io/results_io.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using Edge = std::array<std::size_t, 2>;

fuelsim::UnstructuredHex20Mesh make_mesh() {
    constexpr std::size_t normal_elements = 4;
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::vector<std::array<std::size_t, 8>> corners;
    std::vector<std::int64_t> blocks;
    std::array<std::vector<fuelsim::ElementSide>, 4> faces;
    const auto append_block = [&](double x_lower, double x_upper, double transverse_lower, double transverse_upper,
                                  std::int64_t block) {
        const std::size_t offset = nodes.size();
        std::array<double, normal_elements + 1> x{};
        for (std::size_t index = 0; index <= normal_elements; ++index)
            x[index] =
                x_lower + (x_upper - x_lower) * static_cast<double>(index) / static_cast<double>(normal_elements);
        for (double z : {transverse_lower, transverse_upper})
            for (double y : {transverse_lower, transverse_upper})
                for (double x_value : x) nodes.push_back({x_value, y, z});
        const auto node = [offset, normal_elements](std::size_t ix, std::size_t iy, std::size_t iz) {
            return offset + iz * 2U * (normal_elements + 1U) + iy * (normal_elements + 1U) + ix;
        };
        for (std::size_t ix = 0; ix < normal_elements; ++ix) {
            const std::size_t element = corners.size();
            corners.push_back({node(ix, 0, 0), node(ix + 1, 0, 0), node(ix + 1, 1, 0), node(ix, 1, 0), node(ix, 0, 1),
                node(ix + 1, 0, 1), node(ix + 1, 1, 1), node(ix, 1, 1)});
            blocks.push_back(block);
            if (block == 1 && ix == 0) faces[0].push_back({element, 3});
            if (block == 1 && ix + 1 == normal_elements) faces[1].push_back({element, 1});
            if (block == 2 && ix == 0) faces[2].push_back({element, 3});
            if (block == 2 && ix + 1 == normal_elements) faces[3].push_back({element, 1});
        }
    };
    append_block(0.0, 1.0, -0.1, 1.1, 1);
    append_block(1.0, 2.0, 0.0, 1.0, 2);

    static constexpr std::array<Edge, 12> edge_corners = {{{{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}, {{0, 4}}, {{1, 5}},
        {{2, 6}}, {{3, 7}}, {{4, 5}}, {{5, 6}}, {{6, 7}}, {{7, 4}}}};
    std::map<Edge, std::size_t> edge_nodes;
    std::vector<fuelsim::Hex20Element> elements;
    for (const auto& element_corners : corners) {
        fuelsim::Hex20Element element{};
        std::copy(element_corners.begin(), element_corners.end(), element.nodes.begin());
        for (std::size_t edge = 0; edge < edge_corners.size(); ++edge) {
            const std::size_t first = element_corners[edge_corners[edge][0]];
            const std::size_t second = element_corners[edge_corners[edge][1]];
            const Edge key = {std::min(first, second), std::max(first, second)};
            const auto inserted = edge_nodes.emplace(key, nodes.size());
            if (inserted.second) {
                const auto& a = nodes[first];
                const auto& b = nodes[second];
                nodes.push_back({0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z)});
            }
            element.nodes[8 + edge] = inserted.first->second;
        }
        elements.push_back(element);
    }
    return fuelsim::UnstructuredHex20Mesh(std::move(nodes), std::move(elements), std::move(blocks),
        {{1, "primary"}, {2, "secondary"}}, {},
        {{21, "primary_outer", std::move(faces[0])}, {22, "primary_contact", std::move(faces[1])},
            {23, "secondary_contact", std::move(faces[2])}, {24, "secondary_outer", std::move(faces[3])}});
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_b549_small_c3d20t_mesh <output-hex20.e>\n";
        return 2;
    }
    try {
        const fuelsim::UnstructuredHex20Mesh mesh = make_mesh();
        fuelsim::write_exodus_hex20(argv[1], mesh);
        std::cout << "b549_hex20_elements=" << mesh.elements().size() << '\n'
                  << "b549_hex20_nodes=" << mesh.nodes().size() << '\n'
                  << "b549_corner_temperature_nodes=40\n"
                  << "b549_coupled_dofs=" << 40U + 3U * mesh.nodes().size() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "B5.49 mesh generation failed: " << error.what() << '\n';
        return 1;
    }
}
