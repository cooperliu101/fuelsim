#include "fuelsim/io/results_io.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
fuelsim::UnstructuredHex8Mesh mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double z : {0.0, 1.0})
        for (double y : {0.0, 1.0})
            for (double x : {0.0, 1.0, 2.0}) nodes.push_back({x, y, z});
    const auto node = [](std::size_t x, std::size_t y, std::size_t z) { return z * 6 + y * 3 + x; };
    const std::vector<fuelsim::Hex8Element> elements = {
        {{{node(0, 0, 0), node(1, 0, 0), node(1, 1, 0), node(0, 1, 0), node(0, 0, 1), node(1, 0, 1), node(1, 1, 1),
            node(0, 1, 1)}}},
        {{{node(1, 0, 0), node(2, 0, 0), node(2, 1, 0), node(1, 1, 0), node(1, 0, 1), node(2, 0, 1), node(2, 1, 1),
            node(1, 1, 1)}}},
    };
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), elements, {1, 1}, {{1, "solid"}}, {},
        {{11, "x0", {{0, 3}}}, {12, "x2", {{1, 1}}}, {13, "y0", {{0, 0}, {1, 0}}}, {14, "z0", {{0, 4}, {1, 4}}},
            {15, "y1", {{0, 2}, {1, 2}}}, {16, "z1", {{0, 5}, {1, 5}}}});
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3 || std::string(argv[1]) != "b55") return 2;
    try {
        fuelsim::write_exodus_hex8(argv[2], mesh());
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
