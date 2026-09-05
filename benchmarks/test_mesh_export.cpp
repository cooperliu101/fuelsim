#include "fuelsim/io/results_io.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace b59 {
constexpr double length = 2.0, half_thickness = 0.1, width = 0.25, step_time = 1.0e7;

struct CaseSpec final {
    const char* name;
    std::size_t nx, subdivisions_per_layer, nz;
    bool distorted;
};

constexpr std::array<CaseSpec, 3> cases = {
    CaseSpec{"coarse", 2, 1, 1, false}, CaseSpec{"refined", 4, 2, 2, false}, CaseSpec{"distorted", 4, 2, 2, true}};

struct MeshData final {
    fuelsim::UnstructuredHex8Mesh mesh;
    std::vector<std::size_t> interface_nodes, bottom_right_nodes, top_right_nodes;
    double maximum_distortion;
};

MeshData make_mesh(const CaseSpec& spec) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    const std::size_t ny = 2 * spec.subdivisions_per_layer;
    std::vector<fuelsim::CartesianPoint3> nodes;
    double maximum_distortion = 0.0;
    for (std::size_t iz = 0; iz <= spec.nz; ++iz) {
        const double z = width * static_cast<double>(iz) / static_cast<double>(spec.nz);
        for (std::size_t iy = 0; iy <= ny; ++iy) {
            const double y = -half_thickness + 2.0 * half_thickness * static_cast<double>(iy) / static_cast<double>(ny);
            for (std::size_t ix = 0; ix <= spec.nx; ++ix) {
                const double x = length * static_cast<double>(ix) / static_cast<double>(spec.nx);
                fuelsim::CartesianPoint3 point{x, y, z};
                if (spec.distorted) {
                    const double sx = std::sin(pi * x / length);
                    const double sy = std::sin(pi * (y + half_thickness) / (2.0 * half_thickness));
                    const double sz = std::sin(pi * z / width);
                    point.x += 0.040 * sx * sy * (2.0 * z / width - 1.0);
                    point.y += 0.006 * sx * sy * sz;
                    point.z += 0.010 * sx * sy * sz * (y / half_thickness);
                }
                maximum_distortion = std::max(maximum_distortion,
                    std::sqrt(std::pow(point.x - x, 2) + std::pow(point.y - y, 2) + std::pow(point.z - z, 2)));
                nodes.push_back(point);
            }
        }
    }
    const auto node = [&spec, ny](std::size_t ix, std::size_t iy, std::size_t iz) {
        return iz * (ny + 1) * (spec.nx + 1) + iy * (spec.nx + 1) + ix;
    };
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> block_ids;
    std::vector<fuelsim::ElementSide> bottom, top, left_lower;
    for (std::size_t region = 0; region < 2; ++region) {
        const std::size_t first_y = region == 0 ? 0 : spec.subdivisions_per_layer;
        const std::size_t last_y = region == 0 ? spec.subdivisions_per_layer : ny;
        for (std::size_t iz = 0; iz < spec.nz; ++iz)
            for (std::size_t iy = first_y; iy < last_y; ++iy)
                for (std::size_t ix = 0; ix < spec.nx; ++ix) {
                    const std::size_t element = elements.size();
                    elements.push_back({{{node(ix, iy, iz), node(ix + 1, iy, iz), node(ix + 1, iy + 1, iz),
                        node(ix, iy + 1, iz), node(ix, iy, iz + 1), node(ix + 1, iy, iz + 1),
                        node(ix + 1, iy + 1, iz + 1), node(ix, iy + 1, iz + 1)}}});
                    block_ids.push_back(static_cast<std::int64_t>(region + 1));
                    if (region == 0 && iy == 0) bottom.push_back({element, 0});
                    if (region == 1 && iy + 1 == ny) top.push_back({element, 2});
                    if (region == 0 && ix == 0) left_lower.push_back({element, 3});
                }
    }
    std::vector<std::size_t> interface_nodes, bottom_right_nodes, top_right_nodes;
    for (std::size_t iz = 0; iz <= spec.nz; ++iz) {
        for (std::size_t ix = 0; ix <= spec.nx; ++ix)
            interface_nodes.push_back(node(ix, spec.subdivisions_per_layer, iz));
        bottom_right_nodes.push_back(node(spec.nx, 0, iz));
        top_right_nodes.push_back(node(spec.nx, ny, iz));
    }
    fuelsim::UnstructuredHex8Mesh mesh(std::move(nodes), std::move(elements), std::move(block_ids),
        {{1, "lower"}, {2, "upper"}}, {},
        {{11, "bottom", std::move(bottom)}, {12, "top", std::move(top)}, {13, "left_lower", std::move(left_lower)}});
    return {std::move(mesh), std::move(interface_nodes), std::move(bottom_right_nodes), std::move(top_right_nodes),
        maximum_distortion};
}

} // namespace b59

namespace b544 {
constexpr std::size_t x_nodes = 5, y_nodes = 2, z_nodes = 3;
constexpr std::size_t increments = 10;
constexpr double time_step = 1.0e5;

std::size_t node_index(std::size_t x, std::size_t y, std::size_t z) { return z * x_nodes * y_nodes + y * x_nodes + x; }

fuelsim::CartesianPoint3 node_coordinate(std::size_t x, std::size_t y, std::size_t z) {
    const double x_coordinate = static_cast<double>(x);
    const double y_coordinate = static_cast<double>(y);
    const double z_coordinate = 0.5 * static_cast<double>(z);
    if (x == 0 || x + 1 == x_nodes) return {x_coordinate, y_coordinate, z_coordinate};
    const double alternating = x % 2 == 0 ? -1.0 : 1.0;
    return {x_coordinate + 0.06 * (2.0 * static_cast<double>(y) - 1.0) * (static_cast<double>(z) - 1.0),
        y_coordinate + 0.04 * alternating * (static_cast<double>(z) - 1.0),
        z_coordinate + 0.05 * alternating * (static_cast<double>(y) - 0.5)};
}

fuelsim::UnstructuredHex8Mesh mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (std::size_t z = 0; z < z_nodes; ++z)
        for (std::size_t y = 0; y < y_nodes; ++y)
            for (std::size_t x = 0; x < x_nodes; ++x) nodes.push_back(node_coordinate(x, y, z));
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<fuelsim::ElementSide> left, right, y_high, z_high;
    for (std::size_t z = 0; z + 1 < z_nodes; ++z) {
        for (std::size_t x = 0; x + 1 < x_nodes; ++x) {
            const std::size_t element = elements.size();
            elements.push_back({{{node_index(x, 0, z), node_index(x + 1, 0, z), node_index(x + 1, 1, z),
                node_index(x, 1, z), node_index(x, 0, z + 1), node_index(x + 1, 0, z + 1), node_index(x + 1, 1, z + 1),
                node_index(x, 1, z + 1)}}});
            if (x == 0) left.push_back({element, 3});
            if (x + 2 == x_nodes) right.push_back({element, 1});
            y_high.push_back({element, 2});
            if (z + 2 == z_nodes) z_high.push_back({element, 5});
        }
    }
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), elements, std::vector<std::int64_t>(elements.size(), 1),
        {{1, "solid"}}, {},
        {{11, "left", std::move(left)}, {12, "right", std::move(right)}, {13, "y_high", std::move(y_high)},
            {14, "z_high", std::move(z_high)}});
}
} // namespace b544

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
    if (argc != 3) return 2;
    try {
        if (std::string(argv[1]) == "b55")
            fuelsim::write_exodus_hex8(argv[2], mesh());
        else if (std::string(argv[1]) == "b59_coarse")
            fuelsim::write_exodus_hex8(argv[2], b59::make_mesh(b59::cases[0]).mesh);
        else if (std::string(argv[1]) == "b59_refined")
            fuelsim::write_exodus_hex8(argv[2], b59::make_mesh(b59::cases[1]).mesh);
        else if (std::string(argv[1]) == "b59_distorted")
            fuelsim::write_exodus_hex8(argv[2], b59::make_mesh(b59::cases[2]).mesh);
        else if (std::string(argv[1]) == "b544")
            fuelsim::write_exodus_hex8(argv[2], b544::mesh());
        else if (std::string(argv[1]) == "b510")
            fuelsim::write_exodus_hex8(argv[2],
                fuelsim::UnstructuredHex8Mesh(
                    {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
                    {{{{0, 1, 2, 3, 4, 5, 6, 7}}}}, {1}, {{1, "solid"}}, {},
                    {{11, "left", {{0, 3}}}, {12, "right", {{0, 1}}}, {13, "y0", {{0, 0}}}, {14, "z0", {{0, 4}}}}));
        else
            return 2;
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
