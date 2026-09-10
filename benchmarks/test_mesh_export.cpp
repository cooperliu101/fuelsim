#include "io/results_io.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace b523 {
constexpr double pi = 3.141592653589793238462643383279502884;

fuelsim::UnstructuredHex8Mesh mesh(std::size_t through_thickness_elements = 1,
    std::size_t tangential_elements = 2,
    double distortion = 0.0,
    double initial_gap = 0.0,
    bool nonmatching = false) {
    if (through_thickness_elements == 0 || tangential_elements == 0)
        throw std::invalid_argument("B5.23 mesh divisions must be positive");
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> element_blocks;
    std::array<std::vector<fuelsim::ElementSide>, 10> faces;
    const std::array<std::size_t, 2> tangential_divisions =
        nonmatching ? std::array<std::size_t, 2>{2, 1}
                    : std::array<std::size_t, 2>{tangential_elements, tangential_elements};
    const std::array<double, 2> y_lower =
                                    nonmatching ? std::array<double, 2>{0.0, 0.1} : std::array<double, 2>{0.0, 0.0},
                                y_upper =
                                    nonmatching ? std::array<double, 2>{2.0, 0.9} : std::array<double, 2>{1.0, 1.0},
                                z_lower =
                                    nonmatching ? std::array<double, 2>{-1.0, 0.1} : std::array<double, 2>{0.0, 0.0},
                                z_upper =
                                    nonmatching ? std::array<double, 2>{2.0, 0.9} : std::array<double, 2>{1.0, 1.0};
    const std::array<std::size_t, 2> block_offsets = {0,
        2 * (tangential_divisions[0] + 1) * (through_thickness_elements + 1)};
    const auto node = [&](std::size_t block, std::size_t z, std::size_t y, std::size_t x) {
        return block_offsets[block] + z * (tangential_divisions[block] + 1) * (through_thickness_elements + 1)
               + y * (through_thickness_elements + 1) + x;
    };
    for (std::size_t block = 0; block < 2; ++block)
        for (std::size_t z = 0; z < 2; ++z)
            for (std::size_t y = 0; y <= tangential_divisions[block]; ++y)
                for (std::size_t x = 0; x <= through_thickness_elements; ++x)
                    nodes.push_back({static_cast<double>(block) + (block == 1 ? initial_gap : 0.0)
                                         + static_cast<double>(x) / static_cast<double>(through_thickness_elements)
                                         + distortion
                                               * std::sin(pi * static_cast<double>(y)
                                                          / static_cast<double>(tangential_divisions[block])),
                        y_lower[block]
                            + (y_upper[block] - y_lower[block]) * static_cast<double>(y)
                                  / static_cast<double>(tangential_divisions[block]),
                        z == 0 ? z_lower[block] : z_upper[block]});
    for (std::size_t block = 0; block < 2; ++block)
        for (std::size_t y = 0; y < tangential_divisions[block]; ++y)
            for (std::size_t x = 0; x < through_thickness_elements; ++x) {
                const std::size_t element_index = elements.size();
                elements.push_back({{{node(block, 0, y, x),
                    node(block, 0, y, x + 1),
                    node(block, 0, y + 1, x + 1),
                    node(block, 0, y + 1, x),
                    node(block, 1, y, x),
                    node(block, 1, y, x + 1),
                    node(block, 1, y + 1, x + 1),
                    node(block, 1, y + 1, x)}}});
                element_blocks.push_back(block == 0 ? 1 : 2);
                if (block == 0 && x == 0)
                    faces[0].push_back({element_index, 3});
                if (block == 0 && x + 1 == through_thickness_elements)
                    faces[1].push_back({element_index, 1});
                if (block == 1 && x == 0)
                    faces[2].push_back({element_index, 3});
                if (block == 1 && x + 1 == through_thickness_elements)
                    faces[3].push_back({element_index, 1});
                if (block == 1 && x + 1 == through_thickness_elements && 2 * y < tangential_divisions[block])
                    faces[6].push_back({element_index, 1});
                if (block == 1 && x + 1 == through_thickness_elements && 2 * y >= tangential_divisions[block])
                    faces[7].push_back({element_index, 1});
                if (block == 1 && y == 0)
                    faces[4].push_back({element_index, 0});
                if (block == 1 && y + 1 == tangential_divisions[block])
                    faces[8].push_back({element_index, 2});
                if (block == 1)
                    faces[5].push_back({element_index, 4});
                if (block == 1)
                    faces[9].push_back({element_index, 5});
            }

    std::vector<std::size_t> primary_outer, secondary_outer, secondary_y0, secondary_z0;
    for (std::size_t z = 0; z < 2; ++z)
        for (std::size_t y = 0; y <= tangential_divisions[0]; ++y)
            primary_outer.push_back(node(0, z, y, 0));
    for (std::size_t z = 0; z < 2; ++z)
        for (std::size_t y = 0; y <= tangential_divisions[1]; ++y)
            secondary_outer.push_back(node(1, z, y, through_thickness_elements));
    for (std::size_t z = 0; z < 2; ++z)
        for (std::size_t x = 0; x <= through_thickness_elements; ++x)
            secondary_y0.push_back(node(1, z, 0, x));
    for (std::size_t y = 0; y <= tangential_divisions[1]; ++y)
        for (std::size_t x = 0; x <= through_thickness_elements; ++x)
            secondary_z0.push_back(node(1, 0, y, x));
    std::vector<fuelsim::ElementSide> secondary_all_surface;
    for (const std::size_t face : {2U, 3U, 4U, 5U, 8U, 9U})
        secondary_all_surface.insert(secondary_all_surface.end(), faces[face].begin(), faces[face].end());
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        std::move(elements),
        std::move(element_blocks),
        {{1, "primary"}, {2, "secondary"}},
        {{11, "primary_outer_nodes", primary_outer},
            {12, "secondary_outer_nodes", secondary_outer},
            {13, "secondary_y0", secondary_y0},
            {14, "secondary_z0", secondary_z0}},
        {{21, "primary_outer", faces[0]},
            {22, "primary_contact", faces[1]},
            {23, "secondary_contact", faces[2]},
            {24, "secondary_outer", faces[3]},
            {25, "secondary_y0_surface", faces[4]},
            {26, "secondary_z0_surface", faces[5]},
            {27, "secondary_outer_lower", faces[6]},
            {28, "secondary_outer_upper", faces[7]},
            {29, "secondary_all_surface", secondary_all_surface}});
}
} // namespace b523

namespace b527 {
constexpr double pi = 3.141592653589793238462643383279502884;

struct MeshDivisions final {
    std::size_t fuel_radial, clad_radial, angular, axial;
};

struct RegionGrid final {
    std::size_t offset, radial_nodes, angular_nodes, axial_nodes;
};

std::size_t node(const RegionGrid& grid, std::size_t radial, std::size_t angular, std::size_t axial) {
    return grid.offset + axial * grid.angular_nodes * grid.radial_nodes + angular * grid.radial_nodes + radial;
}

void append_region_nodes(std::vector<fuelsim::CartesianPoint3>& nodes,
    double inner_radius,
    double outer_radius,
    double height,
    std::size_t radial_elements,
    std::size_t angular_elements,
    std::size_t axial_elements) {
    for (std::size_t axial = 0; axial <= axial_elements; ++axial)
        for (std::size_t angular = 0; angular <= angular_elements; ++angular)
            for (std::size_t radial = 0; radial <= radial_elements; ++radial) {
                const double radius = inner_radius
                                      + (outer_radius - inner_radius) * static_cast<double>(radial)
                                            / static_cast<double>(radial_elements),
                             angle = 0.5 * pi * static_cast<double>(angular) / static_cast<double>(angular_elements),
                             z = height * static_cast<double>(axial) / static_cast<double>(axial_elements);
                nodes.push_back({radius * std::cos(angle), radius * std::sin(angle), z});
            }
}

fuelsim::UnstructuredHex8Mesh engineering_mesh(const MeshDivisions& divisions) {
    constexpr double fuel_inner_radius = 1.0e-3, fuel_outer_radius = 4.0e-3, clad_inner_radius = 4.005e-3,
                     clad_outer_radius = 4.7e-3, height = 4.0e-2;
    std::vector<fuelsim::CartesianPoint3> nodes;
    append_region_nodes(nodes,
        fuel_inner_radius,
        fuel_outer_radius,
        height,
        divisions.fuel_radial,
        divisions.angular,
        divisions.axial);
    const std::size_t fuel_node_count = nodes.size();
    append_region_nodes(nodes,
        clad_inner_radius,
        clad_outer_radius,
        height,
        divisions.clad_radial,
        divisions.angular,
        divisions.axial);
    const RegionGrid fuel{0, divisions.fuel_radial + 1, divisions.angular + 1, divisions.axial + 1},
        clad{fuel_node_count, divisions.clad_radial + 1, divisions.angular + 1, divisions.axial + 1};
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> blocks;
    std::array<std::vector<fuelsim::ElementSide>, 10> sides;
    const auto append_elements = [&](const RegionGrid& grid, std::size_t radial_elements, std::int64_t block) {
        for (std::size_t axial = 0; axial < divisions.axial; ++axial)
            for (std::size_t angular = 0; angular < divisions.angular; ++angular)
                for (std::size_t radial = 0; radial < radial_elements; ++radial) {
                    const std::size_t index = elements.size();
                    elements.push_back({{{node(grid, radial, angular, axial),
                        node(grid, radial + 1, angular, axial),
                        node(grid, radial + 1, angular + 1, axial),
                        node(grid, radial, angular + 1, axial),
                        node(grid, radial, angular, axial + 1),
                        node(grid, radial + 1, angular, axial + 1),
                        node(grid, radial + 1, angular + 1, axial + 1),
                        node(grid, radial, angular + 1, axial + 1)}}});
                    blocks.push_back(block);
                    if (block == 1 && radial == 0)
                        sides[0].push_back({index, 3});
                    if (block == 1 && radial + 1 == radial_elements)
                        sides[1].push_back({index, 1});
                    if (block == 2 && radial == 0)
                        sides[2].push_back({index, 3});
                    if (block == 2 && radial + 1 == radial_elements)
                        sides[3].push_back({index, 1});
                    const std::size_t symmetry_offset = block == 1 ? 4 : 7;
                    if (angular == 0)
                        sides[symmetry_offset].push_back({index, 0});
                    if (angular + 1 == divisions.angular)
                        sides[symmetry_offset + 1].push_back({index, 2});
                    if (axial == 0)
                        sides[symmetry_offset + 2].push_back({index, 4});
                }
    };
    append_elements(fuel, divisions.fuel_radial, 1);
    append_elements(clad, divisions.clad_radial, 2);
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        std::move(elements),
        std::move(blocks),
        {{1, "fuel"}, {2, "clad"}},
        {},
        {{10, "fuel_inner", sides[0]},
            {11, "fuel_outer", sides[1]},
            {12, "clad_inner", sides[2]},
            {13, "clad_outer", sides[3]},
            {14, "fuel_symmetry_y", sides[4]},
            {15, "fuel_symmetry_x", sides[5]},
            {16, "fuel_bottom", sides[6]},
            {17, "clad_symmetry_y", sides[7]},
            {18, "clad_symmetry_x", sides[8]},
            {19, "clad_bottom", sides[9]}});
}
} // namespace b527

namespace b59 {
constexpr double length = 2.0, half_thickness = 0.1, width = 0.25, step_time = 1.0e7;

struct CaseSpec final {
    const char* name;
    std::size_t nx, subdivisions_per_layer, nz;
    bool distorted;
};

constexpr std::array<CaseSpec, 3> cases = {CaseSpec{"coarse", 2, 1, 1, false},
    CaseSpec{"refined", 4, 2, 2, false},
    CaseSpec{"distorted", 4, 2, 2, true}};

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
                    elements.push_back({{{node(ix, iy, iz),
                        node(ix + 1, iy, iz),
                        node(ix + 1, iy + 1, iz),
                        node(ix, iy + 1, iz),
                        node(ix, iy, iz + 1),
                        node(ix + 1, iy, iz + 1),
                        node(ix + 1, iy + 1, iz + 1),
                        node(ix, iy + 1, iz + 1)}}});
                    block_ids.push_back(static_cast<std::int64_t>(region + 1));
                    if (region == 0 && iy == 0)
                        bottom.push_back({element, 0});
                    if (region == 1 && iy + 1 == ny)
                        top.push_back({element, 2});
                    if (region == 0 && ix == 0)
                        left_lower.push_back({element, 3});
                }
    }
    std::vector<std::size_t> interface_nodes, bottom_right_nodes, top_right_nodes;
    for (std::size_t iz = 0; iz <= spec.nz; ++iz) {
        for (std::size_t ix = 0; ix <= spec.nx; ++ix)
            interface_nodes.push_back(node(ix, spec.subdivisions_per_layer, iz));
        bottom_right_nodes.push_back(node(spec.nx, 0, iz));
        top_right_nodes.push_back(node(spec.nx, ny, iz));
    }
    fuelsim::UnstructuredHex8Mesh mesh(std::move(nodes),
        std::move(elements),
        std::move(block_ids),
        {{1, "lower"}, {2, "upper"}},
        {},
        {{11, "bottom", std::move(bottom)}, {12, "top", std::move(top)}, {13, "left_lower", std::move(left_lower)}});
    return {std::move(mesh),
        std::move(interface_nodes),
        std::move(bottom_right_nodes),
        std::move(top_right_nodes),
        maximum_distortion};
}

} // namespace b59

namespace b544 {
constexpr std::size_t x_nodes = 5, y_nodes = 2, z_nodes = 3;
constexpr std::size_t increments = 10;
constexpr double time_step = 1.0e5;

std::size_t node_index(std::size_t x, std::size_t y, std::size_t z) {
    return z * x_nodes * y_nodes + y * x_nodes + x;
}

fuelsim::CartesianPoint3 node_coordinate(std::size_t x, std::size_t y, std::size_t z) {
    const double x_coordinate = static_cast<double>(x);
    const double y_coordinate = static_cast<double>(y);
    const double z_coordinate = 0.5 * static_cast<double>(z);
    if (x == 0 || x + 1 == x_nodes)
        return {x_coordinate, y_coordinate, z_coordinate};
    const double alternating = x % 2 == 0 ? -1.0 : 1.0;
    return {x_coordinate + 0.06 * (2.0 * static_cast<double>(y) - 1.0) * (static_cast<double>(z) - 1.0),
        y_coordinate + 0.04 * alternating * (static_cast<double>(z) - 1.0),
        z_coordinate + 0.05 * alternating * (static_cast<double>(y) - 0.5)};
}

fuelsim::UnstructuredHex8Mesh mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (std::size_t z = 0; z < z_nodes; ++z)
        for (std::size_t y = 0; y < y_nodes; ++y)
            for (std::size_t x = 0; x < x_nodes; ++x)
                nodes.push_back(node_coordinate(x, y, z));
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<fuelsim::ElementSide> left, right, y_high, z_high;
    for (std::size_t z = 0; z + 1 < z_nodes; ++z) {
        for (std::size_t x = 0; x + 1 < x_nodes; ++x) {
            const std::size_t element = elements.size();
            elements.push_back({{{node_index(x, 0, z),
                node_index(x + 1, 0, z),
                node_index(x + 1, 1, z),
                node_index(x, 1, z),
                node_index(x, 0, z + 1),
                node_index(x + 1, 0, z + 1),
                node_index(x + 1, 1, z + 1),
                node_index(x, 1, z + 1)}}});
            if (x == 0)
                left.push_back({element, 3});
            if (x + 2 == x_nodes)
                right.push_back({element, 1});
            y_high.push_back({element, 2});
            if (z + 2 == z_nodes)
                z_high.push_back({element, 5});
        }
    }
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        elements,
        std::vector<std::int64_t>(elements.size(), 1),
        {{1, "solid"}},
        {},
        {{11, "left", std::move(left)},
            {12, "right", std::move(right)},
            {13, "y_high", std::move(y_high)},
            {14, "z_high", std::move(z_high)}});
}
} // namespace b544

fuelsim::UnstructuredHex8Mesh mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double z : {0.0, 1.0})
        for (double y : {0.0, 1.0})
            for (double x : {0.0, 1.0, 2.0})
                nodes.push_back({x, y, z});
    const auto node = [](std::size_t x, std::size_t y, std::size_t z) {
        return z * 6 + y * 3 + x;
    };
    const std::vector<fuelsim::Hex8Element> elements = {
        {{{node(0, 0, 0),
            node(1, 0, 0),
            node(1, 1, 0),
            node(0, 1, 0),
            node(0, 0, 1),
            node(1, 0, 1),
            node(1, 1, 1),
            node(0, 1, 1)}}},
        {{{node(1, 0, 0),
            node(2, 0, 0),
            node(2, 1, 0),
            node(1, 1, 0),
            node(1, 0, 1),
            node(2, 0, 1),
            node(2, 1, 1),
            node(1, 1, 1)}}},
    };
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        elements,
        {1, 1},
        {{1, "solid"}},
        {},
        {{11, "x0", {{0, 3}}},
            {12, "x2", {{1, 1}}},
            {13, "y0", {{0, 0}, {1, 0}}},
            {14, "z0", {{0, 4}, {1, 4}}},
            {15, "y1", {{0, 2}, {1, 2}}},
            {16, "z1", {{0, 5}, {1, 5}}}});
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3)
        return 2;
    try {
        if (std::string(argv[1]) == "b55")
            fuelsim::write_exodus_hex8(argv[2], mesh());
        else if (std::string(argv[1]) == "b59_coarse")
            fuelsim::write_exodus_hex8(argv[2], b59::make_mesh(b59::cases[0]).mesh);
        else if (std::string(argv[1]) == "b59_refined")
            fuelsim::write_exodus_hex8(argv[2], b59::make_mesh(b59::cases[1]).mesh);
        else if (std::string(argv[1]) == "b59_distorted")
            fuelsim::write_exodus_hex8(argv[2], b59::make_mesh(b59::cases[2]).mesh);
        else if (std::string(argv[1]) == "b526_cycle")
            fuelsim::write_exodus_hex8(argv[2], b523::mesh(2, 1, 0, 1e-4));
        else if (std::string(argv[1]) == "b526_reversal")
            fuelsim::write_exodus_hex8(argv[2], b523::mesh(2, 1));
        else if (std::string(argv[1]) == "b540")
            fuelsim::write_exodus_hex8(argv[2], b523::mesh(1, 1, 0, 5e-4, true));
        else if (std::string(argv[1]) == "b527")
            fuelsim::write_exodus_hex8(argv[2], b527::engineering_mesh({2, 1, 4, 3}));
        else if (std::string(argv[1]) == "b544")
            fuelsim::write_exodus_hex8(argv[2], b544::mesh());
        else if (std::string(argv[1]) == "b510")
            fuelsim::write_exodus_hex8(argv[2],
                fuelsim::UnstructuredHex8Mesh(
                    {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
                    {{{{0, 1, 2, 3, 4, 5, 6, 7}}}},
                    {1},
                    {{1, "solid"}},
                    {},
                    {{11, "left", {{0, 3}}}, {12, "right", {{0, 1}}}, {13, "y0", {{0, 0}}}, {14, "z0", {{0, 4}}}}));
        else
            return 2;
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
