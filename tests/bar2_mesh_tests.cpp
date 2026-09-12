#include "core/mesh.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct MeshInput final {
    std::vector<fuelsim::RzPoint> nodes;
    std::vector<fuelsim::Bar2Element> elements;
    std::vector<std::int64_t> block_ids;
    std::vector<fuelsim::ElementBlockInfo> blocks;
    std::vector<fuelsim::NodeSet> node_sets;
    std::vector<fuelsim::SideSet> side_sets;
};

MeshInput two_layer_mesh() {
    return {{{0.0, 0.5},
                {0.004, 0.5},
                {0.0045, 0.5},
                {0.005, 0.5},
                {0.0, 1.5},
                {0.004, 1.5},
                {0.0045, 1.5},
                {0.005, 1.5},
                {0.0, 0.0},
                {0.0, 1.0},
                {0.0, 2.0},
                {0.0, 0.0},
                {0.0, 1.0},
                {0.0, 2.0},
                {0.002, 0.5}},
        {{{0, 14}, {8, 9}}, {{14, 1}, {8, 9}}, {{2, 3}, {11, 12}}, {{4, 5}, {9, 10}}, {{6, 7}, {12, 13}}},
        {1, 1, 2, 1, 2},
        {{1, "fuel"}, {2, "clad"}},
        {{1, "fuel_bottom", {8}}, {2, "clad_top", {13}}},
        {{1, "fuel_outer", {{1, 1}, {3, 1}}}, {2, "clad_inner", {{2, 0}, {4, 0}}}}};
}

fuelsim::UnstructuredBar2Mesh make_mesh(const MeshInput& input) {
    return {input.nodes, input.elements, input.block_ids, input.blocks, input.node_sets, input.side_sets};
}

void require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

void reject_mutation(const std::string& description, const std::function<void(MeshInput&)>& mutate) {
    auto input = two_layer_mesh();
    mutate(input);
    try {
        const auto mesh = make_mesh(input);
        (void)mesh;
    } catch (const std::invalid_argument&) {
        return;
    } catch (const std::out_of_range&) {
        return;
    }
    throw std::runtime_error("BAR2 accepted " + description);
}

void test_source_node_ownership() {
    const auto mesh = make_mesh(two_layer_mesh());
    require(mesh.radial_source_node_ids() == std::vector<std::size_t>({0, 1, 2, 3, 4, 5, 6, 7, 14}),
        "BAR2 radial field nodes are unique and ordered by source node number");
    require(mesh.axial_source_node_ids() == std::vector<std::size_t>({8, 9, 10, 11, 12, 13}),
        "BAR2 keeps coincident fuel and cladding axial controls independent");
    require(mesh.elements()[0].axial_nodes == mesh.elements()[1].axial_nodes,
        "BAR2 radial subdivision preserves the shared generalized axial motion");
    require(mesh.elements()[0].axial_nodes[1] == mesh.elements()[3].axial_nodes[0],
        "BAR2 adjacent fuel layers share one axial control");
    require(mesh.elements()[2].axial_nodes[1] == mesh.elements()[4].axial_nodes[0],
        "BAR2 adjacent cladding layers share their own axial control");
    require(mesh.nodes()[0].r == 0.0, "BAR2 permits an element endpoint on the symmetry axis");
    require(mesh.node_sets()[0].nodes == std::vector<std::size_t>({8}),
        "BAR2 preserves node sets containing only axial controls");
    require(mesh.side_set_block_id("fuel_outer") == 1 && mesh.side_set_block_id("clad_inner") == 2,
        "BAR2 radial endpoint side sets preserve region ownership");
    require(mesh.side_set("fuel_outer").sides.front().local_side == 1,
        "BAR2 side numbering selects the outer radial endpoint");

    auto input = two_layer_mesh();
    input.nodes.push_back({0.01, 0.5});
    const auto unused = make_mesh(input);
    require(unused.nodes().size() == 16 && unused.radial_source_node_ids().size() == 9
                && unused.axial_source_node_ids().size() == 6,
        "BAR2 does not create field nodes for unreferenced Exodus nodes");
}

void test_geometry_validation() {
    auto input = two_layer_mesh();
    input.nodes[0].z = std::nextafter(input.nodes[0].z, 1.0);
    require(make_mesh(input).elements().size() == 5, "BAR2 tolerates midpoint rounding");
    reject_mutation("negative radius", [](MeshInput& value) { value.nodes[0].r = -1.0e-8; });
    reject_mutation("nonfinite radial coordinate",
        [](MeshInput& value) { value.nodes[1].r = std::numeric_limits<double>::infinity(); });
    reject_mutation("nonfinite axial control",
        [](MeshInput& value) { value.nodes[8].z = std::numeric_limits<double>::quiet_NaN(); });
    reject_mutation("repeated radial node", [](MeshInput& value) { value.elements[0].nodes[1] = 0; });
    reject_mutation("reversed radial orientation", [](MeshInput& value) { value.elements[0].nodes = {{14, 0}}; });
    reject_mutation("repeated axial control", [](MeshInput& value) { value.elements[0].axial_nodes = {{8, 8}}; });
    reject_mutation("reversed axial orientation", [](MeshInput& value) { value.elements[0].axial_nodes = {{9, 8}}; });
    reject_mutation("radial nodes away from the layer midpoint", [](MeshInput& value) { value.nodes[0].z += 0.01; });
    reject_mutation("overflowing axial height", [](MeshInput& value) {
        value.nodes[8].z = -std::numeric_limits<double>::max();
        value.nodes[9].z = std::numeric_limits<double>::max();
    });
    reject_mutation("overlapping radial elements", [](MeshInput& value) {
        value.nodes.push_back({0.001, 0.5});
        value.elements[1].nodes[0] = value.nodes.size() - 1;
    });
    reject_mutation("duplicated radial volume", [](MeshInput& value) { value.elements[1] = value.elements[0]; });
}

void test_topology_validation() {
    reject_mutation("radial connectivity outside the node table",
        [](MeshInput& value) { value.elements[0].nodes[0] = value.nodes.size(); });
    reject_mutation("axial connectivity outside the node table",
        [](MeshInput& value) { value.elements[0].axial_nodes[0] = value.nodes.size(); });
    reject_mutation("one source node carrying radial and axial fields",
        [](MeshInput& value) { value.elements[0].axial_nodes[0] = 0; });
    reject_mutation("an axial chain splitting at one station",
        [](MeshInput& value) { value.elements[2].axial_nodes[0] = 8; });
    reject_mutation("two axial chains merging at one station",
        [](MeshInput& value) { value.elements[2].axial_nodes[1] = 9; });
    reject_mutation("one radial node shared by independent axial layers",
        [](MeshInput& value) { value.elements[2].nodes[0] = 1; });
    reject_mutation("a non-endpoint side number", [](MeshInput& value) { value.side_sets[0].sides[0].local_side = 2; });
    reject_mutation("a node set outside the node table",
        [](MeshInput& value) { value.node_sets[0].nodes[0] = value.nodes.size(); });
    reject_mutation("a missing element block", [](MeshInput& value) { value.block_ids[0] = 99; });
}
} // namespace

int main() {
    try {
        test_source_node_ownership();
        test_geometry_validation();
        test_topology_validation();
        std::cout << "BAR2 mesh checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "BAR2 mesh checks failed: " << error.what() << '\n';
        return 1;
    }
}
