#include "exodus_fixture.hpp"
#include "fuelsim/results_io.hpp"
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <vector>
namespace {
int hex_value(char digit) {
    if (digit >= '0' && digit <= '9') return digit - '0';
    if (digit >= 'a' && digit <= 'f') return digit - 'a' + 10;
    return -1;
}
bool write_fixture(const char* path) {
    const char* hex = fuelsim::test_data::two_quad_exodus_hex;
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) return false;
    for (std::size_t index = 0; hex[index] != '\0'; index += 2) {
        if (hex[index + 1] == '\0') return false;
        const int high = hex_value(hex[index]);
        const int low = hex_value(hex[index + 1]);
        if (high < 0 || low < 0) return false;
        output.put(static_cast<char>((high << 4) | low));
    }
    return output.good();
}
bool meshes_equal(const fuelsim::UnstructuredQuad4Mesh& lhs, const fuelsim::UnstructuredQuad4Mesh& rhs) {
    if (lhs.nodes().size() != rhs.nodes().size() || lhs.elements().size() != rhs.elements().size() ||
        lhs.element_block_ids() != rhs.element_block_ids() ||
        lhs.element_blocks().size() != rhs.element_blocks().size() ||
        lhs.node_sets().size() != rhs.node_sets().size() || lhs.side_sets().size() != rhs.side_sets().size())
        return false;
    for (std::size_t node = 0; node < lhs.nodes().size(); ++node)
        if (lhs.nodes()[node].r != rhs.nodes()[node].r || lhs.nodes()[node].z != rhs.nodes()[node].z) return false;
    for (std::size_t element = 0; element < lhs.elements().size(); ++element)
        if (lhs.elements()[element].nodes != rhs.elements()[element].nodes) return false;
    for (std::size_t block = 0; block < lhs.element_blocks().size(); ++block)
        if (lhs.element_blocks()[block].id != rhs.element_blocks()[block].id ||
            lhs.element_blocks()[block].name != rhs.element_blocks()[block].name)
            return false;
    for (std::size_t set = 0; set < lhs.node_sets().size(); ++set)
        if (lhs.node_sets()[set].id != rhs.node_sets()[set].id ||
            lhs.node_sets()[set].name != rhs.node_sets()[set].name ||
            lhs.node_sets()[set].nodes != rhs.node_sets()[set].nodes)
            return false;
    for (std::size_t set = 0; set < lhs.side_sets().size(); ++set) {
        if (lhs.side_sets()[set].id != rhs.side_sets()[set].id ||
            lhs.side_sets()[set].name != rhs.side_sets()[set].name ||
            lhs.side_sets()[set].sides.size() != rhs.side_sets()[set].sides.size())
            return false;
        for (std::size_t side = 0; side < lhs.side_sets()[set].sides.size(); ++side)
            if (lhs.side_sets()[set].sides[side].element != rhs.side_sets()[set].sides[side].element ||
                lhs.side_sets()[set].sides[side].local_side != rhs.side_sets()[set].sides[side].local_side)
                return false;
    }
    return true;
}
bool run_tests(const char* path) {
    if (!write_fixture(path)) {
        std::cerr << "Could not write the independent Exodus fixture\n";
        return false;
    }
    const fuelsim::UnstructuredQuad4Mesh fixture = fuelsim::read_exodus_quad4(path);
    if (fixture.nodes().size() != 6 || fixture.elements().size() != 2 || fixture.element_block_ids().size() != 2 ||
        fixture.element_block_ids()[0] != 7 || fixture.element_block_ids()[1] != 7) {
        std::cerr << "Independent Exodus Quad4 fixture was read incorrectly\n";
        return false;
    }
    const fuelsim::UnstructuredQuad4Mesh two_block_mesh(fixture.nodes(), fixture.elements(),
        std::vector<std::int64_t>{7, 9}, std::vector<fuelsim::ElementBlockInfo>{{7, "fuel"}, {9, "clad"}},
        std::vector<fuelsim::NodeSet>{{10, "sample_nodes", {0, 3}}},
        std::vector<fuelsim::SideSet>{{20, "sample_sides", {{{0, 0}, {1, 0}}}}});
    fuelsim::write_exodus_quad4(path, two_block_mesh);
    const fuelsim::UnstructuredQuad4Mesh round_trip = fuelsim::read_exodus_quad4(path);
    if (!meshes_equal(two_block_mesh, round_trip)) {
        std::cerr << "Exodus Quad4 write/read round trip changed the mesh\n";
        return false;
    }
    const fuelsim::UnstructuredQuad4Mesh interleaved_blocks(fixture.nodes(), fixture.elements(),
        std::vector<std::int64_t>{9, 7}, std::vector<fuelsim::ElementBlockInfo>{{7, "fuel"}, {9, "clad"}},
        std::vector<fuelsim::NodeSet>{}, std::vector<fuelsim::SideSet>{{20, "clad_side", {{{0, 0}}}}});
    fuelsim::write_exodus_quad4(path, interleaved_blocks);
    const fuelsim::UnstructuredQuad4Mesh reordered = fuelsim::read_exodus_quad4(path);
    if (reordered.elements()[1].nodes != fixture.elements()[0].nodes ||
        reordered.side_set("clad_side").sides[0].element != 1) {
        std::cerr << "Exodus side set did not follow block-grouped element "
                     "numbering\n";
        return false;
    }
    const fuelsim::UnstructuredHex8Mesh hex_mesh(
        {{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
            {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}}},
        {{{{0, 1, 2, 3, 4, 5, 6, 7}}}}, {12}, {{12, "solid"}}, {{31, "fixed", {0, 3, 4, 7}}},
        {{41, "loaded", {{{0, 1}, {0, 5}}}}});
    fuelsim::write_exodus_hex8(path, hex_mesh);
    const fuelsim::UnstructuredHex8Mesh hex_round_trip = fuelsim::read_exodus_hex8(path);
    if (hex_round_trip.nodes().size() != 8 || hex_round_trip.elements().size() != 1 ||
        hex_round_trip.elements()[0].nodes != hex_mesh.elements()[0].nodes ||
        hex_round_trip.element_blocks()[0].name != "solid" ||
        hex_round_trip.node_set("fixed").nodes != hex_mesh.node_set("fixed").nodes ||
        hex_round_trip.side_set("loaded").sides.size() != 2 ||
        hex_round_trip.side_set("loaded").sides[1].local_side != 5) {
        std::cerr << "Exodus HEX8 write/read round trip changed the mesh\n";
        return false;
    }
    std::cout << "Direct Exodus Quad4 and HEX8 I/O passed\n";
    return true;
}
} // namespace
int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: fuelsim_exodus_mesh_io_tests <output.exo>\n";
        return 2;
    }
    bool passed = false;
    try {
        passed = run_tests(argv[1]);
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; }
    const int remove_error = std::remove(argv[1]);
    if (remove_error != 0) {
        std::cerr << "Failed to remove Exodus test file\n";
        return 3;
    }
    return passed ? 0 : 1;
}
