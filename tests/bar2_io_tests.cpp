#include "io/results_io.hpp"
#include <array>
#include <cmath>
#include <exodusII.h>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void change_attribute(const std::string& path, double value) {
    int cpu = 8, disk = 0;
    float version = 0;
    const int file = ex_open(path.c_str(), EX_WRITE, &cpu, &disk, &version);
    require(file >= 0, "Open BAR2 contract fixture for mutation");
    std::array<double, 4> values{value, 8.0, 8.0, 9.0};
    const int status = ex_put_attr(file, EX_ELEM_BLOCK, 1, values.data());
    const int closed = ex_close(file);
    require(status >= 0 && closed >= 0, "Write BAR2 contract fixture attribute");
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "Expected isolated BAR2 test directory");
        std::filesystem::create_directories(argv[1]);
        const std::string path = (std::filesystem::path(argv[1]) / "bar2_mesh.e").string();
        const fuelsim::UnstructuredBar2Mesh source(
            {{0, .5}, {1, .5}, {2, .5}, {0, 1.5}, {1, 1.5}, {2, 1.5}, {0, 0}, {0, 1}, {0, 2}},
            {{{0, 1}, {6, 7}}, {{1, 2}, {6, 7}}, {{3, 4}, {7, 8}}, {{4, 5}, {7, 8}}},
            {2, 1, 2, 1},
            {{1, "outer"}, {2, "inner"}},
            {{5, "axial_bottom", {6}}, {6, "axial_top", {8}}},
            {{7, "outside", {{1, 1}, {3, 1}}}, {8, "axis", {{0, 0}, {2, 0}}}});
        fuelsim::write_exodus_bar2(path, source);
        const auto result = fuelsim::read_exodus_bar2(path);
        require(result.nodes().size() == 9 && result.elements().size() == 4,
            "BAR2 preserves control nodes and element count");
        require(result.radial_source_node_ids() == source.radial_source_node_ids()
                    && result.axial_source_node_ids() == source.axial_source_node_ids(),
            "BAR2 preserves separate radial and axial node roles");
        require(result.element_block_ids() == std::vector<std::int64_t>({1, 1, 2, 2}),
            "BAR2 output orders connectivity by block");
        const std::array<std::size_t, 4> original{1, 3, 0, 2};
        for (std::size_t i = 0; i < original.size(); ++i) {
            require(result.elements()[i].nodes == source.elements()[original[i]].nodes
                        && result.elements()[i].axial_nodes == source.elements()[original[i]].axial_nodes,
                "BAR2 axial attributes follow reordered radial elements");
        }
        require(result.side_set("outside").sides[0].element == 0 && result.side_set("outside").sides[1].element == 1
                    && result.side_set("axis").sides[0].element == 2,
            "BAR2 side sets follow block reordering");
        require(result.node_sets()[0].nodes == std::vector<std::size_t>({6}),
            "BAR2 control node boundary retains source ID");
        for (double invalid :
            {0.0, 7.5, 10.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
            change_attribute(path, invalid);
            bool rejected = false;
            try {
                (void)fuelsim::read_exodus_bar2(path);
            } catch (const std::exception&) {
                rejected = true;
            }
            require(rejected, "BAR2 rejects invalid axial control node attributes");
        }
        fuelsim::write_exodus_bar2(path, source);
        std::cout << "[PASS] BAR2 Exodus node roles, attributes, sets and invalid references\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
