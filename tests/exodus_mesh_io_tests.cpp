#include "fuelsim/exodus_mesh_io.hpp"

#include "exodus_fixture.hpp"

#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <vector>

namespace {

int hex_value(char digit) {
    if (digit >= '0' && digit <= '9')
        return digit - '0';
    if (digit >= 'a' && digit <= 'f')
        return digit - 'a' + 10;
    return -1;
}

bool write_fixture(const char* path) {
    const char* hex = fuelsim::test_data::two_quad_exodus_hex;
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output)
        return false;

    for (std::size_t index = 0; hex[index] != '\0'; index += 2) {
        if (hex[index + 1] == '\0')
            return false;
        const int high = hex_value(hex[index]);
        const int low = hex_value(hex[index + 1]);
        if (high < 0 || low < 0)
            return false;
        output.put(static_cast<char>((high << 4) | low));
    }
    return output.good();
}

bool meshes_equal(const fuelsim::UnstructuredQuad4Mesh& lhs,
                  const fuelsim::UnstructuredQuad4Mesh& rhs) {
    if (lhs.nodes().size() != rhs.nodes().size() ||
        lhs.elements().size() != rhs.elements().size() ||
        lhs.element_block_ids() != rhs.element_block_ids())
        return false;

    for (std::size_t node = 0; node < lhs.nodes().size(); ++node) {
        if (lhs.nodes()[node].r != rhs.nodes()[node].r ||
            lhs.nodes()[node].z != rhs.nodes()[node].z)
            return false;
    }
    for (std::size_t element = 0; element < lhs.elements().size(); ++element) {
        if (lhs.elements()[element].nodes != rhs.elements()[element].nodes)
            return false;
    }
    return true;
}

bool run_tests(const char* path) {
    if (!write_fixture(path)) {
        std::cerr << "Could not write the independent Exodus fixture\n";
        return false;
    }

    const fuelsim::UnstructuredQuad4Mesh fixture =
        fuelsim::ExodusMeshIo::read_quad4(path);
    if (fixture.nodes().size() != 6 || fixture.elements().size() != 2 ||
        fixture.element_block_ids().size() != 2 ||
        fixture.element_block_ids()[0] != 7 ||
        fixture.element_block_ids()[1] != 7) {
        std::cerr << "Independent Exodus Quad4 fixture was read incorrectly\n";
        return false;
    }

    const fuelsim::UnstructuredQuad4Mesh two_block_mesh(
        fixture.nodes(), fixture.elements(), std::vector<std::int64_t>{7, 9});
    fuelsim::ExodusMeshIo::write_quad4(path, two_block_mesh);
    const fuelsim::UnstructuredQuad4Mesh round_trip =
        fuelsim::ExodusMeshIo::read_quad4(path);
    if (!meshes_equal(two_block_mesh, round_trip)) {
        std::cerr << "Exodus Quad4 write/read round trip changed the mesh\n";
        return false;
    }

    std::cout << "Direct Exodus Quad4 I/O: nodes=6 elements=2 blocks=7,9\n";
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
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
    }

    const int remove_error = std::remove(argv[1]);
    if (remove_error != 0) {
        std::cerr << "Failed to remove Exodus test file\n";
        return 3;
    }
    return passed ? 0 : 1;
}
