#ifndef FUELSIM_TEST_MESH_FIXTURE_HPP
#define FUELSIM_TEST_MESH_FIXTURE_HPP

#include "fuelsim/mesh.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace fuelsim::test {

struct AnnularBlockSpec final {
    std::int64_t id;
    std::string name;
    double inner_radius;
    double outer_radius;
    double length;
    std::size_t radial_elements;
    std::size_t axial_elements;
};

UnstructuredQuad4Mesh
make_disconnected_annular_mesh(const std::vector<AnnularBlockSpec>& blocks);

std::size_t annular_node_id(std::size_t radial_elements,
                            std::size_t radial_index, std::size_t axial_index);

} // namespace fuelsim::test

#endif
