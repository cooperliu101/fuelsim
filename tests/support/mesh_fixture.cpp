#include "support/mesh_fixture.hpp"
#include <limits>
#include <stdexcept>
#include <utility>

namespace fuelsim::test {
std::size_t annular_node_id(std::size_t radial_elements, std::size_t radial_index, std::size_t axial_index) {
    if (radial_index > radial_elements) throw std::out_of_range("annular mesh radial node is out of range");
    const std::size_t radial_nodes = radial_elements + 1;
    if (axial_index > (std::numeric_limits<std::size_t>::max() - radial_index) / radial_nodes)
        throw std::length_error("annular mesh node index overflows");
    return axial_index * radial_nodes + radial_index;
}

UnstructuredQuad4Mesh make_disconnected_annular_mesh(const std::vector<AnnularBlockSpec>& blocks) {
    if (blocks.empty()) throw std::invalid_argument("annular mesh fixture requires at least one block");
    std::vector<RzPoint> nodes;
    std::vector<Quad4Element> elements;
    std::vector<std::int64_t> element_block_ids;
    std::vector<ElementBlockInfo> element_blocks;
    std::vector<SideSet> side_sets;
    element_blocks.reserve(blocks.size());
    side_sets.reserve(4 * blocks.size());
    std::int64_t next_side_set_id = 1;
    for (const AnnularBlockSpec& block : blocks) {
        if (block.id < 0 || block.name.empty())
            throw std::invalid_argument("annular mesh fixture requires a named nonnegative block");
        if (!(block.inner_radius >= 0.0) || !(block.outer_radius > block.inner_radius) || !(block.length > 0.0) ||
            block.radial_elements == 0 || block.axial_elements == 0)
            throw std::invalid_argument("annular mesh fixture has invalid geometry or element count");
        const std::size_t radial_nodes = block.radial_elements + 1;
        const std::size_t axial_nodes = block.axial_elements + 1;
        if (radial_nodes > std::numeric_limits<std::size_t>::max() / axial_nodes)
            throw std::length_error("annular mesh node count overflows");
        const std::size_t node_offset = nodes.size();
        const std::size_t element_offset = elements.size();
        for (std::size_t axial = 0; axial < axial_nodes; ++axial) {
            const double z = block.length * static_cast<double>(axial) / static_cast<double>(block.axial_elements);
            for (std::size_t radial = 0; radial < radial_nodes; ++radial) {
                const double fraction = static_cast<double>(radial) / static_cast<double>(block.radial_elements);
                nodes.push_back({block.inner_radius + (block.outer_radius - block.inner_radius) * fraction, z});
            }
        }
        for (std::size_t axial = 0; axial < block.axial_elements; ++axial) {
            for (std::size_t radial = 0; radial < block.radial_elements; ++radial) {
                const std::size_t lower_left = node_offset + annular_node_id(block.radial_elements, radial, axial);
                const std::size_t lower_right = lower_left + 1;
                const std::size_t upper_left = lower_left + radial_nodes;
                const std::size_t upper_right = upper_left + 1;
                elements.push_back({{{lower_left, lower_right, upper_right, upper_left}}});
                element_block_ids.push_back(block.id);
            }
        }
        std::vector<ElementSide> inner;
        std::vector<ElementSide> outer;
        std::vector<ElementSide> bottom;
        std::vector<ElementSide> top;
        inner.reserve(block.axial_elements);
        outer.reserve(block.axial_elements);
        bottom.reserve(block.radial_elements);
        top.reserve(block.radial_elements);
        for (std::size_t axial = 0; axial < block.axial_elements; ++axial) {
            const std::size_t row = element_offset + axial * block.radial_elements;
            inner.push_back({row, 3});
            outer.push_back({row + block.radial_elements - 1, 1});
        }
        for (std::size_t radial = 0; radial < block.radial_elements; ++radial) {
            bottom.push_back({element_offset + radial, 0});
            top.push_back({element_offset + (block.axial_elements - 1) * block.radial_elements + radial, 2});
        }
        element_blocks.push_back({block.id, block.name});
        side_sets.push_back({next_side_set_id++, block.name + "_inner", std::move(inner)});
        side_sets.push_back({next_side_set_id++, block.name + "_outer", std::move(outer)});
        side_sets.push_back({next_side_set_id++, block.name + "_bottom", std::move(bottom)});
        side_sets.push_back({next_side_set_id++, block.name + "_top", std::move(top)});
    }
    return UnstructuredQuad4Mesh(std::move(nodes), std::move(elements), std::move(element_block_ids),
        std::move(element_blocks), {}, std::move(side_sets));
}
} // namespace fuelsim::test
