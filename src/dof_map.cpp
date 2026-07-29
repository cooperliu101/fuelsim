#include "fuelsim/dof_map.hpp"

#include <limits>
#include <stdexcept>

namespace fuelsim {

DofMap::DofMap(std::size_t node_count) : _node_count(node_count) {
    if (node_count == 0)
        throw std::invalid_argument("DofMap node_count must be positive");
    if (node_count > std::numeric_limits<std::size_t>::max() / 3)
        throw std::length_error("DofMap DOF count overflows");
}

std::size_t DofMap::node_count() const noexcept {
    return _node_count;
}

std::size_t DofMap::dof_count() const noexcept {
    return 3 * _node_count;
}

std::size_t DofMap::dof(Field field, std::size_t node) const {
    if (node >= _node_count)
        throw std::out_of_range("DofMap node index is out of range");

    switch (field) {
    case Field::temperature:
        return node;
    case Field::radial_displacement:
        return _node_count + node;
    case Field::axial_displacement:
        return 2 * _node_count + node;
    }
    throw std::invalid_argument("Unknown DofMap field");
}

std::size_t DofMap::temperature(std::size_t node) const {
    return dof(Field::temperature, node);
}

std::size_t DofMap::radial_displacement(std::size_t node) const {
    return dof(Field::radial_displacement, node);
}

std::size_t DofMap::axial_displacement(std::size_t node) const {
    return dof(Field::axial_displacement, node);
}

LocalDofs
DofMap::local_dofs(const std::array<std::size_t, 4>& global_nodes) const {
    LocalDofs result{};
    for (std::size_t node = 0; node < global_nodes.size(); ++node) {
        result[node] = temperature(global_nodes[node]);
        result[global_nodes.size() + node] =
            radial_displacement(global_nodes[node]);
        result[2 * global_nodes.size() + node] =
            axial_displacement(global_nodes[node]);
    }
    return result;
}

LocalDofs DofMap::element_dofs(const Quad4Element& element) const {
    return local_dofs(element.nodes);
}

} // namespace fuelsim
