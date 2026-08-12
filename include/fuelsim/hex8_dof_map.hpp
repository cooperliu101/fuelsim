#ifndef FUELSIM_HEX8_DOF_MAP_HPP
#define FUELSIM_HEX8_DOF_MAP_HPP

#include "fuelsim/dof_map.hpp"
#include "fuelsim/hex8_local_system.hpp"
#include "fuelsim/nonlinear_problem.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace fuelsim {

class Hex8DofMap final {
  public:
    explicit Hex8DofMap(std::size_t node_count);

    std::size_t node_count() const noexcept;
    std::size_t dof_count() const noexcept;
    const std::vector<FieldDescriptor>& field_layout() const noexcept;

    std::size_t dof(Field field, std::size_t node) const;
    std::size_t temperature(std::size_t node) const;
    std::size_t displacement_x(std::size_t node) const;
    std::size_t displacement_y(std::size_t node) const;
    std::size_t displacement_z(std::size_t node) const;

    Hex8LocalDofs local_dofs(const std::array<std::size_t, 8>& global_nodes) const;
    Quad4FaceLocalDofs face_local_dofs(const std::array<std::size_t, 4>& global_nodes) const;

  private:
    std::size_t _node_count;
    std::vector<FieldDescriptor> _field_layout;
};

} // namespace fuelsim

#endif
