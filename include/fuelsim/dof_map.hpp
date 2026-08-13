#pragma once
#include "fuelsim/nonlinear_problem.hpp"
#include <cstddef>
#include <vector>
namespace fuelsim {
enum class Field {
    temperature,
    radial_displacement,
    axial_displacement,
    displacement_x,
    displacement_y,
    displacement_z,
};
enum class DofLayout {
    axisymmetric_rz,
    cartesian_3d,
};
class DofMap final {
  public:
    DofMap(std::size_t node_count, DofLayout layout);
    std::size_t node_count() const noexcept { return _node_count; }
    std::size_t dof_count() const noexcept { return _field_layout.size() * _node_count; }
    const std::vector<FieldDescriptor>& field_layout() const noexcept { return _field_layout; }
    std::size_t dof(Field field, std::size_t node) const;
    std::size_t temperature(std::size_t node) const { return dof(Field::temperature, node); }
    std::size_t radial_displacement(std::size_t node) const { return dof(Field::radial_displacement, node); }
    std::size_t axial_displacement(std::size_t node) const { return dof(Field::axial_displacement, node); }
    std::size_t displacement_x(std::size_t node) const { return dof(Field::displacement_x, node); }
    std::size_t displacement_y(std::size_t node) const { return dof(Field::displacement_y, node); }
    std::size_t displacement_z(std::size_t node) const { return dof(Field::displacement_z, node); }

  private:
    std::size_t _node_count;
    DofLayout _layout;
    std::vector<FieldDescriptor> _field_layout;
};
} // namespace fuelsim
