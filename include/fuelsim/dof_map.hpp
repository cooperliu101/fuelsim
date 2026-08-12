#ifndef FUELSIM_DOF_MAP_HPP
#define FUELSIM_DOF_MAP_HPP

#include <array>
#include <cstddef>
#include <vector>

#include "fuelsim/local_system.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/nonlinear_problem.hpp"

namespace fuelsim {

enum class Field {
    temperature,
    radial_displacement,
    axial_displacement,
};

class DofMap final {
  public:
    explicit DofMap(std::size_t node_count);

    std::size_t node_count() const noexcept;
    std::size_t dof_count() const noexcept;
    const std::vector<FieldDescriptor>& field_layout() const noexcept;

    std::size_t dof(Field field, std::size_t node) const;
    std::size_t temperature(std::size_t node) const;
    std::size_t radial_displacement(std::size_t node) const;
    std::size_t axial_displacement(std::size_t node) const;

    LocalDofs local_dofs(const std::array<std::size_t, 4>& global_nodes) const;

  private:
    std::size_t _node_count;
    std::vector<FieldDescriptor> _field_layout;
};

} // namespace fuelsim

#endif
