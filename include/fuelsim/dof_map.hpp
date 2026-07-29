#ifndef FUELSIM_DOF_MAP_HPP
#define FUELSIM_DOF_MAP_HPP

#include <array>
#include <cstddef>

#include "fuelsim/local_system.hpp"
#include "fuelsim/mesh.hpp"

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

    std::size_t dof(Field field, std::size_t node) const;
    std::size_t temperature(std::size_t node) const;
    std::size_t radial_displacement(std::size_t node) const;
    std::size_t axial_displacement(std::size_t node) const;

    LocalDofs local_dofs(const std::array<std::size_t, 4>& global_nodes) const;
    LocalDofs element_dofs(const Quad4Element& element) const;

  private:
    std::size_t _node_count;
};

} // namespace fuelsim

#endif
