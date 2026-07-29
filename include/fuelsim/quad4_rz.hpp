#ifndef FUELSIM_QUAD4_RZ_HPP
#define FUELSIM_QUAD4_RZ_HPP

#include <array>
#include <cstddef>

#include "fuelsim/local_system.hpp"
#include "fuelsim/material.hpp"
#include "fuelsim/mesh.hpp"

namespace fuelsim {

constexpr std::size_t quad4_node_count = 4;
constexpr std::size_t quad4_local_dof_count = local_dof_count;
constexpr std::size_t quad4_jacobian_size = local_jacobian_size;

using Quad4Coordinates = std::array<RzPoint, quad4_node_count>;

struct RzQuadraturePoint final {
    std::array<double, quad4_node_count> shape;
    std::array<double, quad4_node_count> gradient_r;
    std::array<double, quad4_node_count> gradient_z;
    double radius;
    double weighted_measure;
};

struct Quad4RzGeometry final {
    std::array<RzQuadraturePoint, 4> points;
};

Quad4RzGeometry make_quad4_rz_geometry(const Quad4Coordinates& coordinates);

class Quad4RzThermoelasticKernel final {
  public:
    Quad4RzThermoelasticKernel(IsotropicThermoelasticMaterial material,
                               double volumetric_heat_source);

    LocalResidual residual(const Quad4RzGeometry& geometry,
                           const LocalValues& state) const;

    LocalSystem linearize(const Quad4RzGeometry& geometry,
                          const LocalValues& state) const;

    std::array<AxisymmetricStressValues, 4>
    stress_values(const Quad4RzGeometry& geometry,
                  const LocalValues& state) const;

  private:
    void residual_ad(const Quad4RzGeometry& geometry,
                     const LocalAdValues& state, LocalAdValues& residual) const;

    IsotropicThermoelasticMaterial _material;
    double _volumetric_heat_source;
};

} // namespace fuelsim

#endif
