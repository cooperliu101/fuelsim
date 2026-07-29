#ifndef FUELSIM_INTERFACE_HPP
#define FUELSIM_INTERFACE_HPP

#include <array>
#include <cstddef>

#include "fuelsim/local_system.hpp"
#include "fuelsim/mesh.hpp"

namespace fuelsim {

constexpr std::size_t line2_interface_side_node_count = 2;
constexpr std::size_t line2_interface_node_count = 4;
constexpr std::size_t line2_interface_local_dof_count = local_dof_count;
constexpr std::size_t line2_interface_quadrature_point_count = 2;
constexpr std::size_t line2_interface_jacobian_size = local_jacobian_size;

using Line2InterfaceSideCoordinates =
    std::array<RzPoint, line2_interface_side_node_count>;

struct Line2RzInterfaceQuadraturePoint final {
    std::array<double, line2_interface_side_node_count> shape;
    double fuel_radius;
    double clad_radius;
    double weighted_measure;
};

struct Line2RzInterfaceGeometry final {
    std::array<Line2RzInterfaceQuadraturePoint,
               line2_interface_quadrature_point_count>
        points;
};

struct GapContactProperties final {
    double gap_conductivity;
    double minimum_gap;
    double penalty;
};

struct InterfaceQuadratureValue final {
    double gap;
    // Positive heat flux transfers energy from fuel to cladding.
    double heat_flux;
    double pressure;
};

using InterfaceQuadratureValues =
    std::array<InterfaceQuadratureValue,
               line2_interface_quadrature_point_count>;

Line2RzInterfaceGeometry make_line2_rz_interface_geometry(
    const Line2InterfaceSideCoordinates& fuel_coordinates,
    const Line2InterfaceSideCoordinates& clad_coordinates);

class Line2RzGapContactKernel final {
  public:
    explicit Line2RzGapContactKernel(GapContactProperties properties);

    const GapContactProperties& properties() const noexcept;

    // Fixed ordering:
    // [Tf0, Tf1, Tc0, Tc1, urf0, urf1, urc0, urc1,
    //  uzf0, uzf1, uzc0, uzc1].
    LocalResidual residual(const Line2RzInterfaceGeometry& geometry,
                           const LocalValues& state) const;

    LocalSystem linearize(const Line2RzInterfaceGeometry& geometry,
                          const LocalValues& state) const;

    InterfaceQuadratureValues
    quadrature_values(const Line2RzInterfaceGeometry& geometry,
                      const LocalValues& state) const;

  private:
    void residual_ad(const Line2RzInterfaceGeometry& geometry,
                     const LocalAdValues& state, LocalAdValues& residual) const;

    GapContactProperties properties_;
};

} // namespace fuelsim

#endif
