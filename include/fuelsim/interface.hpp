#ifndef FUELSIM_INTERFACE_HPP
#define FUELSIM_INTERFACE_HPP

#include <array>
#include <cstddef>

#include "fuelsim/local_system.hpp"
#include "fuelsim/mesh.hpp"

namespace fuelsim {

constexpr std::size_t line2_interface_side_node_count = 2;
constexpr std::size_t line2_interface_node_count = 4;
constexpr std::size_t line2_interface_quadrature_point_count = 2;

using Line2InterfaceSideCoordinates =
    std::array<RzPoint, line2_interface_side_node_count>;

struct Line2RzHeatQuadraturePoint final {
    std::array<double, line2_interface_side_node_count> fuel_shape;
    std::array<double, line2_interface_side_node_count> cladding_shape;
    double fuel_reference_radius;
    double cladding_reference_radius;
};

struct Line2RzHeatGeometry final {
    Line2InterfaceSideCoordinates fuel_coordinates;
    Line2InterfaceSideCoordinates cladding_coordinates;
    std::array<Line2RzHeatQuadraturePoint,
               line2_interface_quadrature_point_count>
        points;
};

struct GapHeatProperties final {
    double gap_conductivity;
    double minimum_gap;
};

struct HeatQuadratureValue final {
    double gap;
    // Positive heat flux transfers energy from fuel to cladding.
    double heat_flux;
    double weighted_measure;
};

using HeatQuadratureValues =
    std::array<HeatQuadratureValue, line2_interface_quadrature_point_count>;

Line2RzHeatGeometry make_line2_rz_heat_geometry(
    const Line2InterfaceSideCoordinates& fuel_coordinates,
    const Line2InterfaceSideCoordinates& cladding_coordinates);

class Line2RzGapHeatKernel final {
  public:
    explicit Line2RzGapHeatKernel(GapHeatProperties properties);

    const GapHeatProperties& properties() const noexcept;

    // Fixed ordering:
    // [Tf0, Tf1, Tc0, Tc1, urf0, urf1, urc0, urc1,
    //  uzf0, uzf1, uzc0, uzc1].
    LocalResidual residual(const Line2RzHeatGeometry& geometry,
                           const LocalValues& state) const;
    LocalSystem linearize(const Line2RzHeatGeometry& geometry,
                          const LocalValues& state) const;
    HeatQuadratureValues quadrature_values(const Line2RzHeatGeometry& geometry,
                                           const LocalValues& state) const;

  private:
    void residual_ad(const Line2RzHeatGeometry& geometry,
                     const LocalAdValues& state, LocalAdValues& residual) const;

    GapHeatProperties _properties;
};

struct NodeToLineRzContactGeometry final {
    Line2InterfaceSideCoordinates fuel_edge_coordinates;
    Line2InterfaceSideCoordinates cladding_segment_coordinates;
    std::size_t secondary_local_node;
    bool cladding_segment_includes_upper_endpoint;
};

struct NormalContactProperties final {
    // Pressure per unit penetration, in Pa/m.
    double penalty;
};

struct ContactPointValue final {
    bool projected;
    double gap;
    double pressure;
    double tributary_area;
    double tributary_length;
    double contact_force;
};

NodeToLineRzContactGeometry make_node_to_line_rz_contact_geometry(
    const Line2InterfaceSideCoordinates& fuel_edge_coordinates,
    const Line2InterfaceSideCoordinates& cladding_segment_coordinates,
    std::size_t secondary_local_node,
    bool cladding_segment_includes_upper_endpoint);

class NodeToLineRzContactKernel final {
  public:
    explicit NodeToLineRzContactKernel(NormalContactProperties properties);

    const NormalContactProperties& properties() const noexcept;

    // The local ordering is identical to Line2RzGapHeatKernel. Only the
    // selected fuel node and the two cladding nodes receive radial residuals.
    LocalResidual residual(const NodeToLineRzContactGeometry& geometry,
                           const LocalValues& state) const;
    LocalSystem linearize(const NodeToLineRzContactGeometry& geometry,
                          const LocalValues& state) const;
    ContactPointValue value(const NodeToLineRzContactGeometry& geometry,
                            const LocalValues& state) const;

  private:
    void residual_ad(const NodeToLineRzContactGeometry& geometry,
                     const LocalAdValues& state, LocalAdValues& residual) const;

    NormalContactProperties _properties;
};

} // namespace fuelsim

#endif
