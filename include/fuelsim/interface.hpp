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
    std::array<double, line2_interface_side_node_count> secondary_shape;
    std::array<double, line2_interface_side_node_count> primary_shape;
    double integration_weight;
    double normal_orientation;
};

struct Line2RzHeatGeometry final {
    Line2InterfaceSideCoordinates secondary_coordinates;
    Line2InterfaceSideCoordinates primary_coordinates;
    bool radial_reference_geometry;
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
    // Positive heat flux transfers energy from secondary to primary.
    double heat_flux;
    double weighted_measure;
};

using HeatQuadratureValues =
    std::array<HeatQuadratureValue, line2_interface_quadrature_point_count>;

Line2RzHeatGeometry make_line2_rz_heat_geometry(
    const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates);

Line2RzHeatGeometry make_line2_rz_heat_geometry(
    const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    double secondary_coordinate_lower, double secondary_coordinate_upper);

class Line2RzGapHeatKernel final {
  public:
    explicit Line2RzGapHeatKernel(GapHeatProperties properties);

    const GapHeatProperties& properties() const noexcept;

    // Fixed ordering:
    // [Ts0, Ts1, Tp0, Tp1, urs0, urs1, urp0, urp1,
    //  uzs0, uzs1, uzp0, uzp1].
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
    Line2InterfaceSideCoordinates secondary_edge_coordinates;
    Line2InterfaceSideCoordinates primary_segment_coordinates;
    std::size_t secondary_local_node;
    bool primary_segment_is_first;
    bool primary_segment_includes_second_endpoint;
    double normal_orientation;
    double reference_primary_fraction;
    bool radial_reference_geometry;
};

struct NormalContactProperties final {
    // Pressure per unit penetration, in Pa/m.
    double penalty;
    // Coulomb coefficient. The tangential penalty equals the normal penalty.
    double friction_coefficient = 0.0;
    bool augmented_lagrangian = false;
};

struct ContactPointHistory final {
    double elastic_tangential_slip = 0.0;
    bool sliding = false;
    // Compressive normal traction carried between augmented outer iterations.
    double normal_multiplier = 0.0;
};

struct ContactPointValue final {
    bool projected;
    double gap;
    double pressure;
    double tributary_area;
    double tributary_length;
    double contact_force;
    double tangential_traction;
    double tangential_force;
    double elastic_tangential_slip;
    bool sliding;
};

NodeToLineRzContactGeometry make_node_to_line_rz_contact_geometry(
    const Line2InterfaceSideCoordinates& secondary_edge_coordinates,
    const Line2InterfaceSideCoordinates& primary_segment_coordinates,
    std::size_t secondary_local_node,
    bool primary_segment_is_first,
    bool primary_segment_includes_upper_endpoint);

class NodeToLineRzContactKernel final {
  public:
    explicit NodeToLineRzContactKernel(NormalContactProperties properties);

    const NormalContactProperties& properties() const noexcept;

    // The local ordering is identical to Line2RzGapHeatKernel. The selected
    // secondary node and the two primary nodes receive radial and axial
    // residuals along the current primary-segment normal.
    LocalResidual residual(const NodeToLineRzContactGeometry& geometry,
                           const LocalValues& state,
                           const LocalValues& committed_state,
                           const ContactPointHistory& history) const;
    LocalSystem linearize(const NodeToLineRzContactGeometry& geometry,
                          const LocalValues& state,
                          const LocalValues& committed_state,
                          const ContactPointHistory& history) const;
    ContactPointValue value(const NodeToLineRzContactGeometry& geometry,
                            const LocalValues& state,
                            const LocalValues& committed_state,
                            const ContactPointHistory& history) const;
    ContactPointHistory trial_history(
        const NodeToLineRzContactGeometry& geometry,
        const LocalValues& state, const LocalValues& committed_state,
        const ContactPointHistory& history) const;

  private:
    void residual_ad(const NodeToLineRzContactGeometry& geometry,
                     const LocalAdValues& state,
                     const LocalValues& committed_state,
                     const ContactPointHistory& history,
                     LocalAdValues& residual) const;

    NormalContactProperties _properties;
};

} // namespace fuelsim

#endif
