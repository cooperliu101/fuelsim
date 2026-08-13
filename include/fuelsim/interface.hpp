#ifndef FUELSIM_INTERFACE_HPP
#define FUELSIM_INTERFACE_HPP
#include "fuelsim/local_system.hpp"
#include "fuelsim/mesh.hpp"
#include <array>
#include <cstddef>
namespace fuelsim {
constexpr std::size_t line2_interface_side_node_count = 2;
constexpr std::size_t line2_interface_node_count = 4;
constexpr std::size_t line2_interface_quadrature_point_count = 2;
using Line2InterfaceSideCoordinates = std::array<RzPoint, line2_interface_side_node_count>;
struct Line2RzHeatQuadraturePoint final {
    std::array<double, line2_interface_side_node_count> secondary_shape;
    std::array<double, line2_interface_side_node_count> primary_shape;
    double integration_weight;
    double normal_orientation;
};
struct Line2RzHeatGeometry final {
    Line2InterfaceSideCoordinates secondary_coordinates;
    Line2InterfaceSideCoordinates primary_coordinates;
    std::array<Line2RzHeatQuadraturePoint, line2_interface_quadrature_point_count> points;
};
// One secondary-side integration point paired with one candidate primary
// segment. The primary interval is half open except for the final segment in
// the complete chain, so a point on an internal primary vertex has one owner.
struct Line2RzHeatPointGeometry final {
    Line2InterfaceSideCoordinates secondary_coordinates;
    Line2InterfaceSideCoordinates primary_coordinates;
    Line2RzHeatQuadraturePoint point;
    bool primary_segment_includes_second_endpoint;
};
struct GapHeatProperties final {
    double gap_conductivity;
    double minimum_gap;
};
struct HeatQuadratureValue final {
    bool projected;
    double gap;
    // Positive heat flux transfers energy from secondary to primary.
    double heat_flux;
    double weighted_measure;
};
// zero_gap_orientation_hint is consulted only when a secondary point rides
// exactly on its primary segment (zero reference normal gap). It carries the
// signed distance of the secondary material (the parent-element centroid)
// from the primary line along the base normal (tangent_z, -tangent_r)/length,
// and the normal orientation is chosen so that moving into the secondary
// material opens the gap. A zero hint means no material-side information is
// available; a zero gap on the segment then remains an error.
Line2RzHeatGeometry make_line2_rz_heat_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates, double zero_gap_orientation_hint);
Line2RzHeatGeometry make_line2_rz_heat_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates, double secondary_coordinate_lower,
    double secondary_coordinate_upper, double zero_gap_orientation_hint);
Line2RzHeatPointGeometry make_line2_rz_heat_point_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    const std::array<double, line2_interface_side_node_count>& secondary_shape, double integration_weight,
    bool primary_segment_includes_second_endpoint, double zero_gap_orientation_hint);
class Line2RzGapHeatKernel final {
  public:
    explicit Line2RzGapHeatKernel(GapHeatProperties properties);
    const GapHeatProperties& properties() const noexcept;
    // Fixed ordering:
    // [Ts0, Ts1, Tp0, Tp1, urs0, urs1, urp0, urp1,
    //  uzs0, uzs1, uzp0, uzp1].
    LocalResidual residual(const Line2RzHeatPointGeometry& geometry, const LocalValues& state) const;
    LocalSystem linearize(const Line2RzHeatPointGeometry& geometry, const LocalValues& state) const;
    HeatQuadratureValue quadrature_value(const Line2RzHeatPointGeometry& geometry, const LocalValues& state) const;

  private:
    void residual_ad(
        const Line2RzHeatPointGeometry& geometry, const LocalAdValues& state, LocalAdValues& residual) const;
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
// zero_gap_orientation_hint follows the same contract as the heat-geometry
// constructors above: it is consulted only when the secondary node rides
// exactly on the primary segment, and a zero hint keeps that case an error.
NodeToLineRzContactGeometry make_node_to_line_rz_contact_geometry(
    const Line2InterfaceSideCoordinates& secondary_edge_coordinates,
    const Line2InterfaceSideCoordinates& primary_segment_coordinates, std::size_t secondary_local_node,
    bool primary_segment_is_first, bool primary_segment_includes_upper_endpoint, double zero_gap_orientation_hint);
class NodeToLineRzContactKernel final {
  public:
    explicit NodeToLineRzContactKernel(NormalContactProperties properties);
    const NormalContactProperties& properties() const noexcept;
    // The local ordering is identical to Line2RzGapHeatKernel. The selected
    // secondary node and the two primary nodes receive radial and axial
    // residuals along the current primary-segment normal.
    LocalResidual residual(const NodeToLineRzContactGeometry& geometry, const LocalValues& state,
        const LocalValues& committed_state, const ContactPointHistory& history) const;
    LocalSystem linearize(const NodeToLineRzContactGeometry& geometry, const LocalValues& state,
        const LocalValues& committed_state, const ContactPointHistory& history) const;
    ContactPointValue value(const NodeToLineRzContactGeometry& geometry, const LocalValues& state,
        const LocalValues& committed_state, const ContactPointHistory& history) const;
    ContactPointHistory trial_history(const NodeToLineRzContactGeometry& geometry, const LocalValues& state,
        const LocalValues& committed_state, const ContactPointHistory& history) const;

  private:
    void residual_ad(const NodeToLineRzContactGeometry& geometry, const LocalAdValues& state,
        const LocalValues& committed_state, const ContactPointHistory& history, LocalAdValues& residual) const;
    NormalContactProperties _properties;
};
} // namespace fuelsim
#endif
