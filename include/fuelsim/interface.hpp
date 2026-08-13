#pragma once
#include "fuelsim/mesh.hpp"
#include "fuelsim/quad4_rz.hpp"
#include <array>
#include <cstddef>
namespace fuelsim {
constexpr std::size_t line2_interface_side_node_count = 2;
constexpr std::size_t line2_interface_node_count = 4;
constexpr std::size_t line2_interface_quadrature_point_count = 2;
using Line2InterfaceSideCoordinates = std::array<RzPoint, line2_interface_side_node_count>;
struct Line2RzHeatQuadraturePoint final {
    std::array<double, line2_interface_side_node_count> secondary_shape, primary_shape;
    double integration_weight, normal_orientation;
};
struct Line2RzHeatGeometry final {
    Line2InterfaceSideCoordinates secondary_coordinates, primary_coordinates;
    std::array<Line2RzHeatQuadraturePoint, line2_interface_quadrature_point_count> points;
};
// One secondary point and primary candidate; half-open internal intervals give each primary vertex one owner.
struct Line2RzHeatPointGeometry final {
    Line2InterfaceSideCoordinates secondary_coordinates, primary_coordinates;
    Line2RzHeatQuadraturePoint point;
    bool primary_segment_includes_second_endpoint;
};
struct GapHeatProperties final {
    double gap_conductivity, minimum_gap;
};
struct HeatQuadratureValue final {
    bool projected;
    double gap;
    double heat_flux, weighted_measure;
};
// The zero-gap hint is the secondary parent centroid's signed distance along the primary base normal; a zero hint
// remains an error, otherwise the chosen normal makes motion into secondary material open the gap.
Line2RzHeatGeometry make_line2_rz_heat_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates, double zero_gap_orientation_hint);
Line2RzHeatGeometry make_line2_rz_heat_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates, double secondary_coordinate_lower,
    double secondary_coordinate_upper, double zero_gap_orientation_hint);
Line2RzHeatPointGeometry make_line2_rz_heat_point_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    const std::array<double, line2_interface_side_node_count>& secondary_shape, double integration_weight,
    bool primary_segment_includes_second_endpoint, double zero_gap_orientation_hint);
LocalResidual compute_line2_rz_gap_heat_residual(
    const GapHeatProperties& properties, const Line2RzHeatPointGeometry& geometry, const LocalValues& state);
LocalSystem compute_line2_rz_gap_heat_system(
    const GapHeatProperties& properties, const Line2RzHeatPointGeometry& geometry, const LocalValues& state);
HeatQuadratureValue compute_line2_rz_gap_heat_value(
    const GapHeatProperties& properties, const Line2RzHeatPointGeometry& geometry, const LocalValues& state);
struct NodeToLineRzContactGeometry final {
    Line2InterfaceSideCoordinates secondary_edge_coordinates, primary_segment_coordinates;
    std::size_t secondary_local_node;
    bool primary_segment_is_first, primary_segment_includes_second_endpoint;
    double normal_orientation, reference_primary_fraction;
};
struct NormalContactProperties final {
    double penalty;
    double friction_coefficient = 0.0;
    bool augmented_lagrangian = false;
};
struct ContactPointHistory final {
    double elastic_tangential_slip = 0.0;
    bool sliding = false;
    double normal_multiplier = 0.0;
};
struct ContactPointValue final {
    bool projected;
    double gap, pressure, tributary_area, tributary_length, contact_force, tangential_traction, tangential_force,
        elastic_tangential_slip;
    bool sliding;
};
NodeToLineRzContactGeometry make_node_to_line_rz_contact_geometry(
    const Line2InterfaceSideCoordinates& secondary_edge_coordinates,
    const Line2InterfaceSideCoordinates& primary_segment_coordinates, std::size_t secondary_local_node,
    bool primary_segment_is_first, bool primary_segment_includes_upper_endpoint, double zero_gap_orientation_hint);
LocalResidual compute_node_to_line_rz_contact_residual(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history);
LocalSystem compute_node_to_line_rz_contact_system(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history);
ContactPointValue compute_node_to_line_rz_contact_value(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history);
ContactPointHistory compute_node_to_line_rz_contact_trial_history(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history);
} // namespace fuelsim
