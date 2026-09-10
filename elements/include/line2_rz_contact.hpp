#pragma once
#include "cax4_types.hpp"
#include "contact_types.hpp"
#include "coordinates.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>
#include <vector>

namespace fuelsim {
constexpr std::size_t line2_interface_side_node_count = 2;

constexpr std::size_t line2_interface_quadrature_point_count = 2;

using Line2InterfaceSideCoordinates = std::array<RzPoint, line2_interface_side_node_count>;

struct Line2RzHeatQuadraturePoint final {
    std::array<double, line2_interface_side_node_count> secondary_shape, primary_shape;
    double integration_weight, normal_orientation;
    bool nodal = false;
};

// One secondary point and primary candidate; half-open internal intervals give each primary vertex one owner.
struct Line2RzHeatPointGeometry final {
    Line2InterfaceSideCoordinates secondary_coordinates, primary_coordinates;
    Line2RzHeatQuadraturePoint point;
    bool primary_segment_includes_second_endpoint;
    bool primary_segment_is_first = false;
};

// The zero-gap hint is the secondary parent centroid's signed distance along the primary base normal; a zero hint
// remains an error, otherwise the chosen normal makes motion into secondary material open the gap.
std::array<Line2RzHeatQuadraturePoint, line2_interface_quadrature_point_count> make_line2_rz_heat_quadrature(
    const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    double secondary_coordinate_lower,
    double secondary_coordinate_upper,
    double zero_gap_orientation_hint);

Line2RzHeatPointGeometry make_line2_rz_heat_point_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    const std::array<double, line2_interface_side_node_count>& secondary_shape,
    double integration_weight,
    bool primary_segment_includes_second_endpoint,
    double zero_gap_orientation_hint);

Cax4LocalResidual compute_line2_rz_gap_heat(const GapHeatProperties& properties,
    const Line2RzHeatPointGeometry& geometry,
    const Cax4LocalValues& state,
    Cax4LocalJacobian* jacobian = nullptr);

HeatQuadratureValue compute_line2_rz_gap_heat_value(const GapHeatProperties& properties,
    const Line2RzHeatPointGeometry& geometry,
    const Cax4LocalValues& state);

ContactProjectionValue compute_line2_rz_heat_projection(const Line2RzHeatPointGeometry& geometry,
    const Cax4LocalValues& state);

struct NodeToLineRzContactGeometry final {
    Line2InterfaceSideCoordinates secondary_edge_coordinates, primary_segment_coordinates;
    std::size_t secondary_local_node;
    bool primary_segment_is_first, primary_segment_includes_second_endpoint;
    double normal_orientation;
};

NodeToLineRzContactGeometry make_node_to_line_rz_contact_geometry(
    const Line2InterfaceSideCoordinates& secondary_edge_coordinates,
    const Line2InterfaceSideCoordinates& primary_segment_coordinates,
    std::size_t secondary_local_node,
    bool primary_segment_is_first,
    bool primary_segment_includes_upper_endpoint,
    double zero_gap_orientation_hint);

Cax4LocalResidual compute_node_to_line_rz_contact(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed_state,
    const ContactPointHistory& history,
    Cax4LocalJacobian* jacobian = nullptr);

ContactPointValue compute_node_to_line_rz_contact_value(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry,
    const Cax4LocalValues& state,
    const Cax4LocalValues& committed_state,
    const ContactPointHistory& history);

ContactProjectionValue compute_node_to_line_rz_contact_projection(const NodeToLineRzContactGeometry& geometry,
    const Cax4LocalValues& state);
} // namespace fuelsim
