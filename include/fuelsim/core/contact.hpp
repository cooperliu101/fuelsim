#pragma once
#include "fuelsim/core/contact_area_rule.hpp"
#include "fuelsim/core/mesh.hpp"
#include "fuelsim/core/rz_quad4.hpp"
#include <array>
#include <cstddef>

namespace fuelsim {
constexpr std::size_t line2_interface_side_node_count = 2;
constexpr std::size_t line2_interface_quadrature_point_count = 2;
using Line2InterfaceSideCoordinates = std::array<RzPoint, line2_interface_side_node_count>;

struct Line2RzHeatQuadraturePoint final {
    std::array<double, line2_interface_side_node_count> secondary_shape, primary_shape;
    double integration_weight, normal_orientation;
};

// One secondary point and primary candidate; half-open internal intervals give each primary vertex one owner.
struct Line2RzHeatPointGeometry final {
    Line2InterfaceSideCoordinates secondary_coordinates, primary_coordinates;
    Line2RzHeatQuadraturePoint point;
    bool primary_segment_includes_second_endpoint;
};

enum class GapHeatConductanceLaw { gas_gap, affine };

struct GapHeatProperties final {
    double gap_conductivity, minimum_gap;
    GapHeatConductanceLaw law = GapHeatConductanceLaw::gas_gap;
    double conductance = 0.0, clearance_derivative = 0.0, pressure_derivative = 0.0, temperature_derivative = 0.0,
           reference_temperature = 0.0, contact_penalty = 0.0;
};

struct HeatQuadratureValue final {
    bool projected;
    double gap;
    double heat_flux, weighted_measure;
};

struct ContactProjectionValue final {
    bool projected;
    double gap;
};

// The zero-gap hint is the secondary parent centroid's signed distance along the primary base normal; a zero hint
// remains an error, otherwise the chosen normal makes motion into secondary material open the gap.
std::array<Line2RzHeatQuadraturePoint, line2_interface_quadrature_point_count> make_line2_rz_heat_quadrature(
    const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates, double secondary_coordinate_lower,
    double secondary_coordinate_upper, double zero_gap_orientation_hint);
Line2RzHeatPointGeometry make_line2_rz_heat_point_geometry(const Line2InterfaceSideCoordinates& secondary_coordinates,
    const Line2InterfaceSideCoordinates& primary_coordinates,
    const std::array<double, line2_interface_side_node_count>& secondary_shape, double integration_weight,
    bool primary_segment_includes_second_endpoint, double zero_gap_orientation_hint);
LocalResidual compute_line2_rz_gap_heat(const GapHeatProperties& properties, const Line2RzHeatPointGeometry& geometry,
    const LocalValues& state, LocalJacobian* jacobian = nullptr);
HeatQuadratureValue compute_line2_rz_gap_heat_value(
    const GapHeatProperties& properties, const Line2RzHeatPointGeometry& geometry, const LocalValues& state);
ContactProjectionValue compute_line2_rz_heat_projection(
    const Line2RzHeatPointGeometry& geometry, const LocalValues& state);

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
    double maximum_elastic_slip = 0.0;
};

struct ContactPointHistory final {
    double elastic_tangential_slip = 0.0;
    bool sliding = false;
    double normal_multiplier = 0.0;
    std::array<double, 3> cartesian_elastic_tangential_slip{};
    bool cartesian_tangent_basis_initialized = false;
    std::array<double, 3> cartesian_contact_normal{}, cartesian_contact_tangent_first{};
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
LocalResidual compute_node_to_line_rz_contact(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history, LocalJacobian* jacobian = nullptr);
ContactPointValue compute_node_to_line_rz_contact_value(const NormalContactProperties& properties,
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state, const LocalValues& committed_state,
    const ContactPointHistory& history);
ContactProjectionValue compute_node_to_line_rz_contact_projection(
    const NodeToLineRzContactGeometry& geometry, const LocalValues& state);

inline constexpr std::size_t quad4_surface_contact_node_count = 8;
inline constexpr std::size_t quad4_surface_contact_local_dof_count = 32;
using Quad4SurfaceContactLocalDofs = std::array<std::size_t, quad4_surface_contact_local_dof_count>;
using Quad4SurfaceContactLocalValues = std::array<double, quad4_surface_contact_local_dof_count>;
using Quad4SurfaceContactLocalResidual = std::array<double, quad4_surface_contact_local_dof_count>;
using Quad4SurfaceContactLocalJacobian =
    std::array<double, quad4_surface_contact_local_dof_count * quad4_surface_contact_local_dof_count>;
using Quad4SurfaceContactLocalAdValues = std::array<adlite::Scalar, quad4_surface_contact_local_dof_count>;

struct Quad4ToQuad4HeatGeometry final {
    std::array<CartesianPoint3, 4> secondary_coordinates, primary_coordinates;
    std::array<double, 4> secondary_shape, secondary_derivative_xi, secondary_derivative_eta;
    double normal_orientation;
};

struct Quad4ReferenceProjectionValue final {
    bool projected;
    std::array<double, 4> primary_shape;
    CartesianPoint3 normal;
    double gap;
};

struct Quad4ToQuad4MechanicalGeometry final {
    std::array<CartesianPoint3, 4> secondary_coordinates, primary_coordinates;
    std::array<double, 4> secondary_shape, secondary_derivative_xi, secondary_derivative_eta;
    std::array<double, 4> secondary_normal_derivative_xi, secondary_normal_derivative_eta;
    double quadrature_weight, normal_orientation, secondary_normal_orientation;
};

struct NodeToQuad4ContactGeometry final {
    std::array<CartesianPoint3, 4> secondary_coordinates, primary_coordinates;
    std::array<std::array<double, 4>, 4> secondary_shapes, secondary_derivatives_xi, secondary_derivatives_eta;
    std::size_t secondary_local_node;
    double normal_orientation;
};

struct CartesianHeatQuadratureValue final {
    bool projected;
    double gap, heat_flux, weighted_measure;
};

struct CartesianContactPointValue final {
    bool projected;
    double gap, pressure, tributary_area, contact_force, tangential_traction, tangential_force, friction_dissipation;
    std::array<double, 3> normal, tangent_first, tangential_traction_vector, elastic_tangential_slip;
    bool sliding;
};

struct Quad4AveragedFrictionGeometryValue final {
    bool projected;
    double area;
    std::array<double, 3> area_normal, area_tangent_first, separation;
    std::array<double, 2> tangential_increment;
};

using Quad4AveragedFrictionGeometryJacobian = std::array<double, 12 * quad4_surface_contact_local_dof_count>;

Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_gap_heat(const GapHeatProperties& properties,
    const Quad4ToQuad4HeatGeometry& geometry, const Quad4SurfaceContactLocalValues& state,
    Quad4SurfaceContactLocalJacobian* jacobian = nullptr);
CartesianHeatQuadratureValue compute_quad4_to_quad4_gap_heat_value(const GapHeatProperties& properties,
    const Quad4ToQuad4HeatGeometry& geometry, const Quad4SurfaceContactLocalValues& state);
ContactProjectionValue compute_quad4_to_quad4_heat_projection(
    const Quad4ToQuad4HeatGeometry& geometry, const Quad4SurfaceContactLocalValues& state);
Quad4ReferenceProjectionValue compute_quad4_reference_projection(
    const std::array<CartesianPoint3, 4>& secondary_coordinates,
    const std::array<CartesianPoint3, 4>& primary_coordinates, const std::array<double, 4>& secondary_shape,
    double normal_orientation);
Quad4SurfaceContactLocalResidual compute_node_to_quad4_contact(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry, const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state, const ContactPointHistory& history,
    Quad4SurfaceContactLocalJacobian* jacobian = nullptr);
CartesianContactPointValue compute_node_to_quad4_contact_value(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry, const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state, const ContactPointHistory& history);
ContactProjectionValue compute_node_to_quad4_contact_projection(
    const NodeToQuad4ContactGeometry& geometry, const Quad4SurfaceContactLocalValues& state);
Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_contact(const NormalContactProperties& properties,
    const Quad4ToQuad4MechanicalGeometry& geometry, const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state, const ContactPointHistory& history,
    Quad4SurfaceContactLocalJacobian* jacobian = nullptr);
Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_tangential_force_geometry(
    const Quad4ToQuad4MechanicalGeometry& geometry, const std::array<double, 4>& secondary_distribution,
    double tangent_orientation, double first_traction, double second_traction,
    const Quad4SurfaceContactLocalValues& state, Quad4SurfaceContactLocalJacobian* jacobian = nullptr);
Quad4AveragedFrictionGeometryValue compute_quad4_averaged_friction_geometry_value(
    const Quad4ToQuad4MechanicalGeometry& geometry, const std::array<double, 4>& secondary_distribution,
    double tangent_orientation, const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state, Quad4AveragedFrictionGeometryJacobian* jacobian = nullptr);
CartesianContactPointValue compute_quad4_to_quad4_contact_value(const NormalContactProperties& properties,
    const Quad4ToQuad4MechanicalGeometry& geometry, const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state, const ContactPointHistory& history);
ContactProjectionValue compute_quad4_to_quad4_contact_projection(
    const Quad4ToQuad4MechanicalGeometry& geometry, const Quad4SurfaceContactLocalValues& state);

inline constexpr std::size_t quad8_surface_contact_temperature_node_count = 8;
inline constexpr std::size_t quad8_surface_contact_displacement_node_count = 16;
inline constexpr std::size_t quad8_surface_contact_local_dof_count =
    quad8_surface_contact_temperature_node_count + 3 * quad8_surface_contact_displacement_node_count;
inline constexpr std::size_t quad8_surface_contact_quadrature_point_count = 9;
using Quad8SurfaceContactLocalDofs = std::array<std::size_t, quad8_surface_contact_local_dof_count>;
using Quad8SurfaceContactLocalValues = std::array<double, quad8_surface_contact_local_dof_count>;
using Quad8SurfaceContactLocalResidual = std::array<double, quad8_surface_contact_local_dof_count>;
using Quad8SurfaceContactLocalJacobian =
    std::array<double, quad8_surface_contact_local_dof_count * quad8_surface_contact_local_dof_count>;
using Quad8SurfaceContactLocalAdValues = std::array<adlite::Scalar, quad8_surface_contact_local_dof_count>;

struct Quad8ToQuad8HeatGeometry final {
    std::array<CartesianPoint3, 8> secondary_coordinates, primary_coordinates;
    std::array<double, 4> secondary_temperature_shape;
    std::array<double, 8> secondary_displacement_shape, secondary_derivative_xi, secondary_derivative_eta;
    double quadrature_weight, normal_orientation;
};

struct Quad8ToQuad8MechanicalGeometry final {
    std::array<CartesianPoint3, 8> secondary_coordinates, primary_coordinates;
    std::array<double, 8> secondary_displacement_shape, secondary_derivative_xi, secondary_derivative_eta;
    std::array<double, 8> primary_displacement_shape, primary_derivative_xi, primary_derivative_eta;
    double quadrature_weight, normal_orientation;
    bool finite_sliding = false;
};

struct Quad8ReferenceProjectionValue final {
    bool projected;
    std::array<double, 8> primary_shape, primary_derivative_xi, primary_derivative_eta;
    CartesianPoint3 normal;
    double gap;
};

struct NodeToQuad8ContactGeometry final {
    std::array<CartesianPoint3, 8> secondary_coordinates, primary_coordinates;
    std::array<std::array<double, 8>, quad8_surface_contact_quadrature_point_count> secondary_shapes;
    std::array<std::array<double, 8>, quad8_surface_contact_quadrature_point_count> secondary_derivatives_xi;
    std::array<std::array<double, 8>, quad8_surface_contact_quadrature_point_count> secondary_derivatives_eta;
    std::array<double, quad8_surface_contact_quadrature_point_count> secondary_quadrature_weights;
    std::size_t secondary_local_node;
    double normal_orientation;
    Quad8NodalAreaRule nodal_area_rule = Quad8NodalAreaRule::positive_lumped;
};

Quad8SurfaceContactLocalResidual compute_quad8_to_quad8_gap_heat(const GapHeatProperties& properties,
    const Quad8ToQuad8HeatGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
    Quad8SurfaceContactLocalJacobian* jacobian = nullptr);
CartesianHeatQuadratureValue compute_quad8_to_quad8_gap_heat_value(const GapHeatProperties& properties,
    const Quad8ToQuad8HeatGeometry& geometry, const Quad8SurfaceContactLocalValues& state);
ContactProjectionValue compute_quad8_to_quad8_heat_projection(
    const Quad8ToQuad8HeatGeometry& geometry, const Quad8SurfaceContactLocalValues& state);
Quad8SurfaceContactLocalResidual compute_quad8_to_quad8_contact(const NormalContactProperties& properties,
    const Quad8ToQuad8MechanicalGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history,
    Quad8SurfaceContactLocalJacobian* jacobian = nullptr);
CartesianContactPointValue compute_quad8_to_quad8_contact_value(const NormalContactProperties& properties,
    const Quad8ToQuad8MechanicalGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history);
ContactProjectionValue compute_quad8_to_quad8_contact_projection(
    const Quad8ToQuad8MechanicalGeometry& geometry, const Quad8SurfaceContactLocalValues& state);
Quad8ReferenceProjectionValue compute_quad8_reference_projection(
    const std::array<CartesianPoint3, 8>& secondary_coordinates,
    const std::array<CartesianPoint3, 8>& primary_coordinates, const std::array<double, 8>& secondary_shape,
    double normal_orientation);
Quad8SurfaceContactLocalResidual compute_node_to_quad8_contact(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history,
    Quad8SurfaceContactLocalJacobian* jacobian = nullptr);
CartesianContactPointValue compute_node_to_quad8_contact_value(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry, const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state, const ContactPointHistory& history);
ContactProjectionValue compute_node_to_quad8_contact_projection(
    const NodeToQuad8ContactGeometry& geometry, const Quad8SurfaceContactLocalValues& state);
} // namespace fuelsim
