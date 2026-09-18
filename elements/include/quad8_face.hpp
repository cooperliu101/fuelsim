#pragma once
#include "contact_types.hpp"
#include "element_types.hpp"
#include "quad4_face.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>
#include <vector>

namespace fuelsim {
inline constexpr std::size_t quad8_face_temperature_node_count = 4;

inline constexpr std::size_t quad8_face_displacement_node_count = 8;

inline constexpr std::size_t quad8_face_local_dof_count = 28;

inline constexpr std::size_t quad8_face_thermal_quadrature_point_count = 4;

inline constexpr std::size_t quad8_face_mechanical_quadrature_point_count = 9;

using Quad8FaceCoordinates = std::array<CartesianPoint3, quad8_face_displacement_node_count>;

using Quad8FaceLocalDofs = std::array<std::size_t, quad8_face_local_dof_count>;

using Quad8FaceLocalValues = std::array<double, quad8_face_local_dof_count>;

using Quad8FaceLocalResidual = std::array<double, quad8_face_local_dof_count>;

using Quad8FaceLocalJacobian = std::array<double, quad8_face_local_dof_count * quad8_face_local_dof_count>;

using Quad8FaceLocalAdValues = std::array<adlite::Scalar, quad8_face_local_dof_count>;

struct Quad8FaceThermalQuadraturePoint final {
    std::array<double, quad8_face_temperature_node_count> temperature_shape;
    std::array<double, quad8_face_displacement_node_count> displacement_shape;
    std::array<double, quad8_face_displacement_node_count> derivative_xi, derivative_eta;
    double quadrature_weight;
    double weighted_measure;
};

struct Quad8FaceMechanicalQuadraturePoint final {
    std::array<double, quad8_face_displacement_node_count> displacement_shape;
    std::array<double, quad8_face_displacement_node_count> derivative_xi, derivative_eta;
    CartesianPoint3 tangent_xi, tangent_eta;
    double quadrature_weight;
};

struct Quad8FaceGeometry final {
    std::array<Quad8FaceThermalQuadraturePoint, quad8_face_thermal_quadrature_point_count> thermal_points;
    std::array<Quad8FaceMechanicalQuadraturePoint, quad8_face_mechanical_quadrature_point_count> mechanical_points;
};

void quad8_shape_values(double xi,
    double eta,
    std::array<double, 8>& shape,
    std::array<double, 8>& derivative_xi,
    std::array<double, 8>& derivative_eta);

Quad8FaceGeometry make_quad8_face_geometry(const Quad8FaceCoordinates& coordinates);

Quad8FaceMechanicalQuadraturePoint make_quad8_face_mechanical_point(const Quad8FaceCoordinates& coordinates,
    double xi,
    double eta,
    double quadrature_weight);

Quad8FaceLocalResidual compute_quad8_face_boundary(const Quad4FaceBoundaryData& data,
    const Quad8FaceGeometry& geometry,
    const Quad8FaceLocalValues& state,
    Quad8FaceLocalJacobian* jacobian = nullptr);
} // namespace fuelsim

namespace fuelsim {
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

struct Quad8HeatPatchSample final {
    Quad8ToQuad8HeatGeometry geometry;
    // Indices in this patch's local state, including shared nodes only once.
    Quad8SurfaceContactLocalDofs local_dofs;
    // Supported faces of one finite-sliding sample share a transfer group.
    bool disk_transfer = false;
    std::size_t transfer_group = 0;
    // Clearance uses the corner's mechanical averaging rule; temperature
    // averages also include the edge-rule samples in the same quadrant.
    double gap_weight = 1.0;
    double secondary_normal_orientation = 1.0;
};

std::vector<double> compute_quad8_gap_heat_patch(const GapHeatProperties& properties,
    const std::vector<Quad8HeatPatchSample>& samples,
    const std::vector<double>& state,
    std::vector<double>* jacobian = nullptr,
    CartesianHeatQuadratureValue* value = nullptr);

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
    const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    Quad8SurfaceContactLocalJacobian* jacobian = nullptr);

CartesianHeatQuadratureValue compute_quad8_to_quad8_gap_heat_value(const GapHeatProperties& properties,
    const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state);

ContactProjectionValue compute_quad8_to_quad8_heat_projection(const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state);

std::array<adlite::Scalar, 8> compute_quad8_primary_shape_derivatives(const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    bool allow_extrapolation = false);

std::array<double, 8> compute_quad8_disk_transfer(const Quad8ToQuad8HeatGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    std::array<adlite::Scalar, 8>* derivatives = nullptr);

Quad8SurfaceContactLocalResidual compute_quad8_to_quad8_contact(const NormalContactProperties& properties,
    const Quad8ToQuad8MechanicalGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history,
    Quad8SurfaceContactLocalJacobian* jacobian = nullptr);

CartesianContactPointValue compute_quad8_to_quad8_contact_value(const NormalContactProperties& properties,
    const Quad8ToQuad8MechanicalGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history);

ContactProjectionValue compute_quad8_to_quad8_contact_projection(const Quad8ToQuad8MechanicalGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state);

Quad8ReferenceProjectionValue compute_quad8_reference_projection(
    const std::array<CartesianPoint3, 8>& secondary_coordinates,
    const std::array<CartesianPoint3, 8>& primary_coordinates,
    const std::array<double, 8>& secondary_shape,
    double normal_orientation,
    bool allow_extrapolation = false);

Quad8SurfaceContactLocalResidual compute_node_to_quad8_contact(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history,
    Quad8SurfaceContactLocalJacobian* jacobian = nullptr);

CartesianContactPointValue compute_node_to_quad8_contact_value(const NormalContactProperties& properties,
    const NodeToQuad8ContactGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state,
    const Quad8SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history);

ContactProjectionValue compute_node_to_quad8_contact_projection(const NodeToQuad8ContactGeometry& geometry,
    const Quad8SurfaceContactLocalValues& state);
} // namespace fuelsim
