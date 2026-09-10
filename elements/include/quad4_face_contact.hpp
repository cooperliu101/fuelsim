#pragma once
#include "contact_types.hpp"
#include "coordinates.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>
#include <vector>

namespace fuelsim {
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
    double gap, distance, xi, eta;
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

struct Quad4AveragedFrictionGeometryValue final {
    bool projected;
    double area;
    std::array<double, 3> area_normal, area_tangent_first, separation;
    std::array<double, 2> tangential_increment;
};

using Quad4AveragedFrictionGeometryJacobian = std::array<double, 12 * quad4_surface_contact_local_dof_count>;

struct Quad4NormalForceAreaValue final {
    bool projected;
    double force, area;
};

using Quad4NormalForceAreaJacobian = std::array<double, 2 * quad4_surface_contact_local_dof_count>;

struct Quad4FiniteRegionNormalGeometryValue final {
    bool projected;
    double area, gap_integral;
    Quad4SurfaceContactLocalResidual unit_pressure_residual;
};

inline constexpr std::size_t quad4_finite_region_normal_geometry_output_count =
    2 + quad4_surface_contact_local_dof_count;

using Quad4FiniteRegionNormalGeometryJacobian =
    std::array<double, quad4_finite_region_normal_geometry_output_count * quad4_surface_contact_local_dof_count>;

Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_gap_heat(const GapHeatProperties& properties,
    const Quad4ToQuad4HeatGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    Quad4SurfaceContactLocalJacobian* jacobian = nullptr);

CartesianHeatQuadratureValue compute_quad4_to_quad4_gap_heat_value(const GapHeatProperties& properties,
    const Quad4ToQuad4HeatGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state);

ContactProjectionValue compute_quad4_to_quad4_heat_projection(const Quad4ToQuad4HeatGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state);

Quad4ReferenceProjectionValue compute_quad4_reference_projection(
    const std::array<CartesianPoint3, 4>& secondary_coordinates,
    const std::array<CartesianPoint3, 4>& primary_coordinates,
    const std::array<double, 4>& secondary_shape,
    double normal_orientation);

Quad4ReferenceProjectionValue compute_quad4_reference_closest_projection(
    const std::array<CartesianPoint3, 4>& secondary_coordinates,
    const std::array<CartesianPoint3, 4>& primary_coordinates,
    const std::array<double, 4>& secondary_shape,
    double normal_orientation);

Quad4SurfaceContactLocalResidual compute_node_to_quad4_contact(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history,
    Quad4SurfaceContactLocalJacobian* jacobian = nullptr);

CartesianContactPointValue compute_node_to_quad4_contact_value(const NormalContactProperties& properties,
    const NodeToQuad4ContactGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history);

ContactProjectionValue compute_node_to_quad4_contact_projection(const NodeToQuad4ContactGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state);

Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_contact(const NormalContactProperties& properties,
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history,
    Quad4SurfaceContactLocalJacobian* jacobian = nullptr);

Quad4SurfaceContactLocalResidual compute_quad4_to_quad4_tangential_force_geometry(
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const std::array<double, 4>& secondary_distribution,
    double tangent_orientation,
    double first_traction,
    double second_traction,
    const Quad4SurfaceContactLocalValues& state,
    Quad4SurfaceContactLocalJacobian* jacobian = nullptr);

Quad4AveragedFrictionGeometryValue compute_quad4_averaged_friction_geometry_value(
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const std::array<double, 4>& secondary_distribution,
    double tangent_orientation,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    Quad4AveragedFrictionGeometryJacobian* jacobian = nullptr);

Quad4NormalForceAreaValue compute_quad4_to_quad4_normal_force_area(const NormalContactProperties& properties,
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    Quad4NormalForceAreaJacobian* jacobian = nullptr);

Quad4FiniteRegionNormalGeometryValue compute_quad4_finite_region_normal_geometry(
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    Quad4FiniteRegionNormalGeometryJacobian* jacobian = nullptr);

CartesianContactPointValue compute_quad4_to_quad4_contact_value(const NormalContactProperties& properties,
    const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state,
    const Quad4SurfaceContactLocalValues& committed_state,
    const ContactPointHistory& history);

ContactProjectionValue compute_quad4_to_quad4_contact_projection(const Quad4ToQuad4MechanicalGeometry& geometry,
    const Quad4SurfaceContactLocalValues& state);
} // namespace fuelsim
