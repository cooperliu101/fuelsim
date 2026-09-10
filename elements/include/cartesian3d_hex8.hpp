#pragma once
#include "coordinates.hpp"
#include "kinematics.hpp"
#include "material.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>
#include <vector>

namespace fuelsim {
inline constexpr std::size_t hex8_node_count = 8;
inline constexpr std::size_t hex8_local_dof_count = 32;
inline constexpr std::size_t quad4_face_node_count = 4;
inline constexpr std::size_t quad4_face_local_dof_count = 16;
using Hex8LocalDofs = std::array<std::size_t, hex8_local_dof_count>;
using Hex8LocalValues = std::array<double, hex8_local_dof_count>;
using Hex8LocalResidual = std::array<double, hex8_local_dof_count>;
using Hex8LocalJacobian = std::array<double, hex8_local_dof_count * hex8_local_dof_count>;
using Hex8LocalAdValues = std::array<adlite::Scalar, hex8_local_dof_count>;
using Quad4FaceLocalDofs = std::array<std::size_t, quad4_face_local_dof_count>;
using Quad4FaceLocalValues = std::array<double, quad4_face_local_dof_count>;
using Quad4FaceLocalResidual = std::array<double, quad4_face_local_dof_count>;
using Quad4FaceLocalJacobian = std::array<double, quad4_face_local_dof_count * quad4_face_local_dof_count>;
using Quad4FaceLocalAdValues = std::array<adlite::Scalar, quad4_face_local_dof_count>;
using Hex8Coordinates = std::array<CartesianPoint3, hex8_node_count>;
using Quad4FaceCoordinates = std::array<CartesianPoint3, quad4_face_node_count>;

struct Hex8QuadraturePoint final {
    std::array<double, hex8_node_count> shape;
    std::array<std::array<double, 3>, hex8_node_count> gradient;
    CartesianPoint3 position;
    double weighted_measure;
};

struct Hex8CapacityPoint final {
    CartesianPoint3 position;
    double weighted_measure;
};

struct Hex8Geometry final {
    std::array<Hex8QuadraturePoint, 8> points;
    std::array<Hex8CapacityPoint, hex8_node_count> capacity_points;
    Hex8QuadraturePoint reduced_point;
    std::array<Hex8CapacityPoint, hex8_node_count> reduced_capacity_points;
    std::array<std::array<double, 3>, hex8_node_count> average_shape_gradient;
    std::array<std::array<double, 4>, hex8_node_count> hourglass_shape;
    std::array<double, 4> thermal_hourglass_coefficients;
    std::array<double, 3> mechanical_hourglass_metrics;
    CartesianPoint3 selective_position;
    double reference_volume, reduced_body_source_measure;
};

struct Quad4FaceQuadraturePoint final {
    std::array<double, quad4_face_node_count> shape;
    std::array<double, quad4_face_node_count> derivative_xi, derivative_eta;
    CartesianPoint3 tangent_xi, tangent_eta;
    std::array<double, quad4_face_node_count> normal_derivative_xi, normal_derivative_eta;
    CartesianPoint3 normal_tangent_xi, normal_tangent_eta;
    double weighted_measure;
};

struct Quad4FaceGeometry final {
    std::array<Quad4FaceQuadraturePoint, 4> points;
    std::array<Quad4FaceQuadraturePoint, 4> thermal_points;
};

Hex8Geometry make_hex8_geometry(const Hex8Coordinates& coordinates);
Quad4FaceQuadraturePoint make_quad4_face_quadrature_point(const Quad4FaceCoordinates& coordinates,
    double xi,
    double eta,
    double quadrature_weight);
Quad4FaceGeometry make_quad4_face_geometry(const Quad4FaceCoordinates& coordinates);

struct CartesianThermoelasticData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source, time;
    StrainFormulation strain_formulation = StrainFormulation::small;
    Hex8ElementFormulation hex8_element_formulation = Hex8ElementFormulation::c3d8t;
    double initial_temperature = 0.0;
};

struct CartesianKinematics final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> current_gradient;
    adlite::Scalar current_weighted_measure{0.0};
};

CartesianKinematics evaluate_cartesian_incremental_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalAdValues& current_state,
    const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation);
void validate_cartesian_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state);
using CartesianMaterialHistory = std::vector<CartesianMaterialPointState>;
Hex8LocalResidual compute_hex8_thermoelastic(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues* committed_state = nullptr,
    double time_step = 0.0,
    Hex8LocalJacobian* jacobian = nullptr);
Hex8LocalResidual compute_hex8_transient(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step,
    Hex8LocalJacobian* jacobian = nullptr,
    bool include_thermal_time_term = true);
CartesianMaterialHistory compute_hex8_transient_update(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state,
    const Hex8LocalValues& committed_state,
    const CartesianMaterialHistory& committed_material,
    double time_step);
double compute_hex8_mechanical_hourglass_energy(const CartesianThermoelasticData& data,
    const Hex8Geometry& geometry,
    const Hex8LocalValues& state);
std::array<SymmetricTensor3Values, 8>
compute_hex8_stress(const CartesianThermoelasticData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state);
enum class CartesianTractionComponent { x, y, z };
enum class Quad4FaceBoundaryKind { pressure, traction, surface_heat_flux, convection };

struct Quad4FaceBoundaryData final {
    Quad4FaceBoundaryKind kind;
    CartesianTractionComponent component;
    double load, ambient_temperature;
    bool use_displaced_geometry = false;
};

Quad4FaceLocalResidual compute_quad4_face_boundary(const Quad4FaceBoundaryData& data,
    const Quad4FaceGeometry& geometry,
    const Quad4FaceLocalValues& state,
    Quad4FaceLocalJacobian* jacobian = nullptr);
} // namespace fuelsim
