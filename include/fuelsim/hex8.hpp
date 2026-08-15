#pragma once
#include "fuelsim/material.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/quad4_rz.hpp"
#include <adlite/adlite.hpp>
#include <array>
#include <cstddef>

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

struct Hex8Geometry final {
    std::array<Hex8QuadraturePoint, 8> points;
};

struct Quad4FaceQuadraturePoint final {
    std::array<double, quad4_face_node_count> shape;
    std::array<double, quad4_face_node_count> derivative_xi, derivative_eta;
    CartesianPoint3 tangent_xi, tangent_eta;
    double weighted_measure;
};

struct Quad4FaceGeometry final {
    std::array<Quad4FaceQuadraturePoint, 4> points;
};

Hex8Geometry make_hex8_geometry(const Hex8Coordinates& coordinates);
Quad4FaceGeometry make_quad4_face_geometry(const Quad4FaceCoordinates& coordinates);

struct Hex8ThermoelasticData final {
    IsotropicThermoelasticMaterial material;
    double volumetric_heat_source, time;
    StrainFormulation strain_formulation = StrainFormulation::small;
};

struct CartesianKinematics final {
    SymmetricTensor3 strain_increment;
    CartesianRotation rotation;
    std::array<std::array<adlite::Scalar, 3>, hex8_node_count> current_gradient;
    adlite::Scalar current_weighted_measure{0.0};
};

CartesianKinematics evaluate_cartesian_incremental_kinematics(const Hex8QuadraturePoint& point,
    const Hex8LocalAdValues& current_state, const Hex8LocalValues& committed_state,
    StrainFormulation strain_formulation);
void validate_cartesian_deformation(const Hex8QuadraturePoint& point, const Hex8LocalValues& state);
using Hex8MaterialHistory = std::array<CartesianMaterialPointState, 8>;
Hex8LocalResidual compute_hex8_thermoelastic(const Hex8ThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues* committed_state = nullptr, double time_step = 0.0,
    Hex8LocalJacobian* jacobian = nullptr);
Hex8LocalResidual compute_hex8_transient(const Hex8ThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues& committed_state, const Hex8MaterialHistory& committed_material,
    double time_step, Hex8LocalJacobian* jacobian = nullptr);
Hex8MaterialHistory compute_hex8_transient_update(const Hex8ThermoelasticData& data, const Hex8Geometry& geometry,
    const Hex8LocalValues& state, const Hex8LocalValues& committed_state, const Hex8MaterialHistory& committed_material,
    double time_step);
std::array<SymmetricTensor3Values, 8> compute_hex8_stress(
    const Hex8ThermoelasticData& data, const Hex8Geometry& geometry, const Hex8LocalValues& state);
enum class CartesianTractionComponent { x, y, z };
enum class Quad4FaceBoundaryKind { pressure, traction, convection };

struct Quad4FaceBoundaryData final {
    Quad4FaceBoundaryKind kind;
    CartesianTractionComponent component;
    double load, ambient_temperature;
    bool use_displaced_geometry = false;
};

Quad4FaceLocalResidual compute_quad4_face_boundary(const Quad4FaceBoundaryData& data, const Quad4FaceGeometry& geometry,
    const Quad4FaceLocalValues& state, Quad4FaceLocalJacobian* jacobian = nullptr);
} // namespace fuelsim
