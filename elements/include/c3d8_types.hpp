#pragma once
#include "element_types.hpp"
#include "material.hpp"
#include "material_types.hpp"
#include <array>
#include <vector>

namespace fuelsim {
inline constexpr std::size_t hex8_node_count = 8;

inline constexpr std::size_t hex8_local_dof_count = 32;

using Hex8LocalDofs = std::array<std::size_t, hex8_local_dof_count>;

using Hex8LocalValues = std::array<double, hex8_local_dof_count>;

using Hex8LocalResidual = std::array<double, hex8_local_dof_count>;

using Hex8LocalJacobian = std::array<double, hex8_local_dof_count * hex8_local_dof_count>;

using Hex8LocalAdValues = std::array<adlite::Scalar, hex8_local_dof_count>;

using Hex8Coordinates = std::array<CartesianPoint3, hex8_node_count>;

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

} // namespace fuelsim

namespace fuelsim::elements {
struct C3d8Input final {
    const IsotropicThermoelasticMaterial& material;
    const Hex8Geometry& geometry;
    const Hex8LocalValues& state;
    const Hex8LocalValues& committed_state;
    const CartesianMaterialHistory* committed_history = nullptr;
    double time_step = 0.0;
    double time = 0.0;
    double volumetric_heat_source = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    bool include_thermal_time_term = false;
    double initial_temperature = 600.0;
    std::array<double, 3> body_acceleration{};
};

struct C3d8Result final {
    Hex8LocalValues residual{};
    Hex8LocalJacobian jacobian{};
    CartesianMaterialHistory history;
    // Returned with a material-history update for the host energy accounting.
    double current_volume = 0.0, committed_volume = 0.0;
    std::array<std::array<double, 9>, 8> incremental_rotations{};
    std::array<SymmetricTensor3Values, 8> stress{};
};

} // namespace fuelsim::elements
