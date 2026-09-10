#pragma once
#include "cartesian_types.hpp"
#include "coordinates.hpp"
#include "element_types.hpp"
#include "material.hpp"
#include <array>
#include <vector>

namespace fuelsim {
inline constexpr std::size_t hex20_temperature_node_count = 8;

inline constexpr std::size_t hex20_displacement_node_count = 20;

inline constexpr std::size_t hex20_local_dof_count = 68;

inline constexpr std::size_t hex20_thermal_quadrature_point_count = 8;

inline constexpr std::size_t hex20_mechanical_quadrature_point_count = 27;

using Hex20LocalDofs = std::array<std::size_t, hex20_local_dof_count>;

using Hex20LocalValues = std::array<double, hex20_local_dof_count>;

using Hex20LocalResidual = std::array<double, hex20_local_dof_count>;

using Hex20LocalJacobian = std::array<double, hex20_local_dof_count * hex20_local_dof_count>;

using Hex20LocalAdValues = std::array<adlite::Scalar, hex20_local_dof_count>;

using Hex20Coordinates = std::array<CartesianPoint3, hex20_displacement_node_count>;

struct Hex20ThermalQuadraturePoint final {
    std::array<double, hex20_temperature_node_count> temperature_shape;
    std::array<std::array<double, 3>, hex20_temperature_node_count> temperature_gradient;
    CartesianPoint3 position;
    double weighted_measure;
};

struct Hex20MechanicalQuadraturePoint final {
    std::array<double, hex20_displacement_node_count> displacement_shape;
    std::array<std::array<double, 3>, hex20_displacement_node_count> displacement_gradient;
    std::array<double, hex20_temperature_node_count> temperature_shape;
    std::array<std::array<double, 3>, hex20_temperature_node_count> temperature_gradient;
    std::array<std::array<double, 3>, hex20_temperature_node_count> source_displacement_gradient;
    CartesianPoint3 position;
    double weighted_measure;
    double source_weighted_measure;
};

struct Hex20Geometry final {
    std::array<Hex20ThermalQuadraturePoint, hex20_thermal_quadrature_point_count> thermal_points;
    std::vector<Hex20MechanicalQuadraturePoint> mechanical_points;
};

void validate_hex20_deformation(const Hex20MechanicalQuadraturePoint& point, const Hex20LocalValues& state);

} // namespace fuelsim

namespace fuelsim::elements {
struct C3d20Input final {
    const IsotropicThermoelasticMaterial& material;
    const Hex20Geometry& geometry;
    const Hex20LocalValues& state;
    const Hex20LocalValues& committed_state;
    const CartesianMaterialHistory* committed_history = nullptr;
    double time_step = 0.0;
    double time = 0.0;
    double volumetric_heat_source = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    bool include_thermal_time_term = false;
    double initial_temperature = 600.0;
};

struct C3d20Result final {
    Hex20LocalValues residual{};
    Hex20LocalJacobian jacobian{};
    CartesianMaterialHistory history;
    std::vector<SymmetricTensor3Values> stress{};
};

} // namespace fuelsim::elements
