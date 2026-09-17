#pragma once
#include "element_types.hpp"
#include "material.hpp"
#include <array>

namespace fuelsim::elements {
// Four corner temperatures, eight x and eight y displacements, extension,
// and the two bounding-plane rotation parameters, in that order.
using Cpeg8Values = std::array<double, 23>;
using Cpeg8Coordinates = std::array<std::array<double, 2>, 8>;
using Cpeg8History = std::array<CartesianMaterialPointState, 9>;

struct Cpeg8Point final {
    std::array<double, 8> shape{}, gradient_x{}, gradient_y{};
    std::array<double, 4> temperature_shape{}, temperature_gradient_x{}, temperature_gradient_y{};
    double x = 0.0, y = 0.0, measure = 0.0;
};

struct Cpeg8Geometry final {
    Cpeg8Coordinates coordinates{};
    std::array<Cpeg8Point, 9> points{};
    std::array<double, 2> reference_point{};
    double thickness = 1.0;
};

struct Cpeg8Input final {
    const IsotropicThermoelasticMaterial& material;
    const Cpeg8Geometry& geometry;
    const Cpeg8Values& state;
    const Cpeg8Values& committed_state;
    const Cpeg8History* committed_history = nullptr;
    double time_step = 0.0, time = 0.0, volumetric_heat_source = 0.0;
    StrainFormulation strain_formulation = StrainFormulation::small;
    bool include_thermal_time_term = false;
    double initial_temperature = 300.0;
    std::array<double, 2> body_acceleration{};
};

struct Cpeg8Result final {
    Cpeg8Values residual{};
    std::array<double, 529> jacobian{};
    Cpeg8History history{};
    std::array<SymmetricTensor3Values, 9> stress{};
    double generated_heat_rate = 0.0, stored_heat_rate = 0.0;
    double elastic_energy_change = 0.0, plastic_dissipation_increment = 0.0, creep_dissipation_increment = 0.0;
};

Cpeg8Geometry
make_cpeg8t_geometry(const Cpeg8Coordinates& coordinates, double thickness, std::array<double, 2> reference_point = {});
Cpeg8Result evaluate_cpeg8t(const Cpeg8Input& input, ElementRequest request = {});

enum class Cpeg8BoundaryKind { pressure, traction_x, traction_y, heat_flux, convection };

struct Cpeg8BoundaryInput final {
    const Cpeg8Geometry& geometry;
    const Cpeg8Values& state;
    std::size_t side;
    Cpeg8BoundaryKind kind;
    double value, ambient = 0.0;
    bool current = false;
};

Cpeg8Result evaluate_cpeg8t_boundary(const Cpeg8BoundaryInput& input, bool jacobian);
} // namespace fuelsim::elements
