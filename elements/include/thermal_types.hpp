#pragma once
#include "element_types.hpp"
#include "material_functions.hpp"
#include <array>
#include <vector>

namespace fuelsim {
enum class ThermalElement { dcax4, dcax8, dc3d8, dc3d20 };
}

namespace fuelsim::elements {
struct ThermalPoint final {
    std::vector<double> shape;
    std::vector<std::array<double, 3>> gradient;
    CartesianPoint3 position{};
    double measure = 0.0;
};

struct ThermalGeometry final {
    std::size_t node_count = 0;
    std::vector<ThermalPoint> points;
    bool axisymmetric = false;
    ThermalElement element = ThermalElement::dcax4;
    std::vector<CartesianPoint3> coordinates;
};

struct ThermalInput final {
    const ThermalFunctionInstance& material;
    const ThermalGeometry& geometry;
    const std::vector<double>& temperature;
    const std::vector<double>& previous_temperature;
    double initial_temperature = 300.0;
    double time = 0.0;
    double time_step = 0.0;
    double heat_source = 0.0;
};

struct ThermalResult final {
    std::vector<double> residual, jacobian;
    std::vector<std::array<double, 3>> heat_flux;
    double generated_heat_rate = 0.0, stored_heat_rate = 0.0;
};
} // namespace fuelsim::elements
