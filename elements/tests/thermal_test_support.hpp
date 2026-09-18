#pragma once
#include "thermal_types.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

inline void check_thermal_element(const fuelsim::elements::ThermalGeometry& geometry,
    fuelsim::elements::ThermalResult (*evaluate)(const fuelsim::elements::ThermalInput&, bool),
    const std::vector<fuelsim::CartesianPoint3>& coordinates,
    bool axisymmetric) {
    using namespace fuelsim;
    using namespace fuelsim::elements;
    ThermalFunctionInstance material;
    material.function = [](const ThermoelasticFunctionInput& input, ThermalPropertyOutput& output) {
        output.conductivity = 10.0 + 0.02 * (input.temperature - 300.0);
        output.density = 1000.0;
        output.specific_heat = 100.0 + 0.1 * (input.temperature - 300.0);
    };
    const auto count = coordinates.size();
    std::vector<double> current(count), old(count), direction(count);
    for (std::size_t n = 0; n < count; ++n) {
        current[n] = 310.0 + 23.0 * std::sin(static_cast<double>(n));
        old[n] = 300.0 + 7.0 * std::cos(static_cast<double>(n));
        direction[n] = std::cos(0.7 * static_cast<double>(n));
    }
    ThermalInput input{material, geometry, current, old, 300.0, 1.0, 0.7, 25.0};
    const auto tangent = evaluate(input, true), residual = evaluate(input, false);
    if (tangent.residual != residual.residual)
        throw std::runtime_error("Residual-only and tangent paths differ");
    std::vector<double> plus = current, minus = current;
    const double epsilon = 1e-4;
    for (std::size_t n = 0; n < count; ++n) {
        plus[n] += epsilon * direction[n];
        minus[n] -= epsilon * direction[n];
    }
    const auto rp = evaluate({material, geometry, plus, old, 300, 1, 0.7, 25}, false);
    const auto rm = evaluate({material, geometry, minus, old, 300, 1, 0.7, 25}, false);
    double error = 0, norm = 0;
    for (std::size_t i = 0; i < count; ++i) {
        double analytic = 0;
        for (std::size_t j = 0; j < count; ++j)
            analytic += tangent.jacobian[i * count + j] * direction[j];
        const double numerical = (rp.residual[i] - rm.residual[i]) / (2 * epsilon);
        error += (analytic - numerical) * (analytic - numerical);
        norm += analytic * analytic;
    }
    if (std::sqrt(error / norm) > 2e-8)
        throw std::runtime_error("Nonlinear thermal directional derivative failed");
    double sum = 0;
    for (double value : residual.residual)
        sum += value;
    if (std::abs(sum - residual.stored_heat_rate + residual.generated_heat_rate) > 1e-9 * std::max(1.0, std::abs(sum)))
        throw std::runtime_error("Thermal source and storage are not conservative");
    material.function = [](const ThermoelasticFunctionInput&, ThermalPropertyOutput& output) {
        output.conductivity = 10.0;
        output.density = 1000.0;
        output.specific_heat = 100.0;
    };
    for (std::size_t n = 0; n < count; ++n)
        current[n] =
            300.0 + 2.0 * coordinates[n].x + 3.0 * coordinates[n].y + (axisymmetric ? 0 : 4.0 * coordinates[n].z);
    const auto patch = evaluate({material, geometry, current, {}, 300, 0, 0, 0}, true);
    for (const auto& flux : patch.heat_flux)
        if (std::abs(flux[0] + 20) > 1e-9 || std::abs(flux[1] + 30) > 1e-9
            || std::abs(flux[2] + (axisymmetric ? 0 : 40)) > 1e-9)
            throw std::runtime_error("Thermal linear-field patch failed");
    const auto capacity = evaluate({material, geometry, current, old, 300, 1, 1, 0}, true);
    double mass = 0.0, off_diagonal = 0.0, volume = 0.0;
    for (const auto& point : geometry.points)
        volume += point.measure;
    for (std::size_t i = 0; i < count; ++i)
        for (std::size_t j = 0; j < count; ++j) {
            const double entry = capacity.jacobian[i * count + j] - patch.jacobian[i * count + j];
            mass += entry;
            if (i != j)
                off_diagonal += std::abs(entry);
        }
    const bool lumped = geometry.element == ThermalElement::dcax4 || geometry.element == ThermalElement::dc3d8;
    if (std::abs(mass - 100000 * volume) > 1e-9 * mass || (lumped ? off_diagonal != 0 : off_diagonal == 0))
        throw std::runtime_error("Thermal capacity structure or fixed reference mass is incorrect");
    if (geometry.element == ThermalElement::dcax8 || geometry.element == ThermalElement::dc3d20) {
        for (std::size_t n = 0; n < count; ++n)
            current[n] = 300 + 2 * coordinates[n].x * coordinates[n].x + 3 * coordinates[n].y * coordinates[n].y
                         + (axisymmetric ? 0 : 4 * coordinates[n].z * coordinates[n].z);
        const auto quadratic = evaluate({material, geometry, current, {}, 300, 0, 0, 0}, false);
        for (std::size_t q = 0; q < geometry.points.size(); ++q) {
            const auto& p = geometry.points[q].position;
            if (std::abs(quadratic.heat_flux[q][0] + 40 * p.x) > 1e-9
                || std::abs(quadratic.heat_flux[q][1] + 60 * p.y) > 1e-9
                || std::abs(quadratic.heat_flux[q][2] + (axisymmetric ? 0 : 80 * p.z)) > 1e-9)
                throw std::runtime_error("Quadratic temperature interpolation omitted midside nodes");
        }
    }
    std::cout << "thermal directional derivative and conservative linear patch passed\n";
}
