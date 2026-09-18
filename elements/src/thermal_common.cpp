#include "thermal_common.hpp"
#include "c3d_common.hpp"
#include "quad8_face.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fuelsim::elements::thermal_detail {
ThermalGeometry volume_geometry(const std::vector<CartesianPoint3>& coordinates, ThermalElement element) {
    const bool axisymmetric = element == ThermalElement::dcax4 || element == ThermalElement::dcax8;
    const bool quadratic = element == ThermalElement::dcax8 || element == ThermalElement::dc3d20;
    const std::size_t count = axisymmetric ? (quadratic ? 8U : 4U) : (quadratic ? 20U : 8U);
    if (coordinates.size() != count)
        throw std::invalid_argument("Thermal element coordinate count mismatch");
    const std::vector<double> points = quadratic ? std::vector<double>{-std::sqrt(0.6), 0.0, std::sqrt(0.6)}
                                                 : std::vector<double>{-std::sqrt(1.0 / 3.0), std::sqrt(1.0 / 3.0)};
    const std::vector<double> weights =
        quadratic ? std::vector<double>{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0} : std::vector<double>{1.0, 1.0};
    ThermalGeometry geometry{count, {}, axisymmetric, element, coordinates};
    for (std::size_t k = 0; k < (axisymmetric ? 1U : points.size()); ++k)
        for (std::size_t j = 0; j < points.size(); ++j)
            for (std::size_t i = 0; i < points.size(); ++i) {
                const double x = points[i], y = points[j], z = axisymmetric ? 0.0 : points[k];
                ThermalPoint point;
                point.shape.resize(count);
                point.gradient.resize(count);
                std::vector<std::array<double, 3>> derivative(count);
                if (!axisymmetric && quadratic) {
                    std::array<double, 20> shape{};
                    std::array<std::array<double, 3>, 20> gradients{};
                    cartesian_detail::evaluate_hex20_shapes(x, y, z, shape, gradients);
                    std::copy(shape.begin(), shape.end(), point.shape.begin());
                    std::copy(gradients.begin(), gradients.end(), derivative.begin());
                } else if (!axisymmetric) {
                    std::array<double, 8> shape{};
                    std::array<std::array<double, 3>, 8> gradients{};
                    c3d8_detail::hex8_shape_values(x, y, z, shape, gradients);
                    std::copy(shape.begin(), shape.end(), point.shape.begin());
                    std::copy(gradients.begin(), gradients.end(), derivative.begin());
                } else if (quadratic) {
                    Quad8FaceMechanicalQuadraturePoint shape{};
                    quad8_shape_values(x, y, shape.displacement_shape, shape.derivative_xi, shape.derivative_eta);
                    for (std::size_t n = 0; n < count; ++n) {
                        point.shape[n] = shape.displacement_shape[n];
                        derivative[n] = {shape.derivative_xi[n], shape.derivative_eta[n], 0.0};
                    }
                } else {
                    constexpr std::array<std::array<double, 2>, 4> signs{{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
                    for (std::size_t n = 0; n < count; ++n) {
                        const double a = signs[n][0], b = signs[n][1];
                        point.shape[n] = 0.25 * (1 + a * x) * (1 + b * y);
                        derivative[n] = {0.25 * a * (1 + b * y), 0.25 * b * (1 + a * x), 0.0};
                    }
                }
                cartesian_detail::Matrix3 mapping{};
                if (axisymmetric)
                    mapping[2][2] = 1.0;
                for (std::size_t n = 0; n < count; ++n) {
                    const auto& c = coordinates[n];
                    if (!std::isfinite(c.x) || !std::isfinite(c.y) || !std::isfinite(c.z))
                        throw std::invalid_argument("Thermal coordinates must be finite");
                    const std::array<double, 3> position{c.x, c.y, c.z};
                    point.position.x += point.shape[n] * c.x;
                    point.position.y += point.shape[n] * c.y;
                    point.position.z += point.shape[n] * c.z;
                    for (std::size_t a = 0; a < (axisymmetric ? 2U : 3U); ++a)
                        for (std::size_t b = 0; b < (axisymmetric ? 2U : 3U); ++b)
                            mapping[a][b] += position[a] * derivative[n][b];
                }
                const double determinant = cartesian_detail::determinant(mapping);
                if (!std::isfinite(determinant) || determinant <= 0.0 || (axisymmetric && point.position.x <= 0.0))
                    throw std::domain_error("Thermal geometry requires positive Jacobian and axisymmetric radius");
                const auto inverse = cartesian_detail::inverse(mapping, determinant);
                point.measure = determinant * weights[i] * weights[j]
                                * (axisymmetric ? 2.0 * std::acos(-1.0) * point.position.x : weights[k]);
                for (std::size_t n = 0; n < count; ++n)
                    for (std::size_t a = 0; a < 3; ++a)
                        for (std::size_t b = 0; b < 3; ++b)
                            point.gradient[n][a] += derivative[n][b] * inverse[b][a];
                geometry.points.push_back(std::move(point));
            }
    return geometry;
}

ThermalResult integrate(const ThermalInput& input, bool jacobian) {
    const std::size_t count = input.geometry.node_count;
    if (input.temperature.size() != count || !input.material.function || !std::isfinite(input.time)
        || !std::isfinite(input.heat_source) || !std::isfinite(input.time_step) || input.time_step < 0.0
        || (input.time_step > 0.0 && input.previous_temperature.size() != count))
        throw std::invalid_argument("Invalid thermal element input");
    const bool lumped =
        input.geometry.element == ThermalElement::dcax4 || input.geometry.element == ThermalElement::dc3d8;
    ThermalResult result;
    result.residual.resize(count);
    if (jacobian)
        result.jacobian.resize(count * count);
    if (lumped && input.time_step > 0) {
        if (input.geometry.coordinates.size() != count)
            throw std::invalid_argument("Lumped thermal capacity requires reference nodal coordinates");
        std::vector<double> weights(count, 0.0);
        if (input.geometry.element == ThermalElement::dc3d8) {
            // Quadrature loops are xi-fast; pair each Gauss point with its HEX8 corner.
            constexpr std::array<std::size_t, 8> paired_point{0, 1, 3, 2, 4, 5, 7, 6};
            for (std::size_t n = 0; n < count; ++n)
                weights[n] = input.geometry.points.at(paired_point[n]).measure;
        } else {
            for (const auto& point : input.geometry.points)
                for (std::size_t n = 0; n < count; ++n)
                    weights[n] += point.shape.at(n) * point.measure;
        }
        for (std::size_t n = 0; n < count; ++n) {
            const auto& coordinate = input.geometry.coordinates[n];
            const MaterialFunctionContext context{input.time,
                coordinate.x,
                input.geometry.axisymmetric ? 0.0 : coordinate.y,
                input.geometry.axisymmetric ? coordinate.y : coordinate.z};
            ThermalPropertyOutput properties{}, initial{};
            input.material.function({jacobian ? adlite::Scalar::independent(input.temperature[n], 0, 1)
                                              : adlite::Scalar(input.temperature[n]),
                                        context},
                properties);
            input.material.function({input.initial_temperature, {0.0, context.x, context.y, context.z}}, initial);
            const double density = initial.density.value(), cp = properties.specific_heat.value();
            if (!std::isfinite(density) || density <= 0 || !std::isfinite(cp) || cp <= 0)
                throw std::domain_error("Lumped thermal capacity requires positive density and specific heat");
            const double rate = (input.temperature[n] - input.previous_temperature[n]) / input.time_step;
            const double storage = weights[n] * density * cp * rate;
            result.residual[n] += storage;
            result.stored_heat_rate += storage;
            if (jacobian) {
                double derivative = 0.0;
                properties.specific_heat.copy_derivatives(&derivative, 1);
                result.jacobian[n * count + n] += weights[n] * density * (cp / input.time_step + derivative * rate);
            }
        }
    }
    for (const auto& point : input.geometry.points) {
        if (point.shape.size() != count || point.gradient.size() != count || !std::isfinite(point.measure)
            || point.measure <= 0.0)
            throw std::invalid_argument("Invalid thermal quadrature layout or measure");
        double temperature = 0.0, rate = 0.0;
        std::array<double, 3> gradient{};
        for (std::size_t n = 0; n < count; ++n) {
            if (!std::isfinite(input.temperature[n]))
                throw std::domain_error("Thermal nodal temperature must be finite");
            temperature += point.shape[n] * input.temperature[n];
            if (input.time_step > 0.0 && !std::isfinite(input.previous_temperature[n]))
                throw std::domain_error("Previous thermal temperature must be finite");
            if (!lumped && input.time_step > 0.0)
                rate += point.shape[n] * (input.temperature[n] - input.previous_temperature[n]) / input.time_step;
            for (std::size_t a = 0; a < 3; ++a)
                gradient[a] += point.gradient[n][a] * input.temperature[n];
        }
        const MaterialFunctionContext context{input.time,
            point.position.x,
            input.geometry.axisymmetric ? 0.0 : point.position.y,
            input.geometry.axisymmetric ? point.position.y : point.position.z};
        ThermalPropertyOutput properties{}, initial{};
        input.material.function(
            {jacobian ? adlite::Scalar::independent(temperature, 0, 1) : adlite::Scalar(temperature), context},
            properties);
        input.material.function({input.initial_temperature, {0.0, context.x, context.y, context.z}}, initial);
        const double conductivity = properties.conductivity.value(),
                     capacity = initial.density.value() * properties.specific_heat.value();
        if (!std::isfinite(conductivity) || conductivity <= 0.0 || !std::isfinite(capacity) || capacity <= 0.0
            || initial.density.value() <= 0.0 || properties.specific_heat.value() <= 0.0)
            throw std::domain_error("Thermal conductivity and reference heat capacity must be positive and finite");
        double dk = 0.0, dc = 0.0;
        if (jacobian) {
            properties.conductivity.copy_derivatives(&dk, 1);
            properties.specific_heat.copy_derivatives(&dc, 1);
            dc *= initial.density.value();
        }
        result.heat_flux.push_back(
            {-conductivity * gradient[0], -conductivity * gradient[1], -conductivity * gradient[2]});
        result.generated_heat_rate += point.measure * input.heat_source;
        result.stored_heat_rate += point.measure * capacity * rate;
        for (std::size_t i = 0; i < count; ++i) {
            double conduction = 0.0;
            for (std::size_t a = 0; a < 3; ++a)
                conduction += point.gradient[i][a] * gradient[a];
            result.residual[i] +=
                point.measure * (conductivity * conduction + point.shape[i] * (capacity * rate - input.heat_source));
            if (jacobian)
                for (std::size_t j = 0; j < count; ++j) {
                    double stiffness = 0.0;
                    for (std::size_t a = 0; a < 3; ++a)
                        stiffness += point.gradient[i][a] * point.gradient[j][a];
                    result.jacobian[i * count + j] +=
                        point.measure
                        * (conductivity * stiffness + dk * point.shape[j] * conduction
                            + point.shape[i] * point.shape[j]
                                  * (dc * rate
                                      + (!lumped && input.time_step > 0.0 ? capacity / input.time_step : 0.0)));
                }
        }
    }
    return result;
}
} // namespace fuelsim::elements::thermal_detail
