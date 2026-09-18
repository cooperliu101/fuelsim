#include "thermal_boundary.hpp"
#include "quad8_face.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::elements {
ThermalSurfacePoint
thermal_surface_point(const std::vector<CartesianPoint3>& coordinates, bool axisymmetric, double x, double y) {
    const std::size_t count = coordinates.size();
    if ((axisymmetric && count != 2 && count != 3) || (!axisymmetric && count != 4 && count != 8))
        throw std::invalid_argument("Invalid thermal boundary topology");
    ThermalPoint point;
    point.shape.resize(count);
    std::vector<double> dx(count), dy(count);
    if (axisymmetric) {
        if (count == 2) {
            point.shape = {(1 - x) / 2, (1 + x) / 2};
            dx = {-0.5, 0.5};
        } else {
            point.shape = {x * (x - 1) / 2, x * (x + 1) / 2, 1 - x * x};
            dx = {x - 0.5, x + 0.5, -2 * x};
        }
    } else if (count == 8) {
        Quad8FaceMechanicalQuadraturePoint shape{};
        quad8_shape_values(x, y, shape.displacement_shape, shape.derivative_xi, shape.derivative_eta);
        point.shape.assign(shape.displacement_shape.begin(), shape.displacement_shape.end());
        dx.assign(shape.derivative_xi.begin(), shape.derivative_xi.end());
        dy.assign(shape.derivative_eta.begin(), shape.derivative_eta.end());
    } else {
        constexpr std::array<std::array<double, 2>, 4> signs{{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
        for (std::size_t n = 0; n < 4; ++n) {
            const double a = signs[n][0], b = signs[n][1];
            point.shape[n] = (1 + a * x) * (1 + b * y) / 4;
            dx[n] = a * (1 + b * y) / 4;
            dy[n] = b * (1 + a * x) / 4;
        }
    }
    std::array<double, 3> a{}, b{};
    for (std::size_t n = 0; n < count; ++n) {
        point.position.x += point.shape[n] * coordinates[n].x;
        point.position.y += point.shape[n] * coordinates[n].y;
        point.position.z += point.shape[n] * coordinates[n].z;
        const std::array<double, 3> c{coordinates[n].x, coordinates[n].y, coordinates[n].z};
        for (std::size_t d = 0; d < 3; ++d) {
            a[d] += dx[n] * c[d];
            b[d] += dy[n] * c[d];
        }
    }
    return {std::move(point), a, b};
}

ThermalGeometry make_thermal_boundary_geometry(const std::vector<CartesianPoint3>& coordinates, bool axisymmetric) {
    const std::size_t count = coordinates.size();
    if ((axisymmetric && count != 2 && count != 3) || (!axisymmetric && count != 4 && count != 8))
        throw std::invalid_argument("Invalid thermal boundary topology");
    ThermalGeometry result;
    result.node_count = count;
    result.axisymmetric = axisymmetric;
    result.coordinates = coordinates;
    const std::array<double, 3> gauss{-std::sqrt(0.6), 0.0, std::sqrt(0.6)}, weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    for (std::size_t j = 0; j < (axisymmetric ? 1U : 3U); ++j)
        for (std::size_t i = 0; i < 3; ++i) {
            const double x = gauss[i], y = gauss[j];
            auto surface = thermal_surface_point(coordinates, axisymmetric, x, y);
            auto point = std::move(surface.point);
            const auto& a = surface.tangent_xi;
            const auto& b = surface.tangent_eta;
            if (axisymmetric)
                point.measure = std::hypot(a[0], a[1]) * 2 * std::acos(-1.0) * point.position.x * weights[i];
            else
                point.measure =
                    std::hypot(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
                    * weights[i] * weights[j];
            if (!std::isfinite(point.measure) || point.measure < 0.0)
                throw std::domain_error("Invalid thermal boundary measure");
            result.points.push_back(std::move(point));
        }
    return result;
}

ThermalResult evaluate_thermal_boundary(const ThermalGeometry& geometry,
    const std::vector<double>& temperature,
    double inward_flux,
    double coefficient,
    double ambient,
    bool jacobian) {
    if (temperature.size() != geometry.node_count || !std::isfinite(inward_flux) || !std::isfinite(coefficient)
        || coefficient < 0 || !std::isfinite(ambient))
        throw std::invalid_argument("Invalid thermal boundary input");
    ThermalResult result;
    const auto count = temperature.size();
    result.residual.resize(count);
    if (jacobian)
        result.jacobian.resize(count * count);
    for (const auto& point : geometry.points) {
        double value = 0.0;
        for (std::size_t i = 0; i < count; ++i)
            value += point.shape[i] * temperature[i];
        for (std::size_t i = 0; i < count; ++i) {
            result.residual[i] += point.measure * point.shape[i] * (coefficient * (value - ambient) - inward_flux);
            if (jacobian)
                for (std::size_t j = 0; j < count; ++j)
                    result.jacobian[i * count + j] += point.measure * point.shape[i] * coefficient * point.shape[j];
        }
    }
    return result;
}
} // namespace fuelsim::elements
