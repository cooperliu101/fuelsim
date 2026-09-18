#include "thermal_interface.hpp"
#include "c3d_common.hpp"
#include "contact_common.hpp"
#include "quad8_face.hpp"
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace fuelsim::elements {
namespace {
ThermalInterfacePoint project_thermal_point(const std::vector<CartesianPoint3>& secondary,
    const std::vector<std::vector<CartesianPoint3>>& primary,
    bool axisymmetric,
    double xi,
    double eta,
    double weight) {
    const auto source = thermal_surface_point(secondary, axisymmetric, xi, eta);
    const auto& a = source.tangent_xi;
    const auto& b = source.tangent_eta;
    std::array<double, 3> normal =
        axisymmetric
            ? std::array<double, 3>{a[1], -a[0], 0}
            : std::array<double, 3>{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    const double norm = std::hypot(normal[0], normal[1], normal[2]);
    if (!std::isfinite(norm) || norm <= 0)
        throw std::domain_error("Thermal interface has degenerate secondary geometry");
    for (double& v : normal)
        v /= norm;
    ThermalInterfacePoint selected;
    selected.secondary_shape = source.point.shape;
    selected.measure = norm * weight * (axisymmetric ? 2 * std::acos(-1.0) * source.point.position.x : 1.0);
    double distance = std::numeric_limits<double>::infinity();
    for (std::size_t candidate = 0; candidate < primary.size(); ++candidate) {
        double x = 0, y = 0, gap = 0;
        bool converged = false;
        for (std::size_t iteration = 0; iteration < 30; ++iteration) {
            const auto target = thermal_surface_point(primary[candidate], axisymmetric, x, y);
            const std::array<double, 3> residual{target.point.position.x + gap * normal[0] - source.point.position.x,
                target.point.position.y + gap * normal[1] - source.point.position.y,
                target.point.position.z + gap * normal[2] - source.point.position.z};
            cartesian_detail::Matrix3 mapping{};
            for (std::size_t d = 0; d < 3; ++d) {
                mapping[d][0] = target.tangent_xi[d];
                mapping[d][1] = axisymmetric ? normal[d] : target.tangent_eta[d];
                mapping[d][2] = axisymmetric ? (d == 2 ? 1.0 : 0.0) : normal[d];
            }
            const double det = cartesian_detail::determinant(mapping);
            if (!std::isfinite(det) || det == 0)
                break;
            const auto inverse = cartesian_detail::inverse(mapping, det);
            std::array<double, 3> delta{};
            for (std::size_t d = 0; d < 3; ++d)
                for (std::size_t e = 0; e < 3; ++e)
                    delta[d] += inverse[d][e] * residual[e];
            x -= delta[0];
            if (axisymmetric)
                gap -= delta[1];
            else {
                y -= delta[1];
                gap -= delta[2];
            }
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(gap))
                break;
            if (std::abs(delta[0]) < 1e-12 && (axisymmetric || std::abs(delta[1]) < 1e-12)) {
                converged = true;
                break;
            }
        }
        if (!converged || std::abs(x) > 1 + 1e-10 || (!axisymmetric && std::abs(y) > 1 + 1e-10))
            continue;
        if (std::abs(gap) < distance) {
            distance = std::abs(gap);
            selected.primary_nodes.clear();
            for (std::size_t node = 0; node < primary[candidate].size(); ++node)
                selected.primary_nodes.push_back({candidate, node});
            selected.primary_shape = thermal_surface_point(primary[candidate], axisymmetric, x, y).point.shape;
            selected.gap = -gap;
        }
    }
    if (!std::isfinite(distance))
        throw std::invalid_argument("Thermal interface secondary integration point has no primary projection");
    return selected;
}
} // namespace

std::vector<ThermalInterfacePoint> make_thermal_interface_points(const std::vector<CartesianPoint3>& secondary,
    const std::vector<std::vector<CartesianPoint3>>& primary,
    bool axisymmetric) {
    std::vector<ThermalInterfacePoint> result;
    if (!axisymmetric && secondary.size() == 8) {
        // Average both full quadratic temperature traces over each positive
        // nodal neighborhood before applying the interface heat-transfer law.
        for (std::size_t constraint = 0; constraint < 8; ++constraint) {
            ThermalInterfacePoint average;
            average.secondary_shape.resize(secondary.size());
            std::map<std::array<std::size_t, 2>, double> primary_moments;
            const double fraction = constraint < 4 ? 1.0 / 24.0 : 5.0 / 24.0;
            for (const auto& sample : abaqus_quad8_primary_transfer_rule(constraint)) {
                const auto point = project_thermal_point(secondary,
                    primary,
                    false,
                    2.0 * sample.first - 1.0,
                    2.0 * sample.second - 1.0,
                    4.0 * fraction * sample.weight);
                average.measure += point.measure;
                average.gap += point.measure * point.gap;
                for (std::size_t node = 0; node < secondary.size(); ++node)
                    average.secondary_shape[node] += point.measure * point.secondary_shape[node];
                for (std::size_t node = 0; node < point.primary_shape.size(); ++node)
                    primary_moments[point.primary_nodes[node]] += point.measure * point.primary_shape[node];
            }
            if (!std::isfinite(average.measure) || average.measure <= 0.0)
                throw std::domain_error("Quadratic thermal interface has invalid neighborhood measure");
            average.gap /= average.measure;
            for (double& value : average.secondary_shape)
                value /= average.measure;
            for (const auto& moment : primary_moments) {
                average.primary_nodes.push_back(moment.first);
                average.primary_shape.push_back(moment.second / average.measure);
            }
            result.push_back(std::move(average));
        }
        return result;
    }
    const std::array<double, 3> gauss{-std::sqrt(0.6), 0.0, std::sqrt(0.6)}, weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    for (std::size_t j = 0; j < (axisymmetric ? 1U : 3U); ++j)
        for (std::size_t i = 0; i < 3; ++i)
            result.push_back(project_thermal_point(secondary,
                primary,
                axisymmetric,
                gauss[i],
                gauss[j],
                weights[i] * (axisymmetric ? 1.0 : weights[j])));
    return result;
}

ThermalResult evaluate_thermal_interface(const ThermalInterfacePoint& point,
    const GapHeatProperties& material,
    const std::vector<double>& temperatures,
    bool jacobian) {
    const auto secondary_count = point.secondary_shape.size(), count = secondary_count + point.primary_shape.size();
    if (temperatures.size() != count || material.pressure_derivative != 0)
        throw std::invalid_argument("Thermal interface requires temperature data and pressure-independent conductance");
    double secondary = 0, primary = 0;
    for (std::size_t i = 0; i < secondary_count; ++i)
        secondary += point.secondary_shape[i] * temperatures[i];
    for (std::size_t i = 0; i < point.primary_shape.size(); ++i)
        primary += point.primary_shape[i] * temperatures[secondary_count + i];
    const auto ts = jacobian ? adlite::Scalar::independent(secondary, 0, 2) : adlite::Scalar(secondary);
    const auto tp = jacobian ? adlite::Scalar::independent(primary, 1, 2) : adlite::Scalar(primary);
    const auto flux = contact_common::gap_conductance(material, point.gap, ts, tp) * (ts - tp);
    std::array<double, 2> derivative{};
    if (jacobian)
        flux.copy_derivatives(derivative.data(), 2);
    ThermalResult result;
    result.residual.resize(count);
    if (jacobian)
        result.jacobian.resize(count * count);
    for (std::size_t i = 0; i < count; ++i) {
        const double weight =
            point.measure
            * (i < secondary_count ? point.secondary_shape[i] : -point.primary_shape[i - secondary_count]);
        result.residual[i] = weight * flux.value();
        if (jacobian)
            for (std::size_t j = 0; j < count; ++j)
                result.jacobian[i * count + j] =
                    weight
                    * (j < secondary_count ? derivative[0] * point.secondary_shape[j]
                                           : derivative[1] * point.primary_shape[j - secondary_count]);
    }
    return result;
}
} // namespace fuelsim::elements
