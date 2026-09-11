#include "cax_common.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double gauss = 0.577350269189625764509148780501957456;

std::array<double, quad4_node_count> shape_functions(double xi, double eta) {
    return {
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta),
    };
}

std::array<double, quad4_node_count> shape_derivative_xi(double eta) {
    return {
        -0.25 * (1.0 - eta),
        0.25 * (1.0 - eta),
        0.25 * (1.0 + eta),
        -0.25 * (1.0 + eta),
    };
}

std::array<double, quad4_node_count> shape_derivative_eta(double xi) {
    return {
        -0.25 * (1.0 - xi),
        -0.25 * (1.0 + xi),
        0.25 * (1.0 + xi),
        0.25 * (1.0 - xi),
    };
}
} // namespace

Quad4RzGeometry cax4_detail::make_quad4_rz_geometry(const Quad4Coordinates& coordinates) {
    const std::array<std::array<double, 2>, 4> locations = {{
        {{-gauss, -gauss}},
        {{gauss, -gauss}},
        {{gauss, gauss}},
        {{-gauss, gauss}},
    }};
    Quad4RzGeometry geometry{};
    geometry.coordinates = coordinates;
    for (std::size_t q = 0; q < locations.size(); ++q) {
        const double xi = locations[q][0], eta = locations[q][1];
        const std::array<double, quad4_node_count> shape = shape_functions(xi, eta);
        const std::array<double, quad4_node_count> derivative_xi = shape_derivative_xi(eta);
        const std::array<double, quad4_node_count> derivative_eta = shape_derivative_eta(xi);
        double dr_dxi = 0.0, dr_deta = 0.0, dz_dxi = 0.0, dz_deta = 0.0, radius = 0.0, axial_coordinate = 0.0;
        for (std::size_t node = 0; node < quad4_node_count; ++node) {
            dr_dxi += derivative_xi[node] * coordinates[node].r;
            dr_deta += derivative_eta[node] * coordinates[node].r;
            dz_dxi += derivative_xi[node] * coordinates[node].z;
            dz_deta += derivative_eta[node] * coordinates[node].z;
            radius += shape[node] * coordinates[node].r;
            axial_coordinate += shape[node] * coordinates[node].z;
        }
        const double determinant = dr_dxi * dz_deta - dr_deta * dz_dxi;
        if (!(determinant > 0.0))
            throw std::invalid_argument("Quad4RzGeometry requires positive Jacobian determinant");
        if (!(radius > 0.0))
            throw std::invalid_argument("Quad4RzGeometry requires positive quadrature radius");
        RzQuadraturePoint& point = geometry.points[q];
        point.shape = shape;
        point.radius = radius;
        point.axial_coordinate = axial_coordinate;
        point.weighted_measure = 2.0 * pi * radius * determinant;
        for (std::size_t node = 0; node < quad4_node_count; ++node) {
            point.gradient_r[node] = (dz_deta * derivative_xi[node] - dz_dxi * derivative_eta[node]) / determinant;
            point.gradient_z[node] = (-dr_deta * derivative_xi[node] + dr_dxi * derivative_eta[node]) / determinant;
        }
    }
    return geometry;
}

} // namespace fuelsim

namespace fuelsim {
AxisymmetricHughesWinget evaluate_axisymmetric_hughes_winget(const adlite::Scalar& hrr,
    const adlite::Scalar& hrz,
    const adlite::Scalar& hzr,
    const adlite::Scalar& hzz) {
    const adlite::Scalar shear = 0.5 * (hrz + hzr), spin = 0.25 * (hrz - hzr), denominator = 1.0 + spin * spin,
                         cosine = (1.0 - spin * spin) / denominator, sine = 2.0 * spin / denominator;
    return {{cosine, sine, -sine, cosine, adlite::Scalar(1.0)},
        cosine * cosine * hrr - 2.0 * cosine * sine * shear + sine * sine * hzz,
        sine * sine * hrr + 2.0 * cosine * sine * shear + cosine * cosine * hzz,
        cosine * sine * (hrr - hzz) + (cosine * cosine - sine * sine) * shear};
}
} // namespace fuelsim
