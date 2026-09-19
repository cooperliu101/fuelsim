#include "dc3d8.hpp"
#include "pellet_thermal.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

int main() {
    using namespace fuelsim;
    using namespace fuelsim::elements;
    const auto material = make_builtin_material_function_registry().bind_thermal("constant_thermophysical",
        {{"conductivity", 3.0}, {"density", 10000.0}, {"specific_heat", 300.0}});
    std::vector<double> k(27 * 27, 0.0), f(27, 0.0);
    std::vector<ThermalGeometry> geometries;
    std::vector<std::array<std::size_t, 8>> connectivity;
    constexpr std::array<std::array<std::size_t, 3>, 8> corners{
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
    for (std::size_t z = 0; z < 2; ++z)
        for (std::size_t y = 0; y < 2; ++y)
            for (std::size_t x = 0; x < 2; ++x) {
                std::array<CartesianPoint3, 8> points{};
                std::array<std::size_t, 8> nodes{};
                for (std::size_t i = 0; i < 8; ++i) {
                    const auto a = x + corners[i][0], b = y + corners[i][1], c = z + corners[i][2];
                    nodes[i] = a + 3 * b + 9 * c;
                    points[i] = {0.002 * double(a), 0.002 * double(b), 0.003 * double(c)};
                }
                geometries.push_back(make_dc3d8_geometry(points));
                connectivity.push_back(nodes);
                const std::vector<double> t(8, 0.0), old;
                const auto result = evaluate_dc3d8({material, geometries.back(), t, old, 300, 0, 0, 1}, true);
                for (std::size_t i = 0; i < 8; ++i) {
                    f[nodes[i]] -= result.residual[i];
                    for (std::size_t j = 0; j < 8; ++j)
                        k[nodes[i] * 27 + nodes[j]] += result.jacobian[i * 8 + j];
                }
            }
    std::vector<std::size_t> boundary;
    for (std::size_t i = 0; i < 27; ++i)
        if (i != 13)
            boundary.push_back(i);
    ExactCondensedPellet pellet(k, f, boundary);
    double maximum = 0, energy = 0, derivative = 0;
    for (int sample = 0; sample < 8; ++sample) {
        std::vector<double> full(27), surface;
        const double q = 1e8 * (1 + sample);
        for (auto node : boundary) {
            full[node] = 600 + 80 * std::sin(double(node + 1) * 0.3 + sample);
            surface.push_back(full[node]);
        }
        double rhs = q * f[13];
        for (auto node : boundary)
            rhs -= k[13 * 27 + node] * full[node];
        full[13] = rhs / k[13 * 27 + 13];
        std::vector<double> raw(27, 0), r, j;
        for (std::size_t e = 0; e < geometries.size(); ++e) {
            std::vector<double> t, old;
            for (auto node : connectivity[e])
                t.push_back(full[node]);
            const auto result = evaluate_dc3d8({material, geometries[e], t, old, 300, 0, 0, q}, false);
            for (std::size_t i = 0; i < 8; ++i)
                raw[connectivity[e][i]] += result.residual[i];
        }
        pellet.evaluate_with_jacobian(surface, q, r, j);
        for (std::size_t i = 0; i < boundary.size(); ++i)
            maximum = std::max(maximum, std::abs(r[i] - raw[boundary[i]]));
        energy = std::max(energy,
            std::abs(std::accumulate(r.begin(), r.end(), 0.0) + q * std::accumulate(f.begin(), f.end(), 0.0)));
        for (std::size_t col = 0; col < surface.size(); ++col) {
            auto plus = surface, minus = surface;
            plus[col] += 0.001;
            minus[col] -= 0.001;
            std::vector<double> rp, rm;
            pellet.evaluate(plus, q, rp);
            pellet.evaluate(minus, q, rm);
            for (std::size_t row = 0; row < surface.size(); ++row)
                derivative =
                    std::max(derivative, std::abs((rp[row] - rm[row]) / 0.002 - j[row * surface.size() + col]));
        }
    }
    std::cout << "boundary_residual_error_W=" << maximum << " energy_error_W=" << energy
              << " tangent_error_W_per_K=" << derivative << '\n';
    if (maximum > 1e-11 || energy > 1e-10 || derivative > 1e-9)
        throw std::runtime_error("Exact condensation verification failed");
}
