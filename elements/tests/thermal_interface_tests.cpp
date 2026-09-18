#include "thermal_interface.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

int main() {
    using namespace fuelsim;
    using namespace fuelsim::elements;
    try {
        const std::vector<CartesianPoint3> secondary{{1, 0, 0}, {1, 1, 0}, {1, 0.5, 0}};
        const std::vector<std::vector<CartesianPoint3>> primary{{{1.01, 0, 0}, {1.01, 0.5, 0}, {1.01, 0.25, 0}},
            {{1.01, 0.5, 0}, {1.01, 1, 0}, {1.01, 0.75, 0}}};
        const auto points = make_thermal_interface_points(secondary, primary, true);
        if (points.size() != 3)
            throw std::runtime_error("Thermal interface duplicated integration points");
        const GapHeatProperties material{1.0, 1e-6};
        double heat = 0;
        for (const auto& point : points) {
            if (std::abs(point.gap - 0.01) > 1e-12)
                throw std::runtime_error("Thermal interface projection gap mismatch");
            const std::vector<double> state{400, 400, 400, 300, 300, 300};
            const auto result = evaluate_thermal_interface(point, material, state, true);
            if (std::abs(std::accumulate(result.residual.begin(), result.residual.end(), 0.0)) > 1e-10)
                throw std::runtime_error("Thermal interface fails pairwise conservation");
            for (std::size_t i = 0; i < 3; ++i)
                heat += result.residual[i];
            for (std::size_t j = 0; j < 6; ++j) {
                auto plus = state, minus = state;
                plus[j] += 1e-4;
                minus[j] -= 1e-4;
                const auto rp = evaluate_thermal_interface(point, material, plus, false);
                const auto rm = evaluate_thermal_interface(point, material, minus, false);
                for (std::size_t i = 0; i < 6; ++i)
                    if (std::abs((rp.residual[i] - rm.residual[i]) / 2e-4 - result.jacobian[6 * i + j]) > 1e-6)
                        throw std::runtime_error("Thermal interface tangent mismatch");
            }
        }
        if (std::abs(heat - 2 * std::acos(-1.0) * 10000) > 1e-8)
            throw std::runtime_error("Thermal interface analytical heat rate mismatch");
        const std::vector<CartesianPoint3>
            face{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0, 0}, {1, 0.5, 0}, {0.5, 1, 0}, {0, 0.5, 0}};
        auto other = face;
        for (auto& node : other)
            node.z = 0.01;
        const auto face_points = make_thermal_interface_points(face, {other}, false);
        GapHeatProperties temperature_law{1.0, 1e-6};
        temperature_law.law = GapHeatConductanceLaw::affine;
        temperature_law.conductance = 100;
        temperature_law.temperature_derivative = 0.2;
        temperature_law.reference_temperature = 300;
        double face_heat = 0;
        for (const auto& point : face_points) {
            std::vector<double> values(16, 300);
            std::fill(values.begin(), values.begin() + 8, 400);
            const auto evaluated = evaluate_thermal_interface(point, temperature_law, values, true);
            if (std::abs(std::accumulate(evaluated.residual.begin(), evaluated.residual.end(), 0.0)) > 1e-10)
                throw std::runtime_error("Quadratic thermal face fails conservation");
            for (std::size_t i = 0; i < 8; ++i)
                face_heat += evaluated.residual[i];
            auto plus = values, minus = values;
            plus[4] += 1e-4;
            minus[4] -= 1e-4;
            const auto rp = evaluate_thermal_interface(point, temperature_law, plus, false);
            const auto rm = evaluate_thermal_interface(point, temperature_law, minus, false);
            for (std::size_t i = 0; i < 16; ++i)
                if (std::abs((rp.residual[i] - rm.residual[i]) / 2e-4 - evaluated.jacobian[16 * i + 4]) > 1e-6)
                    throw std::runtime_error("Temperature-dependent quadratic interface tangent mismatch");
        }
        if (std::abs(face_heat - 11000) > 1e-8)
            throw std::runtime_error("Quadratic face heat rate is incorrect");
        // A single secondary neighborhood can span multiple primary faces.
        // Check the full temperature chain and zero flux for a continuous trace.
        std::vector<std::vector<CartesianPoint3>> split_primary(2, other);
        for (std::size_t side = 0; side < split_primary.size(); ++side)
            for (auto& node : split_primary[side])
                node.x = 0.5 * (node.x + static_cast<double>(side));
        const auto split_points = make_thermal_interface_points(face, split_primary, false);
        if (split_points.size() != 8)
            throw std::runtime_error("Quadratic thermal face must have eight averaged neighborhoods");
        bool spans_faces = false;
        double measure = 0.0;
        for (const auto& point : split_points) {
            measure += point.measure;
            bool first = false, second = false;
            std::vector<double> values;
            for (const auto& node : face)
                values.push_back(315.0 + 7.0 * node.x + 11.0 * node.y);
            for (const auto& location : point.primary_nodes) {
                first = first || location[0] == 0;
                second = second || location[0] == 1;
                const auto& node = split_primary[location[0]][location[1]];
                values.push_back(315.0 + 7.0 * node.x + 11.0 * node.y);
            }
            spans_faces = spans_faces || (first && second);
            const auto equilibrium = evaluate_thermal_interface(point, temperature_law, values, false);
            for (double residual : equilibrium.residual)
                if (std::abs(residual) > 1e-9)
                    throw std::runtime_error("Nonmatching quadratic interface loses a continuous temperature trace");
            for (std::size_t node = 0; node < face.size(); ++node)
                values[node] += 20.0 + static_cast<double>(node);
            const auto evaluated = evaluate_thermal_interface(point, temperature_law, values, true);
            const auto repeated = evaluate_thermal_interface(point, temperature_law, values, false);
            if (evaluated.residual != repeated.residual)
                throw std::runtime_error("Averaged interface residual depends on tangent evaluation");
            if (std::abs(std::accumulate(evaluated.residual.begin(), evaluated.residual.end(), 0.0)) > 1e-9)
                throw std::runtime_error("Averaged interface does not conserve heat across primary candidates");
            for (std::size_t column = 0; column < values.size(); ++column) {
                auto plus = values, minus = values;
                plus[column] += 1e-4;
                minus[column] -= 1e-4;
                const auto rp = evaluate_thermal_interface(point, temperature_law, plus, false);
                const auto rm = evaluate_thermal_interface(point, temperature_law, minus, false);
                for (std::size_t row = 0; row < values.size(); ++row)
                    if (std::abs((rp.residual[row] - rm.residual[row]) / 2e-4
                                 - evaluated.jacobian[values.size() * row + column])
                        > 1e-6)
                        throw std::runtime_error("Averaged nonmatching interface temperature derivative mismatch");
            }
        }
        if (!spans_faces || std::abs(measure - 1.0) > 1e-12)
            throw std::runtime_error("Quadratic neighborhoods failed primary partition or area conservation");
        bool rejected = false;
        try {
            (void)make_thermal_interface_points(secondary, {}, true);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        if (!rejected)
            throw std::runtime_error("Missing thermal projection was not rejected");
        std::cout << "thermal interface projection, conservation, tangent and missing-candidate checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
