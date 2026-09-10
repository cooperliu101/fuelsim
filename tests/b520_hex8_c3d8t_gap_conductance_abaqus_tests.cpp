#include "core/cartesian3d_hex8.hpp"
#include "core/contact.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t node_count = 16;
constexpr std::size_t step_count = 7;
const std::array<std::string, step_count> step_names =
    {"BASE", "GAP_PLUS", "GAP_MINUS", "SECONDARY_PLUS", "SECONDARY_MINUS", "PRIMARY_PLUS", "PRIMARY_MINUS"};
const std::array<std::size_t, 4> secondary_face_nodes = {1, 2, 6, 5};
const std::array<std::size_t, 4> primary_face_nodes = {8, 11, 15, 12};

struct NodalReference final {
    double temperature = 0.0, reaction_heat_flux = 0.0, displacement_x = 0.0, reaction_force_x = 0.0;
};

using Reference = std::array<std::array<NodalReference, node_count>, step_count>;

struct FieldMetrics final {
    double relative_l2 = 0.0, relative_absolute_peak = 0.0, maximum_pointwise_relative = 0.0,
           zero_reference_maximum_absolute = 0.0;
    std::size_t zero_reference_count = 0;
};

struct KernelProbe final {
    std::array<double, node_count> residual{}, displacement_derivative{}, secondary_temperature_derivative{},
        primary_temperature_derivative{};
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        result.push_back(value);
    return result;
}

std::size_t step_index(const std::string& name) {
    const auto found = std::find(step_names.begin(), step_names.end(), name);
    if (found == step_names.end())
        throw std::invalid_argument("Unknown Abaqus B5.20 step " + name);
    return static_cast<std::size_t>(found - step_names.begin());
}

Reference read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus gap-conductance reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,node,temperature_k,reaction_heat_flux_w,u1_m,reaction_force_x_n")
        throw std::invalid_argument("Unexpected Abaqus B5.20 header in " + path);
    Reference result{};
    std::array<std::array<bool, node_count>, step_count> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 6)
            throw std::invalid_argument("Unexpected Abaqus B5.20 column count in " + path);
        const std::size_t step = step_index(values[0]), node = std::stoul(values[1]);
        if (node < 1 || node > node_count || present[step][node - 1])
            throw std::invalid_argument("Invalid or duplicate Abaqus B5.20 node label in " + path);
        present[step][node - 1] = true;
        result[step][node - 1] = {std::stod(values[2]),
            std::stod(values[3]),
            std::stod(values[4]),
            std::stod(values[5])};
    }
    for (const auto& step : present)
        if (std::find(step.begin(), step.end(), false) != step.end())
            throw std::invalid_argument("Abaqus B5.20 reference must contain every node in every step");
    return result;
}

std::array<bool, node_count> interface_mask() {
    std::array<bool, node_count> result{};
    for (const std::size_t node : secondary_face_nodes)
        result[node] = true;
    for (const std::size_t node : primary_face_nodes)
        result[node] = true;
    return result;
}

FieldMetrics field_metrics(const std::array<double, node_count>& actual,
    const std::array<double, node_count>& reference,
    const std::array<bool, node_count>& nonzero_reference) {
    double error_squared = 0.0, reference_squared = 0.0, maximum_error = 0.0, maximum_reference = 0.0;
    FieldMetrics result;
    for (std::size_t node = 0; node < node_count; ++node) {
        const double error = std::abs(actual[node] - reference[node]);
        if (nonzero_reference[node]) {
            error_squared += error * error;
            reference_squared += reference[node] * reference[node];
            maximum_error = std::max(maximum_error, error);
            maximum_reference = std::max(maximum_reference, std::abs(reference[node]));
            result.maximum_pointwise_relative =
                std::max(result.maximum_pointwise_relative, error / std::abs(reference[node]));
        } else {
            ++result.zero_reference_count;
            result.zero_reference_maximum_absolute = std::max(result.zero_reference_maximum_absolute, error);
        }
    }
    result.relative_l2 = std::sqrt(error_squared / reference_squared);
    result.relative_absolute_peak = maximum_error / maximum_reference;
    return result;
}

std::array<double, node_count> reference_residual(const Reference& reference, std::size_t step) {
    std::array<double, node_count> result{};
    for (std::size_t node = 0; node < node_count; ++node)
        result[node] = reference[step][node].reaction_heat_flux;
    return result;
}

std::array<double, node_count>
reference_derivative(const Reference& reference, std::size_t plus, std::size_t minus, double denominator) {
    std::array<double, node_count> result{};
    for (std::size_t node = 0; node < node_count; ++node)
        result[node] =
            (reference[plus][node].reaction_heat_flux - reference[minus][node].reaction_heat_flux) / denominator;
    return result;
}

KernelProbe evaluate_kernel(const fuelsim::GapHeatProperties& properties, double secondary_displacement) {
    const fuelsim::Quad4FaceCoordinates
        secondary = {{{1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}, {1.0, 0.0, 1.0}}},
        primary = {{{1.1, 0.0, 0.0}, {1.1, 1.0, 0.0}, {1.1, 1.0, 1.0}, {1.1, 0.0, 1.0}}};
    const fuelsim::Quad4FaceGeometry face = fuelsim::make_quad4_face_geometry(secondary);
    fuelsim::Quad4SurfaceContactLocalValues state{};
    for (std::size_t node = 0; node < 4; ++node) {
        state[node] = 400.0;
        state[4 + node] = 300.0;
        state[8 + node] = secondary_displacement;
    }
    fuelsim::Quad4SurfaceContactLocalResidual local_residual{};
    fuelsim::Quad4SurfaceContactLocalJacobian local_jacobian{};
    for (const fuelsim::Quad4FaceQuadraturePoint& point : face.points) {
        const fuelsim::Quad4ToQuad4HeatGeometry geometry{secondary,
            primary,
            point.shape,
            point.derivative_xi,
            point.derivative_eta,
            1.0};
        fuelsim::Quad4SurfaceContactLocalJacobian point_jacobian{};
        const fuelsim::Quad4SurfaceContactLocalResidual point_residual =
            fuelsim::compute_quad4_to_quad4_gap_heat(properties, geometry, state, &point_jacobian);
        for (std::size_t row = 0; row < local_residual.size(); ++row) {
            local_residual[row] += point_residual[row];
            for (std::size_t column = 0; column < state.size(); ++column)
                local_jacobian[row * state.size() + column] += point_jacobian[row * state.size() + column];
        }
    }
    KernelProbe result;
    for (std::size_t local_node = 0; local_node < 8; ++local_node) {
        const std::size_t global_node =
            local_node < 4 ? secondary_face_nodes[local_node] : primary_face_nodes[local_node - 4];
        result.residual[global_node] = local_residual[local_node];
        for (std::size_t direction_node = 0; direction_node < 4; ++direction_node) {
            result.displacement_derivative[global_node] += local_jacobian[local_node * 32 + 8 + direction_node];
            result.secondary_temperature_derivative[global_node] += local_jacobian[local_node * 32 + direction_node];
            result.primary_temperature_derivative[global_node] += local_jacobian[local_node * 32 + 4 + direction_node];
        }
    }
    return result;
}

bool metrics_pass(const std::string& prefix,
    const FieldMetrics& metrics,
    double relative_tolerance,
    double zero_absolute_tolerance) {
    std::cout << prefix << "_relative_l2=" << metrics.relative_l2 << '\n'
              << prefix << "_relative_absolute_peak=" << metrics.relative_absolute_peak << '\n'
              << prefix << "_maximum_pointwise_relative=" << metrics.maximum_pointwise_relative << '\n'
              << prefix << "_zero_reference_count=" << metrics.zero_reference_count << '\n'
              << prefix << "_zero_reference_maximum_absolute=" << metrics.zero_reference_maximum_absolute << '\n';
    return check(metrics.relative_l2 < relative_tolerance && metrics.relative_absolute_peak < relative_tolerance
                     && metrics.maximum_pointwise_relative < relative_tolerance
                     && metrics.zero_reference_maximum_absolute < zero_absolute_tolerance,
        prefix + " matches the Abaqus field using all three nonzero-reference metrics and separate zero accounting");
}

bool verify_law(const std::string& name,
    const Reference& reference,
    const fuelsim::GapHeatProperties& properties,
    double displacement,
    double expected_pressure) {
    const KernelProbe kernel = evaluate_kernel(properties, displacement);
    const std::array<bool, node_count> mask = interface_mask();
    const std::array<double, node_count> abaqus_residual = reference_residual(reference, 0);
    const std::array<double, node_count> abaqus_displacement_derivative = reference_derivative(reference, 2, 1, 2.0e-5);
    const std::array<double, node_count> abaqus_secondary_temperature_derivative =
        reference_derivative(reference, 3, 4, 2.0);
    const std::array<double, node_count> abaqus_primary_temperature_derivative =
        reference_derivative(reference, 5, 6, 2.0);
    bool passed = true;
    passed = metrics_pass("b520_" + name + "_residual",
                 field_metrics(kernel.residual, abaqus_residual, mask),
                 1.0e-12,
                 1.0e-10)
             && passed;
    passed = metrics_pass("b520_" + name + "_displacement_derivative",
                 field_metrics(kernel.displacement_derivative, abaqus_displacement_derivative, mask),
                 1.0e-10,
                 1.0e-8)
             && passed;
    passed = metrics_pass("b520_" + name + "_secondary_temperature_derivative",
                 field_metrics(kernel.secondary_temperature_derivative, abaqus_secondary_temperature_derivative, mask),
                 1.0e-10,
                 1.0e-10)
             && passed;
    passed = metrics_pass("b520_" + name + "_primary_temperature_derivative",
                 field_metrics(kernel.primary_temperature_derivative, abaqus_primary_temperature_derivative, mask),
                 1.0e-10,
                 1.0e-10)
             && passed;
    double abaqus_secondary_reaction_force = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        abaqus_secondary_reaction_force += reference[0][node].reaction_force_x;
    const double pressure_error =
        expected_pressure == 0.0 ? std::abs(abaqus_secondary_reaction_force)
                                 : std::abs(abaqus_secondary_reaction_force - expected_pressure) / expected_pressure;
    std::cout << "b520_" << name << "_abaqus_secondary_reaction_force=" << abaqus_secondary_reaction_force << '\n'
              << "b520_" << name << "_pressure_error=" << pressure_error << '\n';
    passed = check(pressure_error < 1.0e-9,
                 name + " Abaqus mechanical reaction confirms the pressure supplied to gap conductance")
             && passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: fuelsim_b520_hex8_c3d8t_gap_conductance_abaqus_tests <clearance.csv> <pressure.csv>\n";
        return 2;
    }
    try {
        fuelsim::GapHeatProperties clearance{};
        clearance.law = fuelsim::GapHeatConductanceLaw::affine;
        clearance.conductance = 100.0;
        clearance.clearance_derivative = -1000.0;
        clearance.temperature_derivative = 0.1;
        clearance.reference_temperature = 350.0;
        fuelsim::GapHeatProperties pressure{};
        pressure.law = fuelsim::GapHeatConductanceLaw::affine;
        pressure.conductance = 50.0;
        pressure.pressure_derivative = 0.02;
        pressure.temperature_derivative = 0.1;
        pressure.reference_temperature = 350.0;
        pressure.contact_penalty = 1.0e5;
        bool passed = verify_law("clearance", read_reference(argv[1]), clearance, 0.05, 0.0);
        passed = verify_law("pressure", read_reference(argv[2]), pressure, 0.11, 1000.0) && passed;
        if (passed)
            std::cout << "[PASS] B5.20 Abaqus C3D8T clearance-, pressure-, and temperature-dependent gap conductance\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
