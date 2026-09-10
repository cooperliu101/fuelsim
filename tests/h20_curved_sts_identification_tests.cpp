#include "c3d20_types.hpp"
#include "io/results_io.hpp"
#include "quad8_face.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t node_count = 8;
constexpr std::size_t component_count = 3;
constexpr double perturbation = 1.0e-6;

struct Record final {
    std::string step, component;
    std::size_t input, output;
    double delta, copen, cpress;
    std::array<double, 3> force;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ','))
        result.push_back(value);
    return result;
}

std::vector<Record> read_records(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.31 operator reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "step,input_local_node,input_component,output_local_node,displacement_delta_m,copen_m,cpress_pa,"
           "cnormf_x_n,cnormf_y_n,cnormf_z_n")
        throw std::invalid_argument("Unexpected H20.31 operator header in " + path);
    std::vector<Record> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 10)
            throw std::invalid_argument("Incomplete H20.31 operator row in " + path);
        result.push_back({values[0],
            values[2],
            static_cast<std::size_t>(std::stoul(values[1])),
            static_cast<std::size_t>(std::stoul(values[3])),
            std::stod(values[4]),
            std::stod(values[5]),
            std::stod(values[6]),
            {std::stod(values[7]), std::stod(values[8]), std::stod(values[9])}});
    }
    return result;
}

const Record& find_record(const std::vector<Record>& records, const std::string& step, std::size_t output) {
    const auto found = std::find_if(records.begin(), records.end(), [&](const Record& record) {
        return record.step == step && record.output == output;
    });
    if (found == records.end())
        throw std::invalid_argument("Missing H20.31 step/output record");
    return *found;
}

fuelsim::CartesianPoint3 cross(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return {first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x};
}

double dot(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

fuelsim::CartesianPoint3 unit(const fuelsim::CartesianPoint3& value) {
    const double measure = std::sqrt(dot(value, value));
    if (!(measure > 0.0))
        throw std::invalid_argument("H20.31 has an undefined effective normal");
    return {value.x / measure, value.y / measure, value.z / measure};
}

fuelsim::Quad8FaceCoordinates first_secondary_face(const fuelsim::UnstructuredHex20Mesh& mesh) {
    static constexpr std::array<std::array<std::size_t, 8>, 6> face_nodes = {{{{0, 1, 5, 4, 8, 13, 16, 12}},
        {{1, 2, 6, 5, 9, 14, 17, 13}},
        {{2, 3, 7, 6, 10, 15, 18, 14}},
        {{3, 0, 4, 7, 11, 12, 19, 15}},
        {{0, 3, 2, 1, 11, 10, 9, 8}},
        {{4, 5, 6, 7, 16, 17, 18, 19}}}};
    const fuelsim::ElementSide& side = mesh.side_set("secondary_contact").sides.front();
    fuelsim::Quad8FaceCoordinates result{};
    for (std::size_t local = 0; local < result.size(); ++local)
        result[local] = mesh.nodes()[mesh.elements()[side.element].nodes[face_nodes[side.local_side][local]]];
    return result;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_h20_31_curved_sts_identification_tests <operator.csv> <quadratic_mesh.e>\n";
        return 2;
    }
    try {
        const std::vector<Record> records = read_records(argv[1]);
        const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(argv[2]);
        const fuelsim::Quad8FaceCoordinates coordinates = first_secondary_face(mesh);
        // The Abaqus S6 output order is the reverse circumferential orientation
        // of Fuelsim's local side 3 order for this generated element.
        const std::array<std::array<double, 2>, 8> locations = {{{0.75, -0.75},
            {-0.75, -0.75},
            {-0.75, 0.75},
            {0.75, 0.75},
            {0.0, -0.5},
            {-0.5, 0.0},
            {0.0, 0.5},
            {0.5, 0.0}}};
        std::array<fuelsim::CartesianPoint3, 8> expected_normals{};
        for (std::size_t output = 0; output < node_count; ++output) {
            const fuelsim::Quad8FaceMechanicalQuadraturePoint point =
                fuelsim::make_quad8_face_mechanical_point(coordinates, locations[output][0], locations[output][1], 1.0);
            expected_normals[output] = unit(cross(point.tangent_xi, point.tangent_eta));
        }

        std::array<double, node_count * node_count * component_count> gap_gradient{};
        std::array<double, node_count * component_count * node_count * component_count> force_tangent{};
        double maximum_delta_error = 0.0, minimum_base_pressure = std::numeric_limits<double>::max();
        for (std::size_t output = 0; output < node_count; ++output)
            minimum_base_pressure = std::min(minimum_base_pressure, find_record(records, "BASE", output + 1).cpress);
        const std::array<std::string, 3> component_names = {"X", "Y", "Z"};
        for (std::size_t input = 0; input < node_count; ++input)
            for (std::size_t component = 0; component < component_count; ++component) {
                const std::string prefix = "S" + std::to_string(input + 1) + '_' + component_names[component] + '_';
                const std::string plus_step = prefix + "PLUS", minus_step = prefix + "MINUS";
                const Record& plus = find_record(records, plus_step, 1);
                const Record& minus = find_record(records, minus_step, 1);
                maximum_delta_error = std::max(maximum_delta_error,
                    std::max(std::abs(plus.delta - perturbation), std::abs(minus.delta + perturbation)));
                for (std::size_t output = 0; output < node_count; ++output) {
                    const Record& output_plus = find_record(records, plus_step, output + 1);
                    const Record& output_minus = find_record(records, minus_step, output + 1);
                    gap_gradient[(output * node_count + input) * component_count + component] =
                        (output_plus.copen - output_minus.copen) / (2.0 * perturbation);
                    for (std::size_t output_component = 0; output_component < component_count; ++output_component) {
                        const std::size_t row = component_count * output + output_component;
                        const std::size_t column = component_count * input + component;
                        force_tangent[row * node_count * component_count + column] =
                            (output_plus.force[output_component] - output_minus.force[output_component])
                            / (2.0 * perturbation);
                    }
                }
            }

        double maximum_normal_alignment_error = 0.0, maximum_axial_gap_gradient = 0.0,
               minimum_effective_normal_dot = 1.0;
        for (std::size_t output = 0; output < node_count; ++output) {
            fuelsim::CartesianPoint3 gradient{};
            for (std::size_t input = 0; input < node_count; ++input) {
                gradient.x += gap_gradient[(output * node_count + input) * component_count];
                gradient.y += gap_gradient[(output * node_count + input) * component_count + 1];
                maximum_axial_gap_gradient = std::max(maximum_axial_gap_gradient,
                    std::abs(gap_gradient[(output * node_count + input) * component_count + 2]));
            }
            maximum_normal_alignment_error =
                std::max(maximum_normal_alignment_error, 1.0 - std::abs(dot(unit(gradient), expected_normals[output])));
            for (std::size_t other = output + 1; other < node_count; ++other)
                minimum_effective_normal_dot = std::min(minimum_effective_normal_dot,
                    std::abs(dot(expected_normals[output], expected_normals[other])));
        }
        double tangent_difference_squared = 0.0, tangent_reference_squared = 0.0;
        const std::size_t tangent_size = node_count * component_count;
        for (std::size_t row = 0; row < tangent_size; ++row)
            for (std::size_t column = 0; column < tangent_size; ++column) {
                const double value = force_tangent[row * tangent_size + column];
                const double transpose = force_tangent[column * tangent_size + row];
                tangent_difference_squared += (value - transpose) * (value - transpose);
                tangent_reference_squared += value * value;
            }
        const double tangent_symmetry_error = std::sqrt(tangent_difference_squared / tangent_reference_squared);
        std::cout << std::scientific << std::setprecision(12)
                  << "h20_31_effective_normal_alignment_error=" << maximum_normal_alignment_error << '\n'
                  << "h20_31_minimum_effective_normal_dot=" << minimum_effective_normal_dot << '\n'
                  << "h20_31_maximum_axial_gap_gradient=" << maximum_axial_gap_gradient << '\n'
                  << "h20_31_force_tangent_symmetry_error=" << tangent_symmetry_error << '\n';
        bool passed = true;
        passed = check(records.size() == 49 * node_count && maximum_delta_error < 1.0e-18,
                     "H20.31 contains the base state and all forty-eight full-vector perturbation steps")
                 && passed;
        passed =
            check(minimum_base_pressure > 0.0, "H20.31 keeps all eight curved secondary constraints active") && passed;
        passed = check(maximum_normal_alignment_error < 1.0e-5 && minimum_effective_normal_dot < 0.99,
                     "H20.31 Abaqus gap gradients follow distinct normals at the identified effective centers")
                 && passed;
        passed =
            check(maximum_axial_gap_gradient < 1.0e-8, "H20.31 cylindrical contact normals have no axial component")
            && passed;
        passed = check(tangent_symmetry_error < 1.0e-10,
                     "H20.31 Abaqus curved-contact force tangent is work-conjugate and symmetric")
                 && passed;
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.31 curved Abaqus identification raised: " << error.what() << '\n';
        return 1;
    }
}
