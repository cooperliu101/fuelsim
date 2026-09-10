#include "c3d20_types.hpp"
#include "contact_types.hpp"
#include "quad8_face_boundary.hpp"
#include "quad8_face_contact.hpp"

#include "io/results_io.hpp"
#include "support/field_error_metrics.hpp"
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
struct AbaqusStep final {
    std::array<double, 3> normal{}, tangential{};
    double slip_1 = 0.0, slip_2 = 0.0;
};

struct FuelsimStep final {
    std::array<double, 3> normal{}, tangential{};
    std::vector<fuelsim::ContactPointHistory> histories;
    bool all_sliding = true;
    std::size_t sliding_points = 0;
    double residual_balance = 0.0, maximum_coulomb_error = 0.0;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double number(const std::string& text, const std::string& path) {
    std::size_t consumed = 0;
    const double result = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(result))
        throw std::invalid_argument("H20.38 CSV contains an invalid number: " + path);
    return result;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        result.push_back(value);
    return result;
}

std::array<AbaqusStep, 2> read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.38 Abaqus reference: " + path);
    std::string line;
    if (!std::getline(input, line)
        || line != "step,normal_x,normal_y,normal_z,tangential_x,tangential_y,tangential_z,slip_1,slip_2")
        throw std::invalid_argument("Unexpected H20.38 CSV header: " + path);
    std::array<AbaqusStep, 2> result{};
    for (std::size_t step = 0; step < result.size(); ++step) {
        if (!std::getline(input, line))
            throw std::invalid_argument("Incomplete H20.38 CSV: " + path);
        const std::vector<std::string> values = split(line);
        if (values.size() != 9 || number(values[0], path) != static_cast<double>(step + 1))
            throw std::invalid_argument("Invalid H20.38 CSV row: " + path);
        for (std::size_t component = 0; component < 3; ++component) {
            result[step].normal[component] = number(values[1 + component], path);
            result[step].tangential[component] = number(values[4 + component], path);
        }
        result[step].slip_1 = number(values[7], path);
        result[step].slip_2 = number(values[8], path);
    }
    if (std::getline(input, line))
        throw std::invalid_argument("H20.38 CSV has extra rows: " + path);
    return result;
}

std::vector<fuelsim::CartesianPoint3> read_input_coordinates(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.38 Abaqus input: " + path);
    std::vector<fuelsim::CartesianPoint3> result;
    bool reading_nodes = false;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.front() == '*') {
            reading_nodes = line == "*Node";
            continue;
        }
        if (!reading_nodes || line.empty())
            continue;
        const std::vector<std::string> values = split(line);
        if (values.size() != 4 || number(values[0], path) != static_cast<double>(result.size() + 1))
            throw std::invalid_argument("Invalid H20.38 Abaqus node row: " + path);
        result.push_back({number(values[1], path), number(values[2], path), number(values[3], path)});
    }
    return result;
}

fuelsim::Quad8FaceCoordinates face_coordinates(const fuelsim::UnstructuredHex20Mesh& mesh,
    const std::string& side_set) {
    static constexpr std::array<std::array<std::size_t, 8>, 6> face_nodes = {{{{0, 1, 5, 4, 8, 13, 16, 12}},
        {{1, 2, 6, 5, 9, 14, 17, 13}},
        {{2, 3, 7, 6, 10, 15, 18, 14}},
        {{3, 0, 4, 7, 11, 12, 19, 15}},
        {{0, 3, 2, 1, 11, 10, 9, 8}},
        {{4, 5, 6, 7, 16, 17, 18, 19}}}};
    const fuelsim::SideSet& set = mesh.side_set(side_set);
    if (set.sides.size() != 1)
        throw std::invalid_argument("H20.38 requires one face in " + side_set);
    const fuelsim::ElementSide& side = set.sides.front();
    fuelsim::Quad8FaceCoordinates result{};
    for (std::size_t node = 0; node < result.size(); ++node)
        result[node] = mesh.nodes().at(mesh.elements().at(side.element).nodes[face_nodes.at(side.local_side)[node]]);
    return result;
}

fuelsim::Quad8SurfaceContactLocalValues first_state() {
    fuelsim::Quad8SurfaceContactLocalValues result{};
    for (std::size_t node = 0; node < 8; ++node) {
        result[8 + node] = -0.01;
        result[24 + node] = 0.0012;
        result[40 + node] = 0.0016;
    }
    return result;
}

fuelsim::Quad8SurfaceContactLocalValues rotated_state(const fuelsim::Quad8FaceCoordinates& secondary,
    const fuelsim::Quad8FaceCoordinates& primary) {
    constexpr double angle = 0.35, qx = -0.01, qy = 0.0012, qz = 0.0016;
    const double cosine = std::cos(angle), sine = std::sin(angle);
    fuelsim::Quad8SurfaceContactLocalValues result{};
    for (std::size_t node = 0; node < 16; ++node) {
        const bool on_secondary = node < 8;
        const fuelsim::CartesianPoint3& reference = on_secondary ? secondary[node] : primary[node - 8];
        const double x = reference.x + (on_secondary ? qx : 0.0), y = reference.y + (on_secondary ? qy : 0.0);
        result[8 + node] = cosine * x - sine * y - reference.x;
        result[24 + node] = sine * x + cosine * y - reference.y;
        result[40 + node] = on_secondary ? qz : 0.0;
    }
    return result;
}

FuelsimStep evaluate(const fuelsim::NormalContactProperties& properties,
    const fuelsim::Quad8FaceCoordinates& secondary,
    const fuelsim::Quad8FaceCoordinates& primary,
    const fuelsim::Quad8SurfaceContactLocalValues& state,
    const fuelsim::Quad8SurfaceContactLocalValues& committed_state,
    const std::vector<fuelsim::ContactPointHistory>& committed_histories) {
    const fuelsim::Quad8FaceGeometry geometry = fuelsim::make_quad8_face_geometry(secondary);
    FuelsimStep result;
    result.histories.resize(geometry.mechanical_points.size());
    for (std::size_t point = 0; point < geometry.mechanical_points.size(); ++point) {
        const fuelsim::Quad8FaceMechanicalQuadraturePoint& quadrature = geometry.mechanical_points[point];
        const fuelsim::Quad8ReferenceProjectionValue reference =
            fuelsim::compute_quad8_reference_projection(secondary, primary, quadrature.displacement_shape, -1.0);
        if (!reference.projected)
            throw std::invalid_argument("H20.38 reference point did not project");
        const fuelsim::Quad8ToQuad8MechanicalGeometry contact_geometry{secondary,
            primary,
            quadrature.displacement_shape,
            quadrature.derivative_xi,
            quadrature.derivative_eta,
            reference.primary_shape,
            reference.primary_derivative_xi,
            reference.primary_derivative_eta,
            quadrature.quadrature_weight,
            -1.0};
        const fuelsim::ContactPointHistory history =
            committed_histories.empty() ? fuelsim::ContactPointHistory{} : committed_histories.at(point);
        const fuelsim::CartesianContactPointValue value = fuelsim::compute_quad8_to_quad8_contact_value(properties,
            contact_geometry,
            state,
            committed_state,
            history);
        const fuelsim::Quad8SurfaceContactLocalResidual residual =
            fuelsim::compute_quad8_to_quad8_contact(properties, contact_geometry, state, committed_state, history);
        for (std::size_t component = 0; component < 3; ++component) {
            result.normal[component] -= value.contact_force * value.normal[component];
            result.tangential[component] -= value.tributary_area * value.tangential_traction_vector[component];
            double balance = 0.0;
            const std::size_t offset = 8 + 16 * component;
            for (std::size_t node = 0; node < 16; ++node)
                balance += residual[offset + node];
            result.residual_balance = std::max(result.residual_balance, std::abs(balance));
        }
        result.all_sliding = result.all_sliding && value.sliding;
        if (value.sliding)
            ++result.sliding_points;
        result.maximum_coulomb_error = std::max(result.maximum_coulomb_error,
            std::abs(value.tangential_traction - properties.friction_coefficient * value.pressure));
        result.histories[point].sliding = value.sliding;
        result.histories[point].cartesian_elastic_tangential_slip = value.elastic_tangential_slip;
        result.histories[point].cartesian_tangent_basis_initialized = true;
        result.histories[point].cartesian_contact_normal = value.normal;
        result.histories[point].cartesian_contact_tangent_first = value.tangent_first;
    }
    return result;
}

std::array<double, 3> rotate(const std::array<double, 3>& value) {
    constexpr double angle = 0.35;
    return {std::cos(angle) * value[0] - std::sin(angle) * value[1],
        std::sin(angle) * value[0] + std::cos(angle) * value[1],
        value[2]};
}

bool compare(const std::string& mesh_path, const std::string& input_path, const std::string& reference_path) {
    const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(mesh_path);
    const std::vector<fuelsim::CartesianPoint3> input_coordinates = read_input_coordinates(input_path);
    if (input_coordinates.size() != mesh.nodes().size())
        throw std::invalid_argument("H20.38 Exodus and Abaqus node counts differ");
    double maximum_coordinate_difference = 0.0;
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node)
        maximum_coordinate_difference = std::max({maximum_coordinate_difference,
            std::abs(mesh.nodes()[node].x - input_coordinates[node].x),
            std::abs(mesh.nodes()[node].y - input_coordinates[node].y),
            std::abs(mesh.nodes()[node].z - input_coordinates[node].z)});
    const fuelsim::Quad8FaceCoordinates primary = face_coordinates(mesh, "primary_contact"),
                                        secondary = face_coordinates(mesh, "secondary_contact");
    const fuelsim::Quad8SurfaceContactLocalValues initial{}, slide = first_state(),
                                                             rotated = rotated_state(secondary, primary);
    const fuelsim::NormalContactProperties properties{1.0e5, 0.5, false, 1.0e-4};
    const FuelsimStep first = evaluate(properties, secondary, primary, slide, initial, {}),
                      second = evaluate(properties, secondary, primary, rotated, slide, first.histories);
    const std::array<AbaqusStep, 2> reference = read_reference(reference_path);

    fuelsim::test::FieldErrorMetrics normal_error, tangential_error, abaqus_normal_rotation, abaqus_tangential_rotation;
    const std::array<double, 3> rotated_reference_normal = rotate(reference[0].normal),
                                rotated_reference_tangential = rotate(reference[0].tangential),
                                rotated_fuelsim_normal = rotate(first.normal),
                                rotated_fuelsim_tangential = rotate(first.tangential);
    for (std::size_t component = 0; component < 3; ++component) {
        normal_error.add(first.normal[component], reference[0].normal[component]);
        normal_error.add(second.normal[component], reference[1].normal[component]);
        tangential_error.add(first.tangential[component], reference[0].tangential[component]);
        tangential_error.add(second.tangential[component], reference[1].tangential[component]);
        abaqus_normal_rotation.add(reference[1].normal[component], rotated_reference_normal[component]);
        abaqus_tangential_rotation.add(reference[1].tangential[component], rotated_reference_tangential[component]);
    }
    fuelsim::test::print_relative_metrics("h20_38_normal_force", normal_error);
    fuelsim::test::print_relative_metrics("h20_38_tangential_force", tangential_error);
    fuelsim::test::print_relative_metrics("h20_38_abaqus_rotated_normal_force", abaqus_normal_rotation);
    fuelsim::test::print_relative_metrics("h20_38_abaqus_rotated_tangential_force", abaqus_tangential_rotation);

    double fuelsim_rotation_error = 0.0, maximum_history_magnitude_error = 0.0;
    bool biaxial_history = true;
    for (std::size_t component = 0; component < 3; ++component) {
        fuelsim_rotation_error = std::max({fuelsim_rotation_error,
            std::abs(second.normal[component] - rotated_fuelsim_normal[component]),
            std::abs(second.tangential[component] - rotated_fuelsim_tangential[component])});
    }
    for (std::size_t point = 0; point < first.histories.size(); ++point) {
        const std::array<double, 3> expected = rotate(first.histories[point].cartesian_elastic_tangential_slip),
                                    actual = second.histories[point].cartesian_elastic_tangential_slip;
        biaxial_history = biaxial_history && std::abs(first.histories[point].cartesian_elastic_tangential_slip[1]) > 0.0
                          && std::abs(first.histories[point].cartesian_elastic_tangential_slip[2]) > 0.0;
        for (std::size_t component = 0; component < 3; ++component)
            maximum_history_magnitude_error =
                std::max(maximum_history_magnitude_error, std::abs(actual[component] - expected[component]));
    }
    const double abaqus_slip_change = std::max(std::abs(reference[1].slip_1 - reference[0].slip_1),
        std::abs(reference[1].slip_2 - reference[0].slip_2));
    std::cout << "h20_38_fuelsim_rotation_maximum_absolute_error=" << fuelsim_rotation_error << '\n'
              << "h20_38_history_rotation_maximum_absolute_error=" << maximum_history_magnitude_error << '\n'
              << "h20_38_abaqus_slip_change=" << abaqus_slip_change << '\n'
              << "h20_38_residual_balance=" << std::max(first.residual_balance, second.residual_balance) << '\n'
              << "h20_38_sliding_points=" << first.sliding_points << ',' << second.sliding_points << '\n'
              << "h20_38_maximum_coulomb_error=" << std::max(first.maximum_coulomb_error, second.maximum_coulomb_error)
              << '\n'
              << "h20_38_maximum_mesh_coordinate_difference=" << maximum_coordinate_difference << '\n';

    bool passed =
        check(maximum_coordinate_difference < 1.0e-14, "H20.38 Abaqus input uses the tracked Exodus node coordinates")
        && check(first.all_sliding && first.sliding_points == 9,
            "H20.38 enters sliding at all nine integration points with two tangential components")
        && check(std::max(first.maximum_coulomb_error, second.maximum_coulomb_error) < 1.0e-8,
            "H20.38 stays on the Coulomb circle before and after rigid rotation")
        && check(biaxial_history, "H20.38 stores two nonzero tangential elastic-slip components")
        && check(fuelsim_rotation_error < 5.0e-10 && maximum_history_magnitude_error < 1.0e-14,
            "H20.38 Fuelsim forces and elastic-slip histories rotate objectively")
        && check(reference[0].slip_1 > 0.0 && reference[0].slip_2 > 0.0 && abaqus_slip_change < 3.0e-8,
            "H20.38 Abaqus retains two nonzero slip components during rigid rotation")
        && check(std::max(first.residual_balance, second.residual_balance) < 1.0e-12,
            "H20.38 surface contact remains action-reaction conservative")
        && check(fuelsim::test::relative_metrics_below(normal_error, 0.01)
                     && normal_error.maximum_zero_reference_difference < 1.0e-10,
            "H20.38 Fuelsim normal-force metrics agree with Abaqus below 1 percent")
        && check(fuelsim::test::relative_metrics_below(tangential_error, 0.01)
                     && tangential_error.maximum_zero_reference_difference < 1.0e-10,
            "H20.38 Fuelsim biaxial tangential-force metrics agree with Abaqus below 1 percent")
        && check(fuelsim::test::relative_metrics_below(abaqus_normal_rotation, 0.001)
                     && fuelsim::test::relative_metrics_below(abaqus_tangential_rotation, 0.001),
            "H20.38 Abaqus contact resultants follow the imposed rigid rotation below 0.1 percent");
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_h20_38_friction_objectivity_abaqus_tests <mesh.e> <input.inp> <reference.csv>\n";
        return 2;
    }
    try {
        const bool passed = compare(argv[1], argv[2], argv[3]);
        if (passed)
            std::cout << "[PASS] H20.38 HEX20 friction objectivity Abaqus comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.38 raised: " << error.what() << '\n';
        return 1;
    }
}
