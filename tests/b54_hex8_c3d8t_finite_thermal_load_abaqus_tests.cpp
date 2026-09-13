#include "c3d8_types.hpp"
#include "c3d8t.hpp"
#include "core/element_evaluation.hpp"
#include "core/element_region_data.hpp"
#include "quad4_face.hpp"
#include "support/test_support.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using NodalValues = std::array<double, fuelsim::hex8_node_count>;
using CapacityMatrix = std::array<NodalValues, fuelsim::hex8_node_count>;
using LoadValues = std::map<std::string, NodalValues>;

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

CapacityMatrix read_capacity(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus finite C3D8T capacity reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "input_local_node,output_local_node,node,temperature_k,reaction_heat_flux_w")
        throw std::invalid_argument("Unexpected Abaqus finite C3D8T capacity header in " + path);
    CapacityMatrix result{};
    std::array<std::array<bool, 8>, 8> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 5)
            throw std::invalid_argument("Unexpected Abaqus finite C3D8T capacity column count in " + path);
        const std::size_t input_node = static_cast<std::size_t>(std::stoul(values[0]));
        const std::size_t output_node = static_cast<std::size_t>(std::stoul(values[1]));
        if (input_node < 1 || input_node > 8 || output_node < 1 || output_node > 8
            || present[input_node - 1][output_node - 1])
            throw std::invalid_argument("Invalid Abaqus finite C3D8T capacity row");
        present[input_node - 1][output_node - 1] = true;
        result[input_node - 1][output_node - 1] = std::stod(values[4]);
    }
    for (const auto& row : present)
        if (std::find(row.begin(), row.end(), false) != row.end())
            throw std::invalid_argument("Incomplete Abaqus finite C3D8T capacity reference");
    return result;
}

LoadValues read_loads(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus finite C3D8T thermal-load reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,element,local_node,reaction_heat_flux_w")
        throw std::invalid_argument("Unexpected Abaqus finite C3D8T thermal-load header in " + path);
    LoadValues result;
    std::map<std::string, std::array<bool, 8>> present;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 4)
            throw std::invalid_argument("Unexpected Abaqus finite C3D8T thermal-load column count in " + path);
        const std::size_t node = static_cast<std::size_t>(std::stoul(values[2]));
        if (node < 1 || node > 8 || present[values[0]][node - 1])
            throw std::invalid_argument("Invalid Abaqus finite C3D8T thermal-load row");
        present[values[0]][node - 1] = true;
        result[values[0]][node - 1] = std::stod(values[3]);
    }
    for (const char* step : {"BODY", "SURFACE", "FILM_BASE", "FILM_ACTIVE"})
        if (result.find(step) == result.end()
            || std::find(present[step].begin(), present[step].end(), false) != present[step].end())
            throw std::invalid_argument(std::string("Incomplete Abaqus finite C3D8T thermal-load step ") + step);
    return result;
}

fuelsim::Hex8Coordinates distorted_coordinates() {
    return {{{0.00, 0.00, 0.00},
        {1.20, 0.10, -0.05},
        {1.10, 1.00, 0.10},
        {-0.10, 0.90, 0.00},
        {0.05, -0.10, 1.00},
        {1.30, 0.00, 1.10},
        {1.00, 1.20, 0.90},
        {-0.20, 1.00, 1.20}}};
}

std::array<double, 3> displacement(const fuelsim::CartesianPoint3& point) {
    return {0.5 * point.x + 0.1 * point.y, 0.25 * point.y + 0.05 * point.z, 0.08 * point.x - 0.2 * point.z};
}

fuelsim::Hex8LocalValues volume_state(const fuelsim::Hex8Coordinates& coordinates, const NodalValues& temperature) {
    fuelsim::Hex8LocalValues result{};
    for (std::size_t node = 0; node < 8; ++node) {
        result[node] = temperature[node];
        const std::array<double, 3> value = displacement(coordinates[node]);
        for (std::size_t component = 0; component < 3; ++component)
            result[8 * (component + 1) + node] = value[component];
    }
    return result;
}

fuelsim::Quad4FaceLocalValues face_state(const fuelsim::Hex8Coordinates& coordinates, const NodalValues& temperature) {
    fuelsim::Quad4FaceLocalValues result{};
    for (std::size_t face_node = 0; face_node < 4; ++face_node) {
        const std::size_t volume_node = face_node + 4;
        result[face_node] = temperature[volume_node];
        const std::array<double, 3> value = displacement(coordinates[volume_node]);
        for (std::size_t component = 0; component < 3; ++component)
            result[4 * (component + 1) + face_node] = value[component];
    }
    return result;
}

double relative_error(const NodalValues& actual, const NodalValues& expected) {
    double difference_squared = 0.0, scale_squared = 0.0;
    for (std::size_t node = 0; node < actual.size(); ++node) {
        const double difference = actual[node] - expected[node];
        difference_squared += difference * difference;
        scale_squared += expected[node] * expected[node];
    }
    return std::sqrt(difference_squared / scale_squared);
}

double relative_error(const CapacityMatrix& actual, const CapacityMatrix& expected) {
    double difference_squared = 0.0, scale_squared = 0.0;
    for (std::size_t input = 0; input < 8; ++input)
        for (std::size_t output = 0; output < 8; ++output) {
            const double difference = actual[input][output] - expected[input][output];
            difference_squared += difference * difference;
            scale_squared += expected[input][output] * expected[input][output];
        }
    return std::sqrt(difference_squared / scale_squared);
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: fuelsim_b54_hex8_c3d8t_finite_thermal_load_abaqus_tests <capacity.csv> <loads.csv>\n";
        return 2;
    }
    try {
        const CapacityMatrix capacity_reference = read_capacity(argv[1]);
        const LoadValues load_reference = read_loads(argv[2]);
        const fuelsim::Hex8Coordinates coordinates = distorted_coordinates();
        const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(coordinates);
        const fuelsim::ThermoelasticProperties properties =
            fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
        const fuelsim::IsotropicThermoelasticMaterial material(properties);
        const NodalValues uniform = {300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0};
        const fuelsim::Hex8LocalValues committed = volume_state(coordinates, uniform);
        const fuelsim::CartesianMaterialHistory history(8);

        CapacityMatrix capacity{}, reference_capacity{};
        for (std::size_t input = 0; input < 8; ++input) {
            NodalValues temperature = uniform;
            temperature[input] = 301.0;
            const fuelsim::Hex8LocalValues state = volume_state(coordinates, temperature);
            const fuelsim::CartesianRegionData current_data{material, 0.0, 1.0, fuelsim::StrainFormulation::finite};
            const fuelsim::CartesianRegionData reference_data{material, 0.0, 1.0, fuelsim::StrainFormulation::small};
            const fuelsim::Hex8LocalResidual current_residual =
                fuelsim::compute_c3d8_transient(current_data, geometry, state, committed, history, 1.0);
            const fuelsim::Hex8LocalResidual reference_residual =
                fuelsim::compute_c3d8_transient(reference_data, geometry, state, committed, history, 1.0);
            std::copy(current_residual.begin(), current_residual.begin() + 8, capacity[input].begin());
            std::copy(reference_residual.begin(), reference_residual.begin() + 8, reference_capacity[input].begin());
        }

        const fuelsim::CartesianRegionData current_body_data{material, 80.0, 1.0, fuelsim::StrainFormulation::finite};
        const fuelsim::CartesianRegionData reference_body_data{material, 80.0, 1.0, fuelsim::StrainFormulation::small};
        const fuelsim::Hex8LocalValues uniform_state = volume_state(coordinates, uniform);
        const fuelsim::Hex8LocalResidual current_body_residual =
            fuelsim::compute_c3d8_thermoelastic(current_body_data, geometry, uniform_state);
        const fuelsim::Hex8LocalResidual reference_body_residual =
            fuelsim::compute_c3d8_thermoelastic(reference_body_data, geometry, uniform_state);
        NodalValues body{}, reference_body{};
        std::copy(current_body_residual.begin(), current_body_residual.begin() + 8, body.begin());
        std::copy(reference_body_residual.begin(), reference_body_residual.begin() + 8, reference_body.begin());

        const fuelsim::Quad4FaceCoordinates face_coordinates = {
            {coordinates[4], coordinates[5], coordinates[6], coordinates[7]}};
        const fuelsim::Quad4FaceGeometry face = fuelsim::make_quad4_face_geometry(face_coordinates);
        const fuelsim::Quad4FaceBoundaryData current_flux = {fuelsim::Quad4FaceBoundaryKind::surface_heat_flux,
            fuelsim::CartesianTractionComponent::x,
            40.0,
            0.0,
            true};
        fuelsim::Quad4FaceBoundaryData reference_flux = current_flux;
        reference_flux.use_displaced_geometry = false;
        const fuelsim::Quad4FaceLocalValues uniform_face_state = face_state(coordinates, uniform);
        const fuelsim::Quad4FaceLocalResidual current_surface_residual =
            fuelsim::compute_quad4_face_boundary(current_flux, face, uniform_face_state);
        const fuelsim::Quad4FaceLocalResidual reference_surface_residual =
            fuelsim::compute_quad4_face_boundary(reference_flux, face, uniform_face_state);
        NodalValues surface{}, reference_surface{};
        std::copy(current_surface_residual.begin(), current_surface_residual.begin() + 4, surface.begin() + 4);
        std::copy(reference_surface_residual.begin(),
            reference_surface_residual.begin() + 4,
            reference_surface.begin() + 4);

        const NodalValues film_temperature = {300.0, 300.0, 300.0, 300.0, 360.0, 410.0, 445.0, 385.0};
        const fuelsim::Quad4FaceLocalValues film_face_state = face_state(coordinates, film_temperature);
        const fuelsim::Quad4FaceBoundaryData current_film = {fuelsim::Quad4FaceBoundaryKind::convection,
            fuelsim::CartesianTractionComponent::x,
            10.0,
            250.0,
            true};
        fuelsim::Quad4FaceBoundaryData reference_film = current_film;
        reference_film.use_displaced_geometry = false;
        const fuelsim::Quad4FaceLocalResidual current_film_residual =
            fuelsim::compute_quad4_face_boundary(current_film, face, film_face_state);
        const fuelsim::Quad4FaceLocalResidual reference_film_residual =
            fuelsim::compute_quad4_face_boundary(reference_film, face, film_face_state);
        NodalValues film{}, reference_film_values{}, abaqus_film{};
        std::copy(current_film_residual.begin(), current_film_residual.begin() + 4, film.begin() + 4);
        std::copy(reference_film_residual.begin(),
            reference_film_residual.begin() + 4,
            reference_film_values.begin() + 4);
        for (std::size_t node = 0; node < 8; ++node)
            abaqus_film[node] = load_reference.at("FILM_ACTIVE")[node] - load_reference.at("FILM_BASE")[node];

        const double capacity_error = relative_error(capacity, capacity_reference);
        const double reference_capacity_error = relative_error(reference_capacity, capacity_reference);
        const double body_error = relative_error(body, load_reference.at("BODY"));
        const double reference_body_error = relative_error(reference_body, load_reference.at("BODY"));
        const double surface_error = relative_error(surface, load_reference.at("SURFACE"));
        const double reference_surface_error = relative_error(reference_surface, load_reference.at("SURFACE"));
        const double film_error = relative_error(film, abaqus_film);
        const double reference_film_error = relative_error(reference_film_values, abaqus_film);
        std::cout << "b54_current_capacity_relative_error=" << capacity_error << '\n'
                  << "b54_reference_capacity_relative_error=" << reference_capacity_error << '\n'
                  << "b54_current_body_source_relative_error=" << body_error << '\n'
                  << "b54_reference_body_source_relative_error=" << reference_body_error << '\n'
                  << "b54_current_surface_heat_flux_relative_error=" << surface_error << '\n'
                  << "b54_reference_surface_heat_flux_relative_error=" << reference_surface_error << '\n'
                  << "b54_current_convection_relative_error=" << film_error << '\n'
                  << "b54_reference_convection_relative_error=" << reference_film_error << '\n';
        const bool passed =
            check(capacity_error < 5.0e-7, "Abaqus and fuelsim finite-deformation capacity matrices agree")
            && check(body_error < 5.0e-7, "Abaqus and fuelsim finite-deformation body-source vectors agree")
            && check(surface_error < 5.0e-7, "Abaqus and fuelsim finite-deformation surface-flux vectors agree")
            && check(film_error < 5.0e-7, "Abaqus and fuelsim finite-deformation convection vectors agree")
            && check(reference_capacity_error > 1.0e-2,
                "the finite-deformation capacity probe distinguishes reference integration")
            && check(reference_body_error > 1.0e-2,
                "the finite-deformation body-source probe distinguishes reference integration")
            && check(reference_surface_error > 1.0e-2,
                "the finite-deformation surface-flux probe distinguishes reference integration")
            && check(reference_film_error > 1.0e-2,
                "the finite-deformation convection probe distinguishes reference integration");
        if (passed)
            std::cout << "[PASS] B5.4 Abaqus C3D8T finite-deformation thermal integration\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.4 Abaqus C3D8T finite-deformation thermal integration raised: " << error.what() << '\n';
        return 1;
    }
}
