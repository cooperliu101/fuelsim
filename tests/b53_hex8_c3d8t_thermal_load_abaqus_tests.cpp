#include "fuelsim/core/cartesian3d_hex8.hpp"
#include "support/material_factory.hpp"
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
using StepValues = std::array<NodalValues, 3>;
using IntegrationVolumes = std::array<double, fuelsim::hex8_node_count>;
using SurfaceGeometryValues = std::map<std::size_t, NodalValues>;

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

std::map<std::string, StepValues> read_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus C3D8T thermal-load reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,element,local_node,temperature_k,reaction_heat_flux_w")
        throw std::invalid_argument("Unexpected Abaqus C3D8T thermal-load header in " + path);
    std::map<std::string, StepValues> result;
    std::map<std::string, std::array<std::array<bool, 8>, 3>> present;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 5)
            throw std::invalid_argument("Unexpected Abaqus C3D8T thermal-load column count in " + path);
        const std::size_t element = static_cast<std::size_t>(std::stoul(values[1]));
        const std::size_t node = static_cast<std::size_t>(std::stoul(values[2]));
        if (element < 1 || element > 3 || node < 1 || node > 8)
            throw std::invalid_argument("Abaqus C3D8T thermal-load label lies outside the probe");
        if (present[values[0]][element - 1][node - 1])
            throw std::invalid_argument("Duplicate Abaqus C3D8T thermal-load row");
        present[values[0]][element - 1][node - 1] = true;
        result[values[0]][element - 1][node - 1] = std::stod(values[4]);
    }
    for (const char* step : {"BODY", "SURFACE", "FILM_BASE", "FILM_ACTIVE"}) {
        if (result.find(step) == result.end())
            throw std::invalid_argument(std::string("Missing Abaqus thermal-load step ") + step);
        for (const auto& element : present[step])
            if (std::find(element.begin(), element.end(), false) != element.end())
                throw std::invalid_argument(std::string("Incomplete Abaqus C3D8T thermal-load step ") + step);
    }
    return result;
}

IntegrationVolumes read_integration_volumes(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus C3D8T integration-volume reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,element,integration_point,temperature_k,integration_volume_m3")
        throw std::invalid_argument("Unexpected Abaqus C3D8T integration-volume header in " + path);
    IntegrationVolumes result{};
    std::array<bool, 8> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 5)
            throw std::invalid_argument("Unexpected Abaqus C3D8T integration-volume column count in " + path);
        if (values[0] != "BODY" || values[1] != "1")
            continue;
        const std::size_t point = static_cast<std::size_t>(std::stoul(values[2]));
        if (point < 1 || point > 8 || present[point - 1])
            throw std::invalid_argument("Invalid Abaqus C3D8T integration-volume point label");
        present[point - 1] = true;
        result[point - 1] = std::stod(values[4]);
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("Incomplete Abaqus C3D8T integration-volume reference");
    return result;
}

SurfaceGeometryValues read_surface_geometry(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus C3D8T surface-geometry reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "element,local_node,reaction_heat_flux_w")
        throw std::invalid_argument("Unexpected Abaqus C3D8T surface-geometry header in " + path);
    SurfaceGeometryValues result;
    std::map<std::size_t, std::array<bool, 8>> present;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 3)
            throw std::invalid_argument("Unexpected Abaqus C3D8T surface-geometry column count in " + path);
        const std::size_t element = static_cast<std::size_t>(std::stoul(values[0]));
        const std::size_t node = static_cast<std::size_t>(std::stoul(values[1]));
        if (node < 1 || node > 8 || present[element][node - 1])
            throw std::invalid_argument("Invalid Abaqus C3D8T surface-geometry row");
        present[element][node - 1] = true;
        result[element][node - 1] = std::stod(values[2]);
    }
    constexpr std::array<std::size_t, 6> expected_elements = {2, 4, 5, 6, 7, 8};
    for (const std::size_t element : expected_elements)
        if (result.find(element) == result.end()
            || std::find(present[element].begin(), present[element].end(), false) != present[element].end())
            throw std::invalid_argument("Incomplete Abaqus C3D8T surface-geometry reference");
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

fuelsim::Quad4FaceCoordinates planar_trapezoid_face() {
    return {{{0.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {1.5, 1.0, 1.0}, {0.0, 1.0, 1.0}}};
}

fuelsim::Quad4FaceCoordinates single_node_warp_face(double height) {
    return {{{0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0 + height}, {0.0, 1.0, 1.0}}};
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

NodalValues legacy_body_load(const fuelsim::Hex8Geometry& geometry, double source) {
    NodalValues result{};
    for (const fuelsim::Hex8QuadraturePoint& point : geometry.points)
        for (std::size_t node = 0; node < result.size(); ++node)
            result[node] -= point.weighted_measure * point.shape[node] * source;
    return result;
}

NodalValues legacy_face_load(const fuelsim::Quad4FaceGeometry& geometry,
    const NodalValues& temperature,
    double load,
    double ambient,
    bool convection) {
    NodalValues result{};
    for (const fuelsim::Quad4FaceQuadraturePoint& point : geometry.points) {
        double point_temperature = 0.0;
        for (std::size_t node = 0; node < 4; ++node)
            point_temperature += point.shape[node] * temperature[4 + node];
        const double flux = convection ? load * (point_temperature - ambient) : -load;
        for (std::size_t node = 0; node < 4; ++node)
            result[4 + node] += point.weighted_measure * point.shape[node] * flux;
    }
    return result;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: fuelsim_b53_hex8_c3d8t_thermal_load_abaqus_tests <nodal.csv> <integration.csv> "
                     "<surface-geometry.csv>\n";
        return 2;
    }
    try {
        const std::map<std::string, StepValues> reference = read_reference(argv[1]);
        const IntegrationVolumes integration_volumes = read_integration_volumes(argv[2]);
        const SurfaceGeometryValues surface_geometry_reference = read_surface_geometry(argv[3]);
        const fuelsim::Hex8Coordinates coordinates = distorted_coordinates();
        const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
        const fuelsim::ThermoelasticProperties properties =
            fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
        const fuelsim::CartesianThermoelasticData data{fuelsim::IsotropicThermoelasticMaterial(properties), 80.0, 1.0};
        fuelsim::Hex8LocalValues volume_state{};
        std::fill(volume_state.begin(), volume_state.begin() + 8, 300.0);
        const fuelsim::Hex8LocalResidual volume_residual =
            fuelsim::compute_hex8_thermoelastic(data, geometry, volume_state);
        NodalValues body{};
        std::copy(volume_residual.begin(), volume_residual.begin() + 8, body.begin());

        const fuelsim::Quad4FaceCoordinates face_coordinates = {
            {coordinates[4], coordinates[5], coordinates[6], coordinates[7]}};
        const fuelsim::Quad4FaceGeometry face = fuelsim::make_quad4_face_geometry(face_coordinates);
        fuelsim::Quad4FaceLocalValues face_state{};
        std::fill(face_state.begin(), face_state.begin() + 4, 300.0);
        const fuelsim::Quad4FaceBoundaryData flux_data = {fuelsim::Quad4FaceBoundaryKind::surface_heat_flux,
            fuelsim::CartesianTractionComponent::x,
            40.0,
            0.0};
        const fuelsim::Quad4FaceLocalResidual flux_residual =
            fuelsim::compute_quad4_face_boundary(flux_data, face, face_state);
        NodalValues surface{};
        std::copy(flux_residual.begin(), flux_residual.begin() + 4, surface.begin() + 4);

        const NodalValues film_temperature = {0.0, 0.0, 0.0, 0.0, 360.0, 410.0, 445.0, 385.0};
        std::copy(film_temperature.begin() + 4, film_temperature.end(), face_state.begin());
        const fuelsim::Quad4FaceBoundaryData film_data = {fuelsim::Quad4FaceBoundaryKind::convection,
            fuelsim::CartesianTractionComponent::x,
            10.0,
            250.0};
        const fuelsim::Quad4FaceLocalResidual film_residual =
            fuelsim::compute_quad4_face_boundary(film_data, face, face_state);
        NodalValues film{};
        std::copy(film_residual.begin(), film_residual.begin() + 4, film.begin() + 4);
        NodalValues abaqus_film{};
        for (std::size_t node = 0; node < 8; ++node)
            abaqus_film[node] = reference.at("FILM_ACTIVE")[2][node] - reference.at("FILM_BASE")[2][node];

        const double body_error = relative_error(body, reference.at("BODY")[0]);
        const double surface_error = relative_error(surface, reference.at("SURFACE")[1]);
        const double film_error = relative_error(film, abaqus_film);
        const double legacy_body_error = relative_error(legacy_body_load(geometry, 80.0), reference.at("BODY")[0]);
        const double legacy_surface_error =
            relative_error(legacy_face_load(face, NodalValues{}, 40.0, 0.0, false), reference.at("SURFACE")[1]);
        const double legacy_film_error =
            relative_error(legacy_face_load(face, film_temperature, 10.0, 250.0, true), abaqus_film);

        double integration_volume_maximum_error = 0.0;
        for (std::size_t point = 0; point < integration_volumes.size(); ++point)
            integration_volume_maximum_error = std::max(integration_volume_maximum_error,
                std::abs(geometry.points[point].weighted_measure - integration_volumes[point]));

        const std::array<std::pair<std::size_t, fuelsim::Quad4FaceCoordinates>, 6> surface_cases = {
            {{2, face_coordinates},
                {4, planar_trapezoid_face()},
                {5, face_coordinates},
                {6, single_node_warp_face(0.1)},
                {7, single_node_warp_face(0.5)},
                {8, single_node_warp_face(1.0)}}};
        double surface_geometry_maximum_error = 0.0;
        for (const auto& [element, surface_coordinates] : surface_cases) {
            const fuelsim::Quad4FaceGeometry surface_geometry = fuelsim::make_quad4_face_geometry(surface_coordinates);
            const fuelsim::Quad4FaceLocalResidual surface_residual =
                fuelsim::compute_quad4_face_boundary(flux_data, surface_geometry, fuelsim::Quad4FaceLocalValues{});
            NodalValues actual{};
            std::copy(surface_residual.begin(), surface_residual.begin() + 4, actual.begin() + 4);
            surface_geometry_maximum_error = std::max(surface_geometry_maximum_error,
                relative_error(actual, surface_geometry_reference.at(element)));
        }
        double identical_face_maximum_difference = 0.0;
        for (std::size_t node = 0; node < 8; ++node)
            identical_face_maximum_difference = std::max(identical_face_maximum_difference,
                std::abs(surface_geometry_reference.at(2)[node] - surface_geometry_reference.at(5)[node]));

        std::cout << "b53_body_source_relative_error=" << body_error << '\n'
                  << "b53_surface_heat_flux_relative_error=" << surface_error << '\n'
                  << "b53_convection_relative_error=" << film_error << '\n'
                  << "b53_integration_volume_maximum_absolute_error=" << integration_volume_maximum_error << '\n'
                  << "b53_surface_geometry_maximum_relative_error=" << surface_geometry_maximum_error << '\n'
                  << "b53_identical_face_maximum_absolute_difference=" << identical_face_maximum_difference << '\n'
                  << "b53_legacy_gauss_body_relative_error=" << legacy_body_error << '\n'
                  << "b53_legacy_gauss_surface_relative_error=" << legacy_surface_error << '\n'
                  << "b53_legacy_gauss_convection_relative_error=" << legacy_film_error << '\n';
        const bool passed =
            check(body_error < 1.0e-12, "Abaqus and fuelsim nodal body-source vectors agree")
            && check(surface_error < 1.0e-12, "Abaqus and fuelsim nodal surface-flux vectors agree")
            && check(film_error < 1.0e-12, "Abaqus and fuelsim nodal convection vectors agree")
            && check(integration_volume_maximum_error < 1.0e-14,
                "Abaqus integration volumes equal the fuelsim standard Gauss-point volume weights")
            && check(surface_geometry_maximum_error < 1.0e-12,
                "Abaqus surface weights agree for planar and warped faces")
            && check(identical_face_maximum_difference < 1.0e-12,
                "Abaqus surface weights depend only on the loaded face geometry")
            && check(legacy_body_error > 1.0e-4, "the distorted body-source probe distinguishes Gauss integration")
            && check(legacy_surface_error > 1.0e-4, "the distorted surface-flux probe distinguishes Gauss integration")
            && check(legacy_film_error > 1.0e-4, "the nonuniform convection probe distinguishes Gauss integration");
        if (passed)
            std::cout << "[PASS] B5.3 Abaqus C3D8T thermal-load integration\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.3 Abaqus C3D8T thermal-load integration raised: " << error.what() << '\n';
        return 1;
    }
}
