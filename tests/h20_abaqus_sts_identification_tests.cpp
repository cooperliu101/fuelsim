#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t node_count = 8;
constexpr double perturbation = 1.0e-6;
constexpr double penalty = 1.0e8;
constexpr double base_closure = 1.0e-4;
using Vector8 = std::array<double, node_count>;
using Matrix8 = std::array<double, node_count * node_count>;

struct OperatorData final {
    Vector8 base_copen{}, base_pressure{}, base_force{};
    Matrix8 plus_copen{}, minus_copen{}, plus_pressure{}, minus_pressure{}, plus_force{}, minus_force{};
    std::array<bool, node_count * 2> observed{};
};

struct BaselineHistory final {
    double area = 0.0;
    std::array<double, 3> force{}, moment{}, center{};
};

struct GeometryHistory final {
    double area = 0.0;
    std::array<double, 3> force{}, moment{}, center{};
};

struct NodalForces final {
    std::array<std::array<double, 3>, node_count> primary{}, secondary{};
    std::size_t primary_count = 0, secondary_count = 0;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> values;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        values.push_back(value);
    return values;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size())
        throw std::invalid_argument("Incomplete H20.26 reference row in " + path);
    return std::stod(values[index]);
}

OperatorData read_operator(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.26 operator reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "step,input_local_node,input_label,output_local_node,output_label,closure_delta_m,copen_m,"
           "cpress_pa,cnormf_x_n,cnormf_y_n,cnormf_z_n")
        throw std::invalid_argument("Unexpected H20.26 operator header in " + path);

    OperatorData result;
    std::size_t base_rows = 0, perturbation_rows = 0;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        const std::string& step = values.at(0);
        const std::size_t output = static_cast<std::size_t>(number(values, 3, path)) - 1;
        if (output >= node_count)
            throw std::invalid_argument("Invalid H20.26 output node in " + path);
        if (step == "BASE") {
            result.base_copen[output] = number(values, 6, path);
            result.base_pressure[output] = number(values, 7, path);
            result.base_force[output] = number(values, 8, path);
            ++base_rows;
            continue;
        }
        const std::size_t input_node = static_cast<std::size_t>(number(values, 1, path)) - 1;
        if (input_node >= node_count)
            throw std::invalid_argument("Invalid H20.26 input node in " + path);
        const bool plus = step.find("_PLUS") != std::string::npos;
        const bool minus = step.find("_MINUS") != std::string::npos;
        if (plus == minus)
            throw std::invalid_argument("Invalid H20.26 perturbation step in " + path);
        const std::size_t entry = output * node_count + input_node;
        Matrix8& copen = plus ? result.plus_copen : result.minus_copen;
        Matrix8& pressure = plus ? result.plus_pressure : result.minus_pressure;
        Matrix8& force = plus ? result.plus_force : result.minus_force;
        copen[entry] = number(values, 6, path);
        pressure[entry] = number(values, 7, path);
        force[entry] = number(values, 8, path);
        result.observed[2 * input_node + (plus ? 0 : 1)] = true;
        ++perturbation_rows;
    }
    if (base_rows != node_count || perturbation_rows != 2 * node_count * node_count)
        throw std::invalid_argument("Unexpected H20.26 operator row count in " + path);
    for (const bool value : result.observed)
        if (!value)
            throw std::invalid_argument("Missing H20.26 perturbation in " + path);
    return result;
}

BaselineHistory read_baseline_history(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.26 history reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "step,contact_area_m2,normal_force_x_n,normal_force_y_n,normal_force_z_n,"
           "normal_moment_x_nm,normal_moment_y_nm,normal_moment_z_nm,center_x_m,center_y_m,center_z_m")
        throw std::invalid_argument("Unexpected H20.26 history header in " + path);
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.at(0) != "BASE")
            continue;
        return {number(values, 1, path),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 8, path), number(values, 9, path), number(values, 10, path)}};
    }
    throw std::invalid_argument("Missing H20.26 BASE history in " + path);
}

GeometryHistory read_geometry_history(const std::string& path, const std::string& requested_step) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.26 geometry history: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "case,step,contact_area_m2,normal_force_x_n,normal_force_y_n,normal_force_z_n,"
           "normal_moment_x_nm,normal_moment_y_nm,normal_moment_z_nm,center_x_m,center_y_m,center_z_m")
        throw std::invalid_argument("Unexpected H20.26 geometry-history header in " + path);
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.at(1) != requested_step)
            continue;
        return {number(values, 2, path),
            {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
            {number(values, 6, path), number(values, 7, path), number(values, 8, path)},
            {number(values, 9, path), number(values, 10, path), number(values, 11, path)}};
    }
    throw std::invalid_argument("Missing H20.26 geometry step " + requested_step + " in " + path);
}

NodalForces read_nodal_forces(const std::string& path, const std::string& requested_step) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.26 geometry nodal reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "case,step,side,local_node,node_label,cnormf_x_n,cnormf_y_n,cnormf_z_n,copen_m,cpress_pa,"
           "coord_x_m,coord_y_m,coord_z_m")
        throw std::invalid_argument("Unexpected H20.26 geometry-nodal header in " + path);
    NodalForces result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.at(1) != requested_step)
            continue;
        const std::size_t node = static_cast<std::size_t>(number(values, 3, path)) - 1;
        if (node >= node_count)
            throw std::invalid_argument("Invalid H20.26 geometry local node in " + path);
        const std::array<double, 3> force{number(values, 5, path), number(values, 6, path), number(values, 7, path)};
        if (values.at(2) == "primary") {
            result.primary[node] = force;
            ++result.primary_count;
        } else if (values.at(2) == "secondary") {
            result.secondary[node] = force;
            ++result.secondary_count;
        } else {
            throw std::invalid_argument("Invalid H20.26 geometry side in " + path);
        }
    }
    if (result.primary_count != node_count || result.secondary_count != node_count)
        throw std::invalid_argument("Incomplete H20.26 geometry nodal step " + requested_step + " in " + path);
    return result;
}

Vector8 solve_linear(Matrix8 matrix, Vector8 right_hand_side) {
    for (std::size_t column = 0; column < node_count; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < node_count; ++row)
            if (std::abs(matrix[row * node_count + column]) > std::abs(matrix[pivot * node_count + column]))
                pivot = row;
        if (!(std::abs(matrix[pivot * node_count + column]) > 1.0e-14))
            throw std::runtime_error("Singular H20.26 identification system");
        if (pivot != column) {
            for (std::size_t entry = column; entry < node_count; ++entry)
                std::swap(matrix[column * node_count + entry], matrix[pivot * node_count + entry]);
            std::swap(right_hand_side[column], right_hand_side[pivot]);
        }
        const double diagonal = matrix[column * node_count + column];
        for (std::size_t row = column + 1; row < node_count; ++row) {
            const double factor = matrix[row * node_count + column] / diagonal;
            for (std::size_t entry = column; entry < node_count; ++entry)
                matrix[row * node_count + entry] -= factor * matrix[column * node_count + entry];
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    Vector8 result{};
    for (std::size_t reverse = 0; reverse < node_count; ++reverse) {
        const std::size_t row = node_count - 1 - reverse;
        double value = right_hand_side[row];
        for (std::size_t column = row + 1; column < node_count; ++column)
            value -= matrix[row * node_count + column] * result[column];
        result[row] = value / matrix[row * node_count + row];
    }
    return result;
}

double relative_frobenius(const Matrix8& value, const Matrix8& reference) {
    double difference_squared = 0.0, reference_squared = 0.0;
    for (std::size_t entry = 0; entry < value.size(); ++entry) {
        difference_squared += (value[entry] - reference[entry]) * (value[entry] - reference[entry]);
        reference_squared += reference[entry] * reference[entry];
    }
    if (!(reference_squared > 0.0))
        throw std::invalid_argument("H20.26 Frobenius reference is zero");
    return std::sqrt(difference_squared / reference_squared);
}

double maximum_pointwise_relative(const Matrix8& value, const Matrix8& reference) {
    double result = 0.0;
    for (std::size_t entry = 0; entry < value.size(); ++entry) {
        if (reference[entry] == 0.0)
            throw std::invalid_argument("H20.26 pointwise reference is zero");
        result = std::max(result, std::abs(value[entry] - reference[entry]) / std::abs(reference[entry]));
    }
    return result;
}

std::array<double, 8> quad8_shape(double xi, double eta) {
    return {0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0),
        0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0),
        0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0),
        0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0),
        0.5 * (1.0 - xi * xi) * (1.0 - eta),
        0.5 * (1.0 + xi) * (1.0 - eta * eta),
        0.5 * (1.0 - xi * xi) * (1.0 + eta),
        0.5 * (1.0 - xi) * (1.0 - eta * eta)};
}

Matrix8 continuum_consistent_operator() {
    const double point = std::sqrt(3.0 / 5.0);
    const std::array<double, 3> points{-point, 0.0, point}, weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    Matrix8 result{};
    for (std::size_t first = 0; first < points.size(); ++first)
        for (std::size_t second = 0; second < points.size(); ++second) {
            const std::array<double, 8> shape = quad8_shape(points[first], points[second]);
            const double factor = penalty * 0.25 * weights[first] * weights[second];
            for (std::size_t row = 0; row < node_count; ++row)
                for (std::size_t column = 0; column < node_count; ++column)
                    result[row * node_count + column] += factor * shape[row] * shape[column];
        }
    return result;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "Usage: fuelsim_h20_26_abaqus_sts_identification_tests <operator.csv> <operator_history.csv> "
                     "<finite_sliding_nodal.csv> <finite_sliding_history.csv> <tilted_nodal.csv> "
                     "<tilted_history.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        const OperatorData data = read_operator(argv[1]);
        const BaselineHistory history = read_baseline_history(argv[2]);
        const NodalForces finite_base = read_nodal_forces(argv[3], "BASE");
        const NodalForces finite_slide = read_nodal_forces(argv[3], "SLIDE");
        const GeometryHistory finite_base_history = read_geometry_history(argv[4], "BASE");
        const GeometryHistory finite_slide_history = read_geometry_history(argv[4], "SLIDE");
        const NodalForces tilted = read_nodal_forces(argv[5], "TILTED");
        const GeometryHistory tilted_history = read_geometry_history(argv[6], "TILTED");
        Matrix8 averaging{}, abaqus_operator{}, pressure_derivative{};
        for (std::size_t entry = 0; entry < averaging.size(); ++entry) {
            averaging[entry] = -(data.plus_copen[entry] - data.minus_copen[entry]) / (2.0 * perturbation);
            abaqus_operator[entry] = (data.plus_force[entry] - data.minus_force[entry]) / (2.0 * perturbation);
            pressure_derivative[entry] =
                (data.plus_pressure[entry] - data.minus_pressure[entry]) / (2.0 * perturbation);
        }

        double row_sum_error = 0.0, averaging_asymmetry = 0.0;
        for (std::size_t row = 0; row < node_count; ++row) {
            double row_sum = 0.0;
            for (std::size_t column = 0; column < node_count; ++column) {
                row_sum += averaging[row * node_count + column];
                averaging_asymmetry = std::max(averaging_asymmetry,
                    std::abs(averaging[row * node_count + column] - averaging[column * node_count + row]));
            }
            row_sum_error = std::max(row_sum_error, std::abs(row_sum - 1.0));
        }

        Matrix8 transpose_system{};
        Vector8 normalized_base_force{};
        for (std::size_t row = 0; row < node_count; ++row) {
            normalized_base_force[row] = data.base_force[row] / (penalty * base_closure);
            for (std::size_t column = 0; column < node_count; ++column)
                transpose_system[row * node_count + column] = averaging[column * node_count + row];
        }
        const Vector8 constraint_areas = solve_linear(transpose_system, normalized_base_force);
        Matrix8 reconstructed{};
        for (std::size_t row = 0; row < node_count; ++row)
            for (std::size_t column = 0; column < node_count; ++column)
                for (std::size_t constraint = 0; constraint < node_count; ++constraint)
                    reconstructed[row * node_count + column] += penalty * averaging[constraint * node_count + row]
                                                                * constraint_areas[constraint]
                                                                * averaging[constraint * node_count + column];

        const Matrix8 continuum = continuum_consistent_operator();
        const double reconstruction_error = relative_frobenius(reconstructed, abaqus_operator);
        const double reconstruction_pointwise = maximum_pointwise_relative(reconstructed, abaqus_operator);
        const double continuum_error = relative_frobenius(continuum, abaqus_operator);
        const double continuum_pointwise = maximum_pointwise_relative(continuum, abaqus_operator);
        Matrix8 penalty_averaging{};
        for (std::size_t entry = 0; entry < penalty_averaging.size(); ++entry)
            penalty_averaging[entry] = penalty * averaging[entry];
        const double displayed_pressure_error = relative_frobenius(pressure_derivative, penalty_averaging);

        double finite_primary_change = 0.0, finite_secondary_change = 0.0, tilted_direction_error = 0.0;
        for (std::size_t node = 0; node < node_count; ++node) {
            for (std::size_t component = 0; component < 3; ++component) {
                finite_primary_change = std::max(finite_primary_change,
                    std::abs(finite_slide.primary[node][component] - finite_base.primary[node][component]));
                finite_secondary_change = std::max(finite_secondary_change,
                    std::abs(finite_slide.secondary[node][component] - finite_base.secondary[node][component]));
            }
            tilted_direction_error = std::max(tilted_direction_error,
                std::abs(tilted.secondary[node][1] + 0.05 * tilted.secondary[node][0]));
        }
        const double finite_center_shift = finite_slide_history.center[1] - finite_base_history.center[1];
        const double tilted_area_reference = std::sqrt(1.0 + 0.05 * 0.05);
        const double tilted_resultant_ratio = tilted_history.force[1] / tilted_history.force[0];

        double area_sum = 0.0, area_error = 0.0, base_copen_error = 0.0, base_pressure_error = 0.0,
               base_force_error = 0.0;
        for (std::size_t node = 0; node < node_count; ++node) {
            const double expected_area = node < 4 ? 1.0 / 24.0 : 5.0 / 24.0;
            const double expected_force = node < 4 ? -penalty * base_closure / 12.0 : penalty * base_closure / 3.0;
            area_sum += constraint_areas[node];
            area_error = std::max(area_error, std::abs(constraint_areas[node] - expected_area));
            base_copen_error = std::max(base_copen_error, std::abs(data.base_copen[node] + base_closure));
            base_pressure_error = std::max(base_pressure_error, std::abs(data.base_pressure[node] - 1.0e4));
            base_force_error = std::max(base_force_error, std::abs(data.base_force[node] - expected_force));
        }

        std::cout << "h20_26_averaging_row_sum_max_error=" << row_sum_error << '\n'
                  << "h20_26_averaging_max_asymmetry=" << averaging_asymmetry << '\n'
                  << "h20_26_constraint_areas=";
        for (std::size_t node = 0; node < node_count; ++node)
            std::cout << (node == 0 ? "" : ",") << constraint_areas[node];
        std::cout << '\n'
                  << "h20_26_constraint_area_sum=" << area_sum << '\n'
                  << "h20_26_abaqus_factorization_relative_frobenius_error=" << reconstruction_error << '\n'
                  << "h20_26_abaqus_factorization_max_pointwise_relative_error=" << reconstruction_pointwise << '\n'
                  << "h20_26_continuum_operator_relative_frobenius_error=" << continuum_error << '\n'
                  << "h20_26_continuum_operator_max_pointwise_relative_error=" << continuum_pointwise << '\n'
                  << "h20_26_displayed_pressure_vs_constraint_relative_frobenius_error=" << displayed_pressure_error
                  << '\n'
                  << "h20_26_finite_sliding_center_y_shift=" << finite_center_shift << '\n'
                  << "h20_26_finite_sliding_primary_force_max_change=" << finite_primary_change << '\n'
                  << "h20_26_finite_sliding_secondary_force_max_change=" << finite_secondary_change << '\n'
                  << "h20_26_tilted_resultant_y_over_x=" << tilted_resultant_ratio << '\n'
                  << "h20_26_tilted_nodal_direction_max_absolute_error=" << tilted_direction_error << '\n';

        bool passed = true;
        passed = check(base_copen_error < 1.0e-16 && base_pressure_error < 1.0e-8 && base_force_error < 1.0e-9,
                     "H20.26 uniform closure gives the exact C3D20 Q8 signed nodal-force pattern")
                 && passed;
        passed =
            check(std::abs(history.area - 1.0) < 1.0e-12 && std::abs(history.force[0] + 1.0e4) < 1.0e-8
                      && std::abs(history.force[1]) < 2.0e-12 && std::abs(history.force[2]) < 2.0e-12
                      && std::abs(history.moment[1] + 5.0e3) < 1.0e-8 && std::abs(history.moment[2] - 5.0e3) < 1.0e-8,
                "H20.26 Abaqus history is consistent with unit area and the applied normal resultant")
            && passed;
        passed = check(row_sum_error < 1.0e-9,
                     "H20.26 node-centered averaged constraints exactly reproduce rigid uniform closure")
                 && passed;
        passed = check(averaging_asymmetry > 0.4,
                     "H20.26 constraint averaging is observably nonsymmetric and is not direct Q8 interpolation")
                 && passed;
        passed = check(area_error < 1.0e-9 && std::abs(area_sum - 1.0) < 1.0e-9,
                     "H20.26 inferred positive constraint areas are 1/24 at corners and 5/24 at edge nodes")
                 && passed;
        passed = check(reconstruction_error < 1.0e-9 && reconstruction_pointwise < 1.0e-8,
                     "H20.26 Abaqus tangent factors as penalty times A-transpose W A")
                 && passed;
        passed = check(continuum_error > 0.19 && continuum_error < 0.20 && continuum_pointwise > 0.84,
                     "H20.26 Abaqus tangent is distinct from direct continuum-consistent Q8 penalty integration")
                 && passed;
        passed = check(displayed_pressure_error > 0.3,
                     "H20.26 displayed CPRESS is a recovered field and does not expose the constraint operator")
                 && passed;
        passed = check(std::abs(finite_base_history.force[0] + 1.0e4) < 1.0e-8
                           && std::abs(finite_slide_history.force[0] + 1.0e4) < 1.0e-8
                           && std::abs(finite_center_shift - 0.1) < 3.0e-8
                           && std::abs(finite_slide_history.moment[2] - 6.0e3) < 1.0e-8,
                     "H20.26 finite-sliding STS transfers the unchanged normal resultant at the current position")
                 && passed;
        passed = check(finite_primary_change > 1.0e3 && finite_secondary_change < 1.0e-8,
                     "H20.26 finite sliding changes primary force transfer while preserving secondary nodal forces")
                 && passed;
        passed = check(std::abs(tilted_history.area - tilted_area_reference) < 1.0e-7
                           && std::abs(tilted_resultant_ratio + 0.05) < 1.0e-7 && tilted_direction_error < 1.0e-8
                           && std::abs(tilted_history.force[2]) < 1.0e-8,
                     "H20.26 small-sliding STS follows the tilted secondary averaged normal")
                 && passed;
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
