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
constexpr std::size_t matching_count = 8;
constexpr std::size_t secondary_count = 13;
constexpr std::size_t primary_count = 23;
constexpr double perturbation = 1.0e-6;
constexpr double penalty = 1.0e8;
constexpr double base_closure = 1.0e-4;
using Vector8 = std::array<double, matching_count>;
using Matrix8 = std::array<double, matching_count * matching_count>;
using Vector13 = std::array<double, secondary_count>;
using Matrix13 = std::array<double, secondary_count * secondary_count>;
using PrimaryVector = std::array<double, primary_count>;
using PrimaryMatrix = std::array<double, primary_count * secondary_count>;
using TransferMatrix = std::array<double, secondary_count * primary_count>;

struct MatchingData final {
    Matrix8 plus_copen{}, minus_copen{};
};

struct NonmatchingData final {
    Vector13 base_copen{}, base_pressure{}, base_force{};
    PrimaryVector base_primary_force{};
    Vector13 open_copen{}, open_pressure{}, open_force{};
    PrimaryVector open_primary_force{};
    Vector13 reclose_copen{}, reclose_pressure{}, reclose_force{};
    PrimaryVector reclose_primary_force{};
    Matrix13 plus_copen{}, minus_copen{}, plus_pressure{}, minus_pressure{}, plus_force{}, minus_force{};
    PrimaryMatrix plus_primary_force{}, minus_primary_force{};
};

struct History final {
    double area = 0.0;
    std::array<double, 3> force{}, moment{}, center{};
};

struct TransferState final {
    Vector8 secondary_copen{}, secondary_pressure{}, secondary_force{};
    PrimaryVector primary_force{};
    std::size_t secondary_rows = 0;
    std::size_t primary_rows = 0;
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
    if (!line.empty() && line.back() == ',')
        values.emplace_back();
    return values;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size() || values[index].empty())
        throw std::invalid_argument("Incomplete H20.27 row in " + path);
    return std::stod(values[index]);
}

MatchingData read_matching(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.26 operator reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "step,input_local_node,input_label,output_local_node,output_label,closure_delta_m,copen_m,"
           "cpress_pa,cnormf_x_n,cnormf_y_n,cnormf_z_n")
        throw std::invalid_argument("Unexpected H20.26 operator header in " + path);
    MatchingData result;
    std::size_t rows = 0;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.at(0) == "BASE")
            continue;
        const std::size_t input_node = static_cast<std::size_t>(number(values, 1, path)) - 1;
        const std::size_t output_node = static_cast<std::size_t>(number(values, 3, path)) - 1;
        if (input_node >= matching_count || output_node >= matching_count)
            throw std::invalid_argument("Invalid H20.26 local node in " + path);
        const bool plus = values.at(0).find("_PLUS") != std::string::npos;
        const bool minus = values.at(0).find("_MINUS") != std::string::npos;
        if (plus == minus)
            throw std::invalid_argument("Invalid H20.26 perturbation step in " + path);
        (plus ? result.plus_copen : result.minus_copen)[output_node * matching_count + input_node] =
            number(values, 6, path);
        ++rows;
    }
    if (rows != 2 * matching_count * matching_count)
        throw std::invalid_argument("Unexpected H20.26 perturbation row count in " + path);
    return result;
}

NonmatchingData read_nonmatching(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.27 operator reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "step,input_local_node,input_label,side,output_local_node,output_label,closure_delta_m,copen_m,"
           "cpress_pa,cnormf_x_n,cnormf_y_n,cnormf_z_n")
        throw std::invalid_argument("Unexpected H20.27 operator header in " + path);

    NonmatchingData result;
    std::size_t base_rows = 0, perturbation_rows = 0, open_rows = 0, reclose_rows = 0;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        const std::string& step = values.at(0);
        const bool secondary = values.at(3) == "secondary";
        if (!secondary && values.at(3) != "primary")
            throw std::invalid_argument("Invalid H20.27 side in " + path);
        const std::size_t output_node = static_cast<std::size_t>(number(values, 4, path)) - 1;
        if ((secondary && output_node >= secondary_count) || (!secondary && output_node >= primary_count))
            throw std::invalid_argument("Invalid H20.27 output node in " + path);
        const double pressure = number(values, 8, path);
        const double force = number(values, 9, path);

        if (step == "BASE" || step == "OPEN" || step == "RECLOSE") {
            Vector13* copen_state = nullptr;
            Vector13* pressure_state = nullptr;
            Vector13* force_state = nullptr;
            PrimaryVector* primary_state = nullptr;
            std::size_t* row_count = nullptr;
            if (step == "BASE") {
                copen_state = &result.base_copen;
                pressure_state = &result.base_pressure;
                force_state = &result.base_force;
                primary_state = &result.base_primary_force;
                row_count = &base_rows;
            } else if (step == "OPEN") {
                copen_state = &result.open_copen;
                pressure_state = &result.open_pressure;
                force_state = &result.open_force;
                primary_state = &result.open_primary_force;
                row_count = &open_rows;
            } else {
                copen_state = &result.reclose_copen;
                pressure_state = &result.reclose_pressure;
                force_state = &result.reclose_force;
                primary_state = &result.reclose_primary_force;
                row_count = &reclose_rows;
            }
            if (secondary) {
                (*copen_state)[output_node] = number(values, 7, path);
                (*pressure_state)[output_node] = pressure;
                (*force_state)[output_node] = force;
            } else {
                (*primary_state)[output_node] = force;
            }
            ++(*row_count);
            continue;
        }

        const std::size_t input_node = static_cast<std::size_t>(number(values, 1, path)) - 1;
        if (input_node >= secondary_count)
            throw std::invalid_argument("Invalid H20.27 input node in " + path);
        const bool plus = step.find("_PLUS") != std::string::npos;
        const bool minus = step.find("_MINUS") != std::string::npos;
        if (plus == minus)
            throw std::invalid_argument("Invalid H20.27 perturbation step in " + path);
        if (secondary) {
            const std::size_t entry = output_node * secondary_count + input_node;
            (plus ? result.plus_copen : result.minus_copen)[entry] = number(values, 7, path);
            (plus ? result.plus_pressure : result.minus_pressure)[entry] = pressure;
            (plus ? result.plus_force : result.minus_force)[entry] = force;
        } else {
            const std::size_t entry = output_node * secondary_count + input_node;
            (plus ? result.plus_primary_force : result.minus_primary_force)[entry] = force;
        }
        ++perturbation_rows;
    }
    const std::size_t state_rows = secondary_count + primary_count;
    if (base_rows != state_rows || open_rows != state_rows || reclose_rows != state_rows
        || perturbation_rows != 2 * secondary_count * state_rows)
        throw std::invalid_argument("Unexpected H20.27 operator row count in " + path);
    return result;
}

History read_history(const std::string& path, const std::string& requested_step) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.27 history reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line
        != "step,contact_area_m2,normal_force_x_n,normal_force_y_n,normal_force_z_n,"
           "normal_moment_x_nm,normal_moment_y_nm,normal_moment_z_nm,center_x_m,center_y_m,center_z_m")
        throw std::invalid_argument("Unexpected H20.27 history header in " + path);
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.at(0) != requested_step)
            continue;
        return {number(values, 1, path),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 8, path), number(values, 9, path), number(values, 10, path)}};
    }
    throw std::invalid_argument("Missing H20.27 history step " + requested_step + " in " + path);
}

TransferState read_transfer_state(const std::string& path, const std::string& requested_step) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.27 transfer reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "step,side,local_node,node_label,copen_m,cpress_pa,cnormf_x_n,cnormf_y_n,cnormf_z_n")
        throw std::invalid_argument("Unexpected H20.27 transfer header in " + path);
    TransferState result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.at(0) != requested_step)
            continue;
        const std::size_t node = static_cast<std::size_t>(number(values, 2, path)) - 1;
        if (values.at(1) == "secondary") {
            if (node >= matching_count)
                throw std::invalid_argument("Invalid H20.27 transfer secondary node");
            result.secondary_copen[node] = number(values, 4, path);
            result.secondary_pressure[node] = number(values, 5, path);
            result.secondary_force[node] = number(values, 6, path);
            ++result.secondary_rows;
        } else if (values.at(1) == "primary") {
            if (node >= primary_count)
                throw std::invalid_argument("Invalid H20.27 transfer primary node");
            result.primary_force[node] = number(values, 6, path);
            ++result.primary_rows;
        } else {
            throw std::invalid_argument("Invalid H20.27 transfer side in " + path);
        }
    }
    if (result.secondary_rows != matching_count || result.primary_rows != primary_count)
        throw std::invalid_argument("Incomplete H20.27 transfer step " + requested_step + " in " + path);
    return result;
}

Vector13 solve13(Matrix13 matrix, Vector13 right_hand_side) {
    for (std::size_t column = 0; column < secondary_count; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < secondary_count; ++row)
            if (std::abs(matrix[row * secondary_count + column]) > std::abs(matrix[pivot * secondary_count + column]))
                pivot = row;
        if (!(std::abs(matrix[pivot * secondary_count + column]) > 1.0e-14))
            throw std::runtime_error("Singular H20.27 identification system");
        if (pivot != column) {
            for (std::size_t entry = column; entry < secondary_count; ++entry)
                std::swap(matrix[column * secondary_count + entry], matrix[pivot * secondary_count + entry]);
            std::swap(right_hand_side[column], right_hand_side[pivot]);
        }
        const double diagonal = matrix[column * secondary_count + column];
        for (std::size_t row = column + 1; row < secondary_count; ++row) {
            const double factor = matrix[row * secondary_count + column] / diagonal;
            for (std::size_t entry = column; entry < secondary_count; ++entry)
                matrix[row * secondary_count + entry] -= factor * matrix[column * secondary_count + entry];
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    Vector13 result{};
    for (std::size_t reverse = 0; reverse < secondary_count; ++reverse) {
        const std::size_t row = secondary_count - 1 - reverse;
        double value = right_hand_side[row];
        for (std::size_t column = row + 1; column < secondary_count; ++column)
            value -= matrix[row * secondary_count + column] * result[column];
        result[row] = value / matrix[row * secondary_count + row];
    }
    return result;
}

double relative_frobenius13(const Matrix13& value, const Matrix13& reference) {
    double difference_squared = 0.0, reference_squared = 0.0;
    for (std::size_t entry = 0; entry < value.size(); ++entry) {
        difference_squared += std::pow(value[entry] - reference[entry], 2);
        reference_squared += reference[entry] * reference[entry];
    }
    return std::sqrt(difference_squared / reference_squared);
}

double relative_frobenius_primary(const PrimaryMatrix& value, const PrimaryMatrix& reference) {
    double difference_squared = 0.0, reference_squared = 0.0;
    for (std::size_t entry = 0; entry < value.size(); ++entry) {
        difference_squared += std::pow(value[entry] - reference[entry], 2);
        reference_squared += reference[entry] * reference[entry];
    }
    return std::sqrt(difference_squared / reference_squared);
}

double relative_frobenius_transfer(const TransferMatrix& value, const TransferMatrix& reference) {
    double difference_squared = 0.0, reference_squared = 0.0;
    for (std::size_t entry = 0; entry < value.size(); ++entry) {
        difference_squared += std::pow(value[entry] - reference[entry], 2);
        reference_squared += reference[entry] * reference[entry];
    }
    return std::sqrt(difference_squared / reference_squared);
}

std::array<double, matching_count> quad8_shape(double xi, double eta) {
    return {0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0),
        0.25 * (1.0 + xi) * (1.0 - eta) * (xi - eta - 1.0),
        0.25 * (1.0 + xi) * (1.0 + eta) * (xi + eta - 1.0),
        0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0),
        0.5 * (1.0 - xi * xi) * (1.0 - eta),
        0.5 * (1.0 + xi) * (1.0 - eta * eta),
        0.5 * (1.0 - xi * xi) * (1.0 + eta),
        0.5 * (1.0 - xi) * (1.0 - eta * eta)};
}

TransferMatrix point_projection_candidate(const Matrix13& averaging) {
    const std::array<double, secondary_count> y{0.0, 1.2, 1.2, 0.0, 0.6, 1.2, 0.6, 0.0, 2.0, 2.0, 1.6, 2.0, 1.6};
    const std::array<double, secondary_count> z{0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0, 1.0, 0.0, 0.5, 1.0};
    const std::array<std::array<std::size_t, matching_count>, 4> face_nodes{{{0, 1, 2, 3, 4, 5, 6, 7},
        {1, 8, 9, 2, 10, 11, 12, 5},
        {8, 13, 14, 9, 15, 16, 17, 11},
        {13, 18, 19, 14, 20, 21, 22, 16}}};
    TransferMatrix projection{};
    for (std::size_t node = 0; node < secondary_count; ++node) {
        const std::size_t face = std::min(static_cast<std::size_t>(y[node] / 0.5), std::size_t{3});
        const double xi = 2.0 * (y[node] - 0.5 * static_cast<double>(face)) / 0.5 - 1.0;
        const std::array<double, matching_count> shape = quad8_shape(xi, 2.0 * z[node] - 1.0);
        for (std::size_t local = 0; local < matching_count; ++local)
            projection[node * primary_count + face_nodes[face][local]] += shape[local];
    }
    TransferMatrix result{};
    for (std::size_t row = 0; row < secondary_count; ++row)
        for (std::size_t column = 0; column < primary_count; ++column)
            for (std::size_t node = 0; node < secondary_count; ++node)
                result[row * primary_count + column] +=
                    averaging[row * secondary_count + node] * projection[node * primary_count + column];
    return result;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: fuelsim_h20_27_abaqus_sts_nonmatching_identification_tests <matching_operator.csv> "
                     "<nonmatching_operator.csv> <nonmatching_history.csv> <transfer_nodal.csv> "
                     "<transfer_history.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        const MatchingData matching_data = read_matching(argv[1]);
        const NonmatchingData data = read_nonmatching(argv[2]);
        Matrix8 matching_averaging{};
        for (std::size_t entry = 0; entry < matching_averaging.size(); ++entry)
            matching_averaging[entry] =
                -(matching_data.plus_copen[entry] - matching_data.minus_copen[entry]) / (2.0 * perturbation);

        Matrix13 averaging{}, secondary_tangent{}, pressure_derivative{};
        PrimaryMatrix primary_tangent{};
        for (std::size_t entry = 0; entry < averaging.size(); ++entry) {
            averaging[entry] = -(data.plus_copen[entry] - data.minus_copen[entry]) / (2.0 * perturbation);
            secondary_tangent[entry] = (data.plus_force[entry] - data.minus_force[entry]) / (2.0 * perturbation);
            pressure_derivative[entry] =
                (data.plus_pressure[entry] - data.minus_pressure[entry]) / (2.0 * perturbation);
        }
        for (std::size_t entry = 0; entry < primary_tangent.size(); ++entry)
            primary_tangent[entry] =
                (data.plus_primary_force[entry] - data.minus_primary_force[entry]) / (2.0 * perturbation);

        double row_sum_error = 0.0, asymmetry = 0.0;
        for (std::size_t row = 0; row < secondary_count; ++row) {
            double row_sum = 0.0;
            for (std::size_t column = 0; column < secondary_count; ++column) {
                row_sum += averaging[row * secondary_count + column];
                asymmetry = std::max(asymmetry,
                    std::abs(averaging[row * secondary_count + column] - averaging[column * secondary_count + row]));
            }
            row_sum_error = std::max(row_sum_error, std::abs(row_sum - 1.0));
        }

        Matrix13 transpose_averaging{};
        Vector13 normalized_base_force{};
        for (std::size_t row = 0; row < secondary_count; ++row) {
            normalized_base_force[row] = data.base_force[row] / (penalty * base_closure);
            for (std::size_t column = 0; column < secondary_count; ++column)
                transpose_averaging[row * secondary_count + column] = averaging[column * secondary_count + row];
        }
        const Vector13 areas = solve13(transpose_averaging, normalized_base_force);
        const Vector13 expected_areas{1.0 / 20.0,
            1.0 / 12.0,
            1.0 / 12.0,
            1.0 / 20.0,
            1.0 / 4.0,
            5.0 / 12.0,
            1.0 / 4.0,
            1.0 / 4.0,
            1.0 / 30.0,
            1.0 / 30.0,
            1.0 / 6.0,
            1.0 / 6.0,
            1.0 / 6.0};
        double area_sum = 0.0, area_error = 0.0;
        for (std::size_t node = 0; node < secondary_count; ++node) {
            area_sum += areas[node];
            area_error = std::max(area_error, std::abs(areas[node] - expected_areas[node]));
        }

        const std::array<std::array<std::size_t, matching_count>, 2> secondary_faces{
            {{0, 1, 2, 3, 4, 5, 6, 7}, {1, 8, 9, 2, 10, 11, 12, 5}}};
        const std::array<double, 2> face_areas{1.2, 0.8};
        Matrix13 assembled_averaging{};
        for (std::size_t face = 0; face < secondary_faces.size(); ++face)
            for (std::size_t local_row = 0; local_row < matching_count; ++local_row) {
                const std::size_t global_row = secondary_faces[face][local_row];
                const double local_area = face_areas[face] * (local_row < 4 ? 1.0 / 24.0 : 5.0 / 24.0);
                for (std::size_t local_column = 0; local_column < matching_count; ++local_column) {
                    const std::size_t global_column = secondary_faces[face][local_column];
                    assembled_averaging[global_row * secondary_count + global_column] +=
                        local_area * matching_averaging[local_row * matching_count + local_column];
                }
            }
        for (std::size_t row = 0; row < secondary_count; ++row)
            for (std::size_t column = 0; column < secondary_count; ++column)
                assembled_averaging[row * secondary_count + column] /= expected_areas[row];
        const double assembled_averaging_error = relative_frobenius13(assembled_averaging, averaging);

        Matrix13 reconstructed_secondary{}, weighted_averaging{}, transpose_weighted_averaging{};
        for (std::size_t row = 0; row < secondary_count; ++row)
            for (std::size_t column = 0; column < secondary_count; ++column) {
                weighted_averaging[row * secondary_count + column] =
                    areas[row] * averaging[row * secondary_count + column];
                transpose_weighted_averaging[column * secondary_count + row] =
                    weighted_averaging[row * secondary_count + column];
                for (std::size_t constraint = 0; constraint < secondary_count; ++constraint)
                    reconstructed_secondary[row * secondary_count + column] +=
                        penalty * averaging[constraint * secondary_count + row] * areas[constraint]
                        * averaging[constraint * secondary_count + column];
            }
        const double secondary_factorization_error = relative_frobenius13(reconstructed_secondary, secondary_tangent);

        TransferMatrix transfer{};
        for (std::size_t primary = 0; primary < primary_count; ++primary) {
            Vector13 right_hand_side{};
            for (std::size_t input_node = 0; input_node < secondary_count; ++input_node)
                right_hand_side[input_node] = -primary_tangent[primary * secondary_count + input_node] / penalty;
            const Vector13 column = solve13(transpose_weighted_averaging, right_hand_side);
            for (std::size_t constraint = 0; constraint < secondary_count; ++constraint)
                transfer[constraint * primary_count + primary] = column[constraint];
        }
        double transfer_row_sum_error = 0.0, transfer_minimum = transfer[0], transfer_maximum = transfer[0];
        for (std::size_t constraint = 0; constraint < secondary_count; ++constraint) {
            double row_sum = 0.0;
            for (std::size_t primary = 0; primary < primary_count; ++primary) {
                const double value = transfer[constraint * primary_count + primary];
                row_sum += value;
                transfer_minimum = std::min(transfer_minimum, value);
                transfer_maximum = std::max(transfer_maximum, value);
            }
            transfer_row_sum_error = std::max(transfer_row_sum_error, std::abs(row_sum - 1.0));
        }

        PrimaryMatrix reconstructed_primary{};
        PrimaryVector reconstructed_base_primary{};
        for (std::size_t primary = 0; primary < primary_count; ++primary)
            for (std::size_t constraint = 0; constraint < secondary_count; ++constraint) {
                reconstructed_base_primary[primary] -=
                    penalty * base_closure * transfer[constraint * primary_count + primary] * areas[constraint];
                for (std::size_t input_node = 0; input_node < secondary_count; ++input_node)
                    reconstructed_primary[primary * secondary_count + input_node] -=
                        penalty * transfer[constraint * primary_count + primary]
                        * weighted_averaging[constraint * secondary_count + input_node];
            }
        const double primary_factorization_error = relative_frobenius_primary(reconstructed_primary, primary_tangent);
        double primary_base_error = 0.0, tangent_conservation_error = 0.0, tangent_scale = 0.0;
        for (std::size_t primary = 0; primary < primary_count; ++primary)
            primary_base_error = std::max(primary_base_error,
                std::abs(reconstructed_base_primary[primary] - data.base_primary_force[primary]));
        for (std::size_t input_node = 0; input_node < secondary_count; ++input_node) {
            double sum = 0.0, scale = 0.0;
            for (std::size_t secondary = 0; secondary < secondary_count; ++secondary) {
                const double value = secondary_tangent[secondary * secondary_count + input_node];
                sum += value;
                scale += std::abs(value);
            }
            for (std::size_t primary = 0; primary < primary_count; ++primary) {
                const double value = primary_tangent[primary * secondary_count + input_node];
                sum += value;
                scale += std::abs(value);
            }
            tangent_conservation_error = std::max(tangent_conservation_error, std::abs(sum));
            tangent_scale = std::max(tangent_scale, scale);
        }
        const double tangent_conservation_relative = tangent_conservation_error / tangent_scale;

        Matrix13 penalty_averaging{};
        for (std::size_t entry = 0; entry < penalty_averaging.size(); ++entry)
            penalty_averaging[entry] = penalty * averaging[entry];
        const double displayed_pressure_error = relative_frobenius13(pressure_derivative, penalty_averaging);
        const TransferMatrix projection = point_projection_candidate(averaging);
        const double point_projection_error = relative_frobenius_transfer(projection, transfer);

        double open_error = 0.0, reclose_error = 0.0;
        for (std::size_t node = 0; node < matching_count; ++node) {
            open_error = std::max(open_error, std::abs(data.open_copen[node] - base_closure));
            open_error = std::max(open_error, std::abs(data.open_pressure[node]));
            open_error = std::max(open_error, std::abs(data.open_force[node]));
            reclose_error = std::max(reclose_error, std::abs(data.reclose_copen[node] - data.base_copen[node]));
            reclose_error = std::max(reclose_error, std::abs(data.reclose_pressure[node] - data.base_pressure[node]));
            reclose_error = std::max(reclose_error, std::abs(data.reclose_force[node] - data.base_force[node]));
        }
        for (std::size_t node = 0; node < primary_count; ++node) {
            open_error = std::max(open_error, std::abs(data.open_primary_force[node]));
            reclose_error =
                std::max(reclose_error, std::abs(data.reclose_primary_force[node] - data.base_primary_force[node]));
        }

        const History base_history = read_history(argv[3], "BASE");
        const History open_history = read_history(argv[3], "OPEN");
        const History reclose_history = read_history(argv[3], "RECLOSE");
        const TransferState transfer_base = read_transfer_state(argv[4], "BASE");
        const TransferState transfer_shift = read_transfer_state(argv[4], "SHIFT");
        const TransferState transfer_open = read_transfer_state(argv[4], "OPEN");
        const TransferState transfer_reclose = read_transfer_state(argv[4], "RECLOSE");
        const TransferState transfer_return = read_transfer_state(argv[4], "RETURN");
        const History transfer_base_history = read_history(argv[5], "BASE");
        const History transfer_shift_history = read_history(argv[5], "SHIFT");
        const History transfer_open_history = read_history(argv[5], "OPEN");
        const History transfer_reclose_history = read_history(argv[5], "RECLOSE");
        const History transfer_return_history = read_history(argv[5], "RETURN");

        double shift_primary_change = 0.0, shift_secondary_change = 0.0, transfer_open_error = 0.0,
               transfer_reclose_error = 0.0, transfer_return_error = 0.0;
        for (std::size_t node = 0; node < matching_count; ++node) {
            shift_secondary_change = std::max(shift_secondary_change,
                std::abs(transfer_shift.secondary_force[node] - transfer_base.secondary_force[node]));
            transfer_open_error =
                std::max(transfer_open_error, std::abs(transfer_open.secondary_copen[node] - base_closure));
            transfer_open_error = std::max(transfer_open_error, std::abs(transfer_open.secondary_pressure[node]));
            transfer_open_error = std::max(transfer_open_error, std::abs(transfer_open.secondary_force[node]));
            transfer_reclose_error = std::max(transfer_reclose_error,
                std::abs(transfer_reclose.secondary_force[node] - transfer_shift.secondary_force[node]));
            transfer_return_error = std::max(transfer_return_error,
                std::abs(transfer_return.secondary_force[node] - transfer_base.secondary_force[node]));
        }
        for (std::size_t node = 0; node < primary_count; ++node) {
            shift_primary_change = std::max(shift_primary_change,
                std::abs(transfer_shift.primary_force[node] - transfer_base.primary_force[node]));
            transfer_open_error = std::max(transfer_open_error, std::abs(transfer_open.primary_force[node]));
            transfer_reclose_error = std::max(transfer_reclose_error,
                std::abs(transfer_reclose.primary_force[node] - transfer_shift.primary_force[node]));
            transfer_return_error = std::max(transfer_return_error,
                std::abs(transfer_return.primary_force[node] - transfer_base.primary_force[node]));
        }

        std::cout << "h20_27_averaging_row_sum_max_error=" << row_sum_error << '\n'
                  << "h20_27_averaging_max_asymmetry=" << asymmetry << '\n'
                  << "h20_27_constraint_area_sum=" << area_sum << '\n'
                  << "h20_27_constraint_area_max_absolute_error=" << area_error << '\n'
                  << "h20_27_facewise_assembled_averaging_relative_frobenius_error=" << assembled_averaging_error
                  << '\n'
                  << "h20_27_secondary_factorization_relative_frobenius_error=" << secondary_factorization_error << '\n'
                  << "h20_27_primary_factorization_relative_frobenius_error=" << primary_factorization_error << '\n'
                  << "h20_27_primary_base_force_max_absolute_error=" << primary_base_error << '\n'
                  << "h20_27_transfer_row_sum_max_error=" << transfer_row_sum_error << '\n'
                  << "h20_27_transfer_minimum=" << transfer_minimum << '\n'
                  << "h20_27_transfer_maximum=" << transfer_maximum << '\n'
                  << "h20_27_tangent_conservation_relative_error=" << tangent_conservation_relative << '\n'
                  << "h20_27_displayed_pressure_vs_constraint_relative_frobenius_error=" << displayed_pressure_error
                  << '\n'
                  << "h20_27_point_projection_transfer_relative_frobenius_error=" << point_projection_error << '\n'
                  << "h20_27_release_max_absolute_error=" << open_error << '\n'
                  << "h20_27_recontact_max_absolute_error=" << reclose_error << '\n'
                  << "h20_27_small_sliding_center_y_shift="
                  << transfer_shift_history.center[1] - transfer_base_history.center[1] << '\n'
                  << "h20_27_small_sliding_primary_force_max_change=" << shift_primary_change << '\n'
                  << "h20_27_small_sliding_secondary_force_max_change=" << shift_secondary_change << '\n'
                  << "h20_27_small_sliding_release_max_absolute_error=" << transfer_open_error << '\n'
                  << "h20_27_small_sliding_recontact_max_absolute_error=" << transfer_reclose_error << '\n'
                  << "h20_27_small_sliding_return_max_absolute_error=" << transfer_return_error << '\n';

        bool passed = true;
        passed = check(row_sum_error < 1.0e-9 && asymmetry > 0.4,
                     "H20.27 nonmatching averaged constraints preserve rigid closure and remain nonsymmetric")
                 && passed;
        passed = check(area_error < 1.0e-9 && std::abs(area_sum - 2.0) < 1.0e-9,
                     "H20.27 constraint areas are the positive facewise assembled C3D20 values")
                 && passed;
        passed = check(assembled_averaging_error < 1.0e-9,
                     "H20.27 nonmatching secondary averaging is exactly assembled from the H20.26 face rule")
                 && passed;
        passed = check(secondary_factorization_error < 1.0e-9,
                     "H20.27 secondary tangent factors as penalty times A-transpose W A")
                 && passed;
        passed = check(primary_factorization_error < 1.0e-12 && primary_base_error < 1.0e-7
                           && transfer_row_sum_error < 1.0e-8 && tangent_conservation_relative < 1.0e-11,
                     "H20.27 inferred primary transfer reproduces tangent, baseline force, and conservation")
                 && passed;
        passed = check(displayed_pressure_error > 0.35,
                     "H20.27 displayed CPRESS remains distinct from the averaged constraint operator")
                 && passed;
        passed = check(point_projection_error > 0.75,
                     "H20.27 primary transfer cannot be replaced by point projection of secondary nodes")
                 && passed;
        passed = check(open_error < 1.0e-10 && reclose_error < 1.0e-8 && open_history.area == 0.0
                           && open_history.force[0] == 0.0 && std::abs(base_history.force[0] + 2.0e4) < 1.0e-8
                           && std::abs(reclose_history.force[0] - base_history.force[0]) < 1.0e-8,
                     "H20.27 nonmatching contact releases to zero and exactly recovers on recontact")
                 && passed;
        passed = check(std::abs(transfer_shift_history.center[1] - transfer_base_history.center[1] - 0.1) < 1.0e-7
                           && shift_primary_change < 1.0e-8 && shift_secondary_change < 1.0e-8
                           && std::abs(transfer_shift_history.force[0] - transfer_base_history.force[0]) < 1.0e-8,
                     "H20.27 small sliding moves the force center but keeps its original primary transfer anchor")
                 && passed;
        passed = check(transfer_open_error < 1.0e-10 && transfer_open_history.area == 0.0
                           && transfer_reclose_error < 1.0e-8 && transfer_return_error < 1.0e-8
                           && std::abs(transfer_reclose_history.force[0] - transfer_shift_history.force[0]) < 1.0e-8
                           && std::abs(transfer_return_history.force[0] - transfer_base_history.force[0]) < 1.0e-8,
                     "H20.27 small sliding preserves fixed anchors through release, recontact, and return")
                 && passed;
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
