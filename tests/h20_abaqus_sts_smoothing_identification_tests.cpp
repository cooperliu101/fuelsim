#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t secondary_count = 8;
constexpr std::size_t primary_count = 225;
constexpr double perturbation = 1.0e-6;
constexpr double penalty = 1.0e8;
using Vector8 = std::array<double, secondary_count>;
using Matrix8 = std::array<double, secondary_count * secondary_count>;

struct Row final {
    std::string step, side;
    std::size_t input_label, output_label;
    double copen, force, y, z;
};

struct IdentifiedOperator final {
    Matrix8 averaging{}, area_matrix{};
    std::vector<double> transfer;
    std::array<double, primary_count> y{}, z{};
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

std::vector<Row> read_rows(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.28 operator reference: " + path);
    std::string line;
    std::getline(input, line);
    if (!line.empty() && line.back() == '\r')
        line.pop_back();
    if (line
        != "step,input_local_node,input_label,side,output_local_node,output_label,closure_delta_m,"
           "copen_m,cnormf_x_n,coord_y_m,coord_z_m")
        throw std::invalid_argument("Unexpected H20.28 operator header in " + path);
    std::vector<Row> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 11)
            throw std::invalid_argument("Incomplete H20.28 operator row in " + path);
        result.push_back({values[0],
            values[3],
            static_cast<std::size_t>(std::stoull(values[2])),
            static_cast<std::size_t>(std::stoull(values[5])),
            values[7].empty() ? 0.0 : std::stod(values[7]),
            std::stod(values[8]),
            std::stod(values[9]),
            std::stod(values[10])});
    }
    return result;
}

Row find_row(const std::vector<Row>& rows, const std::string& step, const std::string& side, std::size_t output_label) {
    const auto found = std::find_if(rows.begin(), rows.end(), [&](const Row& row) {
        return row.step == step && row.side == side && row.output_label == output_label;
    });
    if (found == rows.end())
        throw std::invalid_argument("Missing H20.28 perturbation row");
    return *found;
}

Vector8 solve8(Matrix8 matrix, Vector8 right_hand_side) {
    for (std::size_t column = 0; column < secondary_count; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < secondary_count; ++row)
            if (std::abs(matrix[row * secondary_count + column]) > std::abs(matrix[pivot * secondary_count + column]))
                pivot = row;
        if (!(std::abs(matrix[pivot * secondary_count + column]) > 1.0e-14))
            throw std::runtime_error("Singular H20.28 identification system");
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
    Vector8 result{};
    for (std::size_t reverse = 0; reverse < secondary_count; ++reverse) {
        const std::size_t row = secondary_count - 1 - reverse;
        double value = right_hand_side[row];
        for (std::size_t column = row + 1; column < secondary_count; ++column)
            value -= matrix[row * secondary_count + column] * result[column];
        result[row] = value / matrix[row * secondary_count + row];
    }
    return result;
}

Matrix8 inverse8(const Matrix8& matrix) {
    Matrix8 result{};
    for (std::size_t column = 0; column < secondary_count; ++column) {
        Vector8 right_hand_side{};
        right_hand_side[column] = 1.0;
        const Vector8 solution = solve8(matrix, right_hand_side);
        for (std::size_t row = 0; row < secondary_count; ++row)
            result[row * secondary_count + column] = solution[row];
    }
    return result;
}

Matrix8 multiply8(const Matrix8& first, const Matrix8& second) {
    Matrix8 result{};
    for (std::size_t row = 0; row < secondary_count; ++row)
        for (std::size_t column = 0; column < secondary_count; ++column)
            for (std::size_t inner = 0; inner < secondary_count; ++inner)
                result[row * secondary_count + column] +=
                    first[row * secondary_count + inner] * second[inner * secondary_count + column];
    return result;
}

Matrix8 transpose8(const Matrix8& matrix) {
    Matrix8 result{};
    for (std::size_t row = 0; row < secondary_count; ++row)
        for (std::size_t column = 0; column < secondary_count; ++column)
            result[row * secondary_count + column] = matrix[column * secondary_count + row];
    return result;
}

IdentifiedOperator identify(const std::string& path) {
    const std::vector<Row> rows = read_rows(path);
    std::set<std::size_t> secondary_set, primary_set;
    for (const Row& row : rows) {
        if (row.step != "BASE")
            continue;
        (row.side == "secondary" ? secondary_set : primary_set).insert(row.output_label);
    }
    if (secondary_set.size() != secondary_count || primary_set.size() != primary_count)
        throw std::invalid_argument("Unexpected H20.28 face-node counts in " + path);
    const std::vector<std::size_t> secondary_labels(secondary_set.begin(), secondary_set.end());
    const std::vector<std::size_t> primary_labels(primary_set.begin(), primary_set.end());
    Matrix8 secondary_tangent{};
    std::vector<double> primary_tangent(primary_count * secondary_count);
    IdentifiedOperator result;
    for (std::size_t input_node = 0; input_node < secondary_count; ++input_node) {
        const std::string prefix = "S" + std::to_string(input_node + 1);
        const std::string plus_step = prefix + "_PLUS", minus_step = prefix + "_MINUS";
        const Row input_row = find_row(rows, plus_step, "secondary", secondary_labels.front());
        const auto input_found = std::find(secondary_labels.begin(), secondary_labels.end(), input_row.input_label);
        if (input_found == secondary_labels.end())
            throw std::invalid_argument("Invalid H20.28 input label");
        const std::size_t column = static_cast<std::size_t>(input_found - secondary_labels.begin());
        for (std::size_t output = 0; output < secondary_count; ++output) {
            const Row plus = find_row(rows, plus_step, "secondary", secondary_labels[output]);
            const Row minus = find_row(rows, minus_step, "secondary", secondary_labels[output]);
            result.averaging[output * secondary_count + column] = -(plus.copen - minus.copen) / (2.0 * perturbation);
            secondary_tangent[output * secondary_count + column] = (plus.force - minus.force) / (2.0 * perturbation);
        }
        for (std::size_t output = 0; output < primary_count; ++output) {
            const Row plus = find_row(rows, plus_step, "primary", primary_labels[output]);
            const Row minus = find_row(rows, minus_step, "primary", primary_labels[output]);
            primary_tangent[output * secondary_count + column] = (plus.force - minus.force) / (2.0 * perturbation);
        }
    }
    const Matrix8 inverse = inverse8(result.averaging);
    Matrix8 normalized_tangent = secondary_tangent;
    for (double& value : normalized_tangent)
        value /= penalty;
    result.area_matrix = multiply8(multiply8(transpose8(inverse), normalized_tangent), inverse);
    const Matrix8 weighted_averaging = multiply8(result.area_matrix, result.averaging);
    const Matrix8 transpose_weighted = transpose8(weighted_averaging);
    result.transfer.resize(secondary_count * primary_count);
    for (std::size_t primary = 0; primary < primary_count; ++primary) {
        Vector8 right_hand_side{};
        for (std::size_t column = 0; column < secondary_count; ++column)
            right_hand_side[column] = primary_tangent[primary * secondary_count + column] / penalty;
        const Vector8 transfer = solve8(transpose_weighted, right_hand_side);
        for (std::size_t constraint = 0; constraint < secondary_count; ++constraint)
            result.transfer[constraint * primary_count + primary] = transfer[constraint];
    }
    for (std::size_t primary = 0; primary < primary_count; ++primary) {
        const Row row = find_row(rows, "BASE", "primary", primary_labels[primary]);
        result.y[primary] = row.y;
        result.z[primary] = row.z;
    }
    return result;
}

double relative_difference(const std::vector<double>& first, const std::vector<double>& second) {
    if (first.size() != second.size())
        throw std::logic_error("H20.28 comparison size mismatch");
    double difference_squared = 0.0, reference_squared = 0.0;
    for (std::size_t entry = 0; entry < first.size(); ++entry) {
        difference_squared += std::pow(first[entry] - second[entry], 2);
        reference_squared += second[entry] * second[entry];
    }
    return std::sqrt(difference_squared / reference_squared);
}

std::vector<double> values(const Matrix8& matrix) {
    return {matrix.begin(), matrix.end()};
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_h20_28_abaqus_sts_smoothing_identification_tests "
                     "<default.csv> <linear.csv> <quadratic.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        const IdentifiedOperator default_operator = identify(argv[1]);
        const IdentifiedOperator linear_operator = identify(argv[2]);
        const IdentifiedOperator quadratic_operator = identify(argv[3]);
        const double default_quadratic_averaging =
                         relative_difference(values(default_operator.averaging), values(quadratic_operator.averaging)),
                     default_quadratic_area = relative_difference(values(default_operator.area_matrix),
                         values(quadratic_operator.area_matrix)),
                     default_quadratic_transfer =
                         relative_difference(default_operator.transfer, quadratic_operator.transfer),
                     default_linear_averaging =
                         relative_difference(values(default_operator.averaging), values(linear_operator.averaging)),
                     default_linear_transfer = relative_difference(default_operator.transfer, linear_operator.transfer);
        double default_area_error = 0.0, linear_area_error = 0.0, off_diagonal = 0.0, row_sum_error = 0.0,
               transfer_sum_error = 0.0, center_error = 0.0;
        std::array<std::size_t, secondary_count> default_support{}, linear_support{};
        const std::array<std::array<double, 2>, secondary_count> expected_centers{{{0.125, 0.125},
            {0.875, 0.125},
            {0.125, 0.875},
            {0.875, 0.875},
            {0.5, 0.25},
            {0.5, 0.75},
            {0.25, 0.5},
            {0.75, 0.5}}};
        for (std::size_t row = 0; row < secondary_count; ++row) {
            const double expected_default = row < 4 ? 1.0 / 24.0 : 5.0 / 24.0;
            const double expected_linear = row < 4 ? 1.0 / 48.0 : 11.0 / 48.0;
            default_area_error = std::max(default_area_error,
                std::abs(default_operator.area_matrix[row * secondary_count + row] - expected_default));
            linear_area_error = std::max(linear_area_error,
                std::abs(linear_operator.area_matrix[row * secondary_count + row] - expected_linear));
            double averaging_sum = 0.0, transfer_sum = 0.0, center_y = 0.0, center_z = 0.0;
            for (std::size_t column = 0; column < secondary_count; ++column) {
                averaging_sum += default_operator.averaging[row * secondary_count + column];
                if (row != column)
                    off_diagonal =
                        std::max(off_diagonal, std::abs(default_operator.area_matrix[row * secondary_count + column]));
            }
            for (std::size_t primary = 0; primary < primary_count; ++primary) {
                const double default_value = default_operator.transfer[row * primary_count + primary];
                const double linear_value = linear_operator.transfer[row * primary_count + primary];
                transfer_sum += default_value;
                center_y += default_value * default_operator.y[primary];
                center_z += default_value * default_operator.z[primary];
                if (std::abs(default_value) > 1.0e-8)
                    ++default_support[row];
                if (std::abs(linear_value) > 1.0e-8)
                    ++linear_support[row];
            }
            row_sum_error = std::max(row_sum_error, std::abs(averaging_sum - 1.0));
            transfer_sum_error = std::max(transfer_sum_error, std::abs(transfer_sum - 1.0));
            center_error = std::max({center_error,
                std::abs(center_y - expected_centers[row][0]),
                std::abs(center_z - expected_centers[row][1])});
        }
        std::cout << "h20_28_default_quadratic_averaging_relative_difference=" << default_quadratic_averaging << '\n'
                  << "h20_28_default_quadratic_area_relative_difference=" << default_quadratic_area << '\n'
                  << "h20_28_default_quadratic_transfer_relative_difference=" << default_quadratic_transfer << '\n'
                  << "h20_28_default_linear_averaging_relative_difference=" << default_linear_averaging << '\n'
                  << "h20_28_default_linear_transfer_relative_difference=" << default_linear_transfer << '\n'
                  << "h20_28_default_constraint_area_max_absolute_error=" << default_area_error << '\n'
                  << "h20_28_linear_constraint_area_max_absolute_error=" << linear_area_error << '\n'
                  << "h20_28_area_matrix_off_diagonal_max_absolute=" << off_diagonal << '\n'
                  << "h20_28_averaging_row_sum_max_error=" << row_sum_error << '\n'
                  << "h20_28_transfer_row_sum_max_error=" << transfer_sum_error << '\n'
                  << "h20_28_default_effective_center_max_absolute_error=" << center_error << '\n'
                  << "h20_28_default_transfer_support_counts=";
        for (std::size_t row = 0; row < secondary_count; ++row)
            std::cout << (row == 0 ? "" : ",") << default_support[row];
        std::cout << '\n' << "h20_28_linear_transfer_support_counts=";
        for (std::size_t row = 0; row < secondary_count; ++row)
            std::cout << (row == 0 ? "" : ",") << linear_support[row];
        std::cout << '\n';

        bool passed = true;
        passed = check(default_quadratic_averaging < 1.0e-15 && default_quadratic_area < 1.0e-15
                           && default_quadratic_transfer < 1.0e-15,
                     "H20.28 Abaqus default sliding transition is exactly quadratic smoothing")
                 && passed;
        passed = check(default_linear_averaging > 0.25 && default_linear_transfer > 0.65,
                     "H20.28 linear smoothing changes both secondary averaging and primary transfer")
                 && passed;
        passed = check(default_area_error < 1.0e-10 && linear_area_error < 1.0e-10 && off_diagonal < 1.0e-10,
                     "H20.28 default and linear constraint areas have their identified positive diagonal weights")
                 && passed;
        passed = check(row_sum_error < 1.0e-9 && transfer_sum_error < 1.0e-9 && center_error < 1.0e-10,
                     "H20.28 default constraints preserve translation and the identified centers of action")
                 && passed;
        passed = check(default_support == std::array<std::size_t, secondary_count>{45, 45, 45, 45, 109, 109, 109, 109}
                           && linear_support == std::array<std::size_t, secondary_count>{8, 8, 8, 8, 44, 44, 44, 44},
                     "H20.28 primary redistribution support depends on quadratic versus linear smoothing")
                 && passed;
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
