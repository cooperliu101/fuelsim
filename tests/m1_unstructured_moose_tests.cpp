#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct NodalReference final {
    double temperature;
    double radial_displacement;
    double axial_displacement;
    double r;
    double z;
    bool present;
};

struct PressureReference final {
    double pressure;
    double z;
};

struct ErrorMetrics final {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    double maximum_actual = 0.0;
    double maximum_reference = 0.0;
    double maximum_pointwise_relative = 0.0;
    double maximum_zero_reference_difference = 0.0;
    std::size_t nonzero_reference_count = 0;
    std::size_t zero_reference_count = 0;

    void add(double actual, double reference) {
        const double difference = actual - reference;
        difference_squared += difference * difference;
        reference_squared += reference * reference;
        maximum_actual = std::max(maximum_actual, std::abs(actual));
        maximum_reference = std::max(maximum_reference, std::abs(reference));
        if (reference != 0.0) {
            maximum_pointwise_relative =
                std::max(maximum_pointwise_relative,
                         std::abs(difference) / std::abs(reference));
            ++nonzero_reference_count;
        } else {
            maximum_zero_reference_difference = std::max(
                maximum_zero_reference_difference, std::abs(difference));
            ++zero_reference_count;
        }
    }

    double relative_l2() const {
        if (!(reference_squared > 0.0))
            throw std::domain_error("Reference L2 norm must be positive");
        return std::sqrt(difference_squared / reference_squared);
    }

    double relative_absolute_peak() const {
        if (!(maximum_reference > 0.0))
            throw std::domain_error("Reference maximum norm must be positive");
        return std::abs(maximum_actual - maximum_reference) / maximum_reference;
    }

    double maximum_pointwise_relative_error() const {
        if (nonzero_reference_count == 0)
            throw std::domain_error(
                "Pointwise relative error requires a nonzero reference");
        return maximum_pointwise_relative;
    }
};

void print_metrics(const std::string& name, const ErrorMetrics& metrics) {
    std::cout << name << "_relative_l2=" << metrics.relative_l2() << '\n';
    std::cout << name
              << "_relative_absolute_peak=" << metrics.relative_absolute_peak()
              << '\n';
    std::cout << name << "_maximum_pointwise_relative="
              << metrics.maximum_pointwise_relative_error() << '\n';
    std::cout << name
              << "_zero_reference_count=" << metrics.zero_reference_count
              << '\n';
    std::cout << name << "_maximum_zero_reference_absolute_difference="
              << metrics.maximum_zero_reference_difference << '\n';
}

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool check_additional_metrics(const std::string& name,
                              const ErrorMetrics& metrics, double tolerance) {
    bool passed = check(metrics.relative_absolute_peak() < tolerance,
                        name + " relative absolute-peak error is below 1%");
    passed = check(metrics.maximum_pointwise_relative_error() < tolerance,
                   name + " maximum pointwise relative error is below 1%") &&
             passed;
    return passed;
}

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t comma = line.find(',', begin);
        fields.push_back(line.substr(begin, comma - begin));
        if (comma == std::string::npos)
            break;
        begin = comma + 1;
    }
    return fields;
}

std::size_t column_index(const std::vector<std::string>& header,
                         const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end())
        throw std::invalid_argument("CSV is missing column '" + name + "'");
    return static_cast<std::size_t>(found - header.begin());
}

double csv_value(const std::vector<std::string>& fields, std::size_t column,
                 const std::string& path) {
    if (column >= fields.size())
        throw std::invalid_argument("CSV row is too short in '" + path + "'");
    std::size_t parsed = 0;
    const double value = std::stod(fields[column], &parsed);
    if (parsed != fields[column].size() || !std::isfinite(value))
        throw std::invalid_argument("CSV contains an invalid number in '" +
                                    path + "'");
    return value;
}

std::size_t csv_id(const std::vector<std::string>& fields, std::size_t column,
                   const std::string& path) {
    const double value = csv_value(fields, column, path);
    if (value < 0.0 ||
        value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("CSV contains an invalid node ID in '" +
                                    path + "'");
    return static_cast<std::size_t>(value);
}

std::vector<NodalReference> read_nodal_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MOOSE nodal reference '" +
                                 path + "'");
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("MOOSE nodal reference is empty: " + path);
    const std::vector<std::string> header = split_csv_line(line);
    const std::size_t temperature = column_index(header, "T");
    const std::size_t radial = column_index(header, "disp_x");
    const std::size_t axial = column_index(header, "disp_y");
    const std::size_t id = column_index(header, "id");
    const std::size_t x = column_index(header, "x");
    const std::size_t y = column_index(header, "y");

    std::vector<NodalReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> fields = split_csv_line(line);
        const std::size_t node = csv_id(fields, id, path);
        if (node >= result.size())
            result.resize(node + 1, {0.0, 0.0, 0.0, 0.0, 0.0, false});
        if (result[node].present)
            throw std::invalid_argument("Duplicate node ID in '" + path + "'");
        result[node] = {csv_value(fields, temperature, path),
                        csv_value(fields, radial, path),
                        csv_value(fields, axial, path),
                        csv_value(fields, x, path),
                        csv_value(fields, y, path),
                        true};
    }
    if (result.empty() ||
        std::any_of(result.begin(), result.end(),
                    [](const NodalReference& value) { return !value.present; }))
        throw std::invalid_argument(
            "MOOSE nodal reference IDs must be contiguous");
    return result;
}

std::vector<PressureReference>
read_pressure_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MOOSE pressure reference '" +
                                 path + "'");
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("MOOSE pressure reference is empty: " +
                                    path);
    const std::vector<std::string> header = split_csv_line(line);
    const std::size_t pressure = column_index(header, "contact_pressure");
    const std::size_t y = column_index(header, "y");
    std::vector<PressureReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> fields = split_csv_line(line);
        result.push_back(
            {csv_value(fields, pressure, path), csv_value(fields, y, path)});
    }
    std::sort(result.begin(), result.end(),
              [](const PressureReference& lhs, const PressureReference& rhs) {
                  return lhs.z < rhs.z;
              });
    return result;
}

bool same_coordinate(double lhs, double rhs) {
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= 1.0e-12 * scale;
}

bool structured_conversion_is_rejected(
    const fuelsim::UnstructuredQuad4Mesh& mesh, const std::string& block) {
    try {
        (void)fuelsim::StructuredRzMesh::from_unstructured_block(mesh, block);
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

bool run_comparison(const std::string& input_path,
                    const std::string& nodal_reference_path,
                    const std::string& pressure_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(input_path);
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    bool passed =
        check(structured_conversion_is_rejected(source, "fuel") &&
                  structured_conversion_is_rejected(source, "clad"),
              "MOOSE comparison blocks are genuinely non-tensor Quad4 meshes");

    fuelsim::SteadyProblem problem(definition.steady_definition(), source);
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    const fuelsim::SteadyResult result = fuelsim::solve_steady(
        problem, definition.steady_execution.load_steps, options);
    passed = check(result.completed && result.solve.converged,
                   "non-tensor MOOSE mesh solve converged") &&
             passed;
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "non-tensor load path reuses one PETSc workspace") &&
             passed;

    const std::vector<NodalReference> reference =
        read_nodal_reference(nodal_reference_path);
    passed = check(reference.size() == source.nodes().size(),
                   "MOOSE reference covers every Exodus node") &&
             passed;

    ErrorMetrics temperature;
    ErrorMetrics radial;
    ErrorMetrics axial;
    std::size_t compared_nodes = 0;
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const fuelsim::RegionMesh& mesh = problem.region_mesh(region);
        for (std::size_t local_node = 0; local_node < mesh.nodes().size();
             ++local_node) {
            const std::size_t source_node = mesh.source_node_ids()[local_node];
            const NodalReference& expected = reference.at(source_node);
            const fuelsim::RzPoint& point = mesh.nodes()[local_node];
            if (!same_coordinate(point.r, expected.r) ||
                !same_coordinate(point.z, expected.z))
                throw std::invalid_argument(
                    "MOOSE and fuelsim nodal coordinates do not align");
            const std::size_t global_node =
                problem.region_node_offset(region) + local_node;
            temperature.add(
                result.solve.state[problem.dof_map().temperature(global_node)],
                expected.temperature);
            radial.add(
                result.solve
                    .state[problem.dof_map().radial_displacement(global_node)],
                expected.radial_displacement);
            axial.add(
                result.solve
                    .state[problem.dof_map().axial_displacement(global_node)],
                expected.axial_displacement);
            ++compared_nodes;
        }
    }
    passed = check(compared_nodes == reference.size(),
                   "every MOOSE node is compared exactly once") &&
             passed;

    const std::vector<PressureReference> pressure_reference =
        read_pressure_reference(pressure_reference_path);
    const std::vector<fuelsim::ContactNodeSummary> pressure_values =
        problem.summarize_contact_nodes(0, result.solve.state);
    passed = check(pressure_values.size() == pressure_reference.size(),
                   "MOOSE and fuelsim contact vectors have equal length") &&
             passed;
    ErrorMetrics pressure;
    for (std::size_t node = 0; node < pressure_reference.size(); ++node) {
        if (!same_coordinate(pressure_values.at(node).z,
                             pressure_reference[node].z))
            throw std::invalid_argument(
                "MOOSE and fuelsim contact coordinates do not align");
        pressure.add(pressure_values.at(node).pressure,
                     pressure_reference[node].pressure);
    }

    constexpr double tolerance = 1.0e-2;
    passed = check(temperature.relative_l2() < tolerance,
                   "temperature relative L2 error is below 1%") &&
             passed;
    passed = check(radial.relative_l2() < tolerance,
                   "radial displacement relative L2 error is below 1%") &&
             passed;
    passed = check(axial.relative_l2() < tolerance,
                   "axial displacement relative L2 error is below 1%") &&
             passed;
    passed = check(pressure.relative_l2() < tolerance,
                   "contact pressure relative L2 error is below 1%") &&
             passed;
    passed = check_additional_metrics("temperature", temperature, tolerance) &&
             passed;
    passed =
        check_additional_metrics("radial displacement", radial, tolerance) &&
        passed;
    passed = check_additional_metrics("axial displacement", axial, tolerance) &&
             passed;
    passed =
        check_additional_metrics("contact pressure", pressure, tolerance) &&
        passed;

    print_metrics("unstructured_temperature", temperature);
    print_metrics("unstructured_radial_displacement", radial);
    print_metrics("unstructured_axial_displacement", axial);
    print_metrics("unstructured_contact_pressure", pressure);
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_m1_unstructured_moose_tests "
                     "<case.fsi> <all_nodes.csv> <fuel_surface.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim non-tensor Quad4 MOOSE comparison\n");
        if (!run_comparison(argv[1], argv[2], argv[3]))
            return 1;
        std::cout << "[PASS] non-tensor Quad4 MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
