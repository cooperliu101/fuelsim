#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct NodeReference final {
    fuelsim::CartesianPoint3 point;
    std::array<double, 4> fields;
};

struct ElementReference final {
    std::array<double, 3> values;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

std::size_t column(const std::vector<std::string>& header, const std::string& name, const std::string& path) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("MOOSE CSV is missing '" + name + "': " + path);
    return static_cast<std::size_t>(found - header.begin());
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("MOOSE CSV row is incomplete: " + path);
    const double value = std::stod(values[index]);
    if (!std::isfinite(value)) throw std::invalid_argument("MOOSE CSV contains a non-finite value: " + path);
    return value;
}

std::vector<NodeReference> read_nodes(const std::string& path, std::size_t node_count) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE nodal reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE nodal reference is empty: " + path);
    const auto header = split_csv(line);
    const std::array<std::size_t, 8> columns = {column(header, "id", path), column(header, "x", path),
        column(header, "y", path), column(header, "z", path), column(header, "T", path), column(header, "disp_x", path),
        column(header, "disp_y", path), column(header, "disp_z", path)};
    std::vector<NodeReference> result(node_count);
    std::vector<bool> present(node_count, false);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split_csv(line);
        const std::size_t id = static_cast<std::size_t>(number(values, columns[0], path));
        if (id >= node_count) throw std::invalid_argument("MOOSE node ID is outside the Exodus mesh: " + path);
        const NodeReference candidate = {
            {number(values, columns[1], path), number(values, columns[2], path), number(values, columns[3], path)},
            {number(values, columns[4], path), number(values, columns[5], path), number(values, columns[6], path),
                number(values, columns[7], path)}};
        if (present[id]) {
            const NodeReference& existing = result[id];
            if (std::abs(existing.point.x - candidate.point.x) > 1.0e-12 ||
                std::abs(existing.point.y - candidate.point.y) > 1.0e-12 ||
                std::abs(existing.point.z - candidate.point.z) > 1.0e-12 ||
                !std::equal(existing.fields.begin(), existing.fields.end(), candidate.fields.begin(),
                    [](double left, double right) { return std::abs(left - right) < 1.0e-12; }))
                throw std::invalid_argument("MOOSE shared-node rows disagree: " + path);
            continue;
        }
        result[id] = candidate;
        present[id] = true;
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("MOOSE nodal reference does not cover every Exodus node: " + path);
    return result;
}

std::vector<ElementReference> read_elements(const std::string& path, std::size_t element_count) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE element reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE element reference is empty: " + path);
    const auto header = split_csv(line);
    const std::array<std::size_t, 4> columns = {column(header, "id", path), column(header, "stress_xx", path),
        column(header, "effective_plastic_strain", path), column(header, "effective_creep_strain", path)};
    std::vector<ElementReference> result(element_count);
    std::vector<bool> present(element_count, false);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split_csv(line);
        const std::size_t id = static_cast<std::size_t>(number(values, columns[0], path));
        if (id >= element_count || present[id]) throw std::invalid_argument("Invalid MOOSE element ID: " + path);
        for (std::size_t value = 0; value < result[id].values.size(); ++value)
            result[id].values[value] = number(values, columns[value + 1], path);
        present[id] = true;
    }
    if (std::find(present.begin(), present.end(), false) != present.end())
        throw std::invalid_argument("MOOSE element reference does not cover every Exodus element: " + path);
    return result;
}

bool check_metrics(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics) {
    fuelsim::test::print_relative_metrics(name, metrics);
    return check(fuelsim::test::relative_metrics_below(metrics, 5.0e-3), name + " three errors are below 0.5 percent");
}

bool run(const std::string& input_path, const std::string& node_path, const std::string& element_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::transient ||
        definition.geometry != fuelsim::CaseGeometry::cartesian_3d)
        throw std::invalid_argument("B3.6 requires a transient Cartesian three-dimensional input card");
    const fuelsim::UnstructuredHex8Mesh mesh = fuelsim::read_exodus_hex8(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, mesh);
    fuelsim::SolverOptions options = {definition.solver.absolute_tolerance, definition.solver.relative_tolerance,
        definition.solver.step_tolerance, definition.solver.maximum_iterations};
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    const auto& execution = definition.transient_execution;
    const fuelsim::TransientResult solve = fuelsim::solve_transient(problem,
        {execution.end_time, execution.initial_time_step, execution.minimum_time_step, execution.maximum_time_step,
            execution.growth_factor, execution.cutback_factor, execution.maximum_cutbacks_per_step,
            execution.load_ramp_time},
        options);
    bool passed = check(solve.completed && solve.accepted_steps.size() == 40,
        "the shared-node fuel-plate path accepts forty Backward Euler time steps");
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    passed = check(dofs.node_count() == mesh.nodes().size() && problem.dof_count() == 4 * mesh.nodes().size(),
                 "the fuel plate has one four-field global node per Exodus node") &&
             passed;
    std::vector<std::size_t> occurrences(mesh.nodes().size(), 0U), source_global(mesh.nodes().size(), 0U);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const std::size_t global = dofs.global_node(region, local);
            if (occurrences[source]++ == 0U)
                source_global[source] = global;
            else
                passed =
                    check(source_global[source] == global, "each meat-clad interface node maps to one global node") &&
                    passed;
        }
    }
    const auto nodes = read_nodes(node_path, mesh.nodes().size());
    std::array<fuelsim::test::FieldErrorMetrics, 4> nodal;
    const std::array<fuelsim::Field, 4> nodal_fields = {fuelsim::Field::temperature, fuelsim::Field::displacement_x,
        fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};
    for (std::size_t source = 0; source < mesh.nodes().size(); ++source) {
        if (occurrences[source] == 0U) throw std::invalid_argument("fuelsim did not visit an Exodus node");
        const auto& point = mesh.nodes()[source];
        const auto& reference = nodes[source];
        if (std::max({std::abs(point.x - reference.point.x), std::abs(point.y - reference.point.y),
                std::abs(point.z - reference.point.z)}) > 1.0e-12)
            throw std::invalid_argument("MOOSE and fuelsim node coordinates differ");
        for (std::size_t field = 0; field < 4; ++field)
            nodal[field].add(problem.committed_solution()[dofs.dof(nodal_fields[field], source_global[source])],
                reference.fields[field]);
    }
    const std::array<std::string, 4> nodal_names = {
        "temperature", "displacement_x", "displacement_y", "displacement_z"};
    for (std::size_t field = 0; field < nodal.size(); ++field)
        passed = check_metrics("b36_" + nodal_names[field], nodal[field]) && passed;
    const auto elements = read_elements(element_path, mesh.elements().size());
    std::array<fuelsim::test::FieldErrorMetrics, 3> material;
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        for (std::size_t element = 0; element < region_mesh.elements().size(); ++element) {
            const std::size_t source = region_mesh.source_element_ids()[element];
            std::array<double, 3> average{};
            const auto& history = fuelsim::cartesian::ProblemAccess::material_history(problem, region, element);
            for (const auto& point : history) {
                average[0] += point.stress.xx;
                average[1] += point.equivalent_plastic_strain;
                average[2] += point.equivalent_creep_strain;
            }
            for (std::size_t value = 0; value < average.size(); ++value)
                material[value].add(
                    average[value] / static_cast<double>(history.size()), elements[source].values[value]);
        }
    }
    const std::array<std::string, 3> material_names = {
        "stress_xx", "equivalent_plastic_strain", "equivalent_creep_strain"};
    for (std::size_t value = 0; value < material.size(); ++value)
        passed = check_metrics("b36_" + material_names[value], material[value]) && passed;
    return check(material[1].maximum_reference > 0.0 && material[2].maximum_reference > 0.0,
               "both plasticity and creep are active in the meat-clad plate") &&
           passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_b36_plate_meat_clad_moose_tests <case.fsi> <nodes.csv> <elements.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B3.6 shared-node fuel-plate MOOSE comparison\n");
        if (!run(argv[1], argv[2], argv[3])) return 1;
        std::cout << "[PASS] B3.6 shared-node fuel-plate MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B3.6 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
