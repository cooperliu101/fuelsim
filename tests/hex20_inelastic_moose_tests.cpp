#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
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
struct NodeRow {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 4> values{};
};

struct ElementRow {
    std::size_t id;
    std::array<double, 8> values{};
};

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

std::size_t column(const std::vector<std::string>& header, const std::string& name, const std::string& path) {
    const auto it = std::find(header.begin(), header.end(), name);
    if (it == header.end()) throw std::invalid_argument("Missing column " + name + " in " + path);
    return static_cast<std::size_t>(it - header.begin());
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete row in " + path);
    return std::stod(values[index]);
}

std::vector<NodeRow> read_nodes(const std::string& path, bool temperature) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE reference " + path);
    std::string line;
    std::getline(input, line);
    const auto header = split(line);
    const std::size_t id = column(header, "id", path), x = column(header, "x", path), y = column(header, "y", path),
                      z = column(header, "z", path);
    std::array<std::size_t, 4> field{};
    if (temperature)
        field[0] = column(header, "T", path);
    else {
        field[1] = column(header, "disp_x", path);
        field[2] = column(header, "disp_y", path);
        field[3] = column(header, "disp_z", path);
    }
    std::vector<NodeRow> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        NodeRow row{};
        row.id = static_cast<std::size_t>(number(values, id, path));
        row.point = {number(values, x, path), number(values, y, path), number(values, z, path)};
        if (temperature)
            row.values[0] = number(values, field[0], path);
        else
            for (std::size_t c = 1; c < 4; ++c) row.values[c] = number(values, field[c], path);
        result.push_back(row);
    }
    return result;
}

std::vector<ElementRow> read_elements(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE element reference " + path);
    std::string line;
    std::getline(input, line);
    const auto header = split(line);
    const std::size_t id = column(header, "id", path);
    const std::array<std::string, 8> names = {"stress_xx", "stress_yy", "stress_zz", "stress_xy", "stress_yz",
        "stress_xz", "effective_plastic_strain", "effective_creep_strain"};
    std::array<std::size_t, 8> fields{};
    for (std::size_t i = 0; i < fields.size(); ++i) fields[i] = column(header, names[i], path);
    std::vector<ElementRow> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split(line);
        ElementRow row{};
        row.id = static_cast<std::size_t>(number(values, id, path));
        for (std::size_t i = 0; i < fields.size(); ++i) row.values[i] = number(values, fields[i], path);
        result.push_back(row);
    }
    return result;
}

bool check(bool condition, const std::string& message) {
    if (!condition) std::cerr << "[FAIL] " << message << '\n';
    return condition;
}

bool same(const fuelsim::CartesianMaterialPointState& a, const fuelsim::CartesianMaterialPointState& b) {
    return a.elastic_strain == b.elastic_strain && a.plastic_strain == b.plastic_strain &&
           a.creep_strain == b.creep_strain && a.equivalent_plastic_strain == b.equivalent_plastic_strain &&
           a.equivalent_creep_strain == b.equivalent_creep_strain && a.stress.xx == b.stress.xx &&
           a.stress.yy == b.stress.yy && a.stress.zz == b.stress.zz && a.stress.xy == b.stress.xy &&
           a.stress.yz == b.stress.yz && a.stress.xz == b.stress.xz;
}

bool same_state(const fuelsim::TransientCommittedState& a, const fuelsim::TransientCommittedState& b) {
    if (a.time != b.time || a.load_factor != b.load_factor || a.solution != b.solution ||
        a.cartesian_material_histories.size() != b.cartesian_material_histories.size())
        return false;
    for (std::size_t r = 0; r < a.cartesian_material_histories.size(); ++r)
        for (std::size_t e = 0; e < a.cartesian_material_histories[r].size(); ++e)
            for (std::size_t q = 0; q < 27; ++q)
                if (!same(a.cartesian_material_histories[r][e][q], b.cartesian_material_histories[r][e][q]))
                    return false;
    return true;
}

bool run_branch(const std::string& name, const std::string& case_path, const std::string& temperature_path,
    const std::string& displacement_path, const std::string& element_path) {
    const auto definition = fuelsim::read_case_input(case_path);
    if (definition.geometry != fuelsim::CaseGeometry::cartesian_3d ||
        definition.problem != fuelsim::CaseProblem::transient)
        throw std::invalid_argument("HEX20 inelastic case must be transient Cartesian");
    const auto mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, mesh);
    const auto initial = fuelsim::cartesian::ProblemAccess::committed_state(problem);
    problem.begin_time_step({0.1, 0.1});
    fuelsim::ContributionWorkspace workspace;
    problem.evaluate_contribution(0, problem.committed_solution(), workspace, true);
    problem.evaluate_contribution(0, problem.committed_solution(), workspace, false);
    bool passed = check(same_state(initial, fuelsim::cartesian::ProblemAccess::committed_state(problem)),
        name + " residual/Jacobian preserve committed HEX20 history");
    problem.rollback_time_step();
    passed = check(same_state(initial, fuelsim::cartesian::ProblemAccess::committed_state(problem)),
                 name + " rollback restores committed HEX20 history") &&
             passed;
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    const auto& ex = definition.transient_execution;
    const auto result = fuelsim::solve_transient(problem,
        {ex.end_time, ex.initial_time_step, ex.minimum_time_step, ex.maximum_time_step, ex.growth_factor,
            ex.cutback_factor, ex.maximum_cutbacks_per_step, ex.load_ramp_time},
        options);
    passed = check(result.completed && result.accepted_steps.size() == 10,
                 name + " accepts ten fixed Backward Euler steps") &&
             passed;
    const auto temperature = read_nodes(temperature_path, true), displacement = read_nodes(displacement_path, false);
    if (temperature.size() != 8 || displacement.size() != mesh.nodes().size())
        throw std::invalid_argument("Unexpected HEX20 MOOSE node counts");
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& fields = spatial.field_layout();
    std::array<fuelsim::test::FieldErrorMetrics, 4> nodal;
    double coordinate_error = 0;
    for (const auto& row : temperature) {
        const auto& point = mesh.nodes().at(row.id);
        coordinate_error = std::max(coordinate_error,
            std::max(
                {std::abs(point.x - row.point.x), std::abs(point.y - row.point.y), std::abs(point.z - row.point.z)}));
        nodal[0].add(problem.committed_solution().at(fields[0].begin + spatial.global_temperature_node(0, row.id)),
            row.values[0]);
    }
    for (const auto& row : displacement) {
        const auto& point = mesh.nodes().at(row.id);
        coordinate_error = std::max(coordinate_error,
            std::max(
                {std::abs(point.x - row.point.x), std::abs(point.y - row.point.y), std::abs(point.z - row.point.z)}));
        const auto global = spatial.global_node(0, row.id);
        for (std::size_t c = 0; c < 3; ++c)
            nodal[c + 1].add(problem.committed_solution().at(fields[c + 1].begin + global), row.values[c + 1]);
    }
    for (std::size_t c = 0; c < 4; ++c) {
        fuelsim::test::print_relative_metrics(
            name + "_" +
                std::array<std::string, 4>{"temperature", "displacement_x", "displacement_y", "displacement_z"}[c],
            nodal[c]);
        passed =
            check(fuelsim::test::relative_metrics_below(nodal[c], 5e-3), name + " nodal field is below 0.5 percent") &&
            passed;
    }
    const auto elements = read_elements(element_path);
    if (elements.size() != 1) throw std::invalid_argument("Expected one MOOSE HEX20 element state");
    std::array<fuelsim::test::FieldErrorMetrics, 8> material;
    const auto& history = fuelsim::cartesian::ProblemAccess::material_history(problem, 0, 0);
    for (const auto& point : history)
        for (std::size_t c = 0; c < 8; ++c)
            material[c].add(c < 6 ? std::array<double, 6>{point.stress.xx, point.stress.yy, point.stress.zz,
                                        point.stress.xy, point.stress.yz, point.stress.xz}[c]
                                  : (c == 6 ? point.equivalent_plastic_strain : point.equivalent_creep_strain),
                elements[0].values[c]);
    for (std::size_t c = 0; c < 8; ++c) {
        const bool inactive = material[c].maximum_reference < 1.0e-14;
        if (inactive || (c >= 1 && c <= 5))
            fuelsim::test::print_absolute_metrics(name + "_material_" + std::to_string(c), material[c]);
        else
            fuelsim::test::print_relative_metrics(name + "_material_" + std::to_string(c), material[c]);
        const bool matches = inactive ? material[c].maximum_absolute_difference < 1.0e-10
                                      : (c >= 1 && c <= 5 ? material[c].maximum_absolute_difference < 1.0
                                                          : fuelsim::test::relative_metrics_below(material[c], 5e-3));
        passed = check(matches, name + " HEX20 27-point material state matches MOOSE") && passed;
    }
    return check(coordinate_error < 1e-14, name + " uses the MOOSE HEX20 coordinates") && passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 13) {
        std::cerr << "Usage: fuelsim_hex20_inelastic_moose_tests <plastic case,temp,disp,state> <creep "
                     "case,temp,disp,state> <coupled case,temp,disp,state>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim HEX20 inelastic MOOSE comparison\n");
        bool passed = true;
        for (std::size_t b = 0; b < 3; ++b)
            passed = run_branch(std::array<std::string, 3>{"plastic", "creep", "coupled"}[b], argv[1 + 4 * b],
                         argv[2 + 4 * b], argv[3 + 4 * b], argv[4 + 4 * b]) &&
                     passed;
        if (passed && session.rank() == 0) std::cout << "[PASS] HEX20 merged plastic/creep/coupled MOOSE comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] HEX20 inelastic comparison raised: " << error.what() << '\n';
        return 1;
    }
}
