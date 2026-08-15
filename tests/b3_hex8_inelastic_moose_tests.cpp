#include "fuelsim/case_input.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct NodeReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 4> fields;
};

struct ElementReference final {
    std::size_t id;
    fuelsim::SymmetricTensor3Values stress;
    double equivalent_plastic_strain, equivalent_creep_strain;
};

bool same_point(const fuelsim::CartesianMaterialPointState& left, const fuelsim::CartesianMaterialPointState& right) {
    return left.elastic_strain == right.elastic_strain && left.plastic_strain == right.plastic_strain &&
           left.creep_strain == right.creep_strain &&
           left.equivalent_plastic_strain == right.equivalent_plastic_strain &&
           left.equivalent_creep_strain == right.equivalent_creep_strain && left.stress.xx == right.stress.xx &&
           left.stress.yy == right.stress.yy && left.stress.zz == right.stress.zz &&
           left.stress.xy == right.stress.xy && left.stress.yz == right.stress.yz && left.stress.xz == right.stress.xz;
}

bool same_committed_state(const fuelsim::TransientCommittedState& left, const fuelsim::TransientCommittedState& right) {
    if (left.time != right.time || left.load_factor != right.load_factor || left.solution != right.solution ||
        left.cartesian_material_histories.size() != right.cartesian_material_histories.size())
        return false;
    for (std::size_t region = 0; region < left.cartesian_material_histories.size(); ++region) {
        if (left.cartesian_material_histories[region].size() != right.cartesian_material_histories[region].size())
            return false;
        for (std::size_t element = 0; element < left.cartesian_material_histories[region].size(); ++element)
            for (std::size_t q = 0; q < 8; ++q)
                if (!same_point(left.cartesian_material_histories[region][element][q],
                        right.cartesian_material_histories[region][element][q]))
                    return false;
    }
    return true;
}

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
    const auto position = std::find(header.begin(), header.end(), name);
    if (position == header.end()) throw std::invalid_argument("Missing column '" + name + "' in " + path);
    return static_cast<std::size_t>(position - header.begin());
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete MOOSE row in " + path);
    return std::stod(values[index]);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read three-dimensional MOOSE nodes: " + path);
    std::string line;
    std::getline(input, line);
    const auto header = split_csv(line);
    const std::array<std::size_t, 8> columns = {column(header, "id", path), column(header, "x", path),
        column(header, "y", path), column(header, "z", path), column(header, "T", path), column(header, "disp_x", path),
        column(header, "disp_y", path), column(header, "disp_z", path)};
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, columns[0], path)),
            {number(values, columns[1], path), number(values, columns[2], path), number(values, columns[3], path)},
            {number(values, columns[4], path), number(values, columns[5], path), number(values, columns[6], path),
                number(values, columns[7], path)}});
    }
    return result;
}

std::vector<ElementReference> read_elements(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read three-dimensional MOOSE element states: " + path);
    std::string line;
    std::getline(input, line);
    const auto header = split_csv(line);
    const std::array<std::size_t, 9> columns = {column(header, "id", path), column(header, "stress_xx", path),
        column(header, "stress_yy", path), column(header, "stress_zz", path), column(header, "stress_xy", path),
        column(header, "stress_yz", path), column(header, "stress_xz", path),
        column(header, "effective_plastic_strain", path), column(header, "effective_creep_strain", path)};
    std::vector<ElementReference> result;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, columns[0], path)),
            {number(values, columns[1], path), number(values, columns[2], path), number(values, columns[3], path),
                number(values, columns[4], path), number(values, columns[5], path), number(values, columns[6], path)},
            number(values, columns[7], path), number(values, columns[8], path)});
    }
    return result;
}

bool check_relative(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics) {
    fuelsim::test::print_relative_metrics(name, metrics);
    return check(fuelsim::test::relative_metrics_below(metrics, 5.0e-3), name + " three errors are below 0.5 percent");
}

bool run(const std::string& branch, const std::string& input_path, const std::string& nodal_path,
    const std::string& element_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::transient ||
        definition.geometry != fuelsim::CaseGeometry::cartesian_3d)
        throw std::invalid_argument("Three-dimensional inelastic comparison requires a transient Cartesian input card");
    const auto mesh = fuelsim::read_exodus_hex8(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, mesh);
    const auto initial = fuelsim::cartesian::ProblemAccess::committed_state(problem);
    problem.begin_time_step({0.1, 0.1});
    fuelsim::ContributionWorkspace workspace;
    problem.evaluate_contribution(0, problem.committed_solution(), workspace, true);
    problem.evaluate_contribution(0, problem.committed_solution(), workspace, false);
    bool passed = check(same_committed_state(initial, fuelsim::cartesian::ProblemAccess::committed_state(problem)),
        branch + " repeated residual and Jacobian calls do not mutate committed histories");
    problem.rollback_time_step();
    passed = check(same_committed_state(initial, fuelsim::cartesian::ProblemAccess::committed_state(problem)),
                 branch + " rollback restores the complete committed state") &&
             passed;
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    const auto& execution = definition.transient_execution;
    const auto solve = fuelsim::solve_transient(problem,
        {execution.end_time, execution.initial_time_step, execution.minimum_time_step, execution.maximum_time_step,
            execution.growth_factor, execution.cutback_factor, execution.maximum_cutbacks_per_step,
            execution.load_ramp_time},
        options);
    passed = check(solve.completed && solve.accepted_steps.size() == 10,
                 branch + " accepts the same ten Backward Euler steps as MOOSE") &&
             passed;
    const std::string checkpoint_path = "/tmp/fuelsim_b3_hex8_" + branch + "_checkpoint.bin";
    fuelsim::write_transient_checkpoint(checkpoint_path, problem, 0.1);
    fuelsim::TransientProblem restored(definition.spatial, mesh);
    const double restored_time_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restored);
    passed = check(restored_time_step == 0.1 &&
                       same_committed_state(fuelsim::cartesian::ProblemAccess::committed_state(problem),
                           fuelsim::cartesian::ProblemAccess::committed_state(restored)),
                 branch + " checkpoint restores every active six-component material history exactly") &&
             check(std::remove(checkpoint_path.c_str()) == 0, branch + " checkpoint artifact is removed") && passed;
    const auto references = read_nodes(nodal_path);
    if (references.size() != mesh.nodes().size()) throw std::invalid_argument("MOOSE and fuelsim node counts differ");
    std::array<fuelsim::test::FieldErrorMetrics, 4> nodal;
    std::vector<bool> visited(references.size(), false);
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        const std::size_t offset = fuelsim::cartesian::ProblemAccess::region_node_offset(problem, region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            if (source >= references.size() || references[source].id != source || visited[source])
                throw std::invalid_argument("Three-dimensional source-node mapping is not unique");
            visited[source] = true;
            const auto& expected = references[source];
            const auto& point = mesh.nodes()[source];
            if (std::max({std::abs(point.x - expected.point.x), std::abs(point.y - expected.point.y),
                    std::abs(point.z - expected.point.z)}) >= 1.0e-12)
                throw std::invalid_argument("Three-dimensional source-node coordinates differ");
            nodal[0].add(problem.committed_solution()[dofs.dof(fuelsim::Field::temperature, offset + local)],
                expected.fields[0]);
            nodal[1].add(problem.committed_solution()[dofs.dof(fuelsim::Field::displacement_x, offset + local)],
                expected.fields[1]);
            nodal[2].add(problem.committed_solution()[dofs.dof(fuelsim::Field::displacement_y, offset + local)],
                expected.fields[2]);
            nodal[3].add(problem.committed_solution()[dofs.dof(fuelsim::Field::displacement_z, offset + local)],
                expected.fields[3]);
        }
    }
    for (std::size_t field = 0; field < nodal.size(); ++field)
        passed = check_relative(branch + "_" +
                                    std::array<std::string, 4>{
                                        "temperature", "displacement_x", "displacement_y", "displacement_z"}[field],
                     nodal[field]) &&
                 passed;
    const auto element_references = read_elements(element_path);
    std::array<fuelsim::test::FieldErrorMetrics, 8> material;
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        for (std::size_t element = 0; element < region_mesh.elements().size(); ++element) {
            const std::size_t source = region_mesh.source_element_ids()[element];
            if (source >= element_references.size() || element_references[source].id != source)
                throw std::invalid_argument("Three-dimensional source-element mapping is not unique");
            const auto& expected = element_references[source];
            for (const auto& point : fuelsim::cartesian::ProblemAccess::material_history(problem, region, element)) {
                material[0].add(point.stress.xx, expected.stress.xx);
                material[1].add(point.stress.yy, expected.stress.yy);
                material[2].add(point.stress.zz, expected.stress.zz);
                material[3].add(point.stress.xy, expected.stress.xy);
                material[4].add(point.stress.yz, expected.stress.yz);
                material[5].add(point.stress.xz, expected.stress.xz);
                material[6].add(point.equivalent_plastic_strain, expected.equivalent_plastic_strain);
                material[7].add(point.equivalent_creep_strain, expected.equivalent_creep_strain);
            }
        }
    }
    passed = check_relative(branch + "_stress_xx", material[0]) && passed;
    for (std::size_t component = 1; component < 6; ++component) {
        fuelsim::test::print_absolute_metrics(
            branch + "_near_zero_stress_" + std::to_string(component), material[component]);
        passed = check(material[component].maximum_absolute_difference < 1.0,
                     branch + " physically zero transverse or shear stress differs by less than one pascal") &&
                 passed;
    }
    for (std::size_t mechanism = 0; mechanism < 2; ++mechanism) {
        const std::string name = branch + (mechanism == 0 ? "_equivalent_plastic_strain" : "_equivalent_creep_strain");
        if (material[6 + mechanism].maximum_reference > 0.0)
            passed = check_relative(name, material[6 + mechanism]) && passed;
        else
            passed = check(material[6 + mechanism].maximum_absolute_difference < 1.0e-14,
                         name + " remains exactly inactive") &&
                     passed;
    }
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 10) {
        std::cerr << "Usage: fuelsim_b3_hex8_inelastic_moose_tests "
                     "<plastic.fsi> <plastic_nodes.csv> <plastic_elements.csv> "
                     "<creep.fsi> <creep_nodes.csv> <creep_elements.csv> "
                     "<coupled.fsi> <coupled_nodes.csv> <coupled_elements.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim three-dimensional inelastic MOOSE comparison\n");
        bool passed = true;
        for (std::size_t branch = 0; branch < 3; ++branch)
            passed = run(std::array<std::string, 3>{"plastic", "creep", "coupled"}[branch], argv[1 + 3 * branch],
                         argv[2 + 3 * branch], argv[3 + 3 * branch]) &&
                     passed;
        if (!passed) return 1;
        std::cout << "[PASS] three-dimensional inelastic MOOSE comparison\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] three-dimensional inelastic comparison raised: " << error.what() << '\n';
        return 1;
    }
}
