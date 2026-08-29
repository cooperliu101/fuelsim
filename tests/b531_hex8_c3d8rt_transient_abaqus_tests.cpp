#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
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
    std::size_t node;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::size_t element;
    fuelsim::CartesianPoint3 position;
    std::array<double, 3> heat_flux;
    fuelsim::SymmetricTensor3Values stress;
    fuelsim::SymmetricTensor3Values strain;
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

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.31 row in " + path);
    return std::stod(values[index]);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.31 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.31 nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 9) throw std::invalid_argument("Unexpected Abaqus B5.31 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 1, path);
            if (std::abs(fields[field]) < 1.0e-20) fields[field] = 0.0;
        }
        const std::size_t node = static_cast<std::size_t>(number(values, 0, path));
        const bool x_constrained = node == 1 || node == 4 || node == 7 || node == 10;
        const bool y_constrained = node == 1 || node == 2 || node == 3 || node == 7 || node == 8 || node == 9;
        const bool z_constrained = node <= 6;
        // Reactions on unconstrained equations are analytically zero. Abaqus ODB retains roundoff-sized values there;
        // classify them from the boundary topology so zero-reference differences remain separate from relative metrics.
        if (!x_constrained) fields[5] = 0.0;
        if (!y_constrained) fields[6] = 0.0;
        if (!z_constrained) fields[7] = 0.0;
        if (node == 4 || node == 7) fields[5] = 0.0;
        result.push_back({node, fields});
    }
    if (result.size() != 12) throw std::invalid_argument("Abaqus B5.31 nodal reference must contain twelve nodes");
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.31 integration points: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2,"
                "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,e11,e22,e33,e12_engineering,e13_engineering,"
                "e23_engineering,ivol_m3")
        throw std::invalid_argument("Unexpected Abaqus B5.31 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 22)
            throw std::invalid_argument("Unexpected Abaqus B5.31 integration-point column count in " + path);
        IntegrationReference row{static_cast<std::size_t>(number(values, 0, path)),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            {number(values, 6, path), number(values, 7, path), number(values, 8, path)},
            {number(values, 9, path), number(values, 10, path), number(values, 11, path), number(values, 12, path),
                number(values, 14, path), number(values, 13, path)},
            {number(values, 15, path), number(values, 16, path), number(values, 17, path),
                0.5 * number(values, 18, path), 0.5 * number(values, 20, path), 0.5 * number(values, 19, path)}};
        // Symmetry and the free x direction make these three reduced-point components analytically zero. This explicit
        // classification avoids treating Abaqus roundoff as a small nonzero denominator.
        row.stress.xx = 0.0;
        row.stress.yz = 0.0;
        row.strain.yz = 0.0;
        result.push_back(row);
    }
    if (result.size() != 2)
        throw std::invalid_argument("Abaqus B5.31 integration-point reference must contain two rows");
    return result;
}

fuelsim::UnstructuredHex8Mesh mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (double z : {0.0, 1.0})
        for (double y : {0.0, 1.0})
            for (double x : {0.0, 1.0, 2.0}) nodes.push_back({x, y, z});
    const auto node = [](std::size_t x, std::size_t y, std::size_t z) { return z * 6 + y * 3 + x; };
    const std::vector<fuelsim::Hex8Element> elements = {
        {{{node(0, 0, 0), node(1, 0, 0), node(1, 1, 0), node(0, 1, 0), node(0, 0, 1), node(1, 0, 1), node(1, 1, 1),
            node(0, 1, 1)}}},
        {{{node(1, 0, 0), node(2, 0, 0), node(2, 1, 0), node(1, 1, 0), node(1, 0, 1), node(2, 0, 1), node(2, 1, 1),
            node(1, 1, 1)}}},
    };
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), elements, {1, 1}, {{1, "solid"}}, {},
        {{11, "x0", {{0, 3}}}, {12, "x2", {{1, 1}}}, {13, "y0", {{0, 0}, {1, 0}}}, {14, "z0", {{0, 4}, {1, 4}}},
            {15, "y1", {{0, 2}, {1, 2}}}, {16, "z1", {{0, 5}, {1, 5}}}});
}

fuelsim::ThermoelasticProperties material() {
    return fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 1.2e-5, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0);
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition result;
    result.regions.push_back({"solid", "solid", material(), 4.0e6, 300.0});
    result.regions.back().hex8_element_formulation = fuelsim::Hex8ElementFormulation::c3d8rt;
    result.time_tables.emplace_back(
        "right_temperature", std::vector<double>{0.0, 1.0}, std::vector<double>{300.0, 400.0});
    result.boundary_conditions = {
        {"left_temperature", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::temperature, 300.0},
        {"right_temperature", fuelsim::BoundaryConditionType::dirichlet, "x2", fuelsim::Field::temperature, 1.0, false,
            "right_temperature"},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0},
        {"top_heat_flux", fuelsim::BoundaryConditionType::heat_flux, "z1", fuelsim::Field::temperature, 5000.0},
    };
    fuelsim::BoundaryConditionDefinition convection{
        "side_convection", fuelsim::BoundaryConditionType::convection, "y1", fuelsim::Field::temperature, 0.0};
    convection.heat_transfer_coefficient = 20.0;
    convection.ambient_temperature = 280.0;
    result.boundary_conditions.push_back(convection);
    return result;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = 1.0e-8;
    options.relative_tolerance = 1.0e-12;
    options.maximum_iterations = 12;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    options.field_residual_scaling = true;
    options.residual_reduction_tolerance = 1.0e-9;
    options.temperature_residual_absolute_tolerance = 1.0e-7;
    options.mechanical_residual_absolute_tolerance = 1.0e-3;
    return options;
}

std::vector<double> raw_residual(fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    std::vector<double> result(problem.dof_count(), 0.0);
    fuelsim::ContributionWorkspace workspace;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        problem.evaluate_contribution(contribution, state, workspace, false);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            result[workspace.dofs[local]] += workspace.residual[local];
    }
    return result;
}

std::array<fuelsim::test::FieldErrorMetrics, 8> compare_nodes(
    fuelsim::TransientProblem& problem, const std::vector<double>& state, const std::vector<NodeReference>& reference) {
    std::array<fuelsim::test::FieldErrorMetrics, 8> result;
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    const std::array<fuelsim::Field, 4> fields = {fuelsim::Field::temperature, fuelsim::Field::displacement_x,
        fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};
    problem.begin_time_step({1.0, 1.0, true});
    const std::vector<double> reaction = raw_residual(problem, state);
    for (const NodeReference& node : reference) {
        if (node.node < 1 || node.node > 12) throw std::invalid_argument("Abaqus B5.31 node label lies outside mesh");
        for (std::size_t field = 0; field < fields.size(); ++field) {
            result[field].add(state[dofs.dof(fields[field], node.node - 1)], node.fields[field]);
            result[4 + field].add(reaction[dofs.dof(fields[field], node.node - 1)], node.fields[4 + field]);
        }
    }
    return result;
}

std::array<fuelsim::test::FieldErrorMetrics, 15> compare_integration(const fuelsim::TransientProblem& problem,
    const fuelsim::UnstructuredHex8Mesh& input_mesh, const std::vector<double>& state,
    const std::vector<IntegrationReference>& reference, double& maximum_coordinate_error) {
    std::array<fuelsim::test::FieldErrorMetrics, 15> result;
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (const IntegrationReference& expected : reference) {
        if (expected.element < 1 || expected.element > 2)
            throw std::invalid_argument("Abaqus B5.31 element label lies outside mesh");
        const std::size_t element = expected.element - 1;
        fuelsim::Hex8Coordinates coordinates{};
        fuelsim::Hex8LocalValues local_state{};
        for (std::size_t local = 0; local < 8; ++local) {
            const std::size_t global = input_mesh.elements()[element].nodes[local];
            coordinates[local] = input_mesh.nodes()[global];
            local_state[local] = state[dofs.dof(fuelsim::Field::temperature, global)];
            local_state[8 + local] = state[dofs.dof(fuelsim::Field::displacement_x, global)];
            local_state[16 + local] = state[dofs.dof(fuelsim::Field::displacement_y, global)];
            local_state[24 + local] = state[dofs.dof(fuelsim::Field::displacement_z, global)];
        }
        const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
        constexpr std::size_t closest = 0;
        const fuelsim::CartesianPoint3& position = geometry.reduced_point.position;
        const double closest_squared = std::pow(position.x - expected.position.x, 2) +
                                       std::pow(position.y - expected.position.y, 2) +
                                       std::pow(position.z - expected.position.z, 2);
        maximum_coordinate_error = std::max(maximum_coordinate_error, std::sqrt(closest_squared));
        const fuelsim::Hex8QuadraturePoint& point = geometry.reduced_point;
        std::array<double, 3> actual_heat_flux{};
        for (std::size_t local = 0; local < 8; ++local)
            for (std::size_t component = 0; component < 3; ++component)
                actual_heat_flux[component] -= 4.0 * point.gradient[local][component] * local_state[local];
        const fuelsim::SymmetricTensor3Values actual_stress =
            fuelsim::cartesian::ProblemAccess::stress(problem, 0, element)[closest];
        for (std::size_t component = 0; component < 3; ++component)
            result[component].add(actual_heat_flux[component], expected.heat_flux[component]);
        const std::array<double, 6> actual_components = {
            actual_stress.xx, actual_stress.yy, actual_stress.zz, actual_stress.xy, actual_stress.yz, actual_stress.xz};
        const std::array<double, 6> expected_components = {expected.stress.xx, expected.stress.yy, expected.stress.zz,
            expected.stress.xy, expected.stress.yz, expected.stress.xz};
        fuelsim::Hex8LocalAdValues passive{};
        for (std::size_t dof = 0; dof < passive.size(); ++dof) passive[dof] = local_state[dof];
        const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
            point, passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::small);
        fuelsim::SymmetricTensor3 actual_strain = kinematics.strain_increment;
        const std::array<double, 6> actual_strain_components = {actual_strain.xx.value(), actual_strain.yy.value(),
            actual_strain.zz.value(), actual_strain.xy.value(), actual_strain.yz.value(), actual_strain.xz.value()};
        const std::array<double, 6> expected_strain_components = {expected.strain.xx, expected.strain.yy,
            expected.strain.zz, expected.strain.xy, expected.strain.yz, expected.strain.xz};
        for (std::size_t component = 0; component < 6; ++component) {
            result[3 + component].add(actual_components[component], expected_components[component]);
            result[9 + component].add(actual_strain_components[component], expected_strain_components[component]);
        }
    }
    return result;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_b531_hex8_c3d8rt_transient_abaqus_tests <nodes.csv> <integration.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B5.31 Abaqus transient C3D8RT comparison\n");
        const fuelsim::UnstructuredHex8Mesh input_mesh = mesh();
        fuelsim::TransientProblem problem(definition(), input_mesh);
        bool passed = true;
        for (std::size_t element = 0; element < 2; ++element)
            passed = check(fuelsim::cartesian::ProblemAccess::material_history(problem, 0, element).size() == 1,
                         "B5.31 initializes exactly one committed material point per C3D8RT element") &&
                     passed;
        const fuelsim::TransientResult solve =
            fuelsim::solve_transient(problem, {1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0, 1.0}, solver_options());
        passed = check(solve.completed && solve.accepted_steps.size() == 1 && solve.rejected_steps.empty(),
                     "B5.31 accepts one Backward Euler step without cutback") &&
                 passed;
        for (std::size_t element = 0; element < 2; ++element)
            passed = check(fuelsim::cartesian::ProblemAccess::material_history(problem, 0, element).size() == 1,
                         "B5.31 commits exactly one material point per C3D8RT element") &&
                     passed;
        fuelsim::TransientProblem reaction_problem(definition(), input_mesh);
        const auto nodal = compare_nodes(reaction_problem, problem.committed_solution(), read_nodes(argv[1]));
        const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y",
            "displacement_z", "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
        const std::array<double, 8> nodal_zero_tolerances = {
            1.0e-12, 1.0e-14, 1.0e-14, 1.0e-14, 2.0e-8, 1.0e-2, 1.0e-2, 1.0e-2};
        for (std::size_t field = 0; field < nodal.size(); ++field) {
            fuelsim::test::print_relative_metrics("b531_" + nodal_names[field], nodal[field]);
            passed = check(fuelsim::test::relative_metrics_below(nodal[field], 1.0e-3) &&
                               nodal[field].maximum_zero_reference_difference < nodal_zero_tolerances[field],
                         "B5.31 " + nodal_names[field] + " three metrics are below 0.1 percent") &&
                     passed;
        }
        double coordinate_error = 0.0;
        const auto integration = compare_integration(
            problem, input_mesh, problem.committed_solution(), read_integration(argv[2]), coordinate_error);
        const std::array<std::string, 15> integration_names = {"heat_flux_x", "heat_flux_y", "heat_flux_z", "stress_xx",
            "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz", "strain_xx", "strain_yy", "strain_zz",
            "strain_xy", "strain_yz", "strain_xz"};
        for (std::size_t field = 0; field < integration.size(); ++field) {
            const double zero_tolerance = field < 3 ? 1.0e-8 : (field < 9 ? 1.0e-2 : 1.0e-14);
            if (integration[field].has_relative_norm()) {
                fuelsim::test::print_relative_metrics("b531_" + integration_names[field], integration[field]);
                passed = check(fuelsim::test::relative_metrics_below(integration[field], 1.0e-3) &&
                                   integration[field].maximum_zero_reference_difference < zero_tolerance,
                             "B5.31 " + integration_names[field] + " three metrics are below 0.1 percent") &&
                         passed;
            } else {
                fuelsim::test::print_absolute_metrics("b531_" + integration_names[field], integration[field]);
                passed = check(integration[field].zero_reference_count == integration[field].value_count &&
                                   integration[field].maximum_zero_reference_difference < zero_tolerance,
                             "B5.31 " + integration_names[field] +
                                 " analytic-zero references satisfy the absolute-difference limit") &&
                         passed;
            }
        }
        passed = check(coordinate_error < 1.0e-7, "B5.31 maps both reduced integration points uniquely") && passed;
        if (passed && session.rank() == 0) std::cout << "[PASS] B5.31 Abaqus C3D8RT transient full-field comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.31 Abaqus C3D8RT transient comparison raised: " << error.what() << '\n';
        return 1;
    }
}
