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
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct NodeReference final {
    std::size_t time, node;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::size_t time, element;
    fuelsim::CartesianPoint3 position;
    std::array<double, 3> heat_flux;
    fuelsim::SymmetricTensor3Values stress, strain;
};

struct StepSnapshot final {
    double time;
    std::vector<double> state;
    fuelsim::TransientConservationSummary conservation;
};

class SnapshotObserver final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        _snapshots.push_back({step.time, problem.committed_solution(), step.conservation});
    }

    const std::vector<StepSnapshot>& snapshots() const noexcept { return _snapshots; }

  private:
    std::vector<StepSnapshot> _snapshots;
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
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.8 row in " + path);
    return std::stod(values[index]);
}

std::size_t integer_time(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0 || rounded > 4.0)
        throw std::invalid_argument("Abaqus B5.8 time is not one of 1, 2, 3, or 4 in " + path);
    return static_cast<std::size_t>(rounded);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.8 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.8 nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 10) throw std::invalid_argument("Unexpected Abaqus B5.8 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 2, path);
            if (std::abs(fields[field]) < 1.0e-20) fields[field] = 0.0;
        }
        result.push_back(
            {integer_time(number(values, 0, path), path), static_cast<std::size_t>(number(values, 1, path)), fields});
    }
    if (result.size() != 48) throw std::invalid_argument("Abaqus B5.8 nodal reference must contain 48 rows");
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.8 integration points: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "time_s,element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,hfl_z_w_m2,"
                "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,e11,e22,e33,e12_engineering,e13_engineering,"
                "e23_engineering,ivol_m3")
        throw std::invalid_argument("Unexpected Abaqus B5.8 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 23)
            throw std::invalid_argument("Unexpected Abaqus B5.8 integration-point column count in " + path);
        result.push_back(
            {integer_time(number(values, 0, path), path), static_cast<std::size_t>(number(values, 1, path)),
                {number(values, 3, path), number(values, 4, path), number(values, 5, path)},
                {number(values, 7, path), number(values, 8, path), number(values, 9, path)},
                {number(values, 10, path), number(values, 11, path), number(values, 12, path), number(values, 13, path),
                    number(values, 15, path), number(values, 14, path)},
                {number(values, 16, path), number(values, 17, path), number(values, 18, path),
                    0.5 * number(values, 19, path), 0.5 * number(values, 21, path), 0.5 * number(values, 20, path)}});
    }
    if (result.size() != 64)
        throw std::invalid_argument("Abaqus B5.8 integration-point reference must contain 64 rows");
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

fuelsim::ThermoelasticProperties material_properties() {
    fuelsim::ThermoelasticProperties result =
        fuelsim::test::thermoelastic(0.0, 4.0, 2.0e11, 0.25, 1.2e-5, 300.0, -1.0e8, 1.0e-4, 2.0e-8, 2000.0, 3000.0);
    fuelsim::MaterialFunctionRegistry registry = fuelsim::make_builtin_material_function_registry();
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>(*result.functions);
    functions->thermal = registry.bind_thermal("linear_temperature_thermophysical",
        {{"conductivity", 4.0}, {"density", 2000.0}, {"specific_heat", 3000.0}, {"reference_temperature", 300.0},
            {"conductivity_temperature_coefficient", 0.01}, {"density_temperature_coefficient", -1.0},
            {"specific_heat_temperature_coefficient", 4.0}});
    result.functions = std::move(functions);
    return result;
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition result;
    fuelsim::RegionDefinition region{"solid", "solid", material_properties(), 1.0, 300.0};
    region.heat_source_function = "body_heat";
    result.regions.push_back(std::move(region));
    result.time_tables.emplace_back("right_temperature", std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0},
        std::vector<double>{300.0, 360.0, 420.0, 420.0, 330.0});
    result.time_tables.emplace_back(
        "body_heat", std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0}, std::vector<double>{0.0, 2.0e6, 4.0e6, 0.0, 0.0});
    result.time_tables.emplace_back("surface_heat", std::vector<double>{0.0, 1.0, 2.0, 3.0, 4.0},
        std::vector<double>{0.0, 2.0e3, 5.0e3, 5.0e3, 0.0});
    result.boundary_conditions = {
        {"left_temperature", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::temperature, 300.0},
        {"right_temperature", fuelsim::BoundaryConditionType::dirichlet, "x2", fuelsim::Field::temperature, 1.0, false,
            "right_temperature"},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "x0", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0},
        {"top_heat_flux", fuelsim::BoundaryConditionType::heat_flux, "z1", fuelsim::Field::temperature, 1.0, false,
            "surface_heat"},
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
    options.relative_tolerance = 1.0e-13;
    options.maximum_iterations = 20;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    options.field_residual_scaling = true;
    options.residual_reduction_tolerance = 1.0e-11;
    options.temperature_residual_absolute_tolerance = 1.0e-8;
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

bool metrics_pass(const fuelsim::test::FieldErrorMetrics& metrics, double relative_tolerance, double zero_tolerance) {
    if (metrics.has_relative_norm() && !fuelsim::test::relative_metrics_below(metrics, relative_tolerance))
        return false;
    return metrics.maximum_zero_reference_difference < zero_tolerance;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_b58_hex8_c3d8t_multistep_abaqus_tests <nodes.csv> <integration.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B5.8 Abaqus multi-step C3D8T comparison\n");
        const std::vector<NodeReference> node_reference = read_nodes(argv[1]);
        const std::vector<IntegrationReference> integration_reference = read_integration(argv[2]);
        const fuelsim::UnstructuredHex8Mesh input_mesh = mesh();
        fuelsim::TransientProblem problem(definition(), input_mesh);
        SnapshotObserver observer;
        const fuelsim::TransientResult solve =
            fuelsim::solve_transient(problem, {4.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0, 0.0}, solver_options(), &observer);
        bool passed = check(solve.completed && observer.snapshots().size() == 4 && solve.accepted_steps.size() == 4 &&
                                solve.rejected_steps.empty(),
            "B5.8 accepts the four prescribed one-second increments without cutback");

        std::map<std::size_t, const StepSnapshot*> snapshots;
        for (const StepSnapshot& snapshot : observer.snapshots())
            snapshots.emplace(integer_time(snapshot.time, "Fuelsim B5.8 snapshots"), &snapshot);
        if (snapshots.size() != 4) throw std::invalid_argument("Fuelsim B5.8 snapshots do not cover four unique times");

        std::array<fuelsim::test::FieldErrorMetrics, 8> nodal_metrics;
        fuelsim::TransientProblem reaction_problem(definition(), input_mesh);
        const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(reaction_problem);
        const std::array<fuelsim::Field, 4> fields = {fuelsim::Field::temperature, fuelsim::Field::displacement_x,
            fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};
        for (std::size_t time = 1; time <= 4; ++time) {
            const StepSnapshot& snapshot = *snapshots.at(time);
            reaction_problem.begin_time_step({static_cast<double>(time), 1.0, true});
            const std::vector<double> reaction = raw_residual(reaction_problem, snapshot.state);
            for (const NodeReference& reference : node_reference) {
                if (reference.time != time) continue;
                if (reference.node < 1 || reference.node > 12)
                    throw std::invalid_argument("Abaqus B5.8 node label lies outside the mesh");
                const std::size_t node = reference.node - 1;
                for (std::size_t field = 0; field < fields.size(); ++field) {
                    const std::size_t dof = dofs.dof(fields[field], node);
                    nodal_metrics[field].add(snapshot.state[dof], reference.fields[field]);
                    nodal_metrics[4 + field].add(reaction[dof], reference.fields[4 + field]);
                }
            }
            reaction_problem.commit_time_step(snapshot.state);
            passed = check(snapshot.conservation.relative_thermal_balance < 1.0e-11,
                         "B5.8 Fuelsim thermal conservation closes at time " + std::to_string(time)) &&
                     passed;
        }

        const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y",
            "displacement_z", "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
        const std::array<double, 8> zero_tolerances = {
            1.0e-12, 1.0e-14, 1.0e-14, 1.0e-14, 2.0e-4, 2.0e-2, 2.0e-2, 2.0e-2};
        for (std::size_t field = 0; field < nodal_metrics.size(); ++field) {
            fuelsim::test::print_relative_metrics("b58_" + nodal_names[field], nodal_metrics[field]);
            passed = check(metrics_pass(nodal_metrics[field], 1.0e-3, zero_tolerances[field]),
                         "B5.8 " + nodal_names[field] + " metrics are below the acceptance limits") &&
                     passed;
        }

        std::array<fuelsim::test::FieldErrorMetrics, 15> integration_metrics;
        double maximum_coordinate_error = 0.0;
        constexpr std::array<std::size_t, 8> gauss_to_material_node = {0, 1, 3, 2, 4, 5, 7, 6};
        const fuelsim::IsotropicThermoelasticMaterial constitutive(material_properties());
        for (const IntegrationReference& reference : integration_reference) {
            if (reference.element < 1 || reference.element > 2)
                throw std::invalid_argument("Abaqus B5.8 element label lies outside the mesh");
            const StepSnapshot& snapshot = *snapshots.at(reference.time);
            const std::size_t element = reference.element - 1;
            fuelsim::Hex8Coordinates coordinates{};
            fuelsim::Hex8LocalValues local_state{};
            for (std::size_t local = 0; local < 8; ++local) {
                const std::size_t global = input_mesh.elements()[element].nodes[local];
                coordinates[local] = input_mesh.nodes()[global];
                local_state[local] = snapshot.state[dofs.dof(fuelsim::Field::temperature, global)];
                local_state[8 + local] = snapshot.state[dofs.dof(fuelsim::Field::displacement_x, global)];
                local_state[16 + local] = snapshot.state[dofs.dof(fuelsim::Field::displacement_y, global)];
                local_state[24 + local] = snapshot.state[dofs.dof(fuelsim::Field::displacement_z, global)];
            }
            const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry(coordinates);
            std::size_t closest = 0;
            double closest_squared = std::numeric_limits<double>::max();
            for (std::size_t q = 0; q < 8; ++q) {
                const fuelsim::CartesianPoint3& point = geometry.points[q].position;
                const double distance_squared = std::pow(point.x - reference.position.x, 2) +
                                                std::pow(point.y - reference.position.y, 2) +
                                                std::pow(point.z - reference.position.z, 2);
                if (distance_squared < closest_squared) {
                    closest = q;
                    closest_squared = distance_squared;
                }
            }
            maximum_coordinate_error = std::max(maximum_coordinate_error, std::sqrt(closest_squared));
            const fuelsim::Hex8QuadraturePoint& point = geometry.points[closest];
            const double material_temperature = local_state[gauss_to_material_node[closest]];
            const double conductivity =
                constitutive
                    .conductivity(adlite::Scalar(material_temperature),
                        {static_cast<double>(reference.time), point.position.x, point.position.y, point.position.z})
                    .value();
            std::array<double, 3> actual_heat_flux{};
            for (std::size_t local = 0; local < 8; ++local)
                for (std::size_t component = 0; component < 3; ++component)
                    actual_heat_flux[component] -= conductivity * point.gradient[local][component] * local_state[local];
            fuelsim::CartesianThermoelasticData data{
                constitutive, 0.0, static_cast<double>(reference.time), fuelsim::StrainFormulation::small};
            const fuelsim::SymmetricTensor3Values actual_stress =
                fuelsim::compute_hex8_stress(data, geometry, local_state)[closest];

            fuelsim::Hex8LocalAdValues passive{};
            for (std::size_t dof = 0; dof < passive.size(); ++dof) passive[dof] = local_state[dof];
            double volume = 0.0, average_trace = 0.0;
            for (const fuelsim::Hex8QuadraturePoint& qpoint : geometry.points) {
                const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
                    qpoint, passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::small);
                volume += qpoint.weighted_measure;
                average_trace += qpoint.weighted_measure *
                                 (kinematics.strain_increment.xx.value() + kinematics.strain_increment.yy.value() +
                                     kinematics.strain_increment.zz.value());
            }
            average_trace /= volume;
            const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
                point, passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::small);
            fuelsim::SymmetricTensor3 actual_strain = kinematics.strain_increment;
            const adlite::Scalar correction =
                (average_trace - actual_strain.xx - actual_strain.yy - actual_strain.zz) / 3.0;
            actual_strain.xx += correction;
            actual_strain.yy += correction;
            actual_strain.zz += correction;

            for (std::size_t component = 0; component < 3; ++component)
                integration_metrics[component].add(actual_heat_flux[component], reference.heat_flux[component]);
            const std::array<double, 6> actual_stress_components = {actual_stress.xx, actual_stress.yy,
                actual_stress.zz, actual_stress.xy, actual_stress.yz, actual_stress.xz};
            const std::array<double, 6> reference_stress_components = {reference.stress.xx, reference.stress.yy,
                reference.stress.zz, reference.stress.xy, reference.stress.yz, reference.stress.xz};
            const std::array<double, 6> actual_strain_components = {actual_strain.xx.value(), actual_strain.yy.value(),
                actual_strain.zz.value(), actual_strain.xy.value(), actual_strain.yz.value(), actual_strain.xz.value()};
            const std::array<double, 6> reference_strain_components = {reference.strain.xx, reference.strain.yy,
                reference.strain.zz, reference.strain.xy, reference.strain.yz, reference.strain.xz};
            for (std::size_t component = 0; component < 6; ++component) {
                integration_metrics[3 + component].add(
                    actual_stress_components[component], reference_stress_components[component]);
                integration_metrics[9 + component].add(
                    actual_strain_components[component], reference_strain_components[component]);
            }
        }

        const std::array<std::string, 15> integration_names = {"heat_flux_x", "heat_flux_y", "heat_flux_z", "stress_xx",
            "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz", "strain_xx", "strain_yy", "strain_zz",
            "strain_xy", "strain_yz", "strain_xz"};
        for (std::size_t field = 0; field < integration_metrics.size(); ++field) {
            fuelsim::test::print_relative_metrics("b58_" + integration_names[field], integration_metrics[field]);
            passed = check(metrics_pass(integration_metrics[field], 1.0e-3, 1.0e-10),
                         "B5.8 " + integration_names[field] + " metrics are below 0.1 percent") &&
                     passed;
        }
        passed =
            check(maximum_coordinate_error < 1.0e-12, "B5.8 maps every integration point at all four times uniquely") &&
            passed;
        if (passed && session.rank() == 0)
            std::cout << "[PASS] B5.8 Abaqus C3D8T multi-step temperature-dependent full-field comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.8 Abaqus comparison raised: " << error.what() << '\n';
        return 1;
    }
}
