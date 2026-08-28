#include "fuelsim/solver/solve_workflows.hpp"
#include "support/abaqus_hex8_full_field.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

struct MeshDivisions final {
    std::size_t fuel_radial, clad_radial, angular, axial;
};

struct RegionGrid final {
    std::size_t offset, radial_nodes, angular_nodes, axial_nodes;
};

std::size_t node(const RegionGrid& grid, std::size_t radial, std::size_t angular, std::size_t axial) {
    return grid.offset + axial * grid.angular_nodes * grid.radial_nodes + angular * grid.radial_nodes + radial;
}

void append_region_nodes(std::vector<fuelsim::CartesianPoint3>& nodes, double inner_radius, double outer_radius,
    double height, std::size_t radial_elements, std::size_t angular_elements, std::size_t axial_elements) {
    for (std::size_t axial = 0; axial <= axial_elements; ++axial)
        for (std::size_t angular = 0; angular <= angular_elements; ++angular)
            for (std::size_t radial = 0; radial <= radial_elements; ++radial) {
                const double radius = inner_radius + (outer_radius - inner_radius) * static_cast<double>(radial) /
                                                         static_cast<double>(radial_elements),
                             angle = 0.5 * pi * static_cast<double>(angular) / static_cast<double>(angular_elements),
                             z = height * static_cast<double>(axial) / static_cast<double>(axial_elements);
                nodes.push_back({radius * std::cos(angle), radius * std::sin(angle), z});
            }
}

fuelsim::UnstructuredHex8Mesh engineering_mesh(const MeshDivisions& divisions) {
    constexpr double fuel_inner_radius = 1.0e-3, fuel_outer_radius = 4.0e-3, clad_inner_radius = 4.005e-3,
                     clad_outer_radius = 4.7e-3, height = 4.0e-2;
    std::vector<fuelsim::CartesianPoint3> nodes;
    append_region_nodes(
        nodes, fuel_inner_radius, fuel_outer_radius, height, divisions.fuel_radial, divisions.angular, divisions.axial);
    const std::size_t fuel_node_count = nodes.size();
    append_region_nodes(
        nodes, clad_inner_radius, clad_outer_radius, height, divisions.clad_radial, divisions.angular, divisions.axial);
    const RegionGrid fuel{0, divisions.fuel_radial + 1, divisions.angular + 1, divisions.axial + 1},
        clad{fuel_node_count, divisions.clad_radial + 1, divisions.angular + 1, divisions.axial + 1};
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> blocks;
    std::array<std::vector<fuelsim::ElementSide>, 10> sides;
    const auto append_elements = [&](const RegionGrid& grid, std::size_t radial_elements, std::int64_t block) {
        for (std::size_t axial = 0; axial < divisions.axial; ++axial)
            for (std::size_t angular = 0; angular < divisions.angular; ++angular)
                for (std::size_t radial = 0; radial < radial_elements; ++radial) {
                    const std::size_t index = elements.size();
                    elements.push_back({{{node(grid, radial, angular, axial), node(grid, radial + 1, angular, axial),
                        node(grid, radial + 1, angular + 1, axial), node(grid, radial, angular + 1, axial),
                        node(grid, radial, angular, axial + 1), node(grid, radial + 1, angular, axial + 1),
                        node(grid, radial + 1, angular + 1, axial + 1), node(grid, radial, angular + 1, axial + 1)}}});
                    blocks.push_back(block);
                    if (block == 1 && radial == 0) sides[0].push_back({index, 3});
                    if (block == 1 && radial + 1 == radial_elements) sides[1].push_back({index, 1});
                    if (block == 2 && radial == 0) sides[2].push_back({index, 3});
                    if (block == 2 && radial + 1 == radial_elements) sides[3].push_back({index, 1});
                    const std::size_t symmetry_offset = block == 1 ? 4 : 7;
                    if (angular == 0) sides[symmetry_offset].push_back({index, 0});
                    if (angular + 1 == divisions.angular) sides[symmetry_offset + 1].push_back({index, 2});
                    if (axial == 0) sides[symmetry_offset + 2].push_back({index, 4});
                }
    };
    append_elements(fuel, divisions.fuel_radial, 1);
    append_elements(clad, divisions.clad_radial, 2);
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), std::move(elements), std::move(blocks),
        {{1, "fuel"}, {2, "clad"}}, {},
        {{10, "fuel_inner", sides[0]}, {11, "fuel_outer", sides[1]}, {12, "clad_inner", sides[2]},
            {13, "clad_outer", sides[3]}, {14, "fuel_symmetry_y", sides[4]}, {15, "fuel_symmetry_x", sides[5]},
            {16, "fuel_bottom", sides[6]}, {17, "clad_symmetry_y", sides[7]}, {18, "clad_symmetry_x", sides[8]},
            {19, "clad_bottom", sides[9]}});
}

fuelsim::ThermoelasticProperties material(
    double conductivity, double modulus, double poisson, double expansion, double density, double heat_capacity) {
    return fuelsim::test::thermoelastic(
        0.0, conductivity, modulus, poisson, expansion, 600.0, 0.0, 0.0, 0.0, density, heat_capacity);
}

fuelsim::SpatialDefinition definition(double penalty) {
    fuelsim::SpatialDefinition result;
    result.regions = {{"fuel", "fuel", material(3.0, 2.0e11, 0.30, 1.0e-5, 1.0e4, 300.0), 0.0, 600.0, -1, "",
                          fuelsim::StrainFormulation::finite},
        {"clad", "clad", material(15.0, 1.0e11, 0.32, 5.0e-6, 6.5e3, 330.0), 0.0, 600.0, -1, "",
            fuelsim::StrainFormulation::finite}};
    const std::vector<double> times = {0.0, 10000.0};
    result.time_tables.emplace_back("fuel_temperature", times, std::vector<double>{600.0, 1200.0});
    result.boundary_conditions = {
        {"fuel_temperature", fuelsim::BoundaryConditionType::dirichlet, "fuel_inner", fuelsim::Field::temperature, 1.0,
            false, "fuel_temperature"},
        {"coolant_temperature", fuelsim::BoundaryConditionType::dirichlet, "clad_outer", fuelsim::Field::temperature,
            600.0},
        {"fuel_symmetry_y", fuelsim::BoundaryConditionType::dirichlet, "fuel_symmetry_y",
            fuelsim::Field::displacement_y, 0.0},
        {"fuel_symmetry_x", fuelsim::BoundaryConditionType::dirichlet, "fuel_symmetry_x",
            fuelsim::Field::displacement_x, 0.0},
        {"fuel_bottom", fuelsim::BoundaryConditionType::dirichlet, "fuel_bottom", fuelsim::Field::displacement_z, 0.0},
        {"clad_symmetry_y", fuelsim::BoundaryConditionType::dirichlet, "clad_symmetry_y",
            fuelsim::Field::displacement_y, 0.0},
        {"clad_symmetry_x", fuelsim::BoundaryConditionType::dirichlet, "clad_symmetry_x",
            fuelsim::Field::displacement_x, 0.0},
        {"clad_bottom", fuelsim::BoundaryConditionType::dirichlet, "clad_bottom", fuelsim::Field::displacement_z, 0.0},
    };
    fuelsim::ContactDefinition contact;
    contact.name = "fuel_clad_contact";
    contact.primary = "fuel_outer";
    contact.secondary = "clad_inner";
    contact.thermal = true;
    contact.mechanical = true;
    contact.penalty = penalty;
    contact.friction_coefficient = 0.1;
    contact.friction_slip_tolerance = 0.005;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    contact.gap_heat_conductance_law = fuelsim::GapHeatConductanceLaw::affine;
    contact.gap_conductance = 100.0;
    contact.gap_conductance_pressure_derivative = 1.0e-5;
    result.contacts.push_back(contact);
    return result;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = 1.0e-8;
    result.relative_tolerance = 1.0e-11;
    result.maximum_iterations = 60;
    result.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    result.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    result.field_residual_scaling = true;
    result.temperature_residual_absolute_tolerance = 1.0e-5;
    result.mechanical_residual_absolute_tolerance = 1.0e-3;
    return result;
}

struct Response final {
    bool completed = false, full_field_passed = false;
    std::size_t accepted_steps = 0, rejected_steps = 0, active_contact_nodes = 0;
    int nonlinear_iterations = 0;
    double contact_force = 0.0, contact_heat_rate = 0.0, maximum_pressure = 0.0, maximum_penetration = 0.0,
           fuel_average_temperature = 0.0, clad_average_temperature = 0.0, maximum_displacement = 0.0,
           maximum_equivalent_stress = 0.0, boundary_heat_rate = 0.0, stored_heat_rate = 0.0;
    std::string failure;
};

Response solve_case(const MeshDivisions& divisions, double time_step, double penalty, const std::string& case_name,
    const std::string& reference_directory) {
    const fuelsim::SpatialDefinition case_definition = definition(penalty);
    const fuelsim::UnstructuredHex8Mesh case_mesh = engineering_mesh(divisions);
    fuelsim::TransientProblem problem(case_definition, case_mesh);
    fuelsim::test::AbaqusHex8SnapshotObserver observer;
    const fuelsim::TransientResult solve = fuelsim::solve_transient(
        problem, {10000.0, time_step, time_step, time_step, 1.0, 0.5, 0, 0.0}, solver_options(), &observer);
    Response response;
    response.completed = solve.completed;
    response.accepted_steps = solve.accepted_steps.size();
    response.rejected_steps = solve.rejected_steps.size();
    response.nonlinear_iterations = solve.total_nonlinear_iterations;
    if (!solve.rejected_steps.empty()) response.failure = solve.rejected_steps.back().failure_message;
    if (!solve.completed) return response;
    fuelsim::test::AbaqusHex8FullFieldOptions comparison;
    comparison.case_name = case_name;
    comparison.reference_prefix = reference_directory + "/" + case_name;
    comparison.expected_steps = static_cast<std::size_t>(std::llround(10000.0 / time_step));
    comparison.time_step = time_step;
    comparison.bulk_relative_tolerance = 1.5e-2;
    comparison.displacement_pointwise_relative_tolerance = 5.0e-2;
    comparison.reaction_pointwise_relative_tolerance = 5.0e-2;
    comparison.reaction_pointwise_absolute_tolerance = 1.0e-6;
    comparison.stress_pointwise_relative_tolerance = 2.0e-2;
    comparison.stress_pointwise_absolute_tolerance = 10.0;
    comparison.logarithmic_strain_pointwise_relative_tolerance = 2.0e-2;
    comparison.elastic_strain_pointwise_relative_tolerance = 2.0e-2;
    comparison.elastic_strain_pointwise_absolute_tolerance = 1.0e-10;
    comparison.contact_pointwise_relative_tolerance = 2.5e-2;
    comparison.energy_pointwise_relative_tolerance = 5.0e-2;
    comparison.minimum_contact_state_match_fraction = 0.95;
    comparison.use_contact_summary_total_slip = true;
    response.full_field_passed = fuelsim::test::compare_abaqus_hex8_full_field(
        problem, case_definition, case_mesh, observer.snapshots(), comparison);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.committed_state);
    response.active_contact_nodes = interface.active_contact_nodes;
    response.contact_force = interface.total_contact_force;
    response.contact_heat_rate = interface.total_heat_rate;
    response.maximum_pressure = interface.maximum_contact_pressure;
    response.maximum_penetration = std::max(0.0, -interface.minimum_contact_gap);
    response.boundary_heat_rate = solve.accepted_steps.back().conservation.dirichlet_heat_input_rate;
    response.stored_heat_rate = solve.accepted_steps.back().conservation.stored_heat_rate;
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        double temperature_sum = 0.0;
        const std::size_t region_nodes = spatial.region_mesh(region).nodes().size();
        for (std::size_t local = 0; local < region_nodes; ++local) {
            const std::size_t global = spatial.global_node(region, local);
            temperature_sum += solve.committed_state[spatial.dof(fuelsim::Field::temperature, global)];
            const double ux = solve.committed_state[spatial.dof(fuelsim::Field::displacement_x, global)],
                         uy = solve.committed_state[spatial.dof(fuelsim::Field::displacement_y, global)],
                         uz = solve.committed_state[spatial.dof(fuelsim::Field::displacement_z, global)];
            response.maximum_displacement = std::max(response.maximum_displacement, std::hypot(ux, uy, uz));
        }
        (region == 0 ? response.fuel_average_temperature : response.clad_average_temperature) =
            temperature_sum / static_cast<double>(region_nodes);
        for (std::size_t element = 0; element < spatial.region_mesh(region).elements().size(); ++element)
            for (const fuelsim::CartesianMaterialPointState& point :
                fuelsim::cartesian::ProblemAccess::material_history(problem, region, element)) {
                const double mean = (point.stress.xx + point.stress.yy + point.stress.zz) / 3.0,
                             x = point.stress.xx - mean, y = point.stress.yy - mean, z = point.stress.zz - mean;
                response.maximum_equivalent_stress = std::max(response.maximum_equivalent_stress,
                    std::sqrt(1.5 * (x * x + y * y + z * z +
                                        2.0 * (point.stress.xy * point.stress.xy + point.stress.yz * point.stress.yz +
                                                  point.stress.xz * point.stress.xz))));
            }
    }
    return response;
}

void print_response(const std::string& name, const Response& response) {
    std::cout << "b527_" << name << "_completed=" << response.completed << '\n'
              << "b527_" << name << "_accepted_steps=" << response.accepted_steps << '\n'
              << "b527_" << name << "_rejected_steps=" << response.rejected_steps << '\n'
              << "b527_" << name << "_full_field_passed=" << response.full_field_passed << '\n'
              << "b527_" << name << "_nonlinear_iterations=" << response.nonlinear_iterations << '\n'
              << "b527_" << name << "_active_contact_nodes=" << response.active_contact_nodes << '\n'
              << "b527_" << name << "_contact_force=" << response.contact_force << '\n'
              << "b527_" << name << "_contact_heat_rate=" << response.contact_heat_rate << '\n'
              << "b527_" << name << "_maximum_pressure=" << response.maximum_pressure << '\n'
              << "b527_" << name << "_maximum_penetration=" << response.maximum_penetration << '\n'
              << "b527_" << name << "_fuel_average_temperature=" << response.fuel_average_temperature << '\n'
              << "b527_" << name << "_clad_average_temperature=" << response.clad_average_temperature << '\n'
              << "b527_" << name << "_maximum_displacement=" << response.maximum_displacement << '\n'
              << "b527_" << name << "_maximum_equivalent_stress=" << response.maximum_equivalent_stress << '\n'
              << "b527_" << name << "_boundary_heat_rate=" << response.boundary_heat_rate << '\n'
              << "b527_" << name << "_stored_heat_rate=" << response.stored_heat_rate << '\n';
    if (!response.failure.empty()) std::cout << "b527_" << name << "_failure=" << response.failure << '\n';
}

double relative_change(double left, double right) {
    return std::abs(left - right) / std::max(std::abs(right), std::numeric_limits<double>::min());
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: fuelsim_b527_hex8_engineering_fuel_clad_tests <Abaqus reference directory> [case]\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B5.27 engineering fuel-clad verification\n");
        const std::string reference_directory = argv[1];
        const std::string selected = argc > 2 ? argv[2] : "";
        const std::array<std::string, 6> known_cases = {
            "coarse", "medium", "fine", "half_step", "low_penalty", "high_penalty"};
        if (!selected.empty() && std::find(known_cases.begin(), known_cases.end(), selected) == known_cases.end())
            throw std::invalid_argument("Unknown B5.27 case: " + selected);
        const auto enabled = [&](const std::string& name) { return selected.empty() || selected == name; };
        const Response coarse = enabled("coarse")
                                    ? solve_case({1, 1, 3, 2}, 1000.0, 1.0e14, "b527_coarse", reference_directory)
                                    : Response{};
        const Response medium = enabled("medium")
                                    ? solve_case({2, 1, 4, 3}, 1000.0, 1.0e14, "b527_medium", reference_directory)
                                    : Response{};
        const Response fine =
            enabled("fine") ? solve_case({2, 2, 6, 4}, 1000.0, 1.0e14, "b527_fine", reference_directory) : Response{};
        const Response half_step = enabled("half_step")
                                       ? solve_case({2, 1, 4, 3}, 500.0, 1.0e14, "b527_half_step", reference_directory)
                                       : Response{};
        const Response low_penalty =
            enabled("low_penalty") ? solve_case({2, 1, 4, 3}, 1000.0, 5.0e13, "b527_low_penalty", reference_directory)
                                   : Response{};
        const Response high_penalty =
            enabled("high_penalty") ? solve_case({2, 1, 4, 3}, 1000.0, 2.0e14, "b527_high_penalty", reference_directory)
                                    : Response{};
        std::vector<std::pair<std::string, const Response*>> cases;
        if (enabled("coarse")) cases.emplace_back("coarse", &coarse);
        if (enabled("medium")) cases.emplace_back("medium", &medium);
        if (enabled("fine")) cases.emplace_back("fine", &fine);
        if (enabled("half_step")) cases.emplace_back("half_step", &half_step);
        if (enabled("low_penalty")) cases.emplace_back("low_penalty", &low_penalty);
        if (enabled("high_penalty")) cases.emplace_back("high_penalty", &high_penalty);
        bool completed = true;
        for (const auto& item : cases) {
            print_response(item.first, *item.second);
            completed = completed && item.second->completed && item.second->full_field_passed &&
                        item.second->active_contact_nodes > 0 && item.second->contact_force > 0.0 &&
                        std::abs(item.second->contact_heat_rate) > 0.0;
        }
        bool passed = check(completed, "B5.27 every engineering fuel-clad mesh, time-step, and penalty case completes "
                                       "with active mechanical and thermal contact");
        if (completed && selected.empty()) {
            const double time_temperature_change =
                             relative_change(half_step.clad_average_temperature, medium.clad_average_temperature),
                         time_force_change = relative_change(half_step.contact_force, medium.contact_force),
                         time_stress_change =
                             relative_change(half_step.maximum_equivalent_stress, medium.maximum_equivalent_stress),
                         medium_fine_temperature_change =
                             relative_change(fine.clad_average_temperature, medium.clad_average_temperature),
                         medium_fine_force_change = relative_change(fine.contact_force, medium.contact_force),
                         medium_fine_stress_change =
                             relative_change(fine.maximum_equivalent_stress, medium.maximum_equivalent_stress),
                         penalty_force_change = relative_change(high_penalty.contact_force, low_penalty.contact_force);
            std::cout << "b527_half_step_clad_temperature_relative_change=" << time_temperature_change << '\n'
                      << "b527_half_step_contact_force_relative_change=" << time_force_change << '\n'
                      << "b527_half_step_stress_relative_change=" << time_stress_change << '\n'
                      << "b527_medium_fine_clad_temperature_relative_change=" << medium_fine_temperature_change << '\n'
                      << "b527_medium_fine_contact_force_relative_change=" << medium_fine_force_change << '\n'
                      << "b527_medium_fine_stress_relative_change=" << medium_fine_stress_change << '\n'
                      << "b527_penalty_contact_force_relative_change=" << penalty_force_change << '\n';
            passed =
                check(time_temperature_change < 5.0e-2 && time_force_change < 5.0e-2 && time_stress_change < 5.0e-2,
                    "B5.27 time-step halving changes clad temperature, contact force, and stress by less than five "
                    "percent") &&
                passed;
            passed = check(medium_fine_temperature_change < 1.5e-1 && medium_fine_force_change < 1.5e-1 &&
                               medium_fine_stress_change < 1.5e-1,
                         "B5.27 medium-to-fine engineering mesh changes clad temperature, contact force, and stress by "
                         "less than fifteen percent") &&
                     passed;
            passed =
                check(low_penalty.maximum_penetration > medium.maximum_penetration &&
                          medium.maximum_penetration > high_penalty.maximum_penetration &&
                          penalty_force_change < 5.0e-2,
                    "B5.27 penalty refinement reduces penetration while preserving the engineering contact force") &&
                passed;
        }
        if (passed) std::cout << "[PASS] B5.27 engineering-scale three-dimensional fuel-clad verification\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.27 engineering fuel-clad verification raised: " << error.what() << '\n';
        return 1;
    }
}
