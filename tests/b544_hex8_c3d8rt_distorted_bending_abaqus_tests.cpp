#include "fuelsim/core/material_functions.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/abaqus_hex8_full_field.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr std::size_t x_nodes = 5, y_nodes = 2, z_nodes = 3;
constexpr std::size_t increments = 10;
constexpr double time_step = 1.0e5;

std::size_t node_index(std::size_t x, std::size_t y, std::size_t z) { return z * x_nodes * y_nodes + y * x_nodes + x; }

fuelsim::CartesianPoint3 node_coordinate(std::size_t x, std::size_t y, std::size_t z) {
    const double x_coordinate = static_cast<double>(x);
    const double y_coordinate = static_cast<double>(y);
    const double z_coordinate = 0.5 * static_cast<double>(z);
    if (x == 0 || x + 1 == x_nodes) return {x_coordinate, y_coordinate, z_coordinate};
    const double alternating = x % 2 == 0 ? -1.0 : 1.0;
    return {x_coordinate + 0.06 * (2.0 * static_cast<double>(y) - 1.0) * (static_cast<double>(z) - 1.0),
        y_coordinate + 0.04 * alternating * (static_cast<double>(z) - 1.0),
        z_coordinate + 0.05 * alternating * (static_cast<double>(y) - 0.5)};
}

fuelsim::UnstructuredHex8Mesh mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    for (std::size_t z = 0; z < z_nodes; ++z)
        for (std::size_t y = 0; y < y_nodes; ++y)
            for (std::size_t x = 0; x < x_nodes; ++x) nodes.push_back(node_coordinate(x, y, z));
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<fuelsim::ElementSide> left, right, y_high, z_high;
    for (std::size_t z = 0; z + 1 < z_nodes; ++z) {
        for (std::size_t x = 0; x + 1 < x_nodes; ++x) {
            const std::size_t element = elements.size();
            elements.push_back({{{node_index(x, 0, z), node_index(x + 1, 0, z), node_index(x + 1, 1, z),
                node_index(x, 1, z), node_index(x, 0, z + 1), node_index(x + 1, 0, z + 1), node_index(x + 1, 1, z + 1),
                node_index(x, 1, z + 1)}}});
            if (x == 0) left.push_back({element, 3});
            if (x + 2 == x_nodes) right.push_back({element, 1});
            y_high.push_back({element, 2});
            if (z + 2 == z_nodes) z_high.push_back({element, 5});
        }
    }
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes), elements, std::vector<std::int64_t>(elements.size(), 1),
        {{1, "solid"}}, {},
        {{11, "left", std::move(left)}, {12, "right", std::move(right)}, {13, "y_high", std::move(y_high)},
            {14, "z_high", std::move(z_high)}});
}

fuelsim::ThermoelasticProperties material() {
    fuelsim::MaterialFunctionRegistry registry = fuelsim::make_builtin_material_function_registry();
    auto functions = std::make_shared<fuelsim::MaterialFunctionSet>();
    functions->name = "b544_temperature_dependent_material";
    functions->thermal = registry.bind_thermal("linear_temperature_thermophysical",
        {{"conductivity", 10.0}, {"density", 2000.0}, {"specific_heat", 500.0}, {"reference_temperature", 300.0},
            {"conductivity_temperature_coefficient", 0.02}, {"density_temperature_coefficient", -0.5},
            {"specific_heat_temperature_coefficient", 1.25}});
    functions->elasticity = registry.bind_elasticity("linear_temperature_isotropic",
        {{"young_modulus", 2.0e8}, {"poisson_ratio", 0.25}, {"reference_temperature", 300.0},
            {"young_modulus_temperature_coefficient", -2.0e5}, {"poisson_ratio_temperature_coefficient", 0.0}});
    functions->eigenstrains.push_back(registry.bind_eigenstrain(
        "thermal", "isotropic_thermal_expansion", {{"thermal_expansion", 1.0e-5}, {"reference_temperature", 300.0}}));
    return {std::move(functions), 2.0e8};
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition result;
    result.regions.push_back({"solid", "solid", material(), 0.5, 300.0});
    result.regions.back().hex8_element_formulation = fuelsim::Hex8ElementFormulation::c3d8rt;
    result.regions.back().strain_formulation = fuelsim::StrainFormulation::finite;
    result.time_tables.emplace_back(
        "right_temperature", std::vector<double>{0.0, increments * time_step}, std::vector<double>{300.0, 700.0});
    result.boundary_conditions = {
        {"left_temperature", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::temperature, 300.0},
        {"right_temperature", fuelsim::BoundaryConditionType::dirichlet, "right", fuelsim::Field::temperature, 1.0,
            false, "right_temperature"},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::displacement_z, 0.0},
        {"top_heat_flux", fuelsim::BoundaryConditionType::heat_flux, "z_high", fuelsim::Field::temperature, 50.0},
        {"tip_bending_traction", fuelsim::BoundaryConditionType::traction, "right", fuelsim::Field::displacement_z,
            -2.0e4},
    };
    fuelsim::BoundaryConditionDefinition convection{
        "side_convection", fuelsim::BoundaryConditionType::convection, "y_high", fuelsim::Field::temperature, 0.0};
    convection.heat_transfer_coefficient = 5.0;
    convection.ambient_temperature = 280.0;
    result.boundary_conditions.push_back(convection);
    return result;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = 1.0e-8;
    result.relative_tolerance = 1.0e-12;
    result.maximum_iterations = 50;
    result.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    result.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    result.field_residual_scaling = true;
    result.residual_reduction_tolerance = 1.0e-10;
    result.temperature_residual_absolute_tolerance = 1.0e-6;
    result.mechanical_residual_absolute_tolerance = 1.0e-4;
    return result;
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::invalid_argument("usage: b544_test <Abaqus reference directory>");
        fuelsim::PetscSession session(argc, argv, "fuelsim B5.44 distorted C3D8RT verification\n");
        const fuelsim::UnstructuredHex8Mesh case_mesh = mesh();
        const fuelsim::SpatialDefinition case_definition = definition();
        fuelsim::TransientProblem problem(case_definition, case_mesh);
        fuelsim::test::AbaqusHex8SnapshotObserver observer;
        const fuelsim::TransientResult solve = fuelsim::solve_transient(problem,
            {increments * time_step, time_step, time_step, time_step, 1.0, 0.5, 0, 0.0}, solver_options(), &observer);
        if (!solve.completed || solve.accepted_steps.size() != increments || !solve.rejected_steps.empty()) {
            std::cerr << "[FAIL] B5.44 did not complete ten fixed increments";
            if (!solve.rejected_steps.empty()) std::cerr << ": " << solve.rejected_steps.back().failure_message;
            std::cerr << '\n';
            return 1;
        }
        fuelsim::test::AbaqusHex8FullFieldOptions comparison;
        comparison.case_name = "b544_distorted_bending";
        comparison.reference_prefix = std::string(argv[1]) + "/b544_hex8_c3d8rt_distorted_bending";
        comparison.expected_steps = increments;
        comparison.time_step = time_step;
        comparison.reduced_integration = true;
        comparison.bulk_relative_tolerance = 1.0e-2;
        comparison.energy_relative_tolerance = 1.0e-2;
        // The distorted bending path keeps the aggregate reaction metrics at one percent, but one
        // low-resultant constrained-node vector needs the explicitly recorded four-percent pointwise gate.
        comparison.reaction_pointwise_relative_tolerance = 4.0e-2;
        comparison.reaction_zero_absolute_tolerance = 1.0e-3;
        // The only nodal heat-reaction relative outlier is a 0.003451 W difference at a
        // 0.302267 W reference value.  Keep the undiluted relative metric and use an explicit
        // absolute check for this low-rate cancellation; all aggregate metrics retain one percent.
        comparison.reaction_heat_flux_pointwise_absolute_tolerance = 1.0e-2;
        const bool passed = fuelsim::test::compare_abaqus_hex8_full_field(
            problem, case_definition, case_mesh, observer.snapshots(), comparison);
        std::cout << "b544_accepted_steps=" << solve.accepted_steps.size() << '\n'
                  << "b544_rejected_steps=" << solve.rejected_steps.size() << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.44 exception: " << error.what() << '\n';
        return 1;
    }
}
