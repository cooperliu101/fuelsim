#include "contact_types.hpp"
#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/field_error_metrics.hpp"
#include "support/test_support.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr std::size_t step_count = 5;
constexpr std::size_t node_count = 20;

struct PathStep final {
    const char* name;
    std::array<double, 3> displacement;
    bool active;
};

constexpr std::array<PathStep, step_count> path = {{{"CLOSE", {-0.01, 0.05, 0.04}, true},
    {"SLIDE", {-0.01, 0.95, 0.20}, true},
    {"OPEN", {0.02, 0.95, 0.20}, false},
    {"OPEN_CROSS", {0.02, 0.05, 0.40}, false},
    {"RECONTACT", {-0.01, 0.05, 0.40}, true}}};

struct NodeReference final {
    double temperature = 0.0, reaction_heat_flux = 0.0;
    std::array<double, 3> displacement{}, reaction_force{};
};

using Reference = std::array<std::array<NodeReference, node_count>, step_count>;

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        result.push_back(value);
    return result;
}

Reference read_reference(const std::string& path_name) {
    std::ifstream input(path_name);
    if (!input)
        throw std::runtime_error("Could not read B5.21 Abaqus nodal reference: " + path_name);
    std::string line;
    std::getline(input, line);
    if (line != "step,node,temperature_k,reaction_heat_flux_w,u1_m,u2_m,u3_m,rf1_n,rf2_n,rf3_n")
        throw std::invalid_argument("Unexpected B5.21 Abaqus header");
    Reference result{};
    std::array<std::array<bool, node_count>, step_count> present{};
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 10)
            throw std::invalid_argument("Unexpected B5.21 Abaqus column count");
        const std::size_t step = std::stoul(values[0]), node = std::stoul(values[1]);
        if (step < 1 || step > step_count || node < 1 || node > node_count || present[step - 1][node - 1])
            throw std::invalid_argument("Invalid or duplicate B5.21 Abaqus row");
        present[step - 1][node - 1] = true;
        result[step - 1][node - 1] = {std::stod(values[2]),
            std::stod(values[3]),
            {std::stod(values[4]), std::stod(values[5]), std::stod(values[6])},
            {std::stod(values[7]), std::stod(values[8]), std::stod(values[9])}};
    }
    for (const auto& step : present)
        if (std::find(step.begin(), step.end(), false) != step.end())
            throw std::invalid_argument("B5.21 Abaqus reference is incomplete");
    return result;
}

fuelsim::Hex8Element append_cuboid(std::vector<fuelsim::CartesianPoint3>& nodes,
    std::map<std::array<double, 3>, std::size_t>& node_map,
    double x0,
    double x1,
    double y0,
    double y1,
    double z0,
    double z1) {
    const std::array<fuelsim::CartesianPoint3, 8> points = {{{x0, y0, z0},
        {x1, y0, z0},
        {x1, y1, z0},
        {x0, y1, z0},
        {x0, y0, z1},
        {x1, y0, z1},
        {x1, y1, z1},
        {x0, y1, z1}}};
    fuelsim::Hex8Element element{};
    for (std::size_t local = 0; local < points.size(); ++local) {
        const std::array<double, 3> key = {points[local].x, points[local].y, points[local].z};
        const auto inserted = node_map.emplace(key, nodes.size());
        if (inserted.second)
            nodes.push_back(points[local]);
        element.nodes[local] = inserted.first->second;
    }
    return element;
}

fuelsim::UnstructuredHex8Mesh mesh() {
    std::vector<fuelsim::CartesianPoint3> nodes;
    std::map<std::array<double, 3>, std::size_t> primary_nodes, secondary_nodes;
    const fuelsim::Hex8Element primary_lower = append_cuboid(nodes, primary_nodes, 0.0, 1.0, 0.0, 1.0, -1.0, 2.0),
                               primary_upper = append_cuboid(nodes, primary_nodes, 0.0, 1.0, 1.0, 2.0, -1.0, 2.0),
                               secondary = append_cuboid(nodes, secondary_nodes, 1.0, 2.0, 0.1, 0.9, 0.1, 0.9);
    std::vector<std::size_t> primary_all, secondary_all;
    for (const auto& entry : primary_nodes)
        primary_all.push_back(entry.second);
    for (const auto& entry : secondary_nodes)
        secondary_all.push_back(entry.second);
    std::sort(primary_all.begin(), primary_all.end());
    std::sort(secondary_all.begin(), secondary_all.end());
    return fuelsim::UnstructuredHex8Mesh(std::move(nodes),
        {primary_lower, primary_upper, secondary},
        {1, 1, 2},
        {{1, "primary"}, {2, "secondary"}},
        {{10, "primary_all", primary_all}, {20, "secondary_all", secondary_all}},
        {{40, "primary_contact", {{0, 1}, {1, 1}}}, {50, "secondary_contact", {{2, 3}}}});
}

fuelsim::ThermoelasticProperties material() {
    return fuelsim::test::thermoelastic(0.0, 10.0, 1.0e9, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 1.0, 1.0);
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition result;
    result.regions = {{"primary", "primary", material(), 0.0, 300.0, -1, "", fuelsim::StrainFormulation::finite},
        {"secondary", "secondary", material(), 0.0, 400.0, -1, "", fuelsim::StrainFormulation::finite}};
    fuelsim::ContactDefinition contact;
    contact.name = "thermal_path";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.thermal = true;
    contact.mechanical = true;
    contact.gap_conductivity = 1.0;
    contact.minimum_gap = 1.0;
    contact.penalty = 1.0e5;
    contact.friction_coefficient = 0.3;
    contact.friction_slip_tolerance = 0.025;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::finite;
    contact.gap_heat_conductance_law = fuelsim::GapHeatConductanceLaw::affine;
    contact.gap_conductance = 0.0;
    contact.gap_conductance_pressure_derivative = 0.02;
    result.contacts.push_back(contact);
    return result;
}

std::vector<double> state_for_step(const fuelsim::SteadyProblem& problem, std::size_t step) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> state = problem.initial_state();
    for (std::size_t local = 0; local < spatial.region_mesh(1).nodes().size(); ++local) {
        const std::size_t global = spatial.global_node(1, local);
        state[spatial.dof(fuelsim::Field::temperature, global)] = 400.0;
        state[spatial.dof(fuelsim::Field::displacement_x, global)] = path[step].displacement[0];
        state[spatial.dof(fuelsim::Field::displacement_y, global)] = path[step].displacement[1];
        state[spatial.dof(fuelsim::Field::displacement_z, global)] = path[step].displacement[2];
    }
    return state;
}

std::vector<double> state_for_step(const fuelsim::TransientProblem& problem, std::size_t step) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::vector<double> state = problem.committed_solution();
    for (std::size_t local = 0; local < spatial.region_mesh(1).nodes().size(); ++local) {
        const std::size_t global = spatial.global_node(1, local);
        state[spatial.dof(fuelsim::Field::temperature, global)] = 400.0;
        state[spatial.dof(fuelsim::Field::displacement_x, global)] = path[step].displacement[0];
        state[spatial.dof(fuelsim::Field::displacement_y, global)] = path[step].displacement[1];
        state[spatial.dof(fuelsim::Field::displacement_z, global)] = path[step].displacement[2];
    }
    return state;
}

std::vector<double> flatten_committed_state(const fuelsim::TransientProblem& problem) {
    const fuelsim::TransientCommittedState state = fuelsim::cartesian::ProblemAccess::committed_state(problem);
    std::vector<double> result = {state.time, state.load_factor};
    result.insert(result.end(), state.solution.begin(), state.solution.end());
    for (const auto& region : state.cartesian_material_histories)
        for (const fuelsim::CartesianMaterialHistory& element : region)
            for (const fuelsim::CartesianMaterialPointState& point : element) {
                result.insert(result.end(), point.elastic_strain.begin(), point.elastic_strain.end());
                result.insert(result.end(), point.plastic_strain.begin(), point.plastic_strain.end());
                result.insert(result.end(), point.creep_strain.begin(), point.creep_strain.end());
                result.push_back(point.equivalent_plastic_strain);
                result.push_back(point.equivalent_creep_strain);
                result.insert(result.end(),
                    {point.stress.xx,
                        point.stress.yy,
                        point.stress.zz,
                        point.stress.xy,
                        point.stress.yz,
                        point.stress.xz});
            }
    for (const auto& contact : state.contact_histories)
        for (const fuelsim::ContactPointHistory& history : contact) {
            result.push_back(history.elastic_tangential_slip);
            result.push_back(history.sliding ? 1.0 : 0.0);
            result.push_back(history.normal_multiplier);
            result.insert(result.end(),
                history.cartesian_elastic_tangential_slip.begin(),
                history.cartesian_elastic_tangential_slip.end());
            result.insert(result.end(),
                history.cartesian_total_tangential_slip.begin(),
                history.cartesian_total_tangential_slip.end());
            result.push_back(history.cartesian_tangent_basis_initialized ? 1.0 : 0.0);
            result.insert(result.end(),
                history.cartesian_contact_normal.begin(),
                history.cartesian_contact_normal.end());
            result.insert(result.end(),
                history.cartesian_contact_tangent_first.begin(),
                history.cartesian_contact_tangent_first.end());
        }
    for (const fuelsim::TransientConservationField& field : fuelsim::transient_conservation_fields)
        result.push_back(state.conservation.*(field.member));
    return result;
}

std::vector<double> assembled_contact_residual(const fuelsim::cartesian::SpatialAssembly& spatial,
    const std::vector<double>& state,
    fuelsim::SpatialContributionType type,
    double* conservation_error = nullptr) {
    std::vector<double> global(state.size());
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != type)
            continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index)
            local[index] = state[dofs[index]];
        std::vector<double> residual;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, nullptr);
        for (std::size_t row = 0; row < dofs.size(); ++row)
            global[dofs[row]] += residual[row];
        if (conservation_error != nullptr && type == fuelsim::SpatialContributionType::thermal_contact) {
            double balance = 0.0;
            for (std::size_t row = 0; row < 8; ++row)
                balance += residual[row];
            *conservation_error = std::max(*conservation_error, std::abs(balance));
        }
    }
    return global;
}

bool compare(const Reference& reference) {
    const fuelsim::UnstructuredHex8Mesh source_mesh = mesh();
    fuelsim::SteadyProblem problem(definition(), source_mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    std::map<std::size_t, std::size_t> source_to_global;
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        for (std::size_t local = 0; local < spatial.region_mesh(region).nodes().size(); ++local)
            source_to_global[spatial.region_mesh(region).source_node_ids()[local]] = spatial.global_node(region, local);
    fuelsim::test::FieldErrorMetrics heat_flux, normal_force, heat_rate, contact_force;
    double state_difference = 0.0, conservation_error = 0.0, interface_heat_rate_difference = 0.0;
    double zero_reference_normal_force_difference = 0.0;
    bool state_sequence = true, crossed_primary_face = false;
    std::size_t first_primary_face = 0;
    for (std::size_t step = 0; step < step_count; ++step) {
        std::vector<double> state = state_for_step(problem, step);
        problem.validate_state(state);
        const std::vector<double> thermal = assembled_contact_residual(spatial,
            state,
            fuelsim::SpatialContributionType::thermal_contact,
            &conservation_error);
        const std::vector<double> mechanical =
            assembled_contact_residual(spatial, state, fuelsim::SpatialContributionType::mechanical_contact);
        const fuelsim::InterfaceSummary interface =
            fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, state);
        const std::vector<fuelsim::CartesianContactNodeSummary> nodes =
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        if (step == 0 && !nodes.empty())
            first_primary_face = nodes.front().primary_face;
        if (step == 1 && !nodes.empty() && nodes.front().primary_face != first_primary_face)
            crossed_primary_face = true;
        state_sequence = state_sequence && (interface.active_contact_nodes > 0) == path[step].active;
        double actual_heat_rate = 0.0, reference_heat_rate = 0.0, actual_normal_force = 0.0,
               reference_normal_force = 0.0;
        for (std::size_t source = 0; source < node_count; ++source) {
            const std::size_t global = source_to_global.at(source);
            const std::size_t temperature_dof = spatial.dof(fuelsim::Field::temperature, global);
            heat_flux.add(thermal[temperature_dof], reference[step][source].reaction_heat_flux);
            const std::size_t region = source < 12 ? 0 : 1;
            const double expected_temperature = region == 0 ? 300.0 : 400.0;
            state_difference =
                std::max(state_difference, std::abs(state[temperature_dof] - reference[step][source].temperature));
            state_difference =
                std::max(state_difference, std::abs(expected_temperature - reference[step][source].temperature));
            const std::array<fuelsim::Field, 3> fields = {fuelsim::Field::displacement_x,
                fuelsim::Field::displacement_y,
                fuelsim::Field::displacement_z};
            for (std::size_t component = 0; component < 3; ++component)
                state_difference = std::max(state_difference,
                    std::abs(state[spatial.dof(fields[component], global)]
                             - reference[step][source].displacement[component]));
            if (source >= 12) {
                actual_heat_rate += thermal[temperature_dof];
                reference_heat_rate += reference[step][source].reaction_heat_flux;
                actual_normal_force -= mechanical[spatial.dof(fuelsim::Field::displacement_x, global)];
                reference_normal_force -= reference[step][source].reaction_force[0];
            }
        }
        heat_rate.add(actual_heat_rate, reference_heat_rate);
        if (path[step].active) {
            contact_force.add(actual_normal_force, reference_normal_force);
            normal_force.add(interface.total_contact_force, reference_normal_force);
        } else {
            zero_reference_normal_force_difference = std::max({zero_reference_normal_force_difference,
                std::abs(actual_normal_force - reference_normal_force),
                std::abs(interface.total_contact_force - reference_normal_force)});
        }
        interface_heat_rate_difference =
            std::max(interface_heat_rate_difference, std::abs(interface.total_heat_rate - actual_heat_rate));
    }
    fuelsim::test::print_relative_metrics("b521_nodal_reaction_heat_flux", heat_flux);
    fuelsim::test::print_relative_metrics("b521_total_heat_rate", heat_rate);
    fuelsim::test::print_relative_metrics("b521_assembled_normal_force", contact_force);
    fuelsim::test::print_relative_metrics("b521_reported_normal_force", normal_force);
    std::cout << "b521_nodal_heat_zero_reference_count=" << heat_flux.zero_reference_count << '\n'
              << "b521_nodal_heat_zero_reference_maximum_absolute_difference="
              << heat_flux.maximum_zero_reference_difference << '\n'
              << "b521_state_maximum_absolute_difference=" << state_difference << '\n'
              << "b521_normal_force_zero_reference_count=2\n"
              << "b521_normal_force_zero_reference_maximum_absolute_difference="
              << zero_reference_normal_force_difference << '\n'
              << "b521_thermal_contact_conservation_maximum_absolute=" << conservation_error << '\n'
              << "b521_interface_heat_rate_maximum_absolute_difference=" << interface_heat_rate_difference << '\n';
    constexpr double relative_tolerance = 1.0e-5;
    return check(source_mesh.nodes().size() == node_count, "B5.21 Fuelsim uses the same twenty-node Abaqus geometry")
           && check(state_difference < 3.0e-8, "B5.21 prescribed temperatures and displacements match Abaqus")
           && check(state_sequence && crossed_primary_face,
               "B5.21 reproduces close, cross-face slide, open, open crossing, and recontact states")
           && check(fuelsim::test::relative_metrics_below(heat_flux, relative_tolerance)
                        && fuelsim::test::relative_metrics_below(heat_rate, relative_tolerance),
               "B5.21 nodal and total contact heat rates match Abaqus through the full path")
           && check(fuelsim::test::relative_metrics_below(contact_force, relative_tolerance)
                        && fuelsim::test::relative_metrics_below(normal_force, relative_tolerance)
                        && zero_reference_normal_force_difference < 1.0e-8,
               "B5.21 assembled and reported normal contact forces match Abaqus")
           && check(heat_flux.maximum_zero_reference_difference < 1.0e-10 && conservation_error < 1.0e-12
                        && interface_heat_rate_difference < 1.0e-12,
               "B5.21 open-state heat transfer is zero and every active thermal constraint is exactly conservative");
}

bool check_transaction_and_restart(const std::string& checkpoint_path) {
    const fuelsim::UnstructuredHex8Mesh source_mesh = mesh();
    fuelsim::TransientProblem problem(definition(), source_mesh);
    const std::vector<double> initial = flatten_committed_state(problem);
    problem.begin_time_step({1.0, 1.0, false});
    const std::vector<double> closed = state_for_step(problem, 0);
    problem.validate_state(closed);
    const std::vector<double> trial_heat = assembled_contact_residual(fuelsim::cartesian::ProblemAccess::view(problem),
        closed,
        fuelsim::SpatialContributionType::thermal_contact);
    const double trial_heat_l2 =
        std::sqrt(std::inner_product(trial_heat.begin(), trial_heat.end(), trial_heat.begin(), 0.0));
    problem.rollback_time_step();
    const bool rollback_exact = flatten_committed_state(problem) == initial;

    problem.begin_time_step({1.0, 1.0, false});
    problem.commit_time_step(closed);
    const std::vector<double> committed = flatten_committed_state(problem);
    fuelsim::write_transient_checkpoint(checkpoint_path, problem, 0.25);
    fuelsim::TransientProblem restored(definition(), source_mesh);
    const double restored_time_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restored);
    const bool restart_exact = restored_time_step == 0.25 && flatten_committed_state(restored) == committed;

    restored.begin_time_step({2.0, 1.0, false});
    const std::vector<double> open = state_for_step(restored, 2);
    restored.validate_state(open);
    const std::vector<double> open_heat = assembled_contact_residual(fuelsim::cartesian::ProblemAccess::view(restored),
        open,
        fuelsim::SpatialContributionType::thermal_contact);
    double open_heat_maximum = 0.0;
    for (double value : open_heat)
        open_heat_maximum = std::max(open_heat_maximum, std::abs(value));
    restored.rollback_time_step();
    const bool second_rollback_exact = flatten_committed_state(restored) == committed;
    (void)std::remove(checkpoint_path.c_str());
    std::cout << "b521_trial_contact_heat_l2=" << trial_heat_l2 << '\n'
              << "b521_open_trial_heat_maximum_absolute=" << open_heat_maximum << '\n'
              << "b521_rollback_exact=" << rollback_exact << '\n'
              << "b521_restart_exact=" << restart_exact << '\n';
    return check(trial_heat_l2 > 0.0 && open_heat_maximum == 0.0,
               "B5.21 transaction covers active heat transfer and an exactly open thermal-contact trial")
           && check(rollback_exact && second_rollback_exact,
               "B5.21 rollback restores the exact nodal, material, contact, time, and conservation state")
           && check(restart_exact, "B5.21 checkpoint restart restores the exact coupled thermal-contact state");
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: fuelsim_b521_hex8_c3d8t_thermal_contact_path_abaqus_tests "
                     "<nodal.csv> <checkpoint.bin>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        const bool passed = compare(read_reference(argv[1])) && check_transaction_and_restart(argv[2]);
        if (passed)
            std::cout << "[PASS] B5.21 Abaqus C3D8T thermal-contact path\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.21 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
