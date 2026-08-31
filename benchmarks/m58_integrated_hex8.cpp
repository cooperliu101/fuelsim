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

struct ContactReference final {
    fuelsim::CartesianPoint3 point;
    double pressure;
};

struct ElementReference final {
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> fields;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

void write_dof_reference(const std::string& path, const std::vector<double>& state) {
    std::ofstream output(path, std::ios::out | std::ios::trunc);
    if (!output) throw std::runtime_error("Could not write M5.8 degree-of-freedom reference: " + path);
    output << state.size() << '\n' << std::setprecision(17);
    for (double value : state) output << value << '\n';
    if (!output) throw std::runtime_error("Could not complete M5.8 degree-of-freedom reference: " + path);
}

bool compare_dof_reference(const std::string& path, const std::vector<double>& state) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read M5.8 degree-of-freedom reference: " + path);
    std::size_t count = 0;
    input >> count;
    std::vector<double> reference(count, 0.0);
    for (double& value : reference) input >> value;
    if (!input || reference.size() != state.size())
        throw std::runtime_error("M5.8 degree-of-freedom reference size differs");
    if (state.size() % 4U != 0U) throw std::runtime_error("M5.8 state does not contain four complete fields");
    const std::size_t nodes = state.size() / 4U;
    double maximum_temperature_difference = 0.0, maximum_displacement_difference = 0.0, maximum_tolerance_ratio = 0.0;
    std::size_t maximum_index = 0;
    for (std::size_t index = 0; index < state.size(); ++index) {
        const double difference = std::abs(state[index] - reference[index]);
        const bool temperature = index < nodes;
        const double tolerance =
            temperature ? 1.0e-6 + 1.0e-9 * std::abs(reference[index]) : 1.0e-10 + 1.0e-7 * std::abs(reference[index]);
        if (temperature)
            maximum_temperature_difference = std::max(maximum_temperature_difference, difference);
        else
            maximum_displacement_difference = std::max(maximum_displacement_difference, difference);
        if (difference / tolerance > maximum_tolerance_ratio) {
            maximum_tolerance_ratio = difference / tolerance;
            maximum_index = index;
        }
    }
    std::cout << "m58_mpi_temperature_maximum_absolute_difference=" << maximum_temperature_difference << '\n'
              << "m58_mpi_displacement_maximum_absolute_difference=" << maximum_displacement_difference << '\n'
              << "m58_mpi_maximum_tolerance_ratio=" << maximum_tolerance_ratio << '\n'
              << "m58_mpi_maximum_tolerance_index=" << maximum_index << '\n';
    return check(maximum_tolerance_ratio <= 1.0,
        "M5.8 one-process and multi-process degrees of freedom satisfy absolute and relative tolerances");
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
    if (found == header.end()) throw std::invalid_argument("Missing column '" + name + "' in " + path);
    return static_cast<std::size_t>(found - header.begin());
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete M5.8 MOOSE row in " + path);
    return std::stod(values[index]);
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read M5.8 MOOSE nodes: " + path);
    std::string line;
    std::getline(input, line);
    const auto header = split_csv(line);
    const std::array<std::size_t, 7> columns = {column(header, "x", path), column(header, "y", path),
        column(header, "z", path), column(header, "T", path), column(header, "disp_x", path),
        column(header, "disp_y", path), column(header, "disp_z", path)};
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        result.push_back(
            {{number(values, columns[0], path), number(values, columns[1], path), number(values, columns[2], path)},
                {number(values, columns[3], path), number(values, columns[4], path), number(values, columns[5], path),
                    number(values, columns[6], path)}});
    }
    return result;
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read M5.8 MOOSE contact pressure: " + path);
    std::string line;
    std::getline(input, line);
    const auto header = split_csv(line);
    const std::array<std::size_t, 4> columns = {column(header, "x", path), column(header, "y", path),
        column(header, "z", path), column(header, "contact_pressure", path)};
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        result.push_back(
            {{number(values, columns[0], path), number(values, columns[1], path), number(values, columns[2], path)},
                number(values, columns[3], path)});
    }
    return result;
}

std::vector<ElementReference> read_elements(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read M5.8 MOOSE element states: " + path);
    std::string line;
    std::getline(input, line);
    const auto header = split_csv(line);
    const std::array<std::size_t, 6> columns = {column(header, "x", path), column(header, "y", path),
        column(header, "z", path), column(header, "vonmises_stress", path),
        column(header, "effective_plastic_strain", path), column(header, "effective_creep_strain", path)};
    std::vector<ElementReference> result;
    while (std::getline(input, line)) {
        const auto values = split_csv(line);
        result.push_back({{number(values, columns[0], path), number(values, columns[1], path),
                              number(values, columns[2], path)},
            {number(values, columns[3], path), number(values, columns[4], path), number(values, columns[5], path)}});
    }
    return result;
}

double equivalent_stress(const fuelsim::SymmetricTensor3Values& stress) {
    const double mean = (stress.xx + stress.yy + stress.zz) / 3.0;
    return std::sqrt(1.5 * ((stress.xx - mean) * (stress.xx - mean) + (stress.yy - mean) * (stress.yy - mean) +
                               (stress.zz - mean) * (stress.zz - mean) +
                               2.0 * (stress.xy * stress.xy + stress.yz * stress.yz + stress.xz * stress.xz)));
}

bool same_point(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    constexpr double tolerance = 1.0e-12;
    return std::abs(first.x - second.x) < tolerance && std::abs(first.y - second.y) < tolerance &&
           std::abs(first.z - second.z) < tolerance;
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = definition.solver.absolute_tolerance;
    result.relative_tolerance = definition.solver.relative_tolerance;
    result.step_tolerance = definition.solver.step_tolerance;
    result.maximum_iterations = definition.solver.maximum_iterations;
    result.linear_solver = definition.solver.linear_solver;
    result.preconditioner = definition.solver.preconditioner;
    result.direct_factorization = definition.solver.direct_factorization;
    result.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
    result.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
    result.jacobian_lag = definition.solver.jacobian_lag;
    result.line_search = definition.solver.line_search;
    result.backtracking_fallback = definition.solver.backtracking_fallback;
    result.field_residual_scaling = definition.solver.field_residual_scaling;
    result.residual_reduction_tolerance = definition.solver.residual_reduction_tolerance;
    result.temperature_residual_absolute_tolerance = definition.solver.temperature_residual_absolute_tolerance;
    result.mechanical_residual_absolute_tolerance = definition.solver.mechanical_residual_absolute_tolerance;
    result.temperature_residual_scale = definition.solver.temperature_residual_scale;
    result.mechanical_residual_scale = definition.solver.mechanical_residual_scale;
    return result;
}

fuelsim::TransientTimeOptions time_options(const fuelsim::FuelSimCaseDefinition& definition) {
    const auto& input = definition.transient_execution;
    return {input.end_time, input.initial_time_step, input.minimum_time_step, input.maximum_time_step,
        input.growth_factor, input.cutback_factor, input.maximum_cutbacks_per_step, input.load_ramp_time,
        input.target_nonlinear_iterations, input.iteration_window, input.time_error_relative_tolerance,
        input.temperature_time_absolute_tolerance, input.displacement_time_absolute_tolerance,
        input.time_error_safety_factor, input.strain_history_time_absolute_tolerance,
        input.stress_history_time_absolute_tolerance, input.include_thermal_time_term, input.use_linear_time_predictor};
}

bool below(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics, bool output) {
    if (output) fuelsim::test::print_relative_metrics(name, metrics);
    return check(fuelsim::test::relative_metrics_below(metrics, 5.0e-3),
        name + " relative L2, relative absolute-peak, and maximum pointwise errors are below 0.5 percent");
}

bool below_with_pointwise_tolerance(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics,
    double aggregate_tolerance, double pointwise_tolerance, bool output) {
    if (output) fuelsim::test::print_relative_metrics(name, metrics);
    return check(fuelsim::test::relative_metrics_below_with_pointwise_tolerance(
                     metrics, aggregate_tolerance, pointwise_tolerance),
        name + " relative L2 and relative absolute-peak errors satisfy the aggregate tolerance, and the maximum "
               "pointwise error satisfies its explicit qualified tolerance");
}

bool run(const std::string& input_path, const std::string& nodal_path, const std::string& contact_path,
    const std::string& element_path, const std::string& reference_mode, const std::string& reference_path,
    bool output) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::transient ||
        definition.geometry != fuelsim::CaseGeometry::cartesian_3d || definition.spatial.contacts.size() != 1 ||
        !definition.spatial.contacts[0].thermal || !definition.spatial.contacts[0].mechanical)
        throw std::invalid_argument("M5.8 requires one coupled transient three-dimensional contact pair");
    const fuelsim::UnstructuredHex8Mesh source = fuelsim::read_exodus_hex8(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, source);
    const auto initial_contact =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, problem.committed_solution());
    double initial_minimum_gap = std::numeric_limits<double>::infinity();
    std::size_t initial_projected_nodes = 0;
    for (const fuelsim::CartesianContactNodeSummary& node : initial_contact)
        if (node.projected) {
            initial_minimum_gap = std::min(initial_minimum_gap, node.gap);
            ++initial_projected_nodes;
        }
    const fuelsim::TransientResult solve =
        fuelsim::solve_transient(problem, time_options(definition), solver_options(definition));
    if (output && !solve.completed) {
        std::cerr << "m58_failure_category="
                  << fuelsim::solve_failure_category_name(solve.last_attempt.failure_category)
                  << " m58_failure_message=" << solve.last_attempt.failure_message << '\n';
        for (std::size_t field = 0; field < solve.last_attempt.field_names.size(); ++field)
            std::cerr << "m58_failure_field=" << solve.last_attempt.field_names[field]
                      << " initial=" << solve.last_attempt.initial_field_residual_norms[field]
                      << " final=" << solve.last_attempt.final_field_residual_norms[field]
                      << " scaled_final=" << solve.last_attempt.final_scaled_field_residual_norms[field]
                      << " scale=" << solve.last_attempt.field_residual_scalings[field] << '\n';
    }
    const std::size_t rank = static_cast<std::size_t>(solve.last_attempt.mpi_rank);
    const std::size_t ranks = static_cast<std::size_t>(solve.last_attempt.mpi_size);
    const auto local_partition = problem.contribution_partition(rank, ranks);
    const std::size_t expected_steps = static_cast<std::size_t>(
        std::llround(definition.transient_execution.end_time / definition.transient_execution.initial_time_step));
    bool passed = check(initial_projected_nodes == initial_contact.size() && initial_minimum_gap > 0.0,
                      "M5.8 starts from a positive projected mechanical gap") &&
                  check(solve.completed && solve.accepted_steps.size() == expected_steps,
                      "M5.8 completes every configured equal Backward Euler step") &&
                  check(solve.aggregate_timing.workspace_setups == 1,
                      "M5.8 reuses one PETSc workspace across the complete path") &&
                  check(local_partition.first == solve.last_attempt.local_contribution_begin &&
                            local_partition.second == solve.last_attempt.local_contribution_end,
                      "M5.8 assigns this message-passing rank its declared contribution interval");
    for (std::size_t partition = 0; partition + 1U < ranks; ++partition) {
        const auto left = problem.contribution_partition(partition, ranks);
        const auto right = problem.contribution_partition(partition + 1U, ranks);
        passed = check(left.second == right.first, "M5.8 contribution intervals are contiguous and do not overlap") &&
                 passed;
    }
    passed = check(problem.contribution_partition(0, ranks).first == 0 &&
                       problem.contribution_partition(ranks - 1U, ranks).second == problem.contribution_count(),
                 "M5.8 contribution intervals cover every runtime contribution exactly once") &&
             passed;
    std::ostringstream local_timing;
    local_timing << "m58_local_assembly_rank=" << rank
                 << " residual_seconds=" << solve.aggregate_timing.local_residual_assembly_seconds
                 << " jacobian_seconds=" << solve.aggregate_timing.local_jacobian_assembly_seconds
                 << " contribution_begin=" << local_partition.first << " contribution_end=" << local_partition.second
                 << '\n';
    std::cout << local_timing.str() << std::flush;

    if (reference_mode == "write_dof") {
        if (output && passed) write_dof_reference(reference_path, problem.committed_solution());
        return passed;
    }
    if (reference_mode == "compare_dof") {
        if (output) passed = compare_dof_reference(reference_path, problem.committed_solution()) && passed;
        return passed;
    }

    const auto node_reference = read_nodes(nodal_path);
    if (node_reference.size() != source.nodes().size()) throw std::invalid_argument("M5.8 node counts differ");
    std::array<fuelsim::test::FieldErrorMetrics, 4> nodal;
    fuelsim::test::FieldErrorMetrics radial_displacement, tangential_displacement, tangential_analytic_zero;
    double maximum_coordinate_difference = 0.0;
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        const std::size_t offset = fuelsim::cartesian::ProblemAccess::region_node_offset(problem, region);
        for (std::size_t local = 0; local < mesh.nodes().size(); ++local) {
            const std::size_t source_node = mesh.source_node_ids()[local];
            const fuelsim::CartesianPoint3& point = source.nodes().at(source_node);
            const auto reference = std::find_if(node_reference.begin(), node_reference.end(),
                [&](const NodeReference& value) { return same_point(value.point, point); });
            if (reference == node_reference.end())
                throw std::invalid_argument("M5.8 source-node coordinate is missing");
            const NodeReference& expected = *reference;
            std::array<double, 4> reference_fields = expected.fields;
            // The cylindrical axis has analytical zero transverse displacement. Classify MOOSE roundoff there as a
            // zero reference rather than adding a denominator floor to the pointwise relative metric.
            if (std::hypot(point.x, point.y) < 1.0e-15) {
                reference_fields[1] = 0.0;
                reference_fields[2] = 0.0;
            }
            maximum_coordinate_difference = std::max(maximum_coordinate_difference,
                std::max({std::abs(point.x - expected.point.x), std::abs(point.y - expected.point.y),
                    std::abs(point.z - expected.point.z)}));
            std::array<double, 4> actual{};
            for (std::size_t field = 0; field < nodal.size(); ++field) {
                actual[field] = problem.committed_solution()[dofs.dof(
                    std::array<fuelsim::Field, 4>{fuelsim::Field::temperature, fuelsim::Field::displacement_x,
                        fuelsim::Field::displacement_y, fuelsim::Field::displacement_z}[field],
                    offset + local)];
                nodal[field].add(actual[field], reference_fields[field]);
            }
            const double radius = std::hypot(point.x, point.y);
            if (radius > 1.0e-15) {
                radial_displacement.add((point.x * actual[1] + point.y * actual[2]) / radius,
                    (point.x * reference_fields[1] + point.y * reference_fields[2]) / radius);
                tangential_displacement.add((-point.y * actual[1] + point.x * actual[2]) / radius,
                    (-point.y * reference_fields[1] + point.x * reference_fields[2]) / radius);
                tangential_analytic_zero.add((-point.y * actual[1] + point.x * actual[2]) / radius, 0.0);
            }
        }
    }
    const std::array<std::string, 4> nodal_names = {
        "temperature", "displacement_x", "displacement_y", "displacement_z"};
    const bool abaqus_reference = reference_mode == "abaqus";
    for (std::size_t field = 0; field < nodal.size(); ++field) {
        if (abaqus_reference && (field == 1 || field == 2)) {
            if (output) fuelsim::test::print_relative_metrics("m58_" + nodal_names[field], nodal[field]);
        } else {
            passed = below("m58_" + nodal_names[field], nodal[field], output) && passed;
        }
    }
    if (abaqus_reference) {
        passed =
            below_with_pointwise_tolerance("m58_radial_displacement", radial_displacement, 5.0e-3, 4.0e-2, output) &&
            passed;
        if (output) {
            fuelsim::test::print_relative_metrics("m58_tangential_displacement", tangential_displacement);
            fuelsim::test::print_absolute_metrics("m58_tangential_analytic_zero", tangential_analytic_zero);
        }
        passed = check(tangential_analytic_zero.maximum_absolute_difference < 1.0e-6,
                     "M5.8 analytical-zero tangential displacement remains below one micrometre") &&
                 passed;
    } else if (output) {
        fuelsim::test::print_relative_metrics("m58_radial_displacement", radial_displacement);
        fuelsim::test::print_relative_metrics("m58_tangential_displacement", tangential_displacement);
    }

    const auto contact =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, problem.committed_solution());
    const auto contact_reference = read_contact(contact_path);
    const auto source_nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    fuelsim::test::FieldErrorMetrics pressure;
    std::size_t sliding = 0, crossed_faces = 0;
    for (std::size_t node = 0; node < contact.size(); ++node) {
        const fuelsim::CartesianPoint3& point = source.nodes().at(source_nodes[node]);
        const auto found = std::find_if(contact_reference.begin(), contact_reference.end(),
            [&](const ContactReference& value) { return same_point(value.point, point); });
        if (found == contact_reference.end()) throw std::invalid_argument("M5.8 contact-node coordinate is missing");
        pressure.add(contact[node].pressure, found->pressure);
        if (contact[node].sliding) ++sliding;
        if (contact[node].primary_face != initial_contact[node].primary_face) ++crossed_faces;
    }
    passed = below("m58_contact_pressure", pressure, output) && passed;

    const auto element_reference = read_elements(element_path);
    std::array<fuelsim::test::FieldErrorMetrics, 3> material;
    double maximum_element_coordinate_difference = 0.0;
    const std::size_t cladding = 1;
    const auto& cladding_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, cladding);
    for (std::size_t element = 0; element < cladding_mesh.elements().size(); ++element) {
        fuelsim::CartesianPoint3 centroid{};
        for (std::size_t node : cladding_mesh.elements()[element].nodes) {
            centroid.x += cladding_mesh.nodes()[node].x / 8.0;
            centroid.y += cladding_mesh.nodes()[node].y / 8.0;
            centroid.z += cladding_mesh.nodes()[node].z / 8.0;
        }
        const auto found = std::find_if(element_reference.begin(), element_reference.end(),
            [&](const ElementReference& value) { return same_point(value.point, centroid); });
        if (found == element_reference.end()) throw std::invalid_argument("M5.8 cladding element is missing");
        maximum_element_coordinate_difference = std::max(maximum_element_coordinate_difference,
            std::max({std::abs(centroid.x - found->point.x), std::abs(centroid.y - found->point.y),
                std::abs(centroid.z - found->point.z)}));
        const auto& geometry = fuelsim::cartesian::ProblemAccess::region_element_geometry(problem, cladding, element);
        const auto& history = fuelsim::cartesian::ProblemAccess::material_history(problem, cladding, element);
        std::array<double, 3> average{};
        double measure = 0.0;
        for (std::size_t q = 0; q < history.size(); ++q) {
            const double weight = geometry.points[q].weighted_measure;
            measure += weight;
            average[0] += weight * equivalent_stress(history[q].stress);
            average[1] += weight * history[q].equivalent_plastic_strain;
            average[2] += weight * history[q].equivalent_creep_strain;
        }
        for (std::size_t field = 0; field < material.size(); ++field)
            material[field].add(average[field] / measure, found->fields[field]);
    }
    const std::array<std::string, 3> material_names = {
        "equivalent_stress", "equivalent_plastic_strain", "equivalent_creep_strain"};
    for (std::size_t field = 0; field < material.size(); ++field)
        passed = below("m58_" + material_names[field], material[field], output) && passed;

    const fuelsim::InterfaceSummary interface =
        fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, problem.committed_solution());
    passed = check(maximum_coordinate_difference < 1.0e-12 && maximum_element_coordinate_difference < 1.0e-12,
                 "M5.8 compares MOOSE fields at matching source coordinates") &&
             check(interface.active_contact_nodes > 0 && sliding > 0 && crossed_faces > 0,
                 "M5.8 keeps the contact pair active and exercises sliding across primary faces") &&
             check(interface.total_heat_rate > 0.0 && interface.total_contact_force > 0.0 &&
                       interface.total_tangential_force > 0.0 &&
                       std::abs(problem.last_conservation_summary().interface_heat_imbalance) < 1.0e-12,
                 "M5.8 exercises conservative heat, normal-force, and friction-force transfer") &&
             check(material[1].maximum_actual > 0.0 && material[2].maximum_actual > 0.0,
                 "M5.8 activates cladding plasticity and creep") &&
             passed;
    if (output)
        std::cout << "m58_dofs=" << problem.dof_count() << '\n'
                  << "m58_initial_minimum_gap=" << initial_minimum_gap << '\n'
                  << "m58_accepted_steps=" << solve.accepted_steps.size() << '\n'
                  << "m58_rejected_steps=" << solve.rejected_steps.size() << '\n'
                  << "m58_nonlinear_iterations=" << solve.total_nonlinear_iterations << '\n'
                  << "m58_residual_evaluations=" << solve.aggregate_timing.residual_evaluations << '\n'
                  << "m58_jacobian_evaluations=" << solve.aggregate_timing.jacobian_evaluations << '\n'
                  << "m58_setup_seconds=" << solve.aggregate_timing.setup_seconds << '\n'
                  << "m58_residual_callback_seconds=" << solve.aggregate_timing.residual_callback_seconds << '\n'
                  << "m58_jacobian_callback_seconds=" << solve.aggregate_timing.jacobian_callback_seconds << '\n'
                  << "m58_memory_initial_resident_bytes=" << solve.aggregate_timing.initial_resident_bytes << '\n'
                  << "m58_memory_setup_resident_bytes=" << solve.aggregate_timing.setup_resident_bytes << '\n'
                  << "m58_memory_solve_resident_bytes=" << solve.aggregate_timing.solve_resident_bytes << '\n'
                  << "m58_memory_final_resident_bytes=" << solve.aggregate_timing.final_resident_bytes << '\n'
                  << "m58_memory_minimum_peak_resident_bytes=" << solve.aggregate_timing.minimum_peak_resident_bytes
                  << '\n'
                  << "m58_memory_maximum_peak_resident_bytes=" << solve.aggregate_timing.maximum_peak_resident_bytes
                  << '\n'
                  << "m58_memory_total_peak_resident_bytes=" << solve.aggregate_timing.total_peak_resident_bytes << '\n'
                  << "m58_minimum_residual_assembly_seconds="
                  << solve.aggregate_timing.minimum_residual_assembly_seconds << '\n'
                  << "m58_maximum_residual_assembly_seconds="
                  << solve.aggregate_timing.maximum_residual_assembly_seconds << '\n'
                  << "m58_minimum_jacobian_assembly_seconds="
                  << solve.aggregate_timing.minimum_jacobian_assembly_seconds << '\n'
                  << "m58_maximum_jacobian_assembly_seconds="
                  << solve.aggregate_timing.maximum_jacobian_assembly_seconds << '\n'
                  << "m58_nonlinear_solve_seconds=" << solve.aggregate_timing.nonlinear_solve_seconds << '\n'
                  << "m58_total_seconds=" << solve.aggregate_timing.total_seconds << '\n'
                  << "m58_mpi_rank=" << solve.last_attempt.mpi_rank << '\n'
                  << "m58_mpi_size=" << solve.last_attempt.mpi_size << '\n'
                  << "m58_local_contribution_begin=" << solve.last_attempt.local_contribution_begin << '\n'
                  << "m58_local_contribution_end=" << solve.last_attempt.local_contribution_end << '\n'
                  << "m58_maximum_shadow_state_dofs=" << solve.last_attempt.maximum_shadow_state_dofs << '\n'
                  << "m58_total_shadow_state_dofs=" << solve.last_attempt.total_shadow_state_dofs << '\n'
                  << "m58_total_remote_shadow_state_dofs=" << solve.last_attempt.total_remote_shadow_state_dofs << '\n'
                  << "m58_active_contact_nodes=" << interface.active_contact_nodes << '\n'
                  << "m58_sliding_contact_nodes=" << sliding << '\n'
                  << "m58_nodes_crossing_primary_faces=" << crossed_faces << '\n';
    if (output && passed && reference_mode == "write")
        write_dof_reference(reference_path, problem.committed_solution());
    else if (output && passed && reference_mode == "compare")
        passed = compare_dof_reference(reference_path, problem.committed_solution()) && passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5 && argc != 6 && argc != 7) {
        std::cerr << "Usage: fuelsim_m58_integrated_hex8_benchmark "
                     "<case.fsi> <all-nodes.csv> <contact.csv> <element-state.csv> "
                     "[abaqus|write|compare|write_dof|compare_dof [dof-reference]]\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M5.8 integrated Hex8 benchmark\n");
        const std::string reference_mode = argc >= 6 ? argv[5] : "";
        const std::string reference_path = argc == 7 ? argv[6] : "";
        if (!reference_mode.empty() && reference_mode != "abaqus" && reference_mode != "write" &&
            reference_mode != "compare" && reference_mode != "write_dof" && reference_mode != "compare_dof")
            throw std::invalid_argument(
                "M5.8 reference mode must be abaqus, write, compare, write_dof, or compare_dof");
        if ((reference_mode == "abaqus") != (argc == 6))
            throw std::invalid_argument("M5.8 Abaqus mode takes no degree-of-freedom reference path");
        if (!run(argv[1], argv[2], argv[3], argv[4], reference_mode, reference_path, session.rank() == 0)) return 1;
        if (session.rank() == 0)
            std::cout << (reference_mode == "write_dof" || reference_mode == "compare_dof"
                              ? "[PASS] M5.8 integrated Hex8 degree-of-freedom reference\n"
                          : reference_mode == "abaqus" ? "[PASS] M5.8 integrated Hex8 Abaqus comparison\n"
                                                       : "[PASS] M5.8 integrated Hex8 MOOSE comparison\n");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M5.8 benchmark raised: " << error.what() << '\n';
        return 1;
    }
}
