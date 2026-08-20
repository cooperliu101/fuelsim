#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/moose_field_comparison.hpp"
#include "support/rz_problem_access.hpp"
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
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t separator = line.find(',', begin);
        result.push_back(line.substr(begin, separator - begin));
        if (separator == std::string::npos) return result;
        begin = separator + 1;
    }
}

double final_csv_value(const std::string& path, const std::string& name) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE scalar reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE scalar reference is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("MOOSE scalar reference is missing column: " + name);
    const std::size_t column = static_cast<std::size_t>(found - header.begin());
    std::vector<std::string> final_row;
    while (std::getline(input, line))
        if (!line.empty()) final_row = split_csv(line);
    if (column >= final_row.size())
        throw std::invalid_argument("MOOSE scalar reference has no final value for: " + name);
    return std::stod(final_row[column]);
}

struct CsvTable final {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
};

CsvTable read_csv(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE integration-point reference: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE integration-point reference is empty: " + path);
    CsvTable result;
    result.header = split_csv(line);
    while (std::getline(input, line))
        if (!line.empty()) result.rows.push_back(split_csv(line));
    return result;
}

double csv_value(const CsvTable& table, const std::vector<std::string>& row, const std::string& name) {
    const auto found = std::find(table.header.begin(), table.header.end(), name);
    if (found == table.header.end())
        throw std::invalid_argument("MOOSE integration-point reference is missing column: " + name);
    const std::size_t column = static_cast<std::size_t>(found - table.header.begin());
    if (column >= row.size()) throw std::invalid_argument("MOOSE integration-point reference row is incomplete");
    return std::stod(row[column]);
}

double equivalent_stress(const fuelsim::AxisymmetricStressValues& stress) {
    const double mean = (stress.rr + stress.zz + stress.hoop) / 3.0;
    return std::sqrt(1.5 * ((stress.rr - mean) * (stress.rr - mean) + (stress.zz - mean) * (stress.zz - mean) +
                               (stress.hoop - mean) * (stress.hoop - mean) + 2.0 * stress.rz * stress.rz));
}

struct RegionAverages final {
    double plastic = 0.0;
    double creep = 0.0;
    double stress = 0.0;
    double maximum_plastic = 0.0;
    double maximum_creep = 0.0;

    struct Point final {
        double radius;
        double axial_coordinate;
        double stress;
        double plastic;
        double creep;
    };

    std::vector<Point> points;
};

RegionAverages region_averages(const fuelsim::TransientProblem& problem, std::size_t region) {
    RegionAverages result;
    double measure = 0.0;
    const fuelsim::RegionMesh& mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, region);
    result.points.reserve(mesh.elements().size() * 4);
    for (std::size_t element = 0; element < fuelsim::rz::ProblemAccess::region_mesh(problem, region).elements().size();
        ++element) {
        const fuelsim::Quad4RzGeometry& geometry =
            fuelsim::rz::ProblemAccess::region_element_geometry(problem, region, element);
        const fuelsim::Quad4Element& mesh_element = mesh.elements()[element];
        const fuelsim::Quad4MaterialHistory& history =
            fuelsim::rz::ProblemAccess::material_history(problem, region, element);
        const auto& stresses = fuelsim::rz::ProblemAccess::material_stress(problem, region, element);
        for (std::size_t q = 0; q < history.size(); ++q) {
            const double weight = geometry.points[q].weighted_measure;
            double axial_coordinate = 0.0;
            for (std::size_t node = 0; node < mesh_element.nodes.size(); ++node)
                axial_coordinate += geometry.points[q].shape[node] * mesh.nodes()[mesh_element.nodes[node]].z;
            const double stress = equivalent_stress(stresses[q]);
            measure += weight;
            result.plastic += weight * history[q].equivalent_plastic_strain;
            result.creep += weight * history[q].equivalent_creep_strain;
            result.stress += weight * stress;
            result.maximum_plastic = std::max(result.maximum_plastic, history[q].equivalent_plastic_strain);
            result.maximum_creep = std::max(result.maximum_creep, history[q].equivalent_creep_strain);
            result.points.push_back({geometry.points[q].radius, axial_coordinate, stress,
                history[q].equivalent_plastic_strain, history[q].equivalent_creep_strain});
        }
    }
    result.plastic /= measure;
    result.creep /= measure;
    result.stress /= measure;
    return result;
}

std::vector<RegionAverages::Point> read_cladding_points(
    const std::string& coordinate_path, const std::string& value_path) {
    const CsvTable coordinates = read_csv(coordinate_path);
    const CsvTable values = read_csv(value_path);
    if (values.rows.empty() || coordinates.rows.size() != values.rows.size() * 4)
        throw std::invalid_argument("MOOSE integration-point coordinate and value counts differ");
    std::vector<RegionAverages::Point> result;
    result.reserve(coordinates.rows.size());
    constexpr std::array<std::size_t, 4> moose_qp_from_fuelsim = {0, 1, 3, 2};
    for (const std::vector<std::string>& value_row : values.rows) {
        const std::size_t element = static_cast<std::size_t>(csv_value(values, value_row, "id"));
        for (const std::size_t qp : moose_qp_from_fuelsim) {
            const auto coordinate = std::find_if(
                coordinates.rows.begin(), coordinates.rows.end(), [&coordinates, element, qp](const auto& row) {
                    return static_cast<std::size_t>(csv_value(coordinates, row, "elem_id")) == element &&
                           static_cast<std::size_t>(csv_value(coordinates, row, "qp_id")) == qp;
                });
            if (coordinate == coordinates.rows.end())
                throw std::invalid_argument("MOOSE integration-point coordinate is missing");
            const std::string suffix = std::to_string(qp);
            result.push_back({csv_value(coordinates, *coordinate, "x"), csv_value(coordinates, *coordinate, "y"),
                csv_value(values, value_row, "stress_q" + suffix), csv_value(values, value_row, "plastic_q" + suffix),
                csv_value(values, value_row, "creep_q" + suffix)});
        }
    }
    return result;
}

double relative_error(double actual, double reference) { return std::abs(actual - reference) / std::abs(reference); }

fuelsim::TransientTimeOptions time_options(const fuelsim::FuelSimCaseDefinition& definition) {
    const auto& input = definition.transient_execution;
    return {input.end_time, input.initial_time_step, input.minimum_time_step, input.maximum_time_step,
        input.growth_factor, input.cutback_factor, input.maximum_cutbacks_per_step, input.load_ramp_time,
        input.target_nonlinear_iterations, input.iteration_window, input.time_error_relative_tolerance,
        input.temperature_time_absolute_tolerance, input.displacement_time_absolute_tolerance,
        input.time_error_safety_factor, input.strain_history_time_absolute_tolerance,
        input.stress_history_time_absolute_tolerance, input.include_thermal_time_term};
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = definition.solver.absolute_tolerance;
    result.relative_tolerance = definition.solver.relative_tolerance;
    result.step_tolerance = definition.solver.step_tolerance;
    result.maximum_iterations = definition.solver.maximum_iterations;
    result.backtracking_fallback = definition.solver.backtracking_fallback;
    result.field_residual_scaling = definition.solver.field_residual_scaling;
    result.residual_reduction_tolerance = definition.solver.residual_reduction_tolerance;
    result.temperature_residual_absolute_tolerance = definition.solver.temperature_residual_absolute_tolerance;
    result.mechanical_residual_absolute_tolerance = definition.solver.mechanical_residual_absolute_tolerance;
    result.temperature_residual_scale = definition.solver.temperature_residual_scale;
    result.mechanical_residual_scale = definition.solver.mechanical_residual_scale;
    return result;
}

bool run_case(const std::string& input_path, const std::string& nodal_reference_path,
    const std::string& pressure_reference_path, const std::string& qp_coordinate_path, const std::string& qp_value_path,
    const std::string& scalar_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(input_path);
    if (definition.problem != fuelsim::CaseProblem::transient || definition.spatial.contacts.size() != 1)
        throw std::invalid_argument("M5.7 integrated validation requires one transient contact pair");
    const fuelsim::UnstructuredQuad4Mesh source = fuelsim::read_exodus_quad4(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, source);
    const std::vector<fuelsim::ContactNodeSummary> initial_contact =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, problem.committed_solution());
    const fuelsim::TransientResult solve =
        fuelsim::solve_transient(problem, time_options(definition), solver_options(definition));
    if (!solve.completed) return check(false, "M5.7 integrated adaptive transient completes");
    const std::vector<fuelsim::ContactNodeSummary> contact =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, solve.committed_state);
    const fuelsim::InterfaceSummary interface =
        fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, solve.committed_state);
    const RegionAverages cladding =
        region_averages(problem, fuelsim::rz::ProblemAccess::region_index(problem, "cladding"));
    std::size_t sliding = 0;
    std::size_t crossed_segments = 0;
    double maximum_coulomb_excess = 0.0;
    for (std::size_t node = 0; node < contact.size(); ++node) {
        if (contact[node].sliding) ++sliding;
        if (contact[node].primary_segment >= initial_contact[node].primary_segment + 2) ++crossed_segments;
        maximum_coulomb_excess = std::max(
            maximum_coulomb_excess, std::abs(contact[node].tangential_traction) -
                                        definition.spatial.contacts[0].friction_coefficient * contact[node].pressure);
    }
    double minimum_accepted_step = std::numeric_limits<double>::infinity();
    double maximum_accepted_step = 0.0;
    double maximum_time_error_estimate = 0.0;
    for (const fuelsim::TransientAcceptedStep& step : solve.accepted_steps) {
        minimum_accepted_step = std::min(minimum_accepted_step, step.time_step);
        maximum_accepted_step = std::max(maximum_accepted_step, step.time_step);
        maximum_time_error_estimate = std::max(maximum_time_error_estimate, step.time_error_estimate);
    }
    if (solve.accepted_steps.size() != 18)
        for (std::size_t index = 0; index < solve.accepted_steps.size(); ++index) {
            const fuelsim::TransientAcceptedStep& step = solve.accepted_steps[index];
            std::cerr << "m57_step[" << index << "]=" << std::setprecision(17) << step.time << ',' << step.time_step
                      << ',' << step.next_time_step << ',' << step.time_error_estimate << '\n';
        }
    const std::vector<fuelsim::test::NodalFieldReference> nodal_reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(problem, solve.committed_state, nodal_reference);
    std::vector<double> pressure_coordinates;
    const std::vector<double> pressure_reference =
        fuelsim::test::read_moose_contact_pressure_reference(pressure_reference_path, pressure_coordinates);
    const fuelsim::test::FieldErrorMetrics pressure =
        fuelsim::test::compare_moose_contact_pressure(contact, pressure_reference, pressure_coordinates, 1.0e-12);
    const std::vector<RegionAverages::Point> cladding_reference =
        read_cladding_points(qp_coordinate_path, qp_value_path);
    fuelsim::test::FieldErrorMetrics point_stress;
    fuelsim::test::FieldErrorMetrics point_plastic;
    fuelsim::test::FieldErrorMetrics point_creep;
    double maximum_point_coordinate_error = 0.0;
    if (cladding.points.size() == cladding_reference.size()) {
        for (std::size_t point = 0; point < cladding.points.size(); ++point) {
            maximum_point_coordinate_error = std::max({maximum_point_coordinate_error,
                std::abs(cladding.points[point].radius - cladding_reference[point].radius),
                std::abs(cladding.points[point].axial_coordinate - cladding_reference[point].axial_coordinate)});
            point_stress.add(cladding.points[point].stress, cladding_reference[point].stress);
            point_plastic.add(cladding.points[point].plastic, cladding_reference[point].plastic);
            point_creep.add(cladding.points[point].creep, cladding_reference[point].creep);
        }
    } else {
        maximum_point_coordinate_error = std::numeric_limits<double>::infinity();
    }
    constexpr double tolerance = 5.0e-3;
    bool passed = check(solve.aggregate_timing.workspace_setups == 1,
                      "M5.7 reuses one PETSc workspace across all adaptive steps") &&
                  check(solve.accepted_steps.size() == 18 && solve.time_error_rejections == 0 &&
                            maximum_time_error_estimate > 0.0 && maximum_time_error_estimate < 1.0 &&
                            minimum_accepted_step < maximum_accepted_step,
                      "M5.7 CTest path controls 18 accepted variable-size steps without rejection") &&
                  check(interface.active_contact_nodes == contact.size() && sliding >= contact.size() - 1 &&
                            crossed_segments == contact.size(),
                      "M5.7 keeps every interface node active, at least eight nodes sliding, and every node crossing "
                      "two segments") &&
                  check(std::abs(interface.total_heat_rate) > 0.0 &&
                            std::abs(problem.last_conservation_summary().interface_heat_imbalance) < 1.0e-10,
                      "M5.7 exercises nonzero conservative thermal contact") &&
                  check(maximum_coulomb_excess <= 1.0e-10 * interface.maximum_contact_pressure &&
                            std::abs(interface.total_tangential_force) > 0.0,
                      "M5.7 respects the Coulomb cap and produces nonzero friction force") &&
                  check(cladding.maximum_plastic > 0.0 && cladding.maximum_creep > 0.0,
                      "M5.7 activates cladding plasticity and creep") &&
                  check(fields.node_count == source.nodes().size() && fields.maximum_coordinate_difference < 1.0e-12,
                      "M5.7 compares every MOOSE source node at matching coordinates") &&
                  check(fuelsim::test::relative_metrics_below(fields.temperature, tolerance),
                      "M5.7 temperature field passes all three MOOSE metrics") &&
                  check(fuelsim::test::relative_metrics_below(fields.radial_displacement, tolerance),
                      "M5.7 radial field passes all three MOOSE metrics") &&
                  check(fuelsim::test::relative_metrics_below(fields.axial_displacement, tolerance),
                      "M5.7 axial field passes all three MOOSE metrics") &&
                  check(fuelsim::test::relative_metrics_below(pressure, tolerance),
                      "M5.7 node-to-segment pressure passes all three MOOSE metrics") &&
                  check(cladding.points.size() == cladding_reference.size() && maximum_point_coordinate_error < 1.0e-12,
                      "M5.7 compares every cladding integration point at matching coordinates") &&
                  check(fuelsim::test::relative_metrics_below(point_stress, tolerance),
                      "M5.7 integration-point stress passes all three MOOSE metrics") &&
                  check(fuelsim::test::relative_metrics_below(point_plastic, tolerance),
                      "M5.7 integration-point plastic strain passes all three MOOSE metrics") &&
                  check(fuelsim::test::relative_metrics_below(point_creep, tolerance),
                      "M5.7 integration-point creep strain passes all three MOOSE metrics");
    const double plastic_error =
        relative_error(cladding.plastic, final_csv_value(scalar_reference_path, "average_effective_plastic"));
    const double creep_error =
        relative_error(cladding.creep, final_csv_value(scalar_reference_path, "average_effective_creep"));
    const double stress_error =
        relative_error(cladding.stress, final_csv_value(scalar_reference_path, "average_vonmises_stress"));
    passed = check(plastic_error < tolerance && creep_error < tolerance && stress_error < tolerance,
                 "M5.7 cladding averages pass the MOOSE tolerance") &&
             passed;
    fuelsim::test::print_relative_metrics("m57_temperature", fields.temperature);
    fuelsim::test::print_relative_metrics("m57_radial_displacement", fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m57_axial_displacement", fields.axial_displacement);
    fuelsim::test::print_relative_metrics("m57_contact_pressure", pressure);
    fuelsim::test::print_relative_metrics("m57_qp_equivalent_stress", point_stress);
    fuelsim::test::print_relative_metrics("m57_qp_equivalent_plastic_strain", point_plastic);
    fuelsim::test::print_relative_metrics("m57_qp_equivalent_creep_strain", point_creep);
    std::cout << "m57_accepted_steps=" << solve.accepted_steps.size() << '\n'
              << "m57_time_error_rejections=" << solve.time_error_rejections << '\n'
              << "m57_maximum_time_error_estimate=" << maximum_time_error_estimate << '\n'
              << "m57_minimum_accepted_step=" << minimum_accepted_step << '\n'
              << "m57_maximum_accepted_step=" << maximum_accepted_step << '\n'
              << "m57_active_contact_nodes=" << interface.active_contact_nodes << '\n'
              << "m57_sliding_contact_nodes=" << sliding << '\n'
              << "m57_nodes_crossing_two_segments=" << crossed_segments << '\n'
              << "m57_total_heat_rate=" << interface.total_heat_rate << '\n'
              << "m57_total_tangential_force=" << interface.total_tangential_force << '\n'
              << "m57_average_plastic_relative_error=" << plastic_error << '\n'
              << "m57_average_creep_relative_error=" << creep_error << '\n'
              << "m57_average_stress_relative_error=" << stress_error << '\n';
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "Usage: fuelsim_m57_integrated_validation_tests <case.fsi> <all-nodes.csv> "
                     "<contact.csv> <qp-coordinates.csv> <qp-values.csv> <scalars.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim M5.7 integrated validation\n");
        if (!run_case(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6])) return 1;
        std::cout << "[PASS] M5.7 integrated single-rank verified MOOSE checks\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M5.7 validation raised: " << error.what() << '\n';
        return 1;
    }
}
