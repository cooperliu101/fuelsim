#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double m23_moose_tolerance = 1.0e-3;
constexpr double m41_moose_tolerance = 5.0e-3;

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ','))
        result.push_back(value);
    return result;
}

struct CsvTable final {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
};

CsvTable read_csv(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read MOOSE PCMI CSV: " + path);
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("MOOSE PCMI CSV is empty: " + path);
    CsvTable table;
    table.header = split_csv(line);
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        table.rows.push_back(split_csv(line));
    }
    return table;
}

std::size_t column_index(const CsvTable& table, const std::string& name) {
    const auto found =
        std::find(table.header.begin(), table.header.end(), name);
    if (found == table.header.end())
        throw std::invalid_argument("MOOSE PCMI CSV is missing column: " +
                                    name);
    return static_cast<std::size_t>(found - table.header.begin());
}

double csv_value(const CsvTable& table,
                 const std::vector<std::string>& row,
                 const std::string& name) {
    const std::size_t column = column_index(table, name);
    if (column >= row.size())
        throw std::invalid_argument("MOOSE PCMI CSV row is incomplete");
    return std::stod(row[column]);
}

double relative_error(double actual, double expected) {
    if (!std::isfinite(actual) || !std::isfinite(expected) || expected == 0.0)
        return std::numeric_limits<double>::infinity();
    return std::abs(actual - expected) / std::abs(expected);
}

struct ErrorMetrics final {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    double maximum_actual = 0.0;
    double maximum_reference = 0.0;
    double maximum_pointwise_relative = 0.0;

    void add(double actual, double reference) {
        const double difference = actual - reference;
        difference_squared += difference * difference;
        reference_squared += reference * reference;
        maximum_actual = std::max(maximum_actual, std::abs(actual));
        maximum_reference = std::max(maximum_reference, std::abs(reference));
        maximum_pointwise_relative = std::max(
            maximum_pointwise_relative, relative_error(actual, reference));
    }

    double relative_l2() const {
        return std::sqrt(difference_squared / reference_squared);
    }

    double relative_absolute_peak() const {
        return std::abs(maximum_actual - maximum_reference) / maximum_reference;
    }
};

bool check_metrics(const std::string& name, const ErrorMetrics& metrics,
                   double tolerance) {
    std::cout << name << "_relative_l2=" << metrics.relative_l2() << '\n';
    std::cout << name
              << "_relative_absolute_peak=" << metrics.relative_absolute_peak()
              << '\n';
    std::cout << name << "_maximum_pointwise_relative="
              << metrics.maximum_pointwise_relative << '\n';
    return check(metrics.relative_l2() < tolerance &&
                     metrics.relative_absolute_peak() < tolerance &&
                     metrics.maximum_pointwise_relative < tolerance,
                 name + " three MOOSE error metrics pass");
}

bool check_scalar_metrics(const std::string& name, double actual,
                          double reference, double tolerance) {
    ErrorMetrics metrics;
    metrics.add(actual, reference);
    return check_metrics(name, metrics, tolerance);
}

struct CladdingPointValue final {
    double radius;
    double axial_coordinate;
    double equivalent_stress;
    double equivalent_plastic_strain;
    double equivalent_creep_strain;
};

std::vector<CladdingPointValue>
read_cladding_point_reference(const std::string& coordinate_path,
                              const std::string& value_path) {
    const CsvTable coordinates = read_csv(coordinate_path);
    const CsvTable values = read_csv(value_path);
    if (values.rows.empty() || coordinates.rows.size() != values.rows.size() * 4)
        throw std::invalid_argument(
            "MOOSE PCMI QP coordinate/value counts are inconsistent");
    std::vector<CladdingPointValue> result;
    result.reserve(values.rows.size() * 4);
    constexpr std::array<std::size_t, 4> moose_qp_from_fuelsim = {0, 1, 3, 2};
    for (const std::vector<std::string>& value_row : values.rows) {
        const std::size_t element = static_cast<std::size_t>(
            csv_value(values, value_row, "id"));
        for (const std::size_t qp : moose_qp_from_fuelsim) {
            const auto coordinate = std::find_if(
                coordinates.rows.begin(), coordinates.rows.end(),
                [&coordinates, element, qp](
                    const std::vector<std::string>& row) {
                    return static_cast<std::size_t>(
                               csv_value(coordinates, row, "elem_id")) ==
                               element &&
                           static_cast<std::size_t>(
                               csv_value(coordinates, row, "qp_id")) == qp;
                });
            if (coordinate == coordinates.rows.end())
                throw std::invalid_argument(
                    "MOOSE PCMI QP coordinate is missing");
            const std::string suffix = std::to_string(qp);
            result.push_back({
                csv_value(coordinates, *coordinate, "x"),
                csv_value(coordinates, *coordinate, "y"),
                csv_value(values, value_row, "stress_q" + suffix),
                csv_value(values, value_row, "plastic_q" + suffix),
                csv_value(values, value_row, "creep_q" + suffix),
            });
        }
    }
    return result;
}

double read_final_scalar_reference(const std::string& path,
                                   const std::string& name) {
    const CsvTable table = read_csv(path);
    if (table.rows.empty())
        throw std::invalid_argument("MOOSE PCMI scalar CSV has no rows");
    return csv_value(table, table.rows.back(), name);
}

double read_total_contact_force_reference(const std::string& path) {
    const CsvTable table = read_csv(path);
    double force = 0.0;
    for (const std::vector<std::string>& row : table.rows) {
        force += csv_value(table, row, "contact_pressure") *
                 csv_value(table, row, "nodal_area");
    }
    return force;
}

struct CladdingMetrics final {
    double average_plastic;
    double average_creep;
    double average_equivalent_stress;
    double maximum_plastic;
    double maximum_creep;
    std::vector<CladdingPointValue> points;
};

CladdingMetrics cladding_metrics(const fuelsim::TransientProblem& problem,
                                 std::size_t cladding_region) {
    double measure = 0.0;
    double weighted_plastic = 0.0;
    double weighted_creep = 0.0;
    double weighted_equivalent_stress = 0.0;
    double maximum_plastic = 0.0;
    double maximum_creep = 0.0;
    const fuelsim::RegionMesh& mesh = problem.region_mesh(cladding_region);
    std::vector<CladdingPointValue> points;
    points.reserve(mesh.elements().size() * 4);
    for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
        const fuelsim::Quad4RzGeometry& geometry =
            problem.region_element_geometry(cladding_region, element);
        const fuelsim::Quad4Element& mesh_element = mesh.elements()[element];
        const fuelsim::Quad4MaterialHistory& history =
            problem.material_history(cladding_region, element);
        const std::array<fuelsim::AxisymmetricStressValues, 4>& stresses =
            problem.material_stress(cladding_region, element);
        for (std::size_t q = 0; q < history.size(); ++q) {
            const fuelsim::RzQuadraturePoint& point = geometry.points[q];
            double axial_coordinate = 0.0;
            for (std::size_t node = 0; node < 4; ++node) {
                axial_coordinate += point.shape[node] *
                                    mesh.nodes()[mesh_element.nodes[node]].z;
            }
            const fuelsim::AxisymmetricStressValues& stress = stresses[q];
            const double mean_stress =
                (stress.rr + stress.zz + stress.hoop) / 3.0;
            const double equivalent_stress = std::sqrt(
                1.5 *
                ((stress.rr - mean_stress) * (stress.rr - mean_stress) +
                 (stress.zz - mean_stress) * (stress.zz - mean_stress) +
                 (stress.hoop - mean_stress) * (stress.hoop - mean_stress) +
                 2.0 * stress.rz * stress.rz));
            const double weight = point.weighted_measure;
            measure += weight;
            weighted_plastic += weight * history[q].equivalent_plastic_strain;
            weighted_creep += weight * history[q].equivalent_creep_strain;
            weighted_equivalent_stress += weight * equivalent_stress;
            maximum_plastic =
                std::max(maximum_plastic, history[q].equivalent_plastic_strain);
            maximum_creep =
                std::max(maximum_creep, history[q].equivalent_creep_strain);
            points.push_back({point.radius, axial_coordinate, equivalent_stress,
                              history[q].equivalent_plastic_strain,
                              history[q].equivalent_creep_strain});
        }
    }
    if (!(measure > 0.0))
        throw std::runtime_error(
            "PCMI cladding history requires positive integration measure");
    return {
        weighted_plastic / measure,
        weighted_creep / measure,
        weighted_equivalent_stress / measure,
        maximum_plastic,
        maximum_creep,
        std::move(points),
    };
}

struct PointwiseError final {
    double relative_l2 = 0.0;
    double relative_absolute_peak = 0.0;
    double maximum_absolute = 0.0;
    double maximum_relative = 0.0;
    double maximum_actual = 0.0;
    double maximum_reference = 0.0;
    std::size_t maximum_absolute_point = 0;
    std::size_t maximum_relative_point = 0;
};

void add_pointwise_error(PointwiseError& error, double actual, double expected,
                         std::size_t point, double& difference_squared,
                         double& reference_squared) {
    const double absolute = std::abs(actual - expected);
    const double relative = relative_error(actual, expected);
    difference_squared += absolute * absolute;
    reference_squared += expected * expected;
    error.maximum_actual = std::max(error.maximum_actual, std::abs(actual));
    error.maximum_reference =
        std::max(error.maximum_reference, std::abs(expected));
    if (absolute > error.maximum_absolute) {
        error.maximum_absolute = absolute;
        error.maximum_absolute_point = point;
    }
    if (relative > error.maximum_relative) {
        error.maximum_relative = relative;
        error.maximum_relative_point = point;
    }
}

void finalize_pointwise_error(PointwiseError& error, double difference_squared,
                              double reference_squared) {
    error.relative_l2 = reference_squared > 0.0
                            ? std::sqrt(difference_squared / reference_squared)
                            : std::numeric_limits<double>::infinity();
    error.relative_absolute_peak =
        error.maximum_reference > 0.0
            ? std::abs(error.maximum_actual - error.maximum_reference) /
                  error.maximum_reference
            : std::numeric_limits<double>::infinity();
}

bool fuel_history_is_elastic(const fuelsim::TransientProblem& problem,
                             std::size_t fuel_region) {
    for (std::size_t element = 0;
         element < problem.region_mesh(fuel_region).elements().size();
         ++element) {
        for (const fuelsim::MaterialPointState& point :
             problem.material_history(fuel_region, element)) {
            if (point.equivalent_plastic_strain != 0.0 ||
                point.equivalent_creep_strain != 0.0)
                return false;
        }
    }
    return true;
}

fuelsim::SolverOptions
solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions options{
        definition.solver.absolute_tolerance,
        definition.solver.relative_tolerance,
        definition.solver.step_tolerance,
        definition.solver.maximum_iterations};
    // The tracked MOOSE snapshot uses the fixed 20 x 1 s load path. Keep that
    // comparison path independent of the production default line search,
    // whose contact-onset cutback deliberately changes the time grid.
    options.line_search = fuelsim::SolverOptions::LineSearch::basic;
    options.residual_reduction_tolerance =
        definition.solver.residual_reduction_tolerance;
    options.field_residual_scaling =
        definition.solver.field_residual_scaling;
    options.temperature_residual_absolute_tolerance =
        definition.solver.temperature_residual_absolute_tolerance;
    options.mechanical_residual_absolute_tolerance =
        definition.solver.mechanical_residual_absolute_tolerance;
    return options;
}

fuelsim::TransientTimeOptions
time_options(const fuelsim::FuelSimCaseDefinition& definition) {
    return {definition.transient_execution.end_time,
            definition.transient_execution.initial_time_step,
            definition.transient_execution.minimum_time_step,
            definition.transient_execution.maximum_time_step,
            definition.transient_execution.growth_factor,
            definition.transient_execution.cutback_factor,
            definition.transient_execution.maximum_cutbacks,
            definition.transient_execution.load_ramp_time};
}

bool test_pcmi_coupled_cladding(const std::string& input_path,
                                const std::string& nodal_reference_path,
                                const std::string& pressure_reference_path,
                                const std::string& qp_coordinate_path,
                                const std::string& qp_value_path,
                                const std::string& scalar_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::transient)
        throw std::invalid_argument(
            "PCMI comparison requires a transient input card");
    bool finite_strain = false;
    for (const fuelsim::CaseRegionDefinition& region : definition.regions)
        finite_strain =
            finite_strain || region.spatial.strain_formulation ==
                                 fuelsim::StrainFormulation::finite;
    const double moose_tolerance =
        finite_strain ? m41_moose_tolerance : m23_moose_tolerance;
    const fuelsim::UnstructuredQuad4Mesh imported =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.transient_definition(),
                                      imported);
    const fuelsim::TransientResult result = fuelsim::solve_transient(
        problem, time_options(definition), solver_options(definition));
    if (!result.completed)
        return check(false, "PCMI transient completes all twenty time steps");

    const std::size_t fuel_region = problem.region_index("fuel");
    const std::size_t cladding_region = problem.region_index("cladding");
    const std::vector<double>& state = result.committed_state;
    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(0, state);
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        problem.summarize_contact_nodes(0, state);
    const CladdingMetrics cladding = cladding_metrics(problem, cladding_region);

    const double expected_average_plastic = read_final_scalar_reference(
        scalar_reference_path, "average_effective_plastic");
    const double expected_average_creep = read_final_scalar_reference(
        scalar_reference_path, "average_effective_creep");
    const double expected_average_equivalent_stress =
        read_final_scalar_reference(scalar_reference_path,
                                    "average_vonmises_stress");
    const double expected_total_contact_force =
        read_total_contact_force_reference(pressure_reference_path);
    const std::vector<CladdingPointValue> expected_cladding_points =
        read_cladding_point_reference(qp_coordinate_path, qp_value_path);
    const double average_plastic_error =
        relative_error(cladding.average_plastic, expected_average_plastic);
    const double average_creep_error =
        relative_error(cladding.average_creep, expected_average_creep);
    const double average_equivalent_stress_error = relative_error(
        cladding.average_equivalent_stress, expected_average_equivalent_stress);
    const double total_contact_force_error = relative_error(
        interface.total_contact_force, expected_total_contact_force);
    const std::vector<fuelsim::test::NodalFieldReference> nodal_reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison full_fields =
        fuelsim::test::compare_moose_nodal_fields(
            problem, result.committed_state, nodal_reference);
    std::vector<double> pressure_coordinates;
    const std::vector<double> pressure_reference =
        fuelsim::test::read_moose_contact_pressure_reference(
            pressure_reference_path, pressure_coordinates);
    const fuelsim::test::FieldErrorMetrics pressure_metrics =
        fuelsim::test::compare_moose_contact_pressure(
            contact_nodes, pressure_reference, pressure_coordinates, 1.0e-12);
    PointwiseError point_stress_error;
    PointwiseError point_plastic_error;
    PointwiseError point_creep_error;
    double point_stress_difference_squared = 0.0;
    double point_stress_reference_squared = 0.0;
    double point_plastic_difference_squared = 0.0;
    double point_plastic_reference_squared = 0.0;
    double point_creep_difference_squared = 0.0;
    double point_creep_reference_squared = 0.0;
    double maximum_point_location_error = 0.0;
    if (cladding.points.size() == expected_cladding_points.size()) {
        for (std::size_t point = 0; point < cladding.points.size(); ++point) {
            const CladdingPointValue& actual = cladding.points[point];
            const CladdingPointValue& expected =
                expected_cladding_points[point];
            maximum_point_location_error =
                std::max({maximum_point_location_error,
                          std::abs(actual.radius - expected.radius),
                          std::abs(actual.axial_coordinate -
                                   expected.axial_coordinate)});
            add_pointwise_error(point_stress_error, actual.equivalent_stress,
                                expected.equivalent_stress, point,
                                point_stress_difference_squared,
                                point_stress_reference_squared);
            add_pointwise_error(point_plastic_error,
                                actual.equivalent_plastic_strain,
                                expected.equivalent_plastic_strain, point,
                                point_plastic_difference_squared,
                                point_plastic_reference_squared);
            add_pointwise_error(
                point_creep_error, actual.equivalent_creep_strain,
                expected.equivalent_creep_strain, point,
                point_creep_difference_squared, point_creep_reference_squared);
        }
        finalize_pointwise_error(point_stress_error,
                                 point_stress_difference_squared,
                                 point_stress_reference_squared);
        finalize_pointwise_error(point_plastic_error,
                                 point_plastic_difference_squared,
                                 point_plastic_reference_squared);
        finalize_pointwise_error(point_creep_error,
                                 point_creep_difference_squared,
                                 point_creep_reference_squared);
    } else {
        point_stress_error.relative_l2 =
            std::numeric_limits<double>::infinity();
        point_plastic_error.relative_l2 =
            std::numeric_limits<double>::infinity();
        point_creep_error.relative_l2 = std::numeric_limits<double>::infinity();
        maximum_point_location_error = std::numeric_limits<double>::infinity();
    }

    bool passed =
        check(imported.nodes().size() == 53 &&
                  imported.elements().size() == 34 &&
                  imported.side_set("fuel_right").sides.size() == 4 &&
                  imported.side_set("clad_left").sides.size() == 5,
              "PCMI reads the tracked four-to-five nonmatching interface") &&
        check(result.accepted_steps.size() == 20,
              "PCMI transient commits twenty accepted steps");
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "PCMI transient reuses one PETSc workspace") &&
             passed;
    passed = check(result.total_cutbacks == 0,
                   "PCMI reference path converges without cutbacks") &&
             passed;
    passed = check(fuel_history_is_elastic(problem, fuel_region),
                   "PCMI fuel remains on the elastic material path") &&
             passed;
    passed =
        check(interface.projected_contact_nodes == contact_nodes.size(),
              "PCMI taller cladding contains every fuel-node projection") &&
        passed;
    passed = check(interface.active_contact_nodes > 0 &&
                       interface.minimum_contact_gap < 0.0 &&
                       interface.maximum_contact_pressure > 0.0,
                   "fuel thermal expansion closes the gap and develops "
                   "contact pressure") &&
             passed;
    passed =
        check(cladding.average_plastic > 0.0 && cladding.average_creep > 0.0 &&
                  cladding.maximum_plastic > 0.0 &&
                  cladding.maximum_creep > 0.0,
              "PCMI cladding accumulates plastic and creep history") &&
        passed;
    passed = check(full_fields.node_count == imported.nodes().size() &&
                       full_fields.maximum_coordinate_difference < 1.0e-12,
                   "PCMI compares every MOOSE node at matching coordinates") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(
                       full_fields.temperature, moose_tolerance),
                   "PCMI full-field temperature three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(
                       full_fields.radial_displacement, moose_tolerance),
                   "PCMI full-field radial displacement three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(
                       full_fields.axial_displacement,
                       moose_tolerance),
                   "PCMI full-field axial displacement three errors pass") &&
             passed;
    fuelsim::test::print_relative_metrics("pcmi_temperature",
                                          full_fields.temperature);
    fuelsim::test::print_relative_metrics("pcmi_radial_displacement",
                                          full_fields.radial_displacement);
    fuelsim::test::print_relative_metrics("pcmi_axial_displacement",
                                          full_fields.axial_displacement);
    passed = check_scalar_metrics("pcmi_average_plastic_strain",
                                  cladding.average_plastic,
                                  expected_average_plastic,
                                  moose_tolerance) &&
             passed;
    passed = check_scalar_metrics("pcmi_average_creep_strain",
                                  cladding.average_creep,
                                  expected_average_creep,
                                  moose_tolerance) &&
             passed;
    passed = check_scalar_metrics("pcmi_average_equivalent_stress",
                                  cladding.average_equivalent_stress,
                                  expected_average_equivalent_stress,
                                  moose_tolerance) &&
             passed;
    passed = check(cladding.points.size() == expected_cladding_points.size(),
                   "PCMI and MOOSE cladding integration-point counts match") &&
             passed;
    passed =
        check(maximum_point_location_error < 1.0e-12,
              "PCMI and MOOSE cladding integration-point locations match") &&
        passed;
    passed = check(point_stress_error.relative_l2 < moose_tolerance &&
                       point_stress_error.relative_absolute_peak <
                           moose_tolerance &&
                       point_stress_error.maximum_relative < moose_tolerance,
                   "PCMI pointwise equivalent-stress three errors pass") &&
             passed;
    passed = check(point_plastic_error.relative_l2 < moose_tolerance &&
                       point_plastic_error.relative_absolute_peak <
                           moose_tolerance &&
                       point_plastic_error.maximum_relative < moose_tolerance,
                   "PCMI pointwise plastic-strain three errors pass") &&
             passed;
    passed = check(point_creep_error.relative_l2 < moose_tolerance &&
                       point_creep_error.relative_absolute_peak <
                           moose_tolerance &&
                       point_creep_error.maximum_relative < moose_tolerance,
                   "PCMI pointwise creep-strain three errors pass") &&
             passed;
    passed =
        check(fuelsim::test::relative_metrics_below(
                  pressure_metrics, moose_tolerance),
              "PCMI full-field contact pressure three errors pass") &&
        passed;
    fuelsim::test::print_relative_metrics("pcmi_contact_pressure",
                                          pressure_metrics);
    passed = check_scalar_metrics("pcmi_total_contact_force",
                                  interface.total_contact_force,
                                  expected_total_contact_force,
                                  moose_tolerance) &&
             passed;

    std::cout << "pcmi_minimum_contact_gap=" << interface.minimum_contact_gap
              << '\n';
    std::cout << "pcmi_maximum_contact_pressure="
              << interface.maximum_contact_pressure << '\n';
    std::cout << "pcmi_total_contact_force=" << interface.total_contact_force
              << '\n';
    std::cout << "pcmi_projected_contact_nodes="
              << interface.projected_contact_nodes << '\n';
    std::cout << "pcmi_active_contact_nodes=" << interface.active_contact_nodes
              << '\n';
    std::cout << "pcmi_cladding_average_equivalent_plastic_strain="
              << cladding.average_plastic << '\n';
    std::cout << "pcmi_cladding_average_equivalent_creep_strain="
              << cladding.average_creep << '\n';
    std::cout << "pcmi_cladding_average_equivalent_stress="
              << cladding.average_equivalent_stress << '\n';
    std::cout << "pcmi_cladding_maximum_equivalent_plastic_strain="
              << cladding.maximum_plastic << '\n';
    std::cout << "pcmi_cladding_maximum_equivalent_creep_strain="
              << cladding.maximum_creep << '\n';
    for (std::size_t node = 0; node < contact_nodes.size(); ++node) {
        std::cout << "pcmi_contact_pressure_" << node << '='
                  << contact_nodes[node].pressure << '\n';
    }
    std::cout << "pcmi_moose_average_plastic_relative_error="
              << average_plastic_error << '\n';
    std::cout << "pcmi_moose_average_creep_relative_error="
              << average_creep_error << '\n';
    std::cout << "pcmi_moose_average_equivalent_stress_relative_error="
              << average_equivalent_stress_error << '\n';
    std::cout << "pcmi_moose_qp_maximum_location_absolute_error="
              << maximum_point_location_error << '\n';
    std::cout << "pcmi_moose_qp_equivalent_stress_relative_l2="
              << point_stress_error.relative_l2 << '\n';
    std::cout << "pcmi_moose_qp_equivalent_stress_relative_absolute_peak="
              << point_stress_error.relative_absolute_peak << '\n';
    std::cout << "pcmi_moose_qp_equivalent_stress_maximum_absolute_error="
              << point_stress_error.maximum_absolute << '\n';
    std::cout
        << "pcmi_moose_qp_equivalent_stress_maximum_absolute_local_element="
        << point_stress_error.maximum_absolute_point / 4 << '\n';
    std::cout << "pcmi_moose_qp_equivalent_stress_maximum_absolute_qp="
              << point_stress_error.maximum_absolute_point % 4 << '\n';
    std::cout << "pcmi_moose_qp_equivalent_stress_maximum_relative_error="
              << point_stress_error.maximum_relative << '\n';
    std::cout
        << "pcmi_moose_qp_equivalent_stress_maximum_relative_local_element="
        << point_stress_error.maximum_relative_point / 4 << '\n';
    std::cout << "pcmi_moose_qp_equivalent_stress_maximum_relative_qp="
              << point_stress_error.maximum_relative_point % 4 << '\n';
    std::cout << "pcmi_moose_qp_plastic_strain_relative_l2="
              << point_plastic_error.relative_l2 << '\n';
    std::cout << "pcmi_moose_qp_plastic_strain_relative_absolute_peak="
              << point_plastic_error.relative_absolute_peak << '\n';
    std::cout << "pcmi_moose_qp_plastic_strain_maximum_absolute_error="
              << point_plastic_error.maximum_absolute << '\n';
    std::cout << "pcmi_moose_qp_plastic_strain_maximum_absolute_local_element="
              << point_plastic_error.maximum_absolute_point / 4 << '\n';
    std::cout << "pcmi_moose_qp_plastic_strain_maximum_absolute_qp="
              << point_plastic_error.maximum_absolute_point % 4 << '\n';
    std::cout << "pcmi_moose_qp_plastic_strain_maximum_relative_error="
              << point_plastic_error.maximum_relative << '\n';
    std::cout << "pcmi_moose_qp_plastic_strain_maximum_relative_local_element="
              << point_plastic_error.maximum_relative_point / 4 << '\n';
    std::cout << "pcmi_moose_qp_plastic_strain_maximum_relative_qp="
              << point_plastic_error.maximum_relative_point % 4 << '\n';
    std::cout << "pcmi_moose_qp_creep_strain_relative_l2="
              << point_creep_error.relative_l2 << '\n';
    std::cout << "pcmi_moose_qp_creep_strain_relative_absolute_peak="
              << point_creep_error.relative_absolute_peak << '\n';
    std::cout << "pcmi_moose_qp_creep_strain_maximum_absolute_error="
              << point_creep_error.maximum_absolute << '\n';
    std::cout << "pcmi_moose_qp_creep_strain_maximum_absolute_local_element="
              << point_creep_error.maximum_absolute_point / 4 << '\n';
    std::cout << "pcmi_moose_qp_creep_strain_maximum_absolute_qp="
              << point_creep_error.maximum_absolute_point % 4 << '\n';
    std::cout << "pcmi_moose_qp_creep_strain_maximum_relative_error="
              << point_creep_error.maximum_relative << '\n';
    std::cout << "pcmi_moose_qp_creep_strain_maximum_relative_local_element="
              << point_creep_error.maximum_relative_point / 4 << '\n';
    std::cout << "pcmi_moose_qp_creep_strain_maximum_relative_qp="
              << point_creep_error.maximum_relative_point % 4 << '\n';
    std::cout << "pcmi_moose_total_contact_force_relative_error="
              << total_contact_force_error << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "Usage: fuelsim_m2_pcmi_solver_tests <pcmi.fsi> "
                     "<all-nodes.csv> <contact-pressure.csv> "
                     "<qp-coordinates.csv> <qp-values.csv> <scalars.csv>\n";
        return 2;
    }

    try {
        const std::string input_path = argv[1];
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv,
            "fuelsim M2 PCMI coupled cladding MOOSE comparison tests\n");
        if (!test_pcmi_coupled_cladding(input_path, argv[2], argv[3], argv[4],
                                        argv[5], argv[6]))
            return 1;
        std::cout << "[PASS] fuelsim M2 PCMI coupled cladding tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] PCMI tests raised: " << error.what() << '\n';
        return 1;
    }
}
