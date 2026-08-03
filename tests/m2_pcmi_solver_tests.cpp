#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double nodal_standard_tolerance = 1.0e-3;
constexpr double radial_peak_tolerance = 1.5e-3;
constexpr double radial_pointwise_tolerance = 3.0e-3;
constexpr double contact_pressure_tolerance = 1.2e-3;
constexpr double stress_l2_tolerance = 2.0e-4;
constexpr double stress_peak_tolerance = 5.0e-4;
constexpr double plastic_l2_tolerance = 1.6e-3;
constexpr double plastic_peak_tolerance = 3.0e-3;
constexpr double plastic_pointwise_tolerance = 3.2e-3;
constexpr double creep_l2_tolerance = 1.0e-3;
constexpr double creep_peak_tolerance = 2.0e-3;

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
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

bool same_coordinate(double lhs, double rhs) {
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= 1.0e-12 * scale;
}

std::size_t find_node(const fuelsim::RegionMesh& mesh, double radius,
                      double axial_coordinate) {
    for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
        if (same_coordinate(mesh.nodes()[node].r, radius) &&
            same_coordinate(mesh.nodes()[node].z, axial_coordinate))
            return node;
    }
    throw std::invalid_argument("PCMI comparison point is not in the mesh");
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
                                const std::string& pressure_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::transient)
        throw std::invalid_argument(
            "PCMI comparison requires a transient input card");
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
    const fuelsim::RegionMesh& fuel_mesh = problem.region_mesh(fuel_region);
    const fuelsim::RegionMesh& cladding_mesh =
        problem.region_mesh(cladding_region);
    const std::size_t fuel_center = find_node(fuel_mesh, 0.0, 0.005);
    const std::size_t fuel_surface = find_node(fuel_mesh, 0.004120, 0.005);
    const std::size_t cladding_inner =
        find_node(cladding_mesh, 0.004121, 0.00501);
    const std::size_t cladding_outer =
        find_node(cladding_mesh, 0.004692, 0.00501);
    const std::size_t fuel_axis_top = find_node(fuel_mesh, 0.0, 0.010);
    const std::size_t fuel_offset = problem.region_node_offset(fuel_region);
    const std::size_t cladding_offset =
        problem.region_node_offset(cladding_region);
    const fuelsim::DofMap& dofs = problem.dof_map();

    const std::vector<double>& state = result.committed_state;
    const double temperature_center =
        state[dofs.temperature(fuel_offset + fuel_center)];
    const double temperature_fuel_surface =
        state[dofs.temperature(fuel_offset + fuel_surface)];
    const double temperature_cladding_inner =
        state[dofs.temperature(cladding_offset + cladding_inner)];
    const double radial_fuel_surface =
        state[dofs.radial_displacement(fuel_offset + fuel_surface)];
    const double radial_cladding_inner =
        state[dofs.radial_displacement(cladding_offset + cladding_inner)];
    const double radial_cladding_outer =
        state[dofs.radial_displacement(cladding_offset + cladding_outer)];
    const double axial_fuel_top =
        state[dofs.axial_displacement(fuel_offset + fuel_axis_top)];
    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(0, state);
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        problem.summarize_contact_nodes(0, state);
    const CladdingMetrics cladding = cladding_metrics(problem, cladding_region);

    // Generated by verification/moose/m23_pcmi_coupled_cladding_rz.i. The
    // fuel is elastic, while the cladding uses coupled Norton creep and J2
    // plasticity in the same material update.
    constexpr double expected_temperature_center = 737.32437038895;
    constexpr double expected_temperature_fuel_surface = 613.74977861649;
    constexpr double expected_temperature_cladding_inner = 612.79216505292;
    constexpr double expected_radial_fuel_surface = 2.9177714961082e-6;
    constexpr double expected_radial_cladding_inner = 1.9098157198910e-6;
    constexpr double expected_radial_cladding_outer = 1.7860488488808e-6;
    constexpr double expected_axial_fuel_top = 9.2609839250304e-6;
    constexpr double expected_average_plastic = 2.4156288403990e-4;
    constexpr double expected_average_creep = 1.2769576281297e-4;
    constexpr double expected_average_equivalent_stress = 5483125.7680798;
    constexpr double expected_total_contact_force = 191.06119157558877;
    // MOOSE QP ordering is [--, +-, -+, ++]. References below are reordered
    // to fuelsim's [--, +-, ++, -+] convention within each element.
    constexpr std::array<CladdingPointValue, 32> expected_cladding_points = {{
        {4.1813332490732002e-3, 5.2936878783998997e-4, 5.5379798327145996e6,
         2.6898991635728998e-4, 1.3147469050565000e-4},
        {4.3461667509268002e-3, 5.2936878783998997e-4, 5.5059675508802002e6,
         2.5298377544010002e-4, 1.2921153876656999e-4},
        {4.3461667509268002e-3, 1.9756312121599999e-3, 5.5033337755910000e6,
         2.5166688779551999e-4, 1.2869950224100001e-4},
        {4.1813332490732002e-3, 1.9756312121599999e-3, 5.5348683751563001e6,
         2.6743418757814002e-4, 1.3111171361262000e-4},
        {4.4668332490732003e-3, 5.2936878783998997e-4, 5.4566641267114002e6,
         2.2833206335568000e-4, 1.2548489788044001e-4},
        {4.6316667509268003e-3, 5.2936878783998997e-4, 5.4302810392854000e6,
         2.1514051964269999e-4, 1.2355554119910999e-4},
        {4.6316667509268003e-3, 1.9756312121599999e-3, 5.4298935069199996e6,
         2.1494675345997999e-4, 1.2341438586164000e-4},
        {4.4668332490732003e-3, 1.9756312121599999e-3, 5.4563087917584004e6,
         2.2815439587921000e-4, 1.2533510884110000e-4},
        {4.1813332490732002e-3, 3.0343687878399998e-3, 5.5275841168417996e6,
         2.6379205842088998e-4, 1.3002439365296999e-4},
        {4.3461667509268002e-3, 3.0343687878399998e-3, 5.4955834560783003e6,
         2.4779172803918002e-4, 1.2752955822262999e-4},
        {4.3461667509268002e-3, 4.4806312121600002e-3, 5.4822938984810999e6,
         2.4114694924055999e-4, 1.2581887569787999e-4},
        {4.1813332490732002e-3, 4.4806312121600002e-3, 5.5137954443515996e6,
         2.5689772217579000e-4, 1.2818256306923001e-4},
        {4.4668332490732003e-3, 3.0343687878399998e-3, 5.4476035984185003e6,
         2.2380179920927001e-4, 1.2403420161417999e-4},
        {4.6316667509268003e-3, 3.0343687878399998e-3, 5.4206822509414004e6,
         2.1034112547070000e-4, 1.2205356963625000e-4},
        {4.6316667509268003e-3, 4.4806312121600002e-3, 5.4051499993891995e6,
         2.0257499969460999e-4, 1.2013173565021000e-4},
        {4.4668332490732003e-3, 4.4806312121600002e-3, 5.4314295954711000e6,
         2.1571479773552999e-4, 1.2195859726964000e-4},
        {4.1813332490732002e-3, 5.5393687878400001e-3, 5.4938910727728996e6,
         2.4694553638645000e-4, 1.2500539810722999e-4},
        {4.3461667509268002e-3, 5.5393687878400001e-3, 5.4688791733983001e6,
         2.3443958669912000e-4, 1.2366093035905999e-4},
        {4.3461667509268002e-3, 6.9856312121599996e-3, 5.4564076657935996e6,
         2.2820383289681999e-4, 1.2351657712684000e-4},
        {4.1813332490732002e-3, 6.9856312121599996e-3, 5.4803377481268002e6,
         2.4016887406340000e-4, 1.2469834667407000e-4},
        {4.4668332490732003e-3, 5.5393687878400001e-3, 5.4206495705773998e6,
         2.1032478528872000e-4, 1.2052486163370000e-4},
        {4.6316667509268003e-3, 5.5393687878400001e-3, 5.4018147255169004e6,
         2.0090736275846000e-4, 1.1980113628869000e-4},
        {4.6316667509268003e-3, 6.9856312121599996e-3, 5.3878462224150999e6,
         1.9392311120754001e-4, 1.2005385380186001e-4},
        {4.4668332490732003e-3, 6.9856312121599996e-3, 5.4053743719853004e6,
         2.0268718599262999e-4, 1.2063673105502000e-4},
        {4.1813332490732002e-3, 8.0443687878399995e-3, 5.5354679527030997e6,
         2.6773397635156999e-4, 1.3245590693081001e-4},
        {4.3461667509268002e-3, 8.0443687878399995e-3, 5.5086660131898001e6,
         2.5433300659488998e-4, 1.3121843774056000e-4},
        {4.3461667509268002e-3, 9.4906312121600007e-3, 5.6562687041458003e6,
         3.2813435207290002e-4, 1.4632416039005001e-4},
        {4.1813332490732002e-3, 9.4906312121600007e-3, 5.6897500752656003e6,
         3.4487503763279003e-4, 1.4842972583189001e-4},
        {4.4668332490732003e-3, 8.0443687878399995e-3, 5.4533520390095003e6,
         2.2667601950473001e-4, 1.2812863288300000e-4},
        {4.6316667509268003e-3, 8.0443687878399995e-3, 5.4349630495820995e6,
         2.1748152479104000e-4, 1.2746111025778000e-4},
        {4.6316667509268003e-3, 9.4906312121600007e-3, 5.5556306027846001e6,
         2.7781530139228000e-4, 1.3918892029350999e-4},
        {4.4668332490732003e-3, 9.4906312121600007e-3, 5.5804350404331004e6,
         2.9021752021653001e-4, 1.4038988537948000e-4},
    }};

    const double temperature_center_error =
        relative_error(temperature_center, expected_temperature_center);
    const double temperature_fuel_surface_error = relative_error(
        temperature_fuel_surface, expected_temperature_fuel_surface);
    const double temperature_cladding_inner_error = relative_error(
        temperature_cladding_inner, expected_temperature_cladding_inner);
    const double radial_fuel_surface_error =
        relative_error(radial_fuel_surface, expected_radial_fuel_surface);
    const double radial_cladding_inner_error =
        relative_error(radial_cladding_inner, expected_radial_cladding_inner);
    const double radial_cladding_outer_error =
        relative_error(radial_cladding_outer, expected_radial_cladding_outer);
    const double axial_fuel_top_error =
        relative_error(axial_fuel_top, expected_axial_fuel_top);
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

    bool passed = check(result.accepted_steps.size() == 20,
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
                       full_fields.temperature, nodal_standard_tolerance),
                   "PCMI full-field temperature three errors pass") &&
             passed;
    passed =
        check(full_fields.radial_displacement.relative_l2() <
                  nodal_standard_tolerance &&
                  full_fields.radial_displacement.relative_absolute_peak() <
                      radial_peak_tolerance &&
                  full_fields.radial_displacement
                          .maximum_pointwise_relative_error() <
                      radial_pointwise_tolerance,
              "PCMI full-field radial displacement three errors pass") &&
        passed;
    passed = check(fuelsim::test::relative_metrics_below(
                       full_fields.axial_displacement,
                       nodal_standard_tolerance),
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
                                  nodal_standard_tolerance) &&
             passed;
    passed = check_scalar_metrics("pcmi_average_creep_strain",
                                  cladding.average_creep,
                                  expected_average_creep,
                                  nodal_standard_tolerance) &&
             passed;
    passed = check_scalar_metrics("pcmi_average_equivalent_stress",
                                  cladding.average_equivalent_stress,
                                  expected_average_equivalent_stress,
                                  nodal_standard_tolerance) &&
             passed;
    passed = check(cladding.points.size() == expected_cladding_points.size(),
                   "PCMI and MOOSE cladding integration-point counts match") &&
             passed;
    passed =
        check(maximum_point_location_error < 1.0e-12,
              "PCMI and MOOSE cladding integration-point locations match") &&
        passed;
    passed = check(point_stress_error.relative_l2 < stress_l2_tolerance &&
                       point_stress_error.relative_absolute_peak <
                           stress_peak_tolerance &&
                       point_stress_error.maximum_relative <
                           stress_peak_tolerance,
                   "PCMI pointwise equivalent-stress three errors pass") &&
             passed;
    passed = check(point_plastic_error.relative_l2 < plastic_l2_tolerance &&
                       point_plastic_error.relative_absolute_peak <
                           plastic_peak_tolerance &&
                       point_plastic_error.maximum_relative <
                           plastic_pointwise_tolerance,
                   "PCMI pointwise plastic-strain three errors pass") &&
             passed;
    passed = check(point_creep_error.relative_l2 < creep_l2_tolerance &&
                       point_creep_error.relative_absolute_peak <
                           creep_peak_tolerance &&
                       point_creep_error.maximum_relative <
                           creep_peak_tolerance,
                   "PCMI pointwise creep-strain three errors pass") &&
             passed;
    passed =
        check(fuelsim::test::relative_metrics_below(
                  pressure_metrics, contact_pressure_tolerance),
              "PCMI full-field contact pressure three errors pass") &&
        passed;
    fuelsim::test::print_relative_metrics("pcmi_contact_pressure",
                                          pressure_metrics);
    passed = check_scalar_metrics("pcmi_total_contact_force",
                                  interface.total_contact_force,
                                  expected_total_contact_force,
                                  nodal_standard_tolerance) &&
             passed;

    std::cout << "pcmi_temperature_center=" << temperature_center << '\n';
    std::cout << "pcmi_temperature_fuel_surface=" << temperature_fuel_surface
              << '\n';
    std::cout << "pcmi_temperature_cladding_inner="
              << temperature_cladding_inner << '\n';
    std::cout << "pcmi_radial_fuel_surface=" << radial_fuel_surface << '\n';
    std::cout << "pcmi_radial_cladding_inner=" << radial_cladding_inner << '\n';
    std::cout << "pcmi_radial_cladding_outer=" << radial_cladding_outer << '\n';
    std::cout << "pcmi_axial_fuel_top=" << axial_fuel_top << '\n';
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
    std::cout << "pcmi_moose_temperature_center_relative_error="
              << temperature_center_error << '\n';
    std::cout << "pcmi_moose_temperature_fuel_surface_relative_error="
              << temperature_fuel_surface_error << '\n';
    std::cout << "pcmi_moose_temperature_cladding_inner_relative_error="
              << temperature_cladding_inner_error << '\n';
    std::cout << "pcmi_moose_radial_fuel_surface_relative_error="
              << radial_fuel_surface_error << '\n';
    std::cout << "pcmi_moose_radial_cladding_inner_relative_error="
              << radial_cladding_inner_error << '\n';
    std::cout << "pcmi_moose_radial_cladding_outer_relative_error="
              << radial_cladding_outer_error << '\n';
    std::cout << "pcmi_moose_axial_fuel_top_relative_error="
              << axial_fuel_top_error << '\n';
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
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_m2_pcmi_solver_tests <pcmi.fsi> "
                     "<all-nodes.csv> <contact-pressure.csv>\n";
        return 2;
    }

    try {
        const std::string input_path = argv[1];
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv,
            "fuelsim M2 PCMI coupled cladding MOOSE comparison tests\n");
        if (!test_pcmi_coupled_cladding(input_path, argv[2], argv[3]))
            return 1;
        std::cout << "[PASS] fuelsim M2 PCMI coupled cladding tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] PCMI tests raised: " << error.what() << '\n';
        return 1;
    }
}
