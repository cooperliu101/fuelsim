#include "fuelsim/case_input.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"
#include "support/moose_field_comparison.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
constexpr double moose_relative_tolerance = 1.0e-3;
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}
double relative_error(double actual, double expected) {
    if (!std::isfinite(actual) || !std::isfinite(expected) || expected == 0.0)
        return std::numeric_limits<double>::infinity();
    return std::abs(actual - expected) / std::abs(expected);
}
bool check_scalar_metrics(const std::string& name, double actual, double reference, double tolerance) {
    const double relative_l2 = relative_error(actual, reference);
    const double relative_absolute_peak = std::abs(std::abs(actual) - std::abs(reference)) / std::abs(reference);
    const double maximum_pointwise_relative = relative_l2;
    std::cout << name << "_relative_l2=" << relative_l2 << '\n';
    std::cout << name << "_relative_absolute_peak=" << relative_absolute_peak << '\n';
    std::cout << name << "_maximum_pointwise_relative=" << maximum_pointwise_relative << '\n';
    return check(
        relative_l2 < tolerance && relative_absolute_peak < tolerance && maximum_pointwise_relative < tolerance,
        name + " three MOOSE error metrics pass");
}
fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    return {definition.solver.absolute_tolerance, definition.solver.relative_tolerance,
        definition.solver.step_tolerance, definition.solver.maximum_iterations};
}
fuelsim::TransientTimeOptions time_options(const fuelsim::FuelSimCaseDefinition& definition) {
    return {definition.transient_execution.end_time, definition.transient_execution.initial_time_step,
        definition.transient_execution.minimum_time_step, definition.transient_execution.maximum_time_step,
        definition.transient_execution.growth_factor, definition.transient_execution.cutback_factor,
        definition.transient_execution.maximum_cutbacks_per_step, definition.transient_execution.load_ramp_time};
}
class TransientCaseRun final {
  public:
    explicit TransientCaseRun(const std::string& input_path, fuelsim::TransientStepObserver* observer = nullptr)
        : _definition(fuelsim::read_case_input(input_path)), _source(fuelsim::read_exodus_quad4(_definition.mesh_file)),
          _problem(_definition.spatial, _source) {
        if (_definition.problem != fuelsim::CaseProblem::transient)
            throw std::invalid_argument("M2.2 comparison requires a transient input card");
        if (fuelsim::rz::ProblemAccess::region_count(_problem) != 1 ||
            fuelsim::rz::ProblemAccess::region_mesh(_problem, 0).elements().size() != 1)
            throw std::invalid_argument("M2.2 comparison input requires one region and one Quad4");
        _result = fuelsim::solve_transient(_problem, time_options(_definition), solver_options(_definition), observer);
    }
    const fuelsim::FuelSimCaseDefinition& definition() const noexcept { return _definition; }
    const fuelsim::UnstructuredQuad4Mesh& source() const noexcept { return _source; }
    const fuelsim::TransientProblem& problem() const noexcept { return _problem; }
    const fuelsim::TransientResult& result() const noexcept { return _result; }

  private:
    fuelsim::FuelSimCaseDefinition _definition;
    fuelsim::UnstructuredQuad4Mesh _source;
    fuelsim::TransientProblem _problem;
    fuelsim::TransientResult _result;
};
struct J2HistoryValue final {
    double time = 0.0;
    double axial_displacement = 0.0;
    double axial_plastic = 0.0;
    double axial_stress = 0.0;
    double effective_plastic = 0.0;
    double hoop_plastic = 0.0;
    double radial_plastic = 0.0;
};
class J2HistoryObserver final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        if (fuelsim::rz::ProblemAccess::region_count(problem) != 1 ||
            fuelsim::rz::ProblemAccess::region_mesh(problem, 0).elements().size() != 1)
            throw std::logic_error("J2 history observer requires one region and one Quad4");
        const fuelsim::RegionMesh& mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, 0);
        const double maximum_z = std::max_element(
            mesh.nodes().begin(), mesh.nodes().end(), [](const fuelsim::RzPoint& left, const fuelsim::RzPoint& right) {
                return left.z < right.z;
            })->z;
        const std::vector<double>& solution = problem.committed_solution();
        double axial_displacement = 0.0;
        std::size_t top_node_count = 0;
        for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
            if (mesh.nodes()[node].z != maximum_z) continue;
            axial_displacement += solution.at(fuelsim::rz::ProblemAccess::dof_map(problem).axial_displacement(node));
            ++top_node_count;
        }
        if (top_node_count == 0) throw std::logic_error("J2 history mesh has no top nodes");
        J2HistoryValue value;
        value.time = step.time;
        value.axial_displacement = axial_displacement / static_cast<double>(top_node_count);
        const std::array<fuelsim::AxisymmetricStressValues, 4>& stresses =
            fuelsim::rz::ProblemAccess::material_stress(problem, 0, 0);
        const fuelsim::Quad4MaterialHistory& history = fuelsim::rz::ProblemAccess::material_history(problem, 0, 0);
        for (std::size_t point = 0; point < history.size(); ++point) {
            value.axial_stress += stresses[point].zz / 4.0;
            value.effective_plastic += history[point].equivalent_plastic_strain / 4.0;
            value.radial_plastic += history[point].plastic_strain[0] / 4.0;
            value.axial_plastic += history[point].plastic_strain[1] / 4.0;
            value.hoop_plastic += history[point].plastic_strain[2] / 4.0;
        }
        _values.push_back(value);
    }
    const std::vector<J2HistoryValue>& values() const noexcept { return _values; }

  private:
    std::vector<J2HistoryValue> _values;
};
std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t separator = line.find(',', begin);
        fields.push_back(line.substr(begin, separator - begin));
        if (separator == std::string::npos) return fields;
        begin = separator + 1;
    }
}
std::size_t csv_column(const std::vector<std::string>& header, const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("MOOSE history CSV is missing column '" + name + "'");
    return static_cast<std::size_t>(found - header.begin());
}
double csv_value(const std::vector<std::string>& fields, std::size_t column, const std::string& path) {
    if (column >= fields.size()) throw std::invalid_argument("MOOSE history CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double value = std::stod(fields[column], &parsed);
    if (parsed != fields[column].size() || !std::isfinite(value))
        throw std::invalid_argument("MOOSE history CSV contains an invalid number: " + path);
    return value;
}
std::vector<J2HistoryValue> read_j2_history(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read MOOSE J2 history: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("MOOSE J2 history is empty: " + path);
    const std::vector<std::string> header = split_csv_line(line);
    const std::size_t time = csv_column(header, "time");
    const std::size_t axial_displacement = csv_column(header, "axial_displacement");
    const std::size_t axial_plastic = csv_column(header, "axial_plastic");
    const std::size_t axial_stress = csv_column(header, "axial_stress");
    const std::size_t effective_plastic = csv_column(header, "effective_plastic");
    const std::size_t hoop_plastic = csv_column(header, "hoop_plastic");
    const std::size_t radial_plastic = csv_column(header, "radial_plastic");
    std::vector<J2HistoryValue> values;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_csv_line(line);
        const J2HistoryValue value{csv_value(fields, time, path), csv_value(fields, axial_displacement, path),
            csv_value(fields, axial_plastic, path), csv_value(fields, axial_stress, path),
            csv_value(fields, effective_plastic, path), csv_value(fields, hoop_plastic, path),
            csv_value(fields, radial_plastic, path)};
        if (value.time > 0.0) values.push_back(value);
    }
    if (values.empty()) throw std::invalid_argument("MOOSE J2 history contains no accepted time steps: " + path);
    return values;
}
bool check_history_metrics(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics, double tolerance,
    double zero_reference_absolute_tolerance) {
    fuelsim::test::print_relative_metrics(name, metrics);
    return check(fuelsim::test::relative_metrics_below(metrics, tolerance) &&
                     metrics.maximum_zero_reference_difference <= zero_reference_absolute_tolerance,
        name + " history three MOOSE error metrics and zero references pass");
}
double average_axial_stress(const TransientCaseRun& run) {
    double value = 0.0;
    for (const fuelsim::AxisymmetricStressValues& stress :
        fuelsim::rz::ProblemAccess::material_stress(run.problem(), 0, 0))
        value += stress.zz / 4.0;
    return value;
}
double average_equivalent_plastic(const TransientCaseRun& run) {
    double value = 0.0;
    for (const fuelsim::MaterialPointState& point : fuelsim::rz::ProblemAccess::material_history(run.problem(), 0, 0))
        value += point.equivalent_plastic_strain / 4.0;
    return value;
}
double average_equivalent_creep(const TransientCaseRun& run) {
    double value = 0.0;
    for (const fuelsim::MaterialPointState& point : fuelsim::rz::ProblemAccess::material_history(run.problem(), 0, 0))
        value += point.equivalent_creep_strain / 4.0;
    return value;
}
double maximum_inelastic_trace(const TransientCaseRun& run, bool plastic_strain) {
    double maximum = 0.0;
    for (const fuelsim::MaterialPointState& point : fuelsim::rz::ProblemAccess::material_history(run.problem(), 0, 0)) {
        const std::array<double, 4>& strain = plastic_strain ? point.plastic_strain : point.creep_strain;
        maximum = std::max(maximum, std::abs(strain[0] + strain[1] + strain[2]));
    }
    return maximum;
}
double average_top_displacement(const TransientCaseRun& run) {
    const fuelsim::RegionBoundary top =
        fuelsim::rz::ProblemAccess::region_mesh(run.problem(), 0).map_side_set(run.source(), "top");
    double value = 0.0;
    for (std::size_t node : top.nodes) {
        value += run.result().committed_state.at(
            fuelsim::rz::ProblemAccess::dof_map(run.problem()).axial_displacement(node));
    }
    return value / static_cast<double>(top.nodes.size());
}
bool temperatures_are_600(const TransientCaseRun& run) {
    for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::dof_map(run.problem()).node_count(); ++node)
        if (std::abs(
                run.result().committed_state.at(fuelsim::rz::ProblemAccess::dof_map(run.problem()).temperature(node)) -
                600.0) > 1.0e-12)
            return false;
    return true;
}
bool common_run_checks(const std::string& name, const TransientCaseRun& run, std::size_t expected_steps,
    const std::string& nodal_reference_path) {
    if (!run.result().completed)
        std::cerr << name << " solve failure: "
                  << fuelsim::solve_failure_category_name(run.result().last_attempt.failure_category) << ": "
                  << run.result().last_attempt.failure_message << '\n';
    bool passed = check(run.result().completed, name + " input-card load path converged");
    passed =
        check(run.result().accepted_steps.size() == expected_steps, name + " commits every configured time step") &&
        passed;
    passed = check(run.result().aggregate_timing.workspace_setups == 1, name + " reuses one PETSc workspace") && passed;
    passed = check(temperatures_are_600(run), name + " committed temperature remains 600 K") && passed;
    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(run.problem(), run.result().committed_state, reference);
    passed = check(fields.node_count == run.source().nodes().size() && fields.maximum_coordinate_difference < 1.0e-12,
                 name + " compares every MOOSE node at matching coordinates") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.temperature, moose_relative_tolerance),
                 name + " full-field temperature three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.radial_displacement, moose_relative_tolerance),
                 name + " full-field radial displacement three errors pass") &&
             passed;
    passed = check(fuelsim::test::relative_metrics_below(fields.axial_displacement, moose_relative_tolerance),
                 name + " full-field axial displacement three errors pass") &&
             passed;
    fuelsim::test::print_relative_metrics("m22_" + name + "_temperature", fields.temperature);
    fuelsim::test::print_relative_metrics("m22_" + name + "_radial_displacement", fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m22_" + name + "_axial_displacement", fields.axial_displacement);
    return passed;
}
bool test_j2_moose_comparison(const std::string& input_path, const std::string& nodal_reference_path) {
    constexpr double moose_axial_stress = 201980198.0198;
    constexpr double moose_equivalent_plastic = 0.000990099009901;
    constexpr double moose_top_displacement = 2.0e-6;
    const TransientCaseRun run(input_path);
    const double axial_stress = average_axial_stress(run);
    const double equivalent_plastic = average_equivalent_plastic(run);
    const double top_displacement = average_top_displacement(run);
    bool passed = common_run_checks("j2", run, 10, nodal_reference_path);
    passed = check_scalar_metrics("m22_j2_axial_stress", axial_stress, moose_axial_stress, moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics(
                 "m22_j2_equivalent_plastic", equivalent_plastic, moose_equivalent_plastic, moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics(
                 "m22_j2_top_displacement", top_displacement, moose_top_displacement, moose_relative_tolerance) &&
             passed;
    passed = check(maximum_inelastic_trace(run, true) < 1.0e-12, "J2 plastic strain is trace-free") && passed;
    return passed;
}
bool test_j2_unload_reload_moose_comparison(
    const std::string& input_path, const std::string& nodal_reference_path, const std::string& history_reference_path) {
    J2HistoryObserver observer;
    const TransientCaseRun run(input_path, &observer);
    const std::vector<J2HistoryValue> reference = read_j2_history(history_reference_path);
    const std::vector<J2HistoryValue>& actual = observer.values();
    bool passed = common_run_checks("j2_unload_reload", run, 30, nodal_reference_path);
    passed = check(actual.size() == reference.size() && actual.size() == 30,
                 "J2 unload-reload compares all 30 accepted MOOSE steps") &&
             passed;
    if (actual.size() != reference.size()) return false;
    fuelsim::test::FieldErrorMetrics axial_displacement;
    fuelsim::test::FieldErrorMetrics axial_plastic;
    fuelsim::test::FieldErrorMetrics axial_stress;
    fuelsim::test::FieldErrorMetrics effective_plastic;
    fuelsim::test::FieldErrorMetrics hoop_plastic;
    fuelsim::test::FieldErrorMetrics radial_plastic;
    double maximum_time_difference = 0.0;
    for (std::size_t step = 0; step < actual.size(); ++step) {
        maximum_time_difference = std::max(maximum_time_difference, std::abs(actual[step].time - reference[step].time));
        axial_displacement.add(actual[step].axial_displacement, reference[step].axial_displacement);
        axial_plastic.add(actual[step].axial_plastic, reference[step].axial_plastic);
        axial_stress.add(actual[step].axial_stress, reference[step].axial_stress);
        effective_plastic.add(actual[step].effective_plastic, reference[step].effective_plastic);
        hoop_plastic.add(actual[step].hoop_plastic, reference[step].hoop_plastic);
        radial_plastic.add(actual[step].radial_plastic, reference[step].radial_plastic);
    }
    std::cout << "m22_j2_unload_reload_maximum_time_difference=" << maximum_time_difference << '\n';
    passed = check(maximum_time_difference < 1.0e-12, "J2 unload-reload time coordinates match MOOSE") && passed;
    passed = check_history_metrics(
                 "m22_j2_unload_reload_axial_displacement", axial_displacement, moose_relative_tolerance, 0.0) &&
             passed;
    passed =
        check_history_metrics("m22_j2_unload_reload_axial_plastic", axial_plastic, moose_relative_tolerance, 1.0e-14) &&
        passed;
    passed =
        check_history_metrics("m22_j2_unload_reload_axial_stress", axial_stress, moose_relative_tolerance, 1.0e-6) &&
        passed;
    passed = check_history_metrics(
                 "m22_j2_unload_reload_effective_plastic", effective_plastic, moose_relative_tolerance, 1.0e-14) &&
             passed;
    passed =
        check_history_metrics("m22_j2_unload_reload_hoop_plastic", hoop_plastic, moose_relative_tolerance, 1.0e-14) &&
        passed;
    passed = check_history_metrics(
                 "m22_j2_unload_reload_radial_plastic", radial_plastic, moose_relative_tolerance, 1.0e-14) &&
             passed;
    const J2HistoryValue& first_peak = actual.at(9);
    const J2HistoryValue& unloaded = actual.at(19);
    const J2HistoryValue& reloaded = actual.at(29);
    std::cout << "m22_j2_unload_reload_first_peak_stress=" << first_peak.axial_stress << '\n';
    std::cout << "m22_j2_unload_reload_unloaded_stress=" << unloaded.axial_stress << '\n';
    std::cout << "m22_j2_unload_reload_first_peak_effective_plastic=" << first_peak.effective_plastic << '\n';
    std::cout << "m22_j2_unload_reload_unloaded_effective_plastic=" << unloaded.effective_plastic << '\n';
    std::cout << "m22_j2_unload_reload_reloaded_effective_plastic=" << reloaded.effective_plastic << '\n';
    passed = check(first_peak.axial_stress > 2.0e8 && unloaded.axial_stress < 0.0,
                 "J2 path reaches tensile plasticity then reverses stress during unload") &&
             passed;
    passed = check(std::abs(unloaded.effective_plastic - first_peak.effective_plastic) < 1.0e-14,
                 "J2 elastic unload preserves committed equivalent plastic strain") &&
             passed;
    passed = check(reloaded.effective_plastic > first_peak.effective_plastic,
                 "J2 reload activates additional plastic strain") &&
             passed;
    passed =
        check(maximum_inelastic_trace(run, true) < 1.0e-12, "J2 unload-reload plastic strain remains trace-free") &&
        passed;
    return passed;
}
bool test_norton_moose_comparison(const std::string& input_path, const std::string& nodal_reference_path) {
    constexpr double moose_axial_stress = 99998007.620195;
    constexpr double moose_equivalent_creep = 9.9991036503199e-5;
    constexpr double moose_top_displacement = 5.9997808603447e-7;
    const TransientCaseRun run(input_path);
    const double axial_stress = average_axial_stress(run);
    const double equivalent_creep = average_equivalent_creep(run);
    const double top_displacement = average_top_displacement(run);
    bool passed = common_run_checks("norton", run, 10, nodal_reference_path);
    passed =
        check_scalar_metrics("m22_norton_axial_stress", axial_stress, moose_axial_stress, moose_relative_tolerance) &&
        passed;
    passed = check_scalar_metrics(
                 "m22_norton_equivalent_creep", equivalent_creep, moose_equivalent_creep, moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics(
                 "m22_norton_top_displacement", top_displacement, moose_top_displacement, moose_relative_tolerance) &&
             passed;
    passed = check(maximum_inelastic_trace(run, false) < 1.0e-12, "Norton creep strain is trace-free") && passed;
    return passed;
}
bool test_coupled_moose_comparison(const std::string& displacement_input, const std::string& traction_input,
    const std::string& displacement_reference, const std::string& traction_reference) {
    constexpr double moose_axial_stress = 200999992.08159;
    constexpr double moose_equivalent_plastic = 4.9999406119027e-4;
    constexpr double moose_equivalent_creep = 2.4564701918764e-4;
    constexpr double moose_top_displacement = 1.7506410289082e-6;
    constexpr double analytic_equivalent_plastic = 5.0e-4;
    constexpr double analytic_equivalent_creep = 2.4564818025e-4;
    constexpr double analytic_top_displacement = 1.75064818025e-6;
    const TransientCaseRun traction(traction_input);
    const double axial_stress = average_axial_stress(traction);
    const double equivalent_plastic = average_equivalent_plastic(traction);
    const double equivalent_creep = average_equivalent_creep(traction);
    const double top_displacement = average_top_displacement(traction);
    bool passed = common_run_checks("coupled_traction", traction, 10, traction_reference);
    passed = check_scalar_metrics(
                 "m22_coupled_traction_axial_stress", axial_stress, moose_axial_stress, moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_traction_equivalent_plastic", equivalent_plastic,
                 moose_equivalent_plastic, moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_traction_equivalent_creep", equivalent_creep, moose_equivalent_creep,
                 moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_traction_top_displacement", top_displacement, moose_top_displacement,
                 moose_relative_tolerance) &&
             passed;
    passed = check(relative_error(equivalent_plastic, analytic_equivalent_plastic) < 1.0e-8 &&
                       relative_error(equivalent_creep, analytic_equivalent_creep) < 1.0e-8 &&
                       relative_error(top_displacement, analytic_top_displacement) < 1.0e-8,
                 "coupled traction matches independent uniaxial history") &&
             passed;
    passed =
        check(maximum_inelastic_trace(traction, true) < 1.0e-12 && maximum_inelastic_trace(traction, false) < 1.0e-12,
            "coupled traction inelastic strains are trace-free") &&
        passed;
    constexpr double displacement_moose_stress = 200963368.63564;
    constexpr double displacement_moose_plastic = 4.8168428343307e-4;
    constexpr double displacement_moose_creep = 5.1349887359505e-4;
    constexpr double displacement_moose_displacement = 2.0e-6;
    const TransientCaseRun displacement(displacement_input);
    passed = common_run_checks("coupled_displacement", displacement, 10, displacement_reference) && passed;
    passed = check_scalar_metrics("m22_coupled_displacement_axial_stress", average_axial_stress(displacement),
                 displacement_moose_stress, moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_displacement_equivalent_plastic",
                 average_equivalent_plastic(displacement), displacement_moose_plastic, moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_displacement_equivalent_creep", average_equivalent_creep(displacement),
                 displacement_moose_creep, moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_displacement_top_displacement", average_top_displacement(displacement),
                 displacement_moose_displacement, moose_relative_tolerance) &&
             passed;
    return passed;
}
} // namespace
int main(int argc, char** argv) {
    if (argc != 12) {
        std::cerr << "Usage: fuelsim_m2_inelastic_solver_tests "
                     "<j2.fsi> <norton.fsi> <coupled_displacement.fsi> "
                     "<coupled_traction.fsi> <j2-nodes.csv> "
                     "<norton-nodes.csv> <coupled-displacement-nodes.csv> "
                     "<coupled-traction-nodes.csv> <j2-unload-reload.fsi> "
                     "<j2-unload-reload-nodes.csv> "
                     "<j2-unload-reload-history.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim input-card M2.2 MOOSE comparison tests\n");
        bool passed = test_j2_moose_comparison(argv[1], argv[5]);
        passed = test_norton_moose_comparison(argv[2], argv[6]) && passed;
        passed = test_coupled_moose_comparison(argv[3], argv[4], argv[7], argv[8]) && passed;
        passed = test_j2_unload_reload_moose_comparison(argv[9], argv[10], argv[11]) && passed;
        if (!passed) return 1;
        std::cout << "[PASS] input-card M2.2 MOOSE comparison tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M2.2 tests raised: " << error.what() << '\n';
        return 1;
    }
}
