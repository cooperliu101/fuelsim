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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double comparison_tolerance = 5.0e-3;
// Local tensor components cross zero during reversal.  These pointwise-only
// gates are paired with the unchanged 0.5% L2 and peak gates.  Stress and
// elastic strain need only a narrow 0.6% allowance at the same low-magnitude
// shear-stress reversal point; combined inelastic strain remains below 0.5%.
constexpr double stress_pointwise_tolerance = 6.0e-3;
constexpr double strain_pointwise_tolerance = 6.0e-3;
constexpr double inelastic_pointwise_tolerance = comparison_tolerance;
constexpr double time_tolerance = 1.0e-12;

enum class ExpectedBehavior {
    elastic,
    plastic,
    creep,
    coupled,
};

struct VariantConfig final {
    std::string name;
    ExpectedBehavior behavior = ExpectedBehavior::elastic;
    std::size_t element_count = 0;
    std::size_t step_count = 0;
    double comparison_tolerance = 0.0;
    double stress_pointwise_tolerance = 0.0;
    double elastic_pointwise_tolerance = 0.0;
    double inelastic_pointwise_tolerance = 0.0;
    double plastic_trace_tolerance = 0.0;
    double creep_trace_tolerance = 0.0;
};

VariantConfig variant_config(const std::string& name) {
    if (name == "production")
        return {name,
                ExpectedBehavior::coupled,
                4,
                100,
                comparison_tolerance,
                stress_pointwise_tolerance,
                strain_pointwise_tolerance,
                inelastic_pointwise_tolerance,
                7.0e-6,
                2.0e-7};
    if (name == "elastic_displacement")
        return {name, ExpectedBehavior::elastic, 4, 100,
                comparison_tolerance, 3.0e-5, 2.5e-4, 0.0,
                1.0e-14, 1.0e-14};
    if (name == "plastic_displacement")
        return {name, ExpectedBehavior::plastic, 4, 100,
                comparison_tolerance, 3.0e-3, 3.0e-3, 4.5e-3,
                7.0e-6, 1.0e-14};
    if (name == "creep_displacement")
        return {name, ExpectedBehavior::creep, 4, 100,
                comparison_tolerance, 4.0e-2, 2.5e-3, 2.0e-4,
                1.0e-14, 1.3e-5};
    if (name == "coupled_displacement")
        return {name, ExpectedBehavior::coupled, 4, 100,
                comparison_tolerance, 5.0e-3, 5.0e-3, 4.0e-3,
                7.0e-6, 2.0e-7};
    if (name == "coupled_pressure")
        return {name, ExpectedBehavior::coupled, 4, 100,
                comparison_tolerance, stress_pointwise_tolerance,
                strain_pointwise_tolerance, inelastic_pointwise_tolerance,
                7.0e-6, 2.0e-7};
    if (name == "material_oracle")
        return {name, ExpectedBehavior::coupled, 1, 100,
                5.0e-6, 1.0e-6, 1.0e-6, 1.0e-6,
                2.0e-6, 1.0e-7};
    throw std::invalid_argument("Unknown M4.3 comparison variant: " + name);
}

bool plastic_active(ExpectedBehavior behavior) noexcept {
    return behavior == ExpectedBehavior::plastic ||
           behavior == ExpectedBehavior::coupled;
}

bool creep_active(ExpectedBehavior behavior) noexcept {
    return behavior == ExpectedBehavior::creep ||
           behavior == ExpectedBehavior::coupled;
}

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t end = line.find(',', begin);
        fields.push_back(line.substr(begin, end - begin));
        if (end == std::string::npos)
            return fields;
        begin = end + 1;
    }
}

std::size_t column_index(const std::vector<std::string>& header,
                         const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end())
        throw std::invalid_argument("M4.3 CSV is missing column '" + name +
                                    "'");
    return static_cast<std::size_t>(found - header.begin());
}

double csv_value(const std::vector<std::string>& fields, std::size_t column,
                 const std::string& path) {
    if (column >= fields.size())
        throw std::invalid_argument("M4.3 CSV row is incomplete: " + path);
    std::size_t parsed = 0;
    const double value = std::stod(fields[column], &parsed);
    if (parsed != fields[column].size() || !std::isfinite(value))
        throw std::invalid_argument("M4.3 CSV value is invalid: " + path);
    return value;
}

std::size_t csv_id(const std::vector<std::string>& fields,
                   std::size_t column, const std::string& path) {
    const double value = csv_value(fields, column, path);
    if (value < 0.0 ||
        value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("M4.3 CSV element ID is invalid: " +
                                    path);
    return static_cast<std::size_t>(value);
}

struct ElementSnapshot final {
    std::size_t element_id = 0;
    double radius = 0.0;
    double axial_coordinate = 0.0;
    fuelsim::AxisymmetricStressValues stress{};
    fuelsim::MaterialPointState state{};
};

struct HistorySnapshot final {
    double time = 0.0;
    double reference_height = 0.0;
    double top_radial_displacement = 0.0;
    double top_axial_displacement = 0.0;
    std::vector<ElementSnapshot> elements;
};

class HistoryObserver final : public fuelsim::TransientStepObserver {
  public:
    explicit HistoryObserver(std::size_t expected_element_count)
        : _expected_element_count(expected_element_count) {}

    void accepted_step(const fuelsim::TransientProblem& problem,
                       const fuelsim::TransientAcceptedStep& step) override {
        if (problem.region_count() != 1 ||
            problem.region_mesh(0).elements().size() !=
                _expected_element_count)
            throw std::logic_error(
                "M4.3 observer element count differs from its variant");
        const fuelsim::RegionMesh& mesh = problem.region_mesh(0);
        const std::vector<double>& solution = problem.committed_solution();
        const double maximum_z = std::max_element(
                                     mesh.nodes().begin(), mesh.nodes().end(),
                                     [](const fuelsim::RzPoint& left,
                                        const fuelsim::RzPoint& right) {
                                         return left.z < right.z;
                                     })
                                     ->z;
        const double minimum_z = std::min_element(
                                     mesh.nodes().begin(), mesh.nodes().end(),
                                     [](const fuelsim::RzPoint& left,
                                        const fuelsim::RzPoint& right) {
                                         return left.z < right.z;
                                     })
                                     ->z;
        double radial = 0.0;
        double axial = 0.0;
        std::size_t count = 0;
        for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
            if (std::abs(mesh.nodes()[node].z - maximum_z) > time_tolerance)
                continue;
            radial += solution.at(
                problem.dof_map().radial_displacement(node));
            axial += solution.at(problem.dof_map().axial_displacement(node));
            ++count;
        }
        if (count == 0)
            throw std::logic_error("M4.3 mesh has no top nodes");

        HistorySnapshot snapshot;
        snapshot.time = step.time;
        snapshot.reference_height = maximum_z - minimum_z;
        snapshot.top_radial_displacement =
            radial / static_cast<double>(count);
        snapshot.top_axial_displacement = axial / static_cast<double>(count);
        snapshot.elements.reserve(mesh.elements().size());
        for (std::size_t element = 0; element < mesh.elements().size();
             ++element) {
            ElementSnapshot value;
            value.element_id = mesh.source_element_ids().at(element);
            for (const std::size_t node : mesh.elements()[element].nodes) {
                value.radius += mesh.nodes()[node].r / 4.0;
                value.axial_coordinate += mesh.nodes()[node].z / 4.0;
            }

            const fuelsim::Quad4RzGeometry& geometry =
                problem.region_element_geometry(0, element);
            const std::array<fuelsim::AxisymmetricStressValues, 4>& stresses =
                problem.material_stress(0, element);
            const fuelsim::Quad4MaterialHistory& history =
                problem.material_history(0, element);
            double total_weight = 0.0;
            for (std::size_t point = 0; point < geometry.points.size();
                 ++point) {
                const double weight = geometry.points[point].weighted_measure;
                total_weight += weight;
                value.stress.rr += weight * stresses[point].rr;
                value.stress.zz += weight * stresses[point].zz;
                value.stress.hoop += weight * stresses[point].hoop;
                value.stress.rz += weight * stresses[point].rz;
                for (std::size_t component = 0; component < 4; ++component) {
                    value.state.elastic_strain[component] +=
                        weight * history[point].elastic_strain[component];
                    value.state.plastic_strain[component] +=
                        weight * history[point].plastic_strain[component];
                    value.state.creep_strain[component] +=
                        weight * history[point].creep_strain[component];
                }
                value.state.equivalent_plastic_strain +=
                    weight * history[point].equivalent_plastic_strain;
                value.state.equivalent_creep_strain +=
                    weight * history[point].equivalent_creep_strain;
                _maximum_plastic_trace = std::max(
                    _maximum_plastic_trace,
                    std::abs(history[point].plastic_strain[0] +
                             history[point].plastic_strain[1] +
                             history[point].plastic_strain[2]));
                _maximum_creep_trace = std::max(
                    _maximum_creep_trace,
                    std::abs(history[point].creep_strain[0] +
                             history[point].creep_strain[1] +
                             history[point].creep_strain[2]));
            }
            if (!(total_weight > 0.0))
                throw std::logic_error(
                    "M4.3 element has zero reference volume");
            value.stress.rr /= total_weight;
            value.stress.zz /= total_weight;
            value.stress.hoop /= total_weight;
            value.stress.rz /= total_weight;
            for (std::size_t component = 0; component < 4; ++component) {
                value.state.elastic_strain[component] /= total_weight;
                value.state.plastic_strain[component] /= total_weight;
                value.state.creep_strain[component] /= total_weight;
            }
            value.state.equivalent_plastic_strain /= total_weight;
            value.state.equivalent_creep_strain /= total_weight;
            snapshot.elements.push_back(value);
        }
        std::sort(snapshot.elements.begin(), snapshot.elements.end(),
                  [](const ElementSnapshot& left,
                     const ElementSnapshot& right) {
                      return left.element_id < right.element_id;
                  });
        _snapshots.push_back(std::move(snapshot));
    }

    const std::vector<HistorySnapshot>& snapshots() const noexcept {
        return _snapshots;
    }

    double maximum_plastic_trace() const noexcept {
        return _maximum_plastic_trace;
    }

    double maximum_creep_trace() const noexcept {
        return _maximum_creep_trace;
    }

  private:
    std::size_t _expected_element_count;
    std::vector<HistorySnapshot> _snapshots;
    double _maximum_plastic_trace = 0.0;
    double _maximum_creep_trace = 0.0;
};

struct ReferenceSnapshot final {
    double time = 0.0;
    std::size_t element_id = 0;
    double radius = 0.0;
    double axial_coordinate = 0.0;
    std::array<double, 4> stress{};
    std::array<double, 4> elastic{};
    std::array<double, 4> combined_inelastic{};
    double equivalent_plastic = 0.0;
    double equivalent_creep = 0.0;
};

std::vector<ReferenceSnapshot>
read_reference_history(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read M4.3 MOOSE history: " + path);
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("M4.3 MOOSE history is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const auto column = [&header](const std::string& name) {
        return column_index(header, name);
    };
    const std::size_t time = column("sample_time");
    const std::size_t id = column("id");
    const std::size_t radius = column("x");
    const std::size_t axial_coordinate = column("y");
    const std::array<std::size_t, 4> stress = {
        column("stress_rr"), column("stress_zz"), column("stress_hoop"),
        column("stress_rz")};
    const std::array<std::size_t, 4> elastic = {
        column("elastic_rr"), column("elastic_zz"),
        column("elastic_hoop"), column("elastic_rz")};
    const std::array<std::size_t, 4> combined_inelastic = {
        column("combined_rr"), column("combined_zz"),
        column("combined_hoop"), column("combined_rz")};
    const std::size_t equivalent_plastic = column("effective_plastic");
    const std::size_t equivalent_creep = column("effective_creep");

    std::vector<ReferenceSnapshot> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> fields = split_csv(line);
        ReferenceSnapshot value;
        value.time = csv_value(fields, time, path);
        value.element_id = csv_id(fields, id, path);
        value.radius = csv_value(fields, radius, path);
        value.axial_coordinate = csv_value(fields, axial_coordinate, path);
        for (std::size_t component = 0; component < 4; ++component) {
            value.stress[component] =
                csv_value(fields, stress[component], path);
            value.elastic[component] =
                csv_value(fields, elastic[component], path);
            value.combined_inelastic[component] =
                csv_value(fields, combined_inelastic[component], path);
        }
        value.equivalent_plastic =
            csv_value(fields, equivalent_plastic, path);
        value.equivalent_creep =
            csv_value(fields, equivalent_creep, path);
        result.push_back(value);
    }
    if (result.empty())
        throw std::invalid_argument(
            "M4.3 MOOSE history has no transient rows: " + path);
    std::sort(result.begin(), result.end(),
              [](const ReferenceSnapshot& left,
                 const ReferenceSnapshot& right) {
                  if (left.time != right.time)
                      return left.time < right.time;
                  return left.element_id < right.element_id;
              });
    return result;
}

std::array<double, 4>
stress_components(const fuelsim::AxisymmetricStressValues& stress) {
    return {stress.rr, stress.zz, stress.hoop, stress.rz};
}

bool check_metrics(const std::string& name,
                   const fuelsim::test::FieldErrorMetrics& metrics,
                   double zero_reference_tolerance,
                   double aggregate_tolerance,
                   double pointwise_tolerance) {
    fuelsim::test::print_relative_metrics(name, metrics);
    std::cout << name << "_maximum_absolute_difference="
              << metrics.maximum_absolute_difference << '\n';
    std::cout << name << "_relative_l2_tolerance="
              << aggregate_tolerance << '\n';
    std::cout << name << "_relative_absolute_peak_tolerance="
              << aggregate_tolerance << '\n';
    std::cout << name << "_maximum_pointwise_relative_tolerance="
              << pointwise_tolerance << '\n';
    return check(metrics.relative_l2() < aggregate_tolerance &&
                     metrics.relative_absolute_peak() <
                         aggregate_tolerance &&
                     metrics.maximum_pointwise_relative_error() <
                         pointwise_tolerance &&
                     metrics.maximum_zero_reference_difference <
                         zero_reference_tolerance,
                 name + " three MOOSE metrics and zero-reference error pass");
}

void print_tensor_metric_locations(
    const std::string& name,
    const fuelsim::test::FieldErrorMetrics& metrics,
    const std::vector<ReferenceSnapshot>& reference) {
    constexpr std::array<const char*, 4> components = {"rr", "zz", "hoop",
                                                        "rz"};
    const auto print_location = [&name, &reference, &components](
                                    const std::string& metric,
                                    std::size_t flat_index) {
        const std::size_t row = flat_index / components.size();
        const std::size_t component = flat_index % components.size();
        if (row >= reference.size())
            throw std::logic_error("M4.3 metric index exceeds history size");
        std::cout << name << '_' << metric
                  << "_time=" << reference[row].time << '\n';
        std::cout << name << '_' << metric
                  << "_element_id=" << reference[row].element_id << '\n';
        std::cout << name << '_' << metric
                  << "_component=" << components[component] << '\n';
    };
    print_location("maximum_pointwise_relative_location",
                   metrics.maximum_pointwise_relative_index);
    print_location("maximum_absolute_difference_location",
                   metrics.maximum_absolute_difference_index);
}

const HistorySnapshot& snapshot_at(const std::vector<HistorySnapshot>& values,
                                   double time) {
    const auto found = std::find_if(
        values.begin(), values.end(), [time](const HistorySnapshot& value) {
            return std::abs(value.time - time) < time_tolerance;
        });
    if (found == values.end())
        throw std::invalid_argument("M4.3 accepted history misses event time");
    return *found;
}

bool check_load_path(const std::vector<HistorySnapshot>& snapshots,
                     ExpectedBehavior behavior) {
    const HistorySnapshot& first_stretch = snapshot_at(snapshots, 1.0);
    const HistorySnapshot& positive_shear = snapshot_at(snapshots, 2.0);
    const HistorySnapshot& axial_reversal = snapshot_at(snapshots, 3.0);
    const HistorySnapshot& shear_reversal = snapshot_at(snapshots, 4.0);
    const HistorySnapshot& final_stretch = snapshot_at(snapshots, 5.0);
    const auto close = [](double actual, double expected) {
        return std::abs(actual - expected) < 1.0e-14;
    };
    bool passed =
        check(close(first_stretch.top_radial_displacement, 0.0) &&
                  close(first_stretch.top_axial_displacement, 2.0e-4) &&
                  close(positive_shear.top_radial_displacement, 2.0e-3) &&
                  close(positive_shear.top_axial_displacement, 2.0e-4) &&
                  close(axial_reversal.top_radial_displacement, 2.0e-3) &&
                  close(axial_reversal.top_axial_displacement, -1.0e-4) &&
                  close(shear_reversal.top_radial_displacement, -1.5e-3) &&
                  close(shear_reversal.top_axial_displacement, -1.0e-4) &&
                  close(final_stretch.top_radial_displacement, -1.5e-3) &&
                  close(final_stretch.top_axial_displacement, 6.0e-5),
              "M4.3 hits every noncoaxial load-path event");
    const double shear = positive_shear.top_radial_displacement /
                         positive_shear.reference_height;
    const double axial_stretch =
        1.0 + positive_shear.top_axial_displacement /
                  positive_shear.reference_height;
    const double positive_polar_rotation =
        std::abs(std::atan2(-shear, 1.0 + axial_stretch));
    std::cout << "m43_positive_polar_rotation=" << positive_polar_rotation
              << '\n';
    passed = check(positive_polar_rotation > 0.44,
                   "M4.3 positive-shear stage exceeds 25 degrees rotation") &&
             passed;
    double maximum_plastic_shear = 0.0;
    double maximum_creep_shear = 0.0;
    for (const ElementSnapshot& element : positive_shear.elements) {
        maximum_plastic_shear =
            std::max(maximum_plastic_shear,
                     std::abs(element.state.plastic_strain[3]));
        maximum_creep_shear =
            std::max(maximum_creep_shear,
                     std::abs(element.state.creep_strain[3]));
    }
    if (plastic_active(behavior))
        passed = check(maximum_plastic_shear > 0.1,
                       "M4.3 activates rotated plastic shear history") &&
                 passed;
    else
        passed = check(maximum_plastic_shear == 0.0,
                       "M4.3 inactive plastic history remains zero") &&
                 passed;
    if (creep_active(behavior))
        passed = check(maximum_creep_shear > 1.0e-5,
                       "M4.3 activates rotated creep shear history") &&
                 passed;
    else
        passed = check(maximum_creep_shear == 0.0,
                       "M4.3 inactive creep history remains zero") &&
                 passed;
    return passed;
}

bool run_test(const VariantConfig& variant, const std::string& input_path,
              const std::string& nodal_reference_path,
              const std::string& history_reference_path) {
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(input_path);
    if (definition.problem != fuelsim::CaseProblem::transient)
        throw std::invalid_argument("M4.3 requires a transient input card");
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::ExodusMeshIo::read_quad4(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.transient_definition(),
                                      source);
    const fuelsim::TransientTimeOptions time_options = {
        definition.transient_execution.end_time,
        definition.transient_execution.initial_time_step,
        definition.transient_execution.minimum_time_step,
        definition.transient_execution.maximum_time_step,
        definition.transient_execution.growth_factor,
        definition.transient_execution.cutback_factor,
        definition.transient_execution.maximum_cutbacks,
        definition.transient_execution.load_ramp_time};
    fuelsim::SolverOptions solver_options = {
        definition.solver.absolute_tolerance,
        definition.solver.relative_tolerance,
        definition.solver.step_tolerance,
        definition.solver.maximum_iterations};
    HistoryObserver observer(variant.element_count);
    const fuelsim::TransientResult result = fuelsim::solve_transient(
        problem, time_options, solver_options, &observer);

    bool passed = check(result.completed &&
                            result.accepted_steps.size() ==
                                variant.step_count &&
                            result.rejected_steps.empty(),
                        "M4.3 completes every fixed step without rejection");
    std::cout << "m43_variant=" << variant.name << '\n';
    std::cout << "m43_completed=" << result.completed << '\n';
    std::cout << "m43_accepted_steps=" << result.accepted_steps.size() << '\n';
    std::cout << "m43_rejected_steps=" << result.rejected_steps.size() << '\n';
    if (!result.completed) {
        std::cout << "m43_last_failure_category="
                  << fuelsim::solve_failure_category_name(
                         result.last_attempt.failure_category)
                  << '\n';
        std::cout << "m43_last_failure_message="
                  << result.last_attempt.failure_message << '\n';
        return false;
    }
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "M4.3 reuses one PETSc workspace") &&
             passed;
    passed = check_load_path(observer.snapshots(), variant.behavior) && passed;

    const std::vector<fuelsim::test::NodalFieldReference> nodal_reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison nodal =
        fuelsim::test::compare_moose_nodal_fields(
            problem, result.committed_state, nodal_reference);
    passed = check(nodal.node_count == source.nodes().size() &&
                       nodal.maximum_coordinate_difference < 1.0e-12,
                   "M4.3 compares every MOOSE node at matching coordinates") &&
             passed;
    passed = check_metrics("m43_temperature", nodal.temperature, 1.0e-12,
                           variant.comparison_tolerance,
                           variant.comparison_tolerance) &&
             passed;
    passed = check_metrics("m43_radial_displacement",
                           nodal.radial_displacement, 1.0e-14,
                           variant.comparison_tolerance,
                           variant.comparison_tolerance) &&
             passed;
    passed = check_metrics("m43_axial_displacement",
                           nodal.axial_displacement, 1.0e-14,
                           variant.comparison_tolerance,
                           variant.comparison_tolerance) &&
             passed;

    const std::vector<ReferenceSnapshot> reference =
        read_reference_history(history_reference_path);
    const std::vector<HistorySnapshot>& actual = observer.snapshots();
    passed = check(reference.size() ==
                           actual.size() * variant.element_count,
                   "M4.3 compares every element at every accepted step") &&
             passed;
    if (reference.size() != actual.size() * variant.element_count)
        return false;

    fuelsim::test::FieldErrorMetrics stress;
    fuelsim::test::FieldErrorMetrics elastic;
    fuelsim::test::FieldErrorMetrics combined_inelastic;
    fuelsim::test::FieldErrorMetrics equivalent_plastic;
    fuelsim::test::FieldErrorMetrics equivalent_creep;
    double maximum_time_difference = 0.0;
    double maximum_coordinate_difference = 0.0;
    std::size_t reference_row = 0;
    for (std::size_t step = 0; step < actual.size(); ++step) {
        passed = check(actual[step].elements.size() == variant.element_count,
                       "M4.3 observer records every element per step") &&
                 passed;
        if (actual[step].elements.size() != variant.element_count)
            return false;
        for (const ElementSnapshot& element : actual[step].elements) {
            const ReferenceSnapshot& expected = reference[reference_row];
            maximum_time_difference = std::max(
                maximum_time_difference,
                std::abs(actual[step].time - expected.time));
            maximum_coordinate_difference = std::max(
                {maximum_coordinate_difference,
                 std::abs(element.radius - expected.radius),
                 std::abs(element.axial_coordinate -
                          expected.axial_coordinate)});
            passed = check(element.element_id == expected.element_id,
                           "M4.3 element IDs match at every time step") &&
                     passed;
            const std::array<double, 4> actual_stress =
                stress_components(element.stress);
            for (std::size_t component = 0; component < 4; ++component) {
                stress.add(actual_stress[component],
                           expected.stress[component]);
                elastic.add(element.state.elastic_strain[component],
                            expected.elastic[component]);
                combined_inelastic.add(
                    element.state.plastic_strain[component] +
                        element.state.creep_strain[component],
                    expected.combined_inelastic[component]);
            }
            equivalent_plastic.add(
                element.state.equivalent_plastic_strain,
                expected.equivalent_plastic);
            equivalent_creep.add(element.state.equivalent_creep_strain,
                                 expected.equivalent_creep);
            ++reference_row;
        }
    }
    passed = check(maximum_time_difference < time_tolerance,
                   "M4.3 MOOSE and fuelsim history times match") &&
             passed;
    passed = check(maximum_coordinate_difference < 1.0e-12,
                   "M4.3 MOOSE and fuelsim element centroids match") &&
             passed;
    std::cout << "m43_compared_element_time_rows=" << reference_row << '\n';
    std::cout << "m43_maximum_element_coordinate_difference="
              << maximum_coordinate_difference << '\n';
    std::cout << "m43_maximum_plastic_trace="
              << observer.maximum_plastic_trace()
              << '\n';
    std::cout << "m43_maximum_creep_trace="
              << observer.maximum_creep_trace() << '\n';
    passed = check(observer.maximum_plastic_trace() <
                           variant.plastic_trace_tolerance &&
                       observer.maximum_creep_trace() <
                           variant.creep_trace_tolerance,
                   "M4.3 default-Rashid accumulated trace drift stays below "
                   "its qualified limits") &&
             passed;
    const std::string prefix = "m43_" + variant.name + "_element_qp_average_";
    print_tensor_metric_locations(prefix + "stress", stress, reference);
    print_tensor_metric_locations(prefix + "elastic_strain", elastic,
                                  reference);
    print_tensor_metric_locations(
        prefix + "combined_inelastic_strain", combined_inelastic, reference);
    passed = check_metrics(prefix + "stress", stress, 1.0e-3,
                           variant.comparison_tolerance,
                           variant.stress_pointwise_tolerance) &&
             passed;
    passed = check_metrics(prefix + "elastic_strain", elastic, 1.0e-12,
                           variant.comparison_tolerance,
                           variant.elastic_pointwise_tolerance) &&
             passed;
    if (variant.behavior == ExpectedBehavior::elastic) {
        fuelsim::test::print_absolute_metrics(prefix + "combined_inelastic_strain",
                                              combined_inelastic);
        passed = check(combined_inelastic.maximum_absolute_difference <
                               1.0e-14 &&
                           combined_inelastic.maximum_actual < 1.0e-14 &&
                           combined_inelastic.maximum_reference < 1.0e-14,
                       "M4.3 elastic combined inelastic history stays zero") &&
                 passed;
    } else
        passed = check_metrics(prefix + "combined_inelastic_strain",
                               combined_inelastic, 1.0e-12,
                               variant.comparison_tolerance,
                               variant.inelastic_pointwise_tolerance) &&
                 passed;
    if (plastic_active(variant.behavior))
        passed = check_metrics(prefix + "equivalent_plastic",
                               equivalent_plastic, 1.0e-12,
                               variant.comparison_tolerance,
                               variant.comparison_tolerance) &&
                 passed;
    else {
        fuelsim::test::print_absolute_metrics(prefix + "equivalent_plastic",
                                              equivalent_plastic);
        passed = check(equivalent_plastic.maximum_absolute_difference <
                               1.0e-14 &&
                           equivalent_plastic.maximum_actual < 1.0e-14 &&
                           equivalent_plastic.maximum_reference < 1.0e-14,
                       "M4.3 inactive equivalent plastic history stays zero") &&
                 passed;
    }
    if (creep_active(variant.behavior))
        passed = check_metrics(prefix + "equivalent_creep", equivalent_creep,
                               1.0e-12, variant.comparison_tolerance,
                               variant.comparison_tolerance) &&
                 passed;
    else {
        fuelsim::test::print_absolute_metrics(prefix + "equivalent_creep",
                                              equivalent_creep);
        passed = check(equivalent_creep.maximum_absolute_difference <
                               1.0e-14 &&
                           equivalent_creep.maximum_actual < 1.0e-14 &&
                           equivalent_creep.maximum_reference < 1.0e-14,
                       "M4.3 inactive equivalent creep history stays zero") &&
                 passed;
    }
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_m43_noncoaxial_finite_strain_moose_tests "
                     "<variant> <case.fsi> <all_nodes.csv> <history.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv,
            "fuelsim M4.3 noncoaxial finite-strain MOOSE comparison\n");
        return run_test(variant_config(argv[1]), argv[2], argv[3], argv[4])
                   ? 0
                   : 1;
    } catch (const std::exception& error) {
        std::cerr << "M4.3 noncoaxial finite-strain test failed: "
                  << error.what() << '\n';
        return 1;
    }
}
