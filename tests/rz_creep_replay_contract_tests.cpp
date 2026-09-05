#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/jacobian_check.hpp"
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
#include <utility>
#include <vector>

namespace {
constexpr double time_tolerance = 1.0e-12;

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t end = line.find(',', begin);
        fields.push_back(line.substr(begin, end - begin));
        if (end == std::string::npos) return fields;
        begin = end + 1;
    }
}

std::size_t column_index(const std::vector<std::string>& header, const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end()) throw std::invalid_argument("M4.3 CSV is missing column '" + name + "'");
    return static_cast<std::size_t>(found - header.begin());
}

double csv_value(const std::vector<std::string>& fields, std::size_t column, const std::string& path) {
    if (column >= fields.size()) throw std::invalid_argument("M4.3 CSV row is incomplete: " + path);
    std::size_t parsed = 0;
    const double value = std::stod(fields[column], &parsed);
    if (parsed != fields[column].size() || !std::isfinite(value))
        throw std::invalid_argument("M4.3 CSV value is invalid: " + path);
    return value;
}

std::size_t csv_id(const std::vector<std::string>& fields, std::size_t column, const std::string& path) {
    const double value = csv_value(fields, column, path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(value) != value)
        throw std::invalid_argument("M4.3 CSV element ID is invalid: " + path);
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
    explicit HistoryObserver(std::size_t expected_element_count) : _expected_element_count(expected_element_count) {}

    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        if (fuelsim::rz::ProblemAccess::region_count(problem) != 1 ||
            fuelsim::rz::ProblemAccess::region_mesh(problem, 0).elements().size() != _expected_element_count)
            throw std::logic_error("M4.3 observer element count differs from its variant");
        const fuelsim::RegionMesh& mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, 0);
        const std::vector<double>& solution = problem.committed_solution();
        const double maximum_z = std::max_element(
            mesh.nodes().begin(), mesh.nodes().end(), [](const fuelsim::RzPoint& left, const fuelsim::RzPoint& right) {
                return left.z < right.z;
            })->z;
        const double minimum_z = std::min_element(
            mesh.nodes().begin(), mesh.nodes().end(), [](const fuelsim::RzPoint& left, const fuelsim::RzPoint& right) {
                return left.z < right.z;
            })->z;
        double radial = 0.0;
        double axial = 0.0;
        std::size_t count = 0;
        for (std::size_t node = 0; node < mesh.nodes().size(); ++node) {
            if (std::abs(mesh.nodes()[node].z - maximum_z) > time_tolerance) continue;
            radial += solution.at(
                fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement, node));
            axial +=
                solution.at(fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement, node));
            ++count;
        }
        if (count == 0) throw std::logic_error("M4.3 mesh has no top nodes");
        HistorySnapshot snapshot;
        snapshot.time = step.time;
        snapshot.reference_height = maximum_z - minimum_z;
        snapshot.top_radial_displacement = radial / static_cast<double>(count);
        snapshot.top_axial_displacement = axial / static_cast<double>(count);
        snapshot.elements.reserve(mesh.elements().size());
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            ElementSnapshot value;
            value.element_id = mesh.source_element_ids().at(element);
            for (const std::size_t node : mesh.elements()[element].nodes) {
                value.radius += mesh.nodes()[node].r / 4.0;
                value.axial_coordinate += mesh.nodes()[node].z / 4.0;
            }
            const fuelsim::Quad4RzGeometry& geometry =
                fuelsim::rz::ProblemAccess::region_element_geometry(problem, 0, element);
            const std::array<fuelsim::AxisymmetricStressValues, 4>& stresses =
                fuelsim::rz::ProblemAccess::material_stress(problem, 0, element);
            const fuelsim::Quad4MaterialHistory& history =
                fuelsim::rz::ProblemAccess::material_history(problem, 0, element);
            double total_weight = 0.0;
            for (std::size_t point = 0; point < geometry.points.size(); ++point) {
                const double weight = geometry.points[point].weighted_measure;
                total_weight += weight;
                value.stress.rr += weight * stresses[point].rr;
                value.stress.zz += weight * stresses[point].zz;
                value.stress.hoop += weight * stresses[point].hoop;
                value.stress.rz += weight * stresses[point].rz;
                for (std::size_t component = 0; component < 4; ++component) {
                    value.state.elastic_strain[component] += weight * history[point].elastic_strain[component];
                    value.state.plastic_strain[component] += weight * history[point].plastic_strain[component];
                    value.state.creep_strain[component] += weight * history[point].creep_strain[component];
                }
                value.state.equivalent_plastic_strain += weight * history[point].equivalent_plastic_strain;
                value.state.equivalent_creep_strain += weight * history[point].equivalent_creep_strain;
                _maximum_plastic_trace = std::max(_maximum_plastic_trace,
                    std::abs(history[point].plastic_strain[0] + history[point].plastic_strain[1] +
                             history[point].plastic_strain[2]));
                _maximum_creep_trace = std::max(
                    _maximum_creep_trace, std::abs(history[point].creep_strain[0] + history[point].creep_strain[1] +
                                                   history[point].creep_strain[2]));
            }
            if (!(total_weight > 0.0)) throw std::logic_error("M4.3 element has zero reference volume");
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
            [](const ElementSnapshot& left, const ElementSnapshot& right) {
                return left.element_id < right.element_id;
            });
        _snapshots.push_back(std::move(snapshot));
    }

    const std::vector<HistorySnapshot>& snapshots() const noexcept { return _snapshots; }

    double maximum_plastic_trace() const noexcept { return _maximum_plastic_trace; }

    double maximum_creep_trace() const noexcept { return _maximum_creep_trace; }

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

struct NodalHistoryEntry final {
    double time = 0.0;
    std::size_t node_id = 0;
    fuelsim::test::NodalFieldReference field{};
};

struct NodalHistorySnapshot final {
    double time = 0.0;
    std::vector<fuelsim::test::NodalFieldReference> nodes;
};

std::vector<ReferenceSnapshot> read_reference_history(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read M4.3 MOOSE history: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("M4.3 MOOSE history is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const auto column = [&header](const std::string& name) { return column_index(header, name); };
    const std::size_t time = column("sample_time");
    const std::size_t id = column("id");
    const std::size_t radius = column("x");
    const std::size_t axial_coordinate = column("y");
    const std::array<std::size_t, 4> stress = {
        column("stress_rr"), column("stress_zz"), column("stress_hoop"), column("stress_rz")};
    const std::array<std::size_t, 4> elastic = {
        column("elastic_rr"), column("elastic_zz"), column("elastic_hoop"), column("elastic_rz")};
    const std::array<std::size_t, 4> combined_inelastic = {
        column("combined_rr"), column("combined_zz"), column("combined_hoop"), column("combined_rz")};
    const std::size_t equivalent_plastic = column("effective_plastic");
    const std::size_t equivalent_creep = column("effective_creep");
    std::vector<ReferenceSnapshot> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_csv(line);
        ReferenceSnapshot value;
        value.time = csv_value(fields, time, path);
        value.element_id = csv_id(fields, id, path);
        value.radius = csv_value(fields, radius, path);
        value.axial_coordinate = csv_value(fields, axial_coordinate, path);
        for (std::size_t component = 0; component < 4; ++component) {
            value.stress[component] = csv_value(fields, stress[component], path);
            value.elastic[component] = csv_value(fields, elastic[component], path);
            value.combined_inelastic[component] = csv_value(fields, combined_inelastic[component], path);
        }
        value.equivalent_plastic = csv_value(fields, equivalent_plastic, path);
        value.equivalent_creep = csv_value(fields, equivalent_creep, path);
        result.push_back(value);
    }
    if (result.empty()) throw std::invalid_argument("M4.3 MOOSE history has no transient rows: " + path);
    std::sort(result.begin(), result.end(), [](const ReferenceSnapshot& left, const ReferenceSnapshot& right) {
        if (left.time != right.time) return left.time < right.time;
        return left.element_id < right.element_id;
    });
    return result;
}

std::vector<NodalHistorySnapshot> read_nodal_history(const std::string& path, std::size_t node_count) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read M4.3 MOOSE nodal history: " + path);
    std::string line;
    if (!std::getline(input, line)) throw std::invalid_argument("M4.3 MOOSE nodal history is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const auto column = [&header](const std::string& name) { return column_index(header, name); };
    const std::size_t time = column("sample_time_nodal");
    const std::size_t id = column("id");
    const std::size_t radius = column("x");
    const std::size_t axial_coordinate = column("y");
    const std::size_t temperature = column("T");
    const std::size_t radial_displacement = column("disp_x");
    const std::size_t axial_displacement = column("disp_y");
    std::vector<NodalHistoryEntry> entries;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split_csv(line);
        NodalHistoryEntry entry;
        entry.time = csv_value(fields, time, path);
        entry.node_id = csv_id(fields, id, path);
        entry.field = {csv_value(fields, radius, path), csv_value(fields, axial_coordinate, path),
            csv_value(fields, temperature, path), csv_value(fields, radial_displacement, path),
            csv_value(fields, axial_displacement, path)};
        entries.push_back(entry);
    }
    if (entries.empty()) throw std::invalid_argument("M4.3 MOOSE nodal history has no transient rows: " + path);
    std::sort(entries.begin(), entries.end(), [](const NodalHistoryEntry& left, const NodalHistoryEntry& right) {
        if (left.time != right.time) return left.time < right.time;
        return left.node_id < right.node_id;
    });
    std::vector<NodalHistorySnapshot> result;
    std::vector<bool> present;
    for (const NodalHistoryEntry& entry : entries) {
        if (result.empty() || entry.time != result.back().time) {
            if (!result.empty() && std::any_of(present.begin(), present.end(), [](bool value) { return !value; }))
                throw std::invalid_argument("M4.3 MOOSE nodal history misses a node: " + path);
            result.push_back({entry.time, std::vector<fuelsim::test::NodalFieldReference>(node_count)});
            present.assign(node_count, false);
        }
        if (entry.node_id >= node_count || present[entry.node_id])
            throw std::invalid_argument("M4.3 MOOSE nodal history has an invalid or duplicate node "
                                        "ID: " +
                                        path);
        result.back().nodes[entry.node_id] = entry.field;
        present[entry.node_id] = true;
    }
    if (std::any_of(present.begin(), present.end(), [](bool value) { return !value; }))
        throw std::invalid_argument("M4.3 MOOSE nodal history misses a node: " + path);
    return result;
}

std::array<double, 4> stress_components(const fuelsim::AxisymmetricStressValues& stress) {
    return {stress.rr, stress.zz, stress.hoop, stress.rz};
}

bool check_metrics(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics,
    double zero_reference_tolerance, double aggregate_tolerance, double pointwise_tolerance) {
    fuelsim::test::print_relative_metrics(name, metrics);
    std::cout << name << "_maximum_absolute_difference=" << metrics.maximum_absolute_difference << '\n';
    std::cout << name << "_relative_l2_tolerance=" << aggregate_tolerance << '\n';
    std::cout << name << "_relative_absolute_peak_tolerance=" << aggregate_tolerance << '\n';
    std::cout << name << "_maximum_pointwise_relative_tolerance=" << pointwise_tolerance << '\n';
    return check(metrics.relative_l2() < aggregate_tolerance &&
                     metrics.relative_absolute_peak() < aggregate_tolerance &&
                     metrics.maximum_pointwise_relative_error() < pointwise_tolerance &&
                     metrics.maximum_zero_reference_difference < zero_reference_tolerance,
        name + " three MOOSE metrics and zero-reference error pass");
}

void print_tensor_metric_locations(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics,
    const std::vector<ReferenceSnapshot>& reference) {
    constexpr std::array<const char*, 4> components = {"rr", "zz", "hoop", "rz"};
    const auto print_location = [&name, &reference, &components](const std::string& metric, std::size_t flat_index) {
        const std::size_t row = flat_index / components.size();
        const std::size_t component = flat_index % components.size();
        if (row >= reference.size()) throw std::logic_error("M4.3 metric index exceeds history size");
        std::cout << name << '_' << metric << "_time=" << reference[row].time << '\n';
        std::cout << name << '_' << metric << "_element_id=" << reference[row].element_id << '\n';
        std::cout << name << '_' << metric << "_component=" << components[component] << '\n';
    };
    print_location("maximum_pointwise_relative_location", metrics.maximum_pointwise_relative_index);
    print_location("maximum_absolute_difference_location", metrics.maximum_absolute_difference_index);
}

bool audit_creep_shared_state(const fuelsim::FuelSimCaseDefinition& definition,
    const fuelsim::UnstructuredQuad4Mesh& source, const std::vector<ReferenceSnapshot>& reference,
    const std::string& nodal_history_path) {
    if (definition.spatial.regions.size() != 1 || !definition.spatial.regions[0].material.functions->has_creep() ||
        definition.spatial.regions[0].material.functions->has_plasticity() || !definition.spatial.contacts.empty() ||
        std::any_of(definition.spatial.boundary_conditions.begin(), definition.spatial.boundary_conditions.end(),
            [](const fuelsim::BoundaryConditionDefinition& boundary) {
                return boundary.type != fuelsim::BoundaryConditionType::dirichlet;
            }))
        throw std::invalid_argument("M4.3 creep shared-state audit requires one Norton-only region "
                                    "with displacement boundary conditions only");
    fuelsim::TransientProblem problem(definition.spatial, source);
    if (fuelsim::rz::ProblemAccess::region_count(problem) != 1)
        throw std::invalid_argument("M4.3 creep shared-state audit requires one region");
    const fuelsim::RegionMesh& mesh = fuelsim::rz::ProblemAccess::region_mesh(problem, 0);
    const std::vector<NodalHistorySnapshot> nodal_history =
        read_nodal_history(nodal_history_path, source.nodes().size());
    if (reference.size() != nodal_history.size() * mesh.elements().size())
        throw std::invalid_argument("M4.3 creep shared-state node and element histories differ in "
                                    "length");
    fuelsim::test::FieldErrorMetrics stress;
    fuelsim::test::FieldErrorMetrics elastic;
    fuelsim::test::FieldErrorMetrics combined_inelastic;
    fuelsim::test::FieldErrorMetrics equivalent_creep;
    double maximum_coordinate_difference = 0.0;
    double maximum_boundary_value_difference = 0.0;
    double maximum_free_residual_l2 = 0.0;
    double maximum_free_residual_l2_ratio = 0.0;
    double maximum_free_residual_infinity = 0.0;
    double maximum_free_residual_infinity_ratio = 0.0;
    double maximum_free_residual_ratio_time = 0.0;
    std::size_t reference_row = 0;
    for (const NodalHistorySnapshot& snapshot : nodal_history) {
        std::vector<double> state = problem.committed_solution();
        const std::size_t offset = fuelsim::rz::ProblemAccess::region_node_offset(problem, 0);
        for (std::size_t local = 0; local < mesh.nodes().size(); ++local) {
            const std::size_t source_node = mesh.source_node_ids().at(local);
            const fuelsim::test::NodalFieldReference& field = snapshot.nodes.at(source_node);
            maximum_coordinate_difference =
                std::max({maximum_coordinate_difference, std::abs(mesh.nodes()[local].r - field.radius),
                    std::abs(mesh.nodes()[local].z - field.axial_coordinate)});
            const std::size_t global_node = offset + local;
            state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature, global_node)] =
                field.temperature;
            state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement, global_node)] =
                field.radial_displacement;
            state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement, global_node)] =
                field.axial_displacement;
        }
        const double ramp_time = definition.transient_execution.load_ramp_time;
        const double load_factor = ramp_time > 0.0 ? std::min(snapshot.time / ramp_time, 1.0) : 1.0;
        problem.begin_time_step({snapshot.time, load_factor});
        std::vector<bool> constrained(problem.dof_count(), false);
        for (const fuelsim::DirichletCondition& condition : problem.dirichlet_conditions()) {
            constrained.at(condition.dof) = true;
            maximum_boundary_value_difference =
                std::max(maximum_boundary_value_difference, std::abs(state.at(condition.dof) - condition.value));
        }
        std::vector<double> residual;
        residual = fuelsim::test::assembled_residual(problem, state);
        double local_force_squared = 0.0;
        double local_force_infinity = 0.0;
        for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
            const fuelsim::LocalResidual local = fuelsim::rz::ProblemAccess::contribution_residual(
                problem, contribution, fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, state));
            for (std::size_t row = 4; row < local.size(); ++row) {
                local_force_squared += local[row] * local[row];
                local_force_infinity = std::max(local_force_infinity, std::abs(local[row]));
            }
        }
        double free_residual_squared = 0.0;
        double free_residual_infinity = 0.0;
        for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::dof_map(problem).node_count(); ++node) {
            const std::array<std::size_t, 2> mechanical_dofs = {
                fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::radial_displacement, node),
                fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement, node)};
            for (const std::size_t dof : mechanical_dofs) {
                if (constrained[dof]) continue;
                free_residual_squared += residual[dof] * residual[dof];
                free_residual_infinity = std::max(free_residual_infinity, std::abs(residual[dof]));
            }
        }
        const double free_residual_l2 = std::sqrt(free_residual_squared);
        const double local_force_l2 = std::sqrt(local_force_squared);
        if (!(local_force_l2 > 0.0) || !(local_force_infinity > 0.0))
            throw std::domain_error("M4.3 creep shared-state force scale is zero");
        const double l2_ratio = free_residual_l2 / local_force_l2;
        const double infinity_ratio = free_residual_infinity / local_force_infinity;
        if (l2_ratio > maximum_free_residual_l2_ratio) maximum_free_residual_ratio_time = snapshot.time;
        maximum_free_residual_l2 = std::max(maximum_free_residual_l2, free_residual_l2);
        maximum_free_residual_l2_ratio = std::max(maximum_free_residual_l2_ratio, l2_ratio);
        maximum_free_residual_infinity = std::max(maximum_free_residual_infinity, free_residual_infinity);
        maximum_free_residual_infinity_ratio = std::max(maximum_free_residual_infinity_ratio, infinity_ratio);
        problem.commit_time_step(state);
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            const ReferenceSnapshot& expected = reference.at(reference_row);
            if (std::abs(expected.time - snapshot.time) >= time_tolerance ||
                expected.element_id != mesh.source_element_ids().at(element))
                throw std::invalid_argument("M4.3 creep shared-state history ordering differs");
            const fuelsim::Quad4RzGeometry& geometry =
                fuelsim::rz::ProblemAccess::region_element_geometry(problem, 0, element);
            const std::array<fuelsim::AxisymmetricStressValues, 4>& stresses =
                fuelsim::rz::ProblemAccess::material_stress(problem, 0, element);
            const fuelsim::Quad4MaterialHistory& history =
                fuelsim::rz::ProblemAccess::material_history(problem, 0, element);
            fuelsim::AxisymmetricStressValues average_stress{};
            fuelsim::MaterialPointState average_state{};
            double total_weight = 0.0;
            for (std::size_t point = 0; point < geometry.points.size(); ++point) {
                const double weight = geometry.points[point].weighted_measure;
                total_weight += weight;
                average_stress.rr += weight * stresses[point].rr;
                average_stress.zz += weight * stresses[point].zz;
                average_stress.hoop += weight * stresses[point].hoop;
                average_stress.rz += weight * stresses[point].rz;
                for (std::size_t component = 0; component < 4; ++component) {
                    average_state.elastic_strain[component] += weight * history[point].elastic_strain[component];
                    average_state.plastic_strain[component] += weight * history[point].plastic_strain[component];
                    average_state.creep_strain[component] += weight * history[point].creep_strain[component];
                }
                average_state.equivalent_creep_strain += weight * history[point].equivalent_creep_strain;
            }
            if (!(total_weight > 0.0)) throw std::domain_error("M4.3 creep shared-state element volume is zero");
            std::array<double, 4> actual_stress = stress_components(average_stress);
            for (std::size_t component = 0; component < 4; ++component) {
                actual_stress[component] /= total_weight;
                average_state.elastic_strain[component] /= total_weight;
                average_state.plastic_strain[component] /= total_weight;
                average_state.creep_strain[component] /= total_weight;
                stress.add(actual_stress[component], expected.stress[component]);
                elastic.add(average_state.elastic_strain[component], expected.elastic[component]);
                combined_inelastic.add(average_state.plastic_strain[component] + average_state.creep_strain[component],
                    expected.combined_inelastic[component]);
            }
            average_state.equivalent_creep_strain /= total_weight;
            equivalent_creep.add(average_state.equivalent_creep_strain, expected.equivalent_creep);
            ++reference_row;
        }
    }
    const std::string prefix = "m43_creep_shared_state_";
    bool passed = check_metrics(prefix + "stress", stress, 1.0e-3, 1.0e-8, 1.0e-5);
    passed = check_metrics(prefix + "elastic_strain", elastic, 1.0e-12, 1.0e-8, 1.0e-5) && passed;
    passed = check_metrics(prefix + "combined_inelastic_strain", combined_inelastic, 1.0e-12, 1.0e-8, 1.0e-5) && passed;
    passed = check_metrics(prefix + "equivalent_creep", equivalent_creep, 1.0e-12, 1.0e-8, 1.0e-5) && passed;
    print_tensor_metric_locations(prefix + "stress", stress, reference);
    std::cout << prefix << "compared_node_time_rows=" << nodal_history.size() * source.nodes().size() << '\n';
    std::cout << prefix << "maximum_coordinate_difference=" << maximum_coordinate_difference << '\n';
    std::cout << prefix << "maximum_boundary_value_difference=" << maximum_boundary_value_difference << '\n';
    std::cout << prefix << "maximum_free_residual_l2=" << maximum_free_residual_l2 << '\n';
    std::cout << prefix << "maximum_free_residual_l2_ratio=" << maximum_free_residual_l2_ratio << '\n';
    std::cout << prefix << "maximum_free_residual_infinity=" << maximum_free_residual_infinity << '\n';
    std::cout << prefix << "maximum_free_residual_infinity_ratio=" << maximum_free_residual_infinity_ratio << '\n';
    std::cout << prefix << "maximum_free_residual_ratio_time=" << maximum_free_residual_ratio_time << '\n';
    passed = check(nodal_history.size() == 100 && reference_row == reference.size(),
                 "M4.3 creep shared-state audit covers all nodes, elements, "
                 "and time steps") &&
             passed;
    passed = check(maximum_coordinate_difference < 1.0e-12 && maximum_boundary_value_difference < 1.0e-12,
                 "M4.3 creep shared-state coordinates and prescribed values "
                 "match") &&
             passed;
    passed = check(maximum_free_residual_l2_ratio < 2.0e-5 && maximum_free_residual_infinity_ratio < 2.0e-5,
                 "M4.3 creep MOOSE states satisfy the fuelsim free-DOF "
                 "weak form within 0.002 percent of the local-force scale") &&
             passed;
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    try {
        const auto definition = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_quad4(definition.mesh_file);
        return audit_creep_shared_state(definition, mesh, read_reference_history(argv[2]), argv[3]) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M4.3 shared-state material contract: " << error.what() << '\n';
        return 1;
    }
}
