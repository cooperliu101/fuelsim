#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double preload_time = 1.0e-9;
constexpr double time_step = 0.1;
constexpr double applied_stress = 3.0e8;

const char* case_id(bool finite_strain) { return finite_strain ? "B5.16" : "B5.11"; }

const char* metric_prefix(bool finite_strain) { return finite_strain ? "b516_" : "b511_"; }

struct NodeReference final {
    std::size_t stage, node;
    double time;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::size_t stage;
    double time;
    fuelsim::CartesianPoint3 position;
    double temperature;
    fuelsim::SymmetricTensor3Values stress, strain, elastic_strain, creep_strain;
    double equivalent_creep_strain, integration_volume;
};

struct EnergyReference final {
    std::size_t stage;
    double time, internal_energy, creep_dissipation, external_work;
};

struct StepSnapshot final {
    double time;
    std::vector<double> state;
    std::array<fuelsim::CartesianMaterialPointState, 8> material;
    fuelsim::TransientConservationSummary conservation;
};

class SnapshotObserver final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        StepSnapshot snapshot{step.time, problem.committed_solution(), {}, step.conservation};
        const fuelsim::CartesianMaterialHistory& history =
            fuelsim::cartesian::ProblemAccess::material_history(problem, 0, 0);
        if (history.size() != snapshot.material.size())
            throw std::invalid_argument("B5.11 material history does not contain eight points");
        std::copy(history.begin(), history.end(), snapshot.material.begin());
        _snapshots.push_back(std::move(snapshot));
    }

    const std::vector<StepSnapshot>& snapshots() const noexcept { return _snapshots; }

  private:
    std::vector<StepSnapshot> _snapshots;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.11 row in " + path);
    return std::stod(values[index]);
}

std::size_t positive_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0)
        throw std::invalid_argument("Abaqus B5.11 index is not a positive integer in " + path);
    return static_cast<std::size_t>(rounded);
}

double zero_noise(double value, double threshold) { return std::abs(value) < threshold ? 0.0 : value; }

fuelsim::SymmetricTensor3Values tensor(
    const std::vector<std::string>& values, std::size_t first, const std::string& path, double zero_threshold) {
    fuelsim::SymmetricTensor3Values result = {number(values, first, path), number(values, first + 1, path),
        number(values, first + 2, path), 0.5 * number(values, first + 3, path), 0.5 * number(values, first + 5, path),
        0.5 * number(values, first + 4, path)};
    for (double* value : {&result.xx, &result.yy, &result.zz, &result.xy, &result.yz, &result.xz})
        *value = zero_noise(*value, zero_threshold);
    return result;
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.11 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "stage,time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.11 nodal header in " + path);
    std::vector<NodeReference> result;
    double maximum_active_reaction = 0.0, maximum_raw_transverse_reaction = 0.0;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 11) throw std::invalid_argument("Unexpected Abaqus B5.11 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 3, path);
            if (field == 5) maximum_active_reaction = std::max(maximum_active_reaction, std::abs(fields[field]));
            if (field == 6 || field == 7) {
                maximum_raw_transverse_reaction = std::max(maximum_raw_transverse_reaction, std::abs(fields[field]));
                fields[field] = 0.0;
            } else if (field != 0)
                fields[field] = zero_noise(fields[field], field >= 5 ? 1.0e-4 : (field == 4 ? 1.0e-10 : 1.0e-16));
        }
        result.push_back({positive_integer(number(values, 0, path), path),
            positive_integer(number(values, 2, path), path), number(values, 1, path), fields});
    }
    if (result.size() != 80) throw std::invalid_argument("Abaqus B5.11 nodal reference must contain 80 rows");
    const double relative_residual = maximum_raw_transverse_reaction / maximum_active_reaction;
    std::cout << "abaqus_raw_transverse_reaction_maximum_absolute=" << maximum_raw_transverse_reaction << '\n'
              << "abaqus_raw_transverse_reaction_relative_active_scale=" << relative_residual << '\n';
    if (!(maximum_active_reaction > 0.0 && relative_residual < 1.0e-7))
        throw std::invalid_argument("Abaqus transverse reaction residual is too large to canonicalize as zero");
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path, bool finite_strain) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.11 integration points: " + path);
    std::string line;
    std::getline(input, line);
    const std::string strain_header = finite_strain
                                          ? "le11,le22,le33,le12_engineering,le13_engineering,le23_engineering,"
                                          : "e11,e22,e33,e12_engineering,e13_engineering,e23_engineering,";
    if (line != "stage,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,"
                "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa," +
                    strain_header +
                    "ee11,ee22,ee33,ee12_engineering,ee13_engineering,ee23_engineering,"
                    "ce11,ce22,ce33,ce12_engineering,ce13_engineering,ce23_engineering,ceeq,ivol_m3")
        throw std::invalid_argument("Unexpected Abaqus B5.11 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    double maximum_active_stress = 0.0, maximum_raw_transverse_stress = 0.0;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 34)
            throw std::invalid_argument("Unexpected Abaqus B5.11 integration-point column count in " + path);
        if (positive_integer(number(values, 2, path), path) != 1)
            throw std::invalid_argument("Abaqus B5.11 reference must contain one element");
        fuelsim::SymmetricTensor3Values stress = {number(values, 8, path), number(values, 9, path),
            number(values, 10, path), number(values, 11, path), number(values, 13, path), number(values, 12, path)};
        maximum_active_stress = std::max(maximum_active_stress, std::abs(stress.xx));
        maximum_raw_transverse_stress =
            std::max({maximum_raw_transverse_stress, std::abs(stress.yy), std::abs(stress.zz)});
        for (double* value : {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz})
            *value = zero_noise(*value, 1.0e-3);
        stress.yy = 0.0;
        stress.zz = 0.0;
        result.push_back({positive_integer(number(values, 0, path), path), number(values, 1, path),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)}, number(values, 7, path),
            stress, tensor(values, 14, path, 1.0e-15), tensor(values, 20, path, 1.0e-15),
            tensor(values, 26, path, 1.0e-15), number(values, 32, path), number(values, 33, path)});
    }
    if (result.size() != 80)
        throw std::invalid_argument("Abaqus B5.11 integration-point reference must contain 80 rows");
    const double relative_residual = maximum_raw_transverse_stress / maximum_active_stress;
    std::cout << "abaqus_raw_transverse_stress_maximum_absolute=" << maximum_raw_transverse_stress << '\n'
              << "abaqus_raw_transverse_stress_relative_active_scale=" << relative_residual << '\n';
    if (!(maximum_active_stress > 0.0 && relative_residual < 1.0e-7))
        throw std::invalid_argument("Abaqus transverse stress residual is too large to canonicalize as zero");
    return result;
}

std::vector<EnergyReference> read_energy(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.11 energy: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "stage,time_s,internal_energy_j,creep_dissipation_j,external_work_j")
        throw std::invalid_argument("Unexpected Abaqus B5.11 energy header in " + path);
    std::vector<EnergyReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 5) throw std::invalid_argument("Unexpected Abaqus B5.11 energy column count in " + path);
        result.push_back({positive_integer(number(values, 0, path), path), number(values, 1, path),
            number(values, 2, path), number(values, 3, path), number(values, 4, path)});
    }
    if (result.size() != 10) throw std::invalid_argument("Abaqus B5.11 energy reference must contain ten rows");
    return result;
}

fuelsim::UnstructuredHex8Mesh mesh() {
    return fuelsim::UnstructuredHex8Mesh({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
                                             {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}},
        {{{{0, 1, 2, 3, 4, 5, 6, 7}}}}, {1}, {{1, "solid"}}, {},
        {{11, "left", {{0, 3}}}, {12, "right", {{0, 1}}}, {13, "y0", {{0, 0}}}, {14, "z0", {{0, 4}}}});
}

fuelsim::ThermoelasticProperties material() {
    return fuelsim::test::with_norton(
        fuelsim::test::thermoelastic(0.0, 1.0, 2.0e11, 0.3, 0.0, 600.0, 0.0, 0.0, 0.0, 1.0, 1.0), 1.0e-4, 1.0e8, 3.0,
        600.0);
}

fuelsim::SpatialDefinition definition(bool finite_strain) {
    fuelsim::SpatialDefinition result;
    result.regions.push_back({"solid", "solid", material(), 0.0, 600.0, -1, "",
        finite_strain ? fuelsim::StrainFormulation::finite : fuelsim::StrainFormulation::small});
    result.boundary_conditions = {
        {"temperature", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::temperature, 600.0},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0},
        {"pull_x", fuelsim::BoundaryConditionType::traction, "right", fuelsim::Field::displacement_x, applied_stress,
            true},
    };
    if (finite_strain) {
        result.boundary_conditions.back().use_displaced_geometry = false;
        result.boundary_conditions.back().configuration_explicit = true;
    }
    return result;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = 1.0e-10;
    options.relative_tolerance = 1.0e-12;
    options.maximum_iterations = 40;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    options.field_residual_scaling = true;
    options.residual_reduction_tolerance = 1.0e-11;
    options.temperature_residual_absolute_tolerance = 1.0e-10;
    options.mechanical_residual_absolute_tolerance = 1.0e-2;
    return options;
}

fuelsim::TransientTimeOptions preload_options() {
    return {preload_time, preload_time, preload_time, preload_time, 1.0, 0.5, 0, preload_time};
}

fuelsim::TransientTimeOptions hold_options() {
    return {preload_time + 10.0 * time_step, time_step, time_step, time_step, 1.0, 0.5, 0, 0.0};
}

std::size_t stage_from_time(double time) {
    const double value = (time - preload_time) / time_step;
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0 || rounded > 10.0)
        throw std::invalid_argument("Fuelsim B5.11 time does not identify one of ten hold stages");
    return static_cast<std::size_t>(rounded);
}

std::vector<double> raw_residual(fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    std::vector<double> result(problem.dof_count(), 0.0);
    fuelsim::ContributionWorkspace workspace;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        problem.evaluate_contribution(contribution, state, workspace, false);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            result[workspace.dofs[local]] += workspace.residual[local];
    }
    return result;
}

bool metrics_pass(const fuelsim::test::FieldErrorMetrics& metrics, double relative_tolerance, double zero_tolerance) {
    if (metrics.has_relative_norm() && !fuelsim::test::relative_metrics_below(metrics, relative_tolerance))
        return false;
    return metrics.maximum_zero_reference_difference < zero_tolerance;
}

void print_metrics(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics) {
    if (metrics.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, metrics);
    else
        fuelsim::test::print_absolute_metrics(name, metrics);
}

std::array<double, 6> components(const fuelsim::SymmetricTensor3Values& value) {
    return {value.xx, value.yy, value.zz, value.xy, value.yz, value.xz};
}

fuelsim::SymmetricTensor3Values logarithmic_strain(
    const fuelsim::Hex8QuadraturePoint& point, const fuelsim::Hex8LocalValues& state) {
    double f00 = 1.0, f01 = 0.0, f10 = 0.0, f11 = 1.0, f22 = 1.0;
    for (std::size_t node = 0; node < 8; ++node) {
        f00 += state[8 + node] * point.gradient[node][0];
        f01 += state[8 + node] * point.gradient[node][1];
        f10 += state[16 + node] * point.gradient[node][0];
        f11 += state[16 + node] * point.gradient[node][1];
        f22 += state[24 + node] * point.gradient[node][2];
    }
    const double b00 = f00 * f00 + f01 * f01;
    const double b01 = f00 * f10 + f01 * f11;
    const double b11 = f10 * f10 + f11 * f11;
    const double mean = 0.5 * (b00 + b11);
    const double radius = std::hypot(0.5 * (b00 - b11), b01);
    const double first = mean + radius, second = mean - radius;
    if (!(first > 0.0) || !(second > 0.0) || !(f22 > 0.0))
        throw std::domain_error("Finite Norton comparison requires a positive left stretch tensor");
    if (radius <= 1.0e-14 * mean) return {0.5 * std::log(mean), 0.5 * std::log(mean), std::log(f22), 0.0, 0.0, 0.0};
    const double beta = (std::log(first) - std::log(second)) / (first - second);
    const double alpha = (first * std::log(second) - second * std::log(first)) / (first - second);
    return {0.5 * (alpha + beta * b00), 0.5 * (alpha + beta * b11), std::log(f22), 0.5 * beta * b01, 0.0, 0.0};
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests "
                     "<small|finite> <nodes.csv> <integration.csv> <energy.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim Abaqus Norton comparison\n");
        const std::string branch = argv[1];
        if (branch != "small" && branch != "finite")
            throw std::invalid_argument("Abaqus Norton comparison branch must be small or finite");
        const bool finite_strain = branch == "finite";
        const std::vector<NodeReference> node_reference = read_nodes(argv[2]);
        const std::vector<IntegrationReference> integration_reference = read_integration(argv[3], finite_strain);
        const std::vector<EnergyReference> energy_reference = read_energy(argv[4]);
        const fuelsim::UnstructuredHex8Mesh input_mesh = mesh();
        fuelsim::TransientProblem problem(definition(finite_strain), input_mesh);
        const fuelsim::TransientResult preload = fuelsim::solve_transient(problem, preload_options(), solver_options());
        SnapshotObserver observer;
        const fuelsim::TransientResult solve =
            fuelsim::solve_transient(problem, hold_options(), solver_options(), &observer);
        bool passed =
            check(preload.completed && preload.accepted_steps.size() == 1 && preload.rejected_steps.empty(),
                std::string(case_id(finite_strain)) + " establishes the constant-force preload without cutback") &&
            check(solve.completed && observer.snapshots().size() == 10 && solve.accepted_steps.size() == 10 &&
                      solve.rejected_steps.empty(),
                std::string(case_id(finite_strain)) +
                    " accepts the ten constant-force creep increments without cutback");
        std::map<std::size_t, const StepSnapshot*> snapshots;
        for (const StepSnapshot& snapshot : observer.snapshots())
            snapshots.emplace(stage_from_time(snapshot.time), &snapshot);

        std::array<fuelsim::test::FieldErrorMetrics, 8> nodal_metrics;
        fuelsim::TransientProblem reaction_problem(definition(finite_strain), input_mesh);
        const fuelsim::TransientResult reaction_preload =
            fuelsim::solve_transient(reaction_problem, preload_options(), solver_options());
        passed = check(reaction_preload.completed,
                     std::string(case_id(finite_strain)) + " reaction replay establishes the same preload") &&
                 passed;
        const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(reaction_problem);
        const std::array<fuelsim::Field, 4> fields = {fuelsim::Field::temperature, fuelsim::Field::displacement_x,
            fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};
        for (std::size_t stage = 1; stage <= 10; ++stage) {
            const StepSnapshot& snapshot = *snapshots.at(stage);
            reaction_problem.begin_time_step({snapshot.time, 1.0, true});
            const std::vector<double> reaction = raw_residual(reaction_problem, snapshot.state);
            for (const NodeReference& reference : node_reference) {
                if (reference.stage != stage) continue;
                if (reference.node < 1 || reference.node > 8 ||
                    std::abs(reference.time - (snapshot.time - preload_time)) > 1.0e-12)
                    throw std::invalid_argument("Abaqus B5.11 nodal label or time lies outside the hold path");
                const std::size_t node = reference.node - 1;
                for (std::size_t field = 0; field < fields.size(); ++field) {
                    const std::size_t dof = dofs.dof(fields[field], node);
                    nodal_metrics[field].add(snapshot.state[dof], reference.fields[field]);
                    nodal_metrics[4 + field].add(reaction[dof], reference.fields[4 + field]);
                }
            }
            reaction_problem.commit_time_step(snapshot.state);
            passed = check(snapshot.conservation.creep_dissipation_increment >= -1.0e-8,
                         std::string(case_id(finite_strain)) + " creep dissipation is nonnegative at stage " +
                             std::to_string(stage)) &&
                     passed;
        }
        const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y",
            "displacement_z", "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
        const std::array<double, 8> nodal_zero_tolerances = {
            1.0e-12, 1.0e-14, 1.0e-14, 1.0e-14, 1.0e-10, 1.0e-2, 1.0e-2, 1.0e-2};
        for (std::size_t field = 0; field < nodal_metrics.size(); ++field) {
            print_metrics(metric_prefix(finite_strain) + nodal_names[field], nodal_metrics[field]);
            passed = check(metrics_pass(nodal_metrics[field], 1.0e-3, nodal_zero_tolerances[field]),
                         std::string(case_id(finite_strain)) + " " + nodal_names[field] +
                             " metrics are below the acceptance limits") &&
                     passed;
        }

        std::array<fuelsim::test::FieldErrorMetrics, 27> integration_metrics;
        double maximum_coordinate_error = 0.0;
        const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry({{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}}});
        for (const IntegrationReference& reference : integration_reference) {
            const StepSnapshot& snapshot = *snapshots.at(reference.stage);
            std::size_t closest = 0;
            double closest_squared = std::numeric_limits<double>::max();
            for (std::size_t q = 0; q < 8; ++q) {
                fuelsim::CartesianPoint3 point = geometry.points[q].position;
                if (finite_strain)
                    for (std::size_t node = 0; node < 8; ++node) {
                        point.x += geometry.points[q].shape[node] *
                                   snapshot.state[dofs.dof(fuelsim::Field::displacement_x, node)];
                        point.y += geometry.points[q].shape[node] *
                                   snapshot.state[dofs.dof(fuelsim::Field::displacement_y, node)];
                        point.z += geometry.points[q].shape[node] *
                                   snapshot.state[dofs.dof(fuelsim::Field::displacement_z, node)];
                    }
                const double distance_squared = std::pow(point.x - reference.position.x, 2) +
                                                std::pow(point.y - reference.position.y, 2) +
                                                std::pow(point.z - reference.position.z, 2);
                if (distance_squared < closest_squared) {
                    closest = q;
                    closest_squared = distance_squared;
                }
            }
            maximum_coordinate_error = std::max(maximum_coordinate_error, std::sqrt(closest_squared));
            fuelsim::Hex8LocalAdValues passive{};
            fuelsim::Hex8LocalValues passive_values{};
            for (std::size_t node = 0; node < 8; ++node) {
                passive[node] = snapshot.state[dofs.dof(fuelsim::Field::temperature, node)];
                passive[8 + node] = snapshot.state[dofs.dof(fuelsim::Field::displacement_x, node)];
                passive[16 + node] = snapshot.state[dofs.dof(fuelsim::Field::displacement_y, node)];
                passive[24 + node] = snapshot.state[dofs.dof(fuelsim::Field::displacement_z, node)];
            }
            for (std::size_t local = 0; local < passive.size(); ++local) passive_values[local] = passive[local].value();
            const fuelsim::CartesianMaterialPointState& actual = snapshot.material[closest];
            const std::array<double, 6> actual_stress = components(actual.stress);
            const std::array<double, 6> expected_stress = components(reference.stress);
            std::array<double, 6> actual_total{};
            if (finite_strain)
                actual_total = components(logarithmic_strain(geometry.points[closest], passive_values));
            else {
                double volume = 0.0, average_trace = 0.0;
                for (const fuelsim::Hex8QuadraturePoint& point : geometry.points) {
                    const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
                        point, passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::small);
                    volume += point.weighted_measure;
                    average_trace += point.weighted_measure *
                                     (kinematics.strain_increment.xx.value() + kinematics.strain_increment.yy.value() +
                                         kinematics.strain_increment.zz.value());
                }
                average_trace /= volume;
                const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
                    geometry.points[closest], passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::small);
                fuelsim::SymmetricTensor3 actual_strain = kinematics.strain_increment;
                const adlite::Scalar correction =
                    (average_trace - actual_strain.xx - actual_strain.yy - actual_strain.zz) / 3.0;
                actual_strain.xx += correction;
                actual_strain.yy += correction;
                actual_strain.zz += correction;
                actual_total = {actual_strain.xx.value(), actual_strain.yy.value(), actual_strain.zz.value(),
                    actual_strain.xy.value(), actual_strain.yz.value(), actual_strain.xz.value()};
            }
            const std::array<double, 6> expected_total = components(reference.strain);
            const std::array<double, 6> expected_elastic = components(reference.elastic_strain);
            const std::array<double, 6> expected_creep = components(reference.creep_strain);
            for (std::size_t component = 0; component < 6; ++component) {
                integration_metrics[component].add(actual_stress[component], expected_stress[component]);
                integration_metrics[6 + component].add(actual_total[component], expected_total[component]);
                integration_metrics[12 + component].add(actual.elastic_strain[component], expected_elastic[component]);
                integration_metrics[18 + component].add(actual.creep_strain[component], expected_creep[component]);
            }
            integration_metrics[24].add(actual.equivalent_creep_strain, reference.equivalent_creep_strain);
            integration_metrics[25].add(passive[0].value(), reference.temperature);
            const double actual_measure =
                finite_strain ? fuelsim::evaluate_cartesian_incremental_kinematics(geometry.points[closest], passive,
                                    passive_values, fuelsim::StrainFormulation::finite)
                                    .current_weighted_measure.value()
                              : geometry.points[closest].weighted_measure;
            integration_metrics[26].add(actual_measure, reference.integration_volume);
        }
        const std::array<std::string, 27> integration_names = {"stress_xx", "stress_yy", "stress_zz", "stress_xy",
            "stress_yz", "stress_xz", "strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz",
            "elastic_strain_xx", "elastic_strain_yy", "elastic_strain_zz", "elastic_strain_xy", "elastic_strain_yz",
            "elastic_strain_xz", "creep_strain_xx", "creep_strain_yy", "creep_strain_zz", "creep_strain_xy",
            "creep_strain_yz", "creep_strain_xz", "equivalent_creep_strain", "material_temperature",
            "integration_volume"};
        for (std::size_t field = 0; field < integration_metrics.size(); ++field) {
            print_metrics(metric_prefix(finite_strain) + integration_names[field], integration_metrics[field]);
            const double zero_tolerance = field < 6 ? 1.0e-2 : (field == 26 ? 1.0e-14 : 1.0e-14);
            passed = check(metrics_pass(integration_metrics[field], 1.0e-3, zero_tolerance),
                         std::string(case_id(finite_strain)) + " " + integration_names[field] +
                             " metrics are below 0.1 percent") &&
                     passed;
        }
        passed = check(maximum_coordinate_error < (finite_strain ? 1.0e-10 : 1.0e-12),
                     std::string(case_id(finite_strain)) + " maps every integration point at all ten accepted times") &&
                 passed;
        std::cout << metric_prefix(finite_strain)
                  << "maximum_integration_coordinate_absolute_difference=" << maximum_coordinate_error << '\n';

        std::array<fuelsim::test::FieldErrorMetrics, 3> energy_metrics;
        double cumulative_elastic = preload.accepted_steps.front().conservation.elastic_energy_change;
        double cumulative_creep = preload.accepted_steps.front().conservation.creep_dissipation_increment;
        for (const EnergyReference& reference : energy_reference) {
            const StepSnapshot& snapshot = *snapshots.at(reference.stage);
            cumulative_elastic += snapshot.conservation.elastic_energy_change;
            cumulative_creep += snapshot.conservation.creep_dissipation_increment;
            energy_metrics[0].add(cumulative_elastic + cumulative_creep, reference.internal_energy);
            energy_metrics[1].add(cumulative_creep, reference.creep_dissipation);
            energy_metrics[2].add(cumulative_elastic + cumulative_creep, reference.external_work);
        }
        const std::array<std::string, 3> energy_names = {"internal_energy", "creep_dissipation", "external_work"};
        for (std::size_t field = 0; field < energy_metrics.size(); ++field) {
            print_metrics(metric_prefix(finite_strain) + energy_names[field], energy_metrics[field]);
            passed = check(metrics_pass(energy_metrics[field], 1.0e-3, 1.0e-8), std::string(case_id(finite_strain)) +
                                                                                    " " + energy_names[field] +
                                                                                    " metrics are below 0.1 percent") &&
                     passed;
        }
        passed = check(snapshots.at(1)->material[0].equivalent_creep_strain > 0.0 &&
                           snapshots.at(10)->material[0].equivalent_creep_strain >
                               9.9 * snapshots.at(1)->material[0].equivalent_creep_strain,
                     std::string(case_id(finite_strain)) +
                         " accumulates positive Norton creep throughout the ten constant-force holds") &&
                 passed;
        if (passed && session.rank() == 0)
            std::cout << "[PASS] " << case_id(finite_strain) << " Abaqus C3D8T constant-force Norton comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.11 Abaqus comparison raised: " << error.what() << '\n';
        return 1;
    }
}
