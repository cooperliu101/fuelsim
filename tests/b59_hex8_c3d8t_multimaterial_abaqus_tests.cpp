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
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr double length = 2.0, half_thickness = 0.1, width = 0.25, step_time = 1.0e7;

struct CaseSpec final {
    const char* name;
    std::size_t nx, subdivisions_per_layer, nz;
    bool distorted;
};

constexpr std::array<CaseSpec, 3> cases = {
    CaseSpec{"coarse", 2, 1, 1, false}, CaseSpec{"refined", 4, 2, 2, false}, CaseSpec{"distorted", 4, 2, 2, true}};

struct NodeReference final {
    std::string case_name;
    std::size_t stage, node;
    double time;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::string case_name;
    std::size_t stage, element;
    double time;
    fuelsim::CartesianPoint3 position;
    double temperature;
    std::array<double, 3> heat_flux;
    fuelsim::SymmetricTensor3Values stress, strain;
    double integration_volume;
};

struct StepSnapshot final {
    double time, load_factor;
    std::vector<double> state;
    fuelsim::TransientConservationSummary conservation;
};

struct MeshData final {
    fuelsim::UnstructuredHex8Mesh mesh;
    std::vector<std::size_t> interface_nodes, bottom_right_nodes, top_right_nodes;
    double maximum_distortion;
};

struct BendingSummary final {
    double peak, cooled;
};

class SnapshotObserver final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        _snapshots.push_back({step.time, step.load_factor, problem.committed_solution(), step.conservation});
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
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.9 row in " + path);
    return std::stod(values[index]);
}

std::size_t positive_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0)
        throw std::invalid_argument("Abaqus B5.9 index is not a positive integer in " + path);
    return static_cast<std::size_t>(rounded);
}

double canonical_zero(double value, double tolerance) { return std::abs(value) < tolerance ? 0.0 : value; }

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.9 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "case,stage,time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.9 nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 12) throw std::invalid_argument("Unexpected Abaqus B5.9 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 4, path);
            const double tolerance = field >= 5 ? 1.0e-5 : (field == 4 ? 1.0e-8 : 1.0e-16);
            if (field != 0) fields[field] = canonical_zero(fields[field], tolerance);
        }
        result.push_back({values[0], positive_integer(number(values, 1, path), path),
            positive_integer(number(values, 3, path), path), number(values, 2, path), fields});
    }
    if (result.size() != 672) throw std::invalid_argument("Abaqus B5.9 nodal reference must contain 672 rows");
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.9 integration points: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "case,stage,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,hfl_x_w_m2,hfl_y_w_m2,"
                "hfl_z_w_m2,s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,e11,e22,e33,e12_engineering,"
                "e13_engineering,e23_engineering,ivol_m3")
        throw std::invalid_argument("Unexpected Abaqus B5.9 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 25)
            throw std::invalid_argument("Unexpected Abaqus B5.9 integration-point column count in " + path);
        std::array<double, 3> heat_flux = {number(values, 9, path), number(values, 10, path), number(values, 11, path)};
        for (double& value : heat_flux) value = canonical_zero(value, 1.0e-8);
        fuelsim::SymmetricTensor3Values stress = {number(values, 12, path), number(values, 13, path),
            number(values, 14, path), number(values, 15, path), number(values, 17, path), number(values, 16, path)};
        fuelsim::SymmetricTensor3Values strain = {number(values, 18, path), number(values, 19, path),
            number(values, 20, path), 0.5 * number(values, 21, path), 0.5 * number(values, 23, path),
            0.5 * number(values, 22, path)};
        for (double* value : {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz})
            *value = canonical_zero(*value, 1.0e-4);
        for (double* value : {&strain.xx, &strain.yy, &strain.zz, &strain.xy, &strain.yz, &strain.xz})
            *value = canonical_zero(*value, 1.0e-15);
        result.push_back({values[0], positive_integer(number(values, 1, path), path),
            positive_integer(number(values, 3, path), path), number(values, 2, path),
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)}, number(values, 8, path),
            heat_flux, stress, strain, number(values, 24, path)});
    }
    if (result.size() != 2176)
        throw std::invalid_argument("Abaqus B5.9 integration-point reference must contain 2176 rows");
    return result;
}

MeshData make_mesh(const CaseSpec& spec) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    const std::size_t ny = 2 * spec.subdivisions_per_layer;
    std::vector<fuelsim::CartesianPoint3> nodes;
    double maximum_distortion = 0.0;
    for (std::size_t iz = 0; iz <= spec.nz; ++iz) {
        const double z = width * static_cast<double>(iz) / static_cast<double>(spec.nz);
        for (std::size_t iy = 0; iy <= ny; ++iy) {
            const double y = -half_thickness + 2.0 * half_thickness * static_cast<double>(iy) / static_cast<double>(ny);
            for (std::size_t ix = 0; ix <= spec.nx; ++ix) {
                const double x = length * static_cast<double>(ix) / static_cast<double>(spec.nx);
                fuelsim::CartesianPoint3 point{x, y, z};
                if (spec.distorted) {
                    const double sx = std::sin(pi * x / length);
                    const double sy = std::sin(pi * (y + half_thickness) / (2.0 * half_thickness));
                    const double sz = std::sin(pi * z / width);
                    point.x += 0.040 * sx * sy * (2.0 * z / width - 1.0);
                    point.y += 0.006 * sx * sy * sz;
                    point.z += 0.010 * sx * sy * sz * (y / half_thickness);
                }
                maximum_distortion = std::max(maximum_distortion,
                    std::sqrt(std::pow(point.x - x, 2) + std::pow(point.y - y, 2) + std::pow(point.z - z, 2)));
                nodes.push_back(point);
            }
        }
    }
    const auto node = [&spec, ny](std::size_t ix, std::size_t iy, std::size_t iz) {
        return iz * (ny + 1) * (spec.nx + 1) + iy * (spec.nx + 1) + ix;
    };
    std::vector<fuelsim::Hex8Element> elements;
    std::vector<std::int64_t> block_ids;
    std::vector<fuelsim::ElementSide> bottom, top, left_lower;
    for (std::size_t region = 0; region < 2; ++region) {
        const std::size_t first_y = region == 0 ? 0 : spec.subdivisions_per_layer;
        const std::size_t last_y = region == 0 ? spec.subdivisions_per_layer : ny;
        for (std::size_t iz = 0; iz < spec.nz; ++iz)
            for (std::size_t iy = first_y; iy < last_y; ++iy)
                for (std::size_t ix = 0; ix < spec.nx; ++ix) {
                    const std::size_t element = elements.size();
                    elements.push_back({{{node(ix, iy, iz), node(ix + 1, iy, iz), node(ix + 1, iy + 1, iz),
                        node(ix, iy + 1, iz), node(ix, iy, iz + 1), node(ix + 1, iy, iz + 1),
                        node(ix + 1, iy + 1, iz + 1), node(ix, iy + 1, iz + 1)}}});
                    block_ids.push_back(static_cast<std::int64_t>(region + 1));
                    if (region == 0 && iy == 0) bottom.push_back({element, 0});
                    if (region == 1 && iy + 1 == ny) top.push_back({element, 2});
                    if (region == 0 && ix == 0) left_lower.push_back({element, 3});
                }
    }
    std::vector<std::size_t> interface_nodes, bottom_right_nodes, top_right_nodes;
    for (std::size_t iz = 0; iz <= spec.nz; ++iz) {
        for (std::size_t ix = 0; ix <= spec.nx; ++ix)
            interface_nodes.push_back(node(ix, spec.subdivisions_per_layer, iz));
        bottom_right_nodes.push_back(node(spec.nx, 0, iz));
        top_right_nodes.push_back(node(spec.nx, ny, iz));
    }
    fuelsim::UnstructuredHex8Mesh mesh(std::move(nodes), std::move(elements), std::move(block_ids),
        {{1, "lower"}, {2, "upper"}}, {},
        {{11, "bottom", std::move(bottom)}, {12, "top", std::move(top)}, {13, "left_lower", std::move(left_lower)}});
    return {std::move(mesh), std::move(interface_nodes), std::move(bottom_right_nodes), std::move(top_right_nodes),
        maximum_distortion};
}

fuelsim::ThermoelasticProperties lower_material() {
    return fuelsim::test::thermoelastic(0.0, 8.0, 1.2e11, 0.28, 1.8e-5, 300.0, 0.0, 0.0, 0.0, 1.0e4, 300.0);
}

fuelsim::ThermoelasticProperties upper_material() {
    return fuelsim::test::thermoelastic(0.0, 24.0, 2.0e11, 0.30, 7.0e-6, 300.0, 0.0, 0.0, 0.0, 6500.0, 500.0);
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition result;
    result.regions.push_back({"lower", "lower", lower_material(), 0.0, 300.0});
    result.regions.push_back({"upper", "upper", upper_material(), 0.0, 300.0});
    result.time_tables.emplace_back("bottom_temperature",
        std::vector<double>{0.0, step_time, 2.0 * step_time, 3.0 * step_time, 4.0 * step_time},
        std::vector<double>{300.0, 480.0, 600.0, 600.0, 330.0});
    result.boundary_conditions = {
        {"bottom_temperature", fuelsim::BoundaryConditionType::dirichlet, "bottom", fuelsim::Field::temperature, 1.0,
            false, "bottom_temperature"},
        {"top_temperature", fuelsim::BoundaryConditionType::dirichlet, "top", fuelsim::Field::temperature, 300.0},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "left_lower", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "left_lower", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "left_lower", fuelsim::Field::displacement_z, 0.0},
    };
    return result;
}

fuelsim::SpatialDefinition steady_cooled_definition() {
    fuelsim::SpatialDefinition result = definition();
    result.time_tables.clear();
    result.boundary_conditions[0].value = 330.0;
    result.boundary_conditions[0].function.clear();
    return result;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = 1.0e-8;
    options.relative_tolerance = 1.0e-13;
    options.maximum_iterations = 20;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    options.field_residual_scaling = true;
    options.residual_reduction_tolerance = 1.0e-11;
    options.temperature_residual_absolute_tolerance = 1.0e-8;
    options.mechanical_residual_absolute_tolerance = 1.0e-2;
    return options;
}

std::size_t stage_from_time(double time, const std::string& source) {
    const double value = time / step_time;
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0 || rounded > 4.0)
        throw std::invalid_argument("B5.9 time does not identify one of four stages in " + source);
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

std::vector<std::size_t> source_global_nodes(const fuelsim::UnstructuredHex8Mesh& mesh,
    const fuelsim::TransientProblem& problem, const std::vector<std::size_t>& interface_nodes) {
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    std::vector<std::size_t> result(mesh.nodes().size(), std::numeric_limits<std::size_t>::max());
    std::vector<std::size_t> occurrences(mesh.nodes().size(), 0);
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local];
            const std::size_t global = dofs.global_node(region, local);
            if (occurrences[source]++ != 0 && result[source] != global)
                throw std::invalid_argument("B5.9 shared source node maps to different global nodes");
            result[source] = global;
        }
    }
    for (std::size_t source = 0; source < mesh.nodes().size(); ++source) {
        const bool interface =
            std::find(interface_nodes.begin(), interface_nodes.end(), source) != interface_nodes.end();
        if (occurrences[source] != (interface ? 2U : 1U))
            throw std::invalid_argument("B5.9 material-interface node ownership is not one-to-two");
    }
    if (dofs.node_count() != mesh.nodes().size())
        throw std::invalid_argument("B5.9 shared-node layout does not contain one global node per source node");
    return result;
}

std::vector<std::pair<std::size_t, std::size_t>> source_element_locations(
    const fuelsim::UnstructuredHex8Mesh& mesh, const fuelsim::TransientProblem& problem) {
    const std::size_t invalid = std::numeric_limits<std::size_t>::max();
    std::vector<std::pair<std::size_t, std::size_t>> result(mesh.elements().size(), {invalid, invalid});
    for (std::size_t region = 0; region < fuelsim::cartesian::ProblemAccess::region_count(problem); ++region) {
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        for (std::size_t local = 0; local < region_mesh.elements().size(); ++local)
            result.at(region_mesh.source_element_ids()[local]) = {region, local};
    }
    if (std::find(result.begin(), result.end(), std::pair<std::size_t, std::size_t>{invalid, invalid}) != result.end())
        throw std::invalid_argument("B5.9 source element did not map into a material region");
    return result;
}

double average_displacement_x(const fuelsim::spatial_detail::SpatialLayout& dofs,
    const std::vector<std::size_t>& source_to_global, const std::vector<std::size_t>& source_nodes,
    const std::vector<double>& state) {
    double sum = 0.0;
    for (const std::size_t source : source_nodes)
        sum += state[dofs.dof(fuelsim::Field::displacement_x, source_to_global.at(source))];
    return sum / static_cast<double>(source_nodes.size());
}

double interface_flux_imbalance(const MeshData& data, const std::vector<std::size_t>& source_to_global, bool& passed,
    const std::string& case_name) {
    fuelsim::SteadyProblem problem(steady_cooled_definition(), data.mesh);
    const fuelsim::SteadyResult solve = fuelsim::solve_steady(problem, {1, 0.5, 0, 1.0e-6}, solver_options());
    passed = check(solve.completed && solve.solve.converged,
                 "B5.9 " + case_name + " cooled steady interface-flux audit converges") &&
             passed;
    const auto& view = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    std::vector<std::size_t> interface_dofs;
    for (const std::size_t source : data.interface_nodes)
        interface_dofs.push_back(dofs.dof(fuelsim::Field::temperature, source_to_global.at(source)));
    std::array<double, 2> side_flux{};
    fuelsim::ContributionWorkspace workspace;
    for (std::size_t contribution = 0; contribution < view.volume_contribution_count(); ++contribution) {
        const std::size_t region = view.element_location(contribution).first;
        problem.evaluate_contribution(contribution, solve.solve.state, workspace, false);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            if (std::find(interface_dofs.begin(), interface_dofs.end(), workspace.dofs[local]) != interface_dofs.end())
                side_flux[region] += workspace.residual[local];
    }
    const double imbalance = std::abs(side_flux[0] + side_flux[1]);
    const double scale = std::max(std::abs(side_flux[0]), std::abs(side_flux[1]));
    const double relative = imbalance / scale;
    std::cout << "b59_interface_lower_heat_rate=" << side_flux[0] << " b59_interface_upper_heat_rate=" << side_flux[1]
              << " b59_interface_heat_rate_relative_imbalance=" << relative << '\n';
    return relative;
}

BendingSummary run_case(const CaseSpec& spec, const std::vector<NodeReference>& node_reference,
    const std::vector<IntegrationReference>& integration_reference, bool& passed) {
    MeshData data = make_mesh(spec);
    fuelsim::TransientProblem problem(definition(), data.mesh);
    const std::vector<std::size_t> source_to_global = source_global_nodes(data.mesh, problem, data.interface_nodes);
    const std::vector<std::pair<std::size_t, std::size_t>> element_locations =
        source_element_locations(data.mesh, problem);
    SnapshotObserver observer;
    const fuelsim::TransientResult solve = fuelsim::solve_transient(
        problem, {4.0 * step_time, step_time, step_time, step_time, 1.0, 0.5, 0, 0.0}, solver_options(), &observer);
    passed = check(solve.completed && observer.snapshots().size() == 4 && solve.accepted_steps.size() == 4 &&
                       solve.rejected_steps.empty(),
                 std::string("B5.9 ") + spec.name + " mesh accepts four prescribed increments without cutback") &&
             passed;
    std::map<std::size_t, const StepSnapshot*> snapshots;
    for (const StepSnapshot& snapshot : observer.snapshots())
        snapshots.emplace(stage_from_time(snapshot.time, "Fuelsim snapshots"), &snapshot);
    if (snapshots.size() != 4) throw std::invalid_argument("Fuelsim B5.9 snapshots do not cover four unique stages");

    std::array<fuelsim::test::FieldErrorMetrics, 8> nodal_metrics;
    fuelsim::TransientProblem reaction_problem(definition(), data.mesh);
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(reaction_problem);
    const std::array<fuelsim::Field, 4> fields = {fuelsim::Field::temperature, fuelsim::Field::displacement_x,
        fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};
    std::size_t compared_nodes = 0;
    for (std::size_t stage = 1; stage <= 4; ++stage) {
        const StepSnapshot& snapshot = *snapshots.at(stage);
        reaction_problem.begin_time_step({snapshot.time, snapshot.load_factor, true});
        const std::vector<double> reaction = raw_residual(reaction_problem, snapshot.state);
        for (const NodeReference& reference : node_reference) {
            if (reference.case_name != spec.name || reference.stage != stage) continue;
            if (std::abs(reference.time - snapshot.time) > 1.0e-6 || reference.node < 1 ||
                reference.node > data.mesh.nodes().size())
                throw std::invalid_argument("Abaqus B5.9 nodal label or time lies outside the Fuelsim case");
            const std::size_t global = source_to_global.at(reference.node - 1);
            for (std::size_t field = 0; field < fields.size(); ++field) {
                const std::size_t dof = dofs.dof(fields[field], global);
                nodal_metrics[field].add(snapshot.state[dof], reference.fields[field]);
                nodal_metrics[4 + field].add(reaction[dof], reference.fields[4 + field]);
            }
            ++compared_nodes;
        }
        reaction_problem.commit_time_step(snapshot.state);
        std::cout << "b59_" << spec.name << "_stage_" << stage
                  << "_relative_thermal_balance=" << snapshot.conservation.relative_thermal_balance
                  << " stored_heat_rate=" << snapshot.conservation.stored_heat_rate
                  << " dirichlet_heat_input_rate=" << snapshot.conservation.dirichlet_heat_input_rate
                  << " global_thermal_balance=" << snapshot.conservation.global_thermal_balance
                  << " unconstrained_thermal_residual_l2=" << snapshot.conservation.unconstrained_thermal_residual_l2
                  << '\n';
        const double thermal_scale = std::abs(snapshot.conservation.stored_heat_rate) +
                                     std::abs(snapshot.conservation.dirichlet_heat_input_rate);
        const bool relative_balance_required = thermal_scale >= 1.0e-3;
        const bool thermal_balance_passed =
            std::abs(snapshot.conservation.global_thermal_balance) < 1.0e-10 &&
            snapshot.conservation.unconstrained_thermal_residual_l2 < 1.0e-10 &&
            (!relative_balance_required || snapshot.conservation.relative_thermal_balance < 1.0e-10);
        passed =
            check(thermal_balance_passed,
                std::string("B5.9 ") + spec.name + " thermal conservation closes at stage " + std::to_string(stage) +
                    (relative_balance_required ? " by relative and absolute metrics"
                                               : " by the separately gated near-zero absolute metric")) &&
            passed;
    }
    passed = check(compared_nodes == 4 * data.mesh.nodes().size(),
                 std::string("B5.9 ") + spec.name + " compares every node at every stage") &&
             passed;

    const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y", "displacement_z",
        "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
    const std::array<double, 8> zero_tolerances = {1.0e-12, 1.0e-14, 1.0e-14, 1.0e-14, 2.0e-4, 5.0e-2, 5.0e-2, 5.0e-2};
    for (std::size_t field = 0; field < nodal_metrics.size(); ++field) {
        print_metrics(std::string("b59_") + spec.name + '_' + nodal_names[field], nodal_metrics[field]);
        passed = check(metrics_pass(nodal_metrics[field], 1.0e-3, zero_tolerances[field]),
                     std::string("B5.9 ") + spec.name + ' ' + nodal_names[field] +
                         " metrics are below the acceptance limits") &&
                 passed;
    }

    std::array<fuelsim::test::FieldErrorMetrics, 17> integration_metrics;
    std::size_t compared_points = 0;
    double maximum_coordinate_error = 0.0;
    constexpr std::array<std::size_t, 8> gauss_to_material_node = {0, 1, 3, 2, 4, 5, 7, 6};
    for (const IntegrationReference& reference : integration_reference) {
        if (reference.case_name != spec.name) continue;
        if (reference.element < 1 || reference.element > data.mesh.elements().size())
            throw std::invalid_argument("Abaqus B5.9 element label lies outside the Fuelsim mesh");
        const StepSnapshot& snapshot = *snapshots.at(reference.stage);
        if (std::abs(reference.time - snapshot.time) > 1.0e-6)
            throw std::invalid_argument("Abaqus B5.9 integration-point time does not match the Fuelsim snapshot");
        const auto location = element_locations.at(reference.element - 1);
        const std::size_t region = location.first, element = location.second;
        const auto& region_mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
        const fuelsim::Hex8Geometry& geometry =
            fuelsim::cartesian::ProblemAccess::region_element_geometry(problem, region, element);
        fuelsim::Hex8LocalValues local_state{};
        for (std::size_t local = 0; local < 8; ++local) {
            const std::size_t global = dofs.global_node(region, region_mesh.elements()[element].nodes[local]);
            local_state[local] = snapshot.state[dofs.dof(fuelsim::Field::temperature, global)];
            local_state[8 + local] = snapshot.state[dofs.dof(fuelsim::Field::displacement_x, global)];
            local_state[16 + local] = snapshot.state[dofs.dof(fuelsim::Field::displacement_y, global)];
            local_state[24 + local] = snapshot.state[dofs.dof(fuelsim::Field::displacement_z, global)];
        }
        std::size_t closest = 0;
        double closest_squared = std::numeric_limits<double>::max();
        for (std::size_t q = 0; q < 8; ++q) {
            const fuelsim::CartesianPoint3& point = geometry.points[q].position;
            const double distance_squared = std::pow(point.x - reference.position.x, 2) +
                                            std::pow(point.y - reference.position.y, 2) +
                                            std::pow(point.z - reference.position.z, 2);
            if (distance_squared < closest_squared) {
                closest = q;
                closest_squared = distance_squared;
            }
        }
        maximum_coordinate_error = std::max(maximum_coordinate_error, std::sqrt(closest_squared));
        const fuelsim::Hex8QuadraturePoint& point = geometry.points[closest];
        const fuelsim::IsotropicThermoelasticMaterial constitutive(
            fuelsim::cartesian::ProblemAccess::region(problem, region).material);
        const double material_temperature = local_state[gauss_to_material_node[closest]];
        const double conductivity = constitutive
                                        .conductivity(adlite::Scalar(material_temperature),
                                            {snapshot.time, point.position.x, point.position.y, point.position.z})
                                        .value();
        std::array<double, 3> actual_heat_flux{};
        for (std::size_t local = 0; local < 8; ++local)
            for (std::size_t component = 0; component < 3; ++component)
                actual_heat_flux[component] -= conductivity * point.gradient[local][component] * local_state[local];
        fuelsim::CartesianThermoelasticData material_data{
            constitutive, 0.0, snapshot.time, fuelsim::StrainFormulation::small};
        const fuelsim::SymmetricTensor3Values actual_stress =
            fuelsim::compute_hex8_stress(material_data, geometry, local_state)[closest];

        fuelsim::Hex8LocalAdValues passive{};
        for (std::size_t dof = 0; dof < passive.size(); ++dof) passive[dof] = local_state[dof];
        double volume = 0.0, average_trace = 0.0;
        for (const fuelsim::Hex8QuadraturePoint& qpoint : geometry.points) {
            const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
                qpoint, passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::small);
            volume += qpoint.weighted_measure;
            average_trace += qpoint.weighted_measure *
                             (kinematics.strain_increment.xx.value() + kinematics.strain_increment.yy.value() +
                                 kinematics.strain_increment.zz.value());
        }
        average_trace /= volume;
        const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
            point, passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::small);
        fuelsim::SymmetricTensor3 actual_strain = kinematics.strain_increment;
        const adlite::Scalar correction =
            (average_trace - actual_strain.xx - actual_strain.yy - actual_strain.zz) / 3.0;
        actual_strain.xx += correction;
        actual_strain.yy += correction;
        actual_strain.zz += correction;

        for (std::size_t component = 0; component < 3; ++component)
            integration_metrics[component].add(actual_heat_flux[component], reference.heat_flux[component]);
        const std::array<double, 6> actual_stress_components = {
            actual_stress.xx, actual_stress.yy, actual_stress.zz, actual_stress.xy, actual_stress.yz, actual_stress.xz};
        const std::array<double, 6> reference_stress_components = {reference.stress.xx, reference.stress.yy,
            reference.stress.zz, reference.stress.xy, reference.stress.yz, reference.stress.xz};
        const std::array<double, 6> actual_strain_components = {actual_strain.xx.value(), actual_strain.yy.value(),
            actual_strain.zz.value(), actual_strain.xy.value(), actual_strain.yz.value(), actual_strain.xz.value()};
        const std::array<double, 6> reference_strain_components = {reference.strain.xx, reference.strain.yy,
            reference.strain.zz, reference.strain.xy, reference.strain.yz, reference.strain.xz};
        for (std::size_t component = 0; component < 6; ++component) {
            integration_metrics[3 + component].add(
                actual_stress_components[component], reference_stress_components[component]);
            integration_metrics[9 + component].add(
                actual_strain_components[component], reference_strain_components[component]);
        }
        integration_metrics[15].add(material_temperature, reference.temperature);
        integration_metrics[16].add(point.weighted_measure, reference.integration_volume);
        ++compared_points;
    }
    passed = check(compared_points == 32 * data.mesh.elements().size(),
                 std::string("B5.9 ") + spec.name + " compares all eight points in every element and stage") &&
             passed;
    const std::array<std::string, 17> integration_names = {"heat_flux_x", "heat_flux_y", "heat_flux_z", "stress_xx",
        "stress_yy", "stress_zz", "stress_xy", "stress_yz", "stress_xz", "strain_xx", "strain_yy", "strain_zz",
        "strain_xy", "strain_yz", "strain_xz", "material_temperature", "integration_volume"};
    const std::array<double, 17> integration_zero_tolerances = {1.0e-8, 1.0e-8, 1.0e-8, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3,
        1.0e-3, 1.0e-3, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-14};
    for (std::size_t field = 0; field < integration_metrics.size(); ++field) {
        print_metrics(std::string("b59_") + spec.name + '_' + integration_names[field], integration_metrics[field]);
        passed =
            check(metrics_pass(integration_metrics[field], 1.0e-3, integration_zero_tolerances[field]),
                std::string("B5.9 ") + spec.name + ' ' + integration_names[field] + " metrics are below 0.1 percent") &&
            passed;
    }
    passed = check(maximum_coordinate_error < 1.0e-12,
                 std::string("B5.9 ") + spec.name + " maps all Abaqus integration points uniquely") &&
             passed;
    passed = check(interface_flux_imbalance(data, source_to_global, passed, spec.name) < 1.0e-10,
                 std::string("B5.9 ") + spec.name + " cooled steady material-interface heat rates balance") &&
             passed;
    if (spec.distorted)
        passed = check(data.maximum_distortion > 5.0e-3,
                     "B5.9 distorted mesh applies a material nonzero controlled nodal perturbation") &&
                 passed;
    else
        passed = check(data.maximum_distortion == 0.0,
                     std::string("B5.9 ") + spec.name + " comparison uses an undistorted mesh") &&
                 passed;

    const auto& solved_dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    const auto bending = [&](std::size_t stage) {
        const std::vector<double>& state = snapshots.at(stage)->state;
        return average_displacement_x(solved_dofs, source_to_global, data.top_right_nodes, state) -
               average_displacement_x(solved_dofs, source_to_global, data.bottom_right_nodes, state);
    };
    const BendingSummary result{bending(2), bending(4)};
    std::cout << "b59_" << spec.name << "_peak_axial_bending_indicator=" << result.peak << " b59_" << spec.name
              << "_cooled_axial_bending_indicator=" << result.cooled << '\n';
    passed =
        check(std::abs(result.peak) > 1.0e-4,
            std::string("B5.9 ") + spec.name + " develops resolved two-layer thermal bending") &&
        check(std::abs(result.cooled) > 1.0e-5 && std::abs(result.cooled) < 0.2 * std::abs(result.peak),
            std::string("B5.9 ") + spec.name + " retains the smaller deformation caused by the cooled 330 K state") &&
        passed;
    return result;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: fuelsim_b59_hex8_c3d8t_multimaterial_abaqus_tests <nodes.csv> <integration.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B5.9 Abaqus multi-material C3D8T comparison\n");
        const std::vector<NodeReference> nodes = read_nodes(argv[1]);
        const std::vector<IntegrationReference> integration = read_integration(argv[2]);
        bool passed = true;
        std::array<BendingSummary, 3> bending{};
        for (std::size_t index = 0; index < cases.size(); ++index)
            bending[index] = run_case(cases[index], nodes, integration, passed);
        const double distortion_sensitivity = std::abs(bending[2].peak - bending[1].peak) / std::abs(bending[1].peak);
        std::cout << "b59_refined_distorted_bending_relative_difference=" << distortion_sensitivity << '\n';
        passed =
            check(distortion_sensitivity < 0.1,
                "B5.9 controlled distortion changes the refined thermal-bending indicator by less than 10 percent") &&
            passed;
        if (passed && session.rank() == 0)
            std::cout << "[PASS] B5.9 Abaqus C3D8T shared-node multi-material refinement comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B5.9 Abaqus comparison raised: " << error.what() << '\n';
        return 1;
    }
}
