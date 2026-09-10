#include "solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/field_error_metrics.hpp"
#include "support/material_factory.hpp"
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

constexpr std::array<CaseSpec, 2> cases = {CaseSpec{"coarse", 2, 1, 1, false}, CaseSpec{"distorted", 4, 2, 2, true}};

struct MeshData final {
    fuelsim::UnstructuredHex8Mesh mesh;
    std::vector<std::size_t> interface_nodes, bottom_right_nodes, top_right_nodes;
    double maximum_distortion;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
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
                    elements.push_back({{{node(ix, iy, iz),
                        node(ix + 1, iy, iz),
                        node(ix + 1, iy + 1, iz),
                        node(ix, iy + 1, iz),
                        node(ix, iy, iz + 1),
                        node(ix + 1, iy, iz + 1),
                        node(ix + 1, iy + 1, iz + 1),
                        node(ix, iy + 1, iz + 1)}}});
                    block_ids.push_back(static_cast<std::int64_t>(region + 1));
                    if (region == 0 && iy == 0)
                        bottom.push_back({element, 0});
                    if (region == 1 && iy + 1 == ny)
                        top.push_back({element, 2});
                    if (region == 0 && ix == 0)
                        left_lower.push_back({element, 3});
                }
    }
    std::vector<std::size_t> interface_nodes, bottom_right_nodes, top_right_nodes;
    for (std::size_t iz = 0; iz <= spec.nz; ++iz) {
        for (std::size_t ix = 0; ix <= spec.nx; ++ix)
            interface_nodes.push_back(node(ix, spec.subdivisions_per_layer, iz));
        bottom_right_nodes.push_back(node(spec.nx, 0, iz));
        top_right_nodes.push_back(node(spec.nx, ny, iz));
    }
    fuelsim::UnstructuredHex8Mesh mesh(std::move(nodes),
        std::move(elements),
        std::move(block_ids),
        {{1, "lower"}, {2, "upper"}},
        {},
        {{11, "bottom", std::move(bottom)}, {12, "top", std::move(top)}, {13, "left_lower", std::move(left_lower)}});
    return {std::move(mesh),
        std::move(interface_nodes),
        std::move(bottom_right_nodes),
        std::move(top_right_nodes),
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
        {"bottom_temperature",
            fuelsim::BoundaryConditionType::dirichlet,
            "bottom",
            fuelsim::Field::temperature,
            1.0,
            false,
            "bottom_temperature"},
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

std::vector<std::size_t> source_global_nodes(const fuelsim::UnstructuredHex8Mesh& mesh,
    const fuelsim::TransientProblem& problem,
    const std::vector<std::size_t>& interface_nodes) {
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

double interface_flux_imbalance(const MeshData& data,
    const std::vector<std::size_t>& source_to_global,
    bool& passed,
    const std::string& case_name) {
    fuelsim::SteadyProblem problem(steady_cooled_definition(), data.mesh);
    const fuelsim::SteadyResult solve = fuelsim::solve_steady(problem, {1, 0.5, 0, 1.0e-6}, solver_options());
    passed = check(solve.completed && solve.solve.converged,
                 "B5.9 " + case_name + " cooled steady interface-flux audit converges")
             && passed;
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

} // namespace

int main(int argc, char** argv) {
    try {
        fuelsim::PetscSession session(argc, argv, "B5.9 internal shared-node and interface conservation checks");
        bool passed = true;
        for (const auto& spec : cases) {
            const auto data = make_mesh(spec);
            fuelsim::TransientProblem problem(definition(), data.mesh);
            const auto source_to_global = source_global_nodes(data.mesh, problem, data.interface_nodes);
            if (spec.distorted)
                passed = check(interface_flux_imbalance(data, source_to_global, passed, spec.name) < 1e-10,
                             std::string("B5.9 ") + spec.name + " cooled interface heat rates balance")
                         && passed;
            passed = check(spec.distorted ? data.maximum_distortion > 5e-3 : data.maximum_distortion == 0,
                         "B5.9 mesh perturbation matches the specified branch")
                     && passed;
        }
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
