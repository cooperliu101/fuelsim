#include "core/cartesian3d_hex20.hpp"
#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/exodus_result_reader.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct GeometryEvidence final {
    double minimum_contact_radius = std::numeric_limits<double>::max();
    double maximum_contact_radius_error = 0.0;
    double maximum_face_planarity_error = 0.0;
    double minimum_face_normal_dot = 1.0;
    double minimum_within_face_normal_dot = 1.0;
};

struct ContactNumericalEvidence final {
    std::array<double, 3> residual_balance{};
    double maximum_jacobian_directional_error = 0.0;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::CartesianPoint3 subtract(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return {first.x - second.x, first.y - second.y, first.z - second.z};
}

fuelsim::CartesianPoint3 cross(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return {first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x};
}

double dot(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

fuelsim::CartesianPoint3 unit(const fuelsim::CartesianPoint3& value) {
    const double measure = std::sqrt(dot(value, value));
    if (!(measure > 0.0))
        throw std::invalid_argument("H20.30 geometry has an undefined surface normal");
    return {value.x / measure, value.y / measure, value.z / measure};
}

GeometryEvidence inspect_contact_geometry(const fuelsim::UnstructuredHex20Mesh& mesh) {
    static constexpr std::array<std::array<std::size_t, 8>, 6> face_nodes = {{{{0, 1, 5, 4, 8, 13, 16, 12}},
        {{1, 2, 6, 5, 9, 14, 17, 13}},
        {{2, 3, 7, 6, 10, 15, 18, 14}},
        {{3, 0, 4, 7, 11, 12, 19, 15}},
        {{0, 3, 2, 1, 11, 10, 9, 8}},
        {{4, 5, 6, 7, 16, 17, 18, 19}}}};
    GeometryEvidence result;
    std::vector<fuelsim::CartesianPoint3> face_normals;
    for (const fuelsim::ElementSide& side : mesh.side_set("secondary_contact").sides) {
        fuelsim::Quad8FaceCoordinates coordinates{};
        for (std::size_t local = 0; local < coordinates.size(); ++local) {
            coordinates[local] = mesh.nodes()[mesh.elements()[side.element].nodes[face_nodes[side.local_side][local]]];
            const double radius = std::hypot(coordinates[local].x, coordinates[local].y);
            result.minimum_contact_radius = std::min(result.minimum_contact_radius, radius);
            result.maximum_contact_radius_error = std::max(result.maximum_contact_radius_error, std::abs(radius - 1.0));
        }
        const fuelsim::CartesianPoint3 plane_normal =
            unit(cross(subtract(coordinates[1], coordinates[0]), subtract(coordinates[3], coordinates[0])));
        for (const fuelsim::CartesianPoint3& point : coordinates)
            result.maximum_face_planarity_error = std::max(result.maximum_face_planarity_error,
                std::abs(dot(subtract(point, coordinates[0]), plane_normal)));
        const fuelsim::Quad8FaceMechanicalQuadraturePoint center =
            fuelsim::make_quad8_face_mechanical_point(coordinates, 0.0, 0.0, 1.0);
        face_normals.push_back(unit(cross(center.tangent_xi, center.tangent_eta)));
        const fuelsim::Quad8FaceMechanicalQuadraturePoint first =
            fuelsim::make_quad8_face_mechanical_point(coordinates, -0.75, -0.75, 1.0);
        const fuelsim::Quad8FaceMechanicalQuadraturePoint second =
            fuelsim::make_quad8_face_mechanical_point(coordinates, 0.75, -0.75, 1.0);
        result.minimum_within_face_normal_dot = std::min(result.minimum_within_face_normal_dot,
            std::abs(dot(unit(cross(first.tangent_xi, first.tangent_eta)),
                unit(cross(second.tangent_xi, second.tangent_eta)))));
    }
    for (std::size_t first = 0; first < face_normals.size(); ++first)
        for (std::size_t second = first + 1; second < face_normals.size(); ++second)
            result.minimum_face_normal_dot =
                std::min(result.minimum_face_normal_dot, std::abs(dot(face_normals[first], face_normals[second])));
    return result;
}

ContactNumericalEvidence inspect_contact_numerics(const fuelsim::SteadyProblem& problem,
    const fuelsim::cartesian::SpatialAssembly& spatial,
    const std::vector<double>& state) {
    ContactNumericalEvidence result;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        problem.contribution_dofs(contribution, dofs);
        std::vector<double> local_state(dofs.size());
        for (std::size_t local = 0; local < dofs.size(); ++local)
            local_state[local] = state[dofs[local]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local_state, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local_state.size()), plus = local_state, minus = local_state;
        constexpr double perturbation = 1.0e-10;
        for (std::size_t local = 0; local < local_state.size(); ++local) {
            direction[local] = std::sin(static_cast<double>(local + 1));
            plus[local] += perturbation * direction[local];
            minus[local] -= perturbation * direction[local];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local_state.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local_state.size(); ++column)
                analytic += jacobian[row * local_state.size() + column] * direction[column];
            const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += std::pow(analytic - finite_difference, 2);
            reference_squared += finite_difference * finite_difference;
        }
        if (!(reference_squared > 0.0))
            throw std::logic_error("H20.30 active contact has a zero Jacobian direction");
        result.maximum_jacobian_directional_error =
            std::max(result.maximum_jacobian_directional_error, std::sqrt(difference_squared / reference_squared));
        for (std::size_t local = 0; local < dofs.size(); ++local)
            for (std::size_t component = 0; component < result.residual_balance.size(); ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (dofs[local] >= field.begin && dofs[local] < field.end)
                    result.residual_balance[component] += residual[local];
            }
    }
    return result;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4)
        return 2;
    try {
        const std::string name = argv[3];
        if (name != "faceted_cylinder" && name != "quadratic_cylinder" && name != "quadratic_friction")
            throw std::invalid_argument("Unknown H20.30 geometry");
        const auto definition = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
        fuelsim::SteadyProblem problem(definition.spatial, mesh);
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        std::vector<double> state(problem.dof_count());
        const auto output = fuelsim::test::read_final_exodus_results(argv[2]);
        if (output.nodes.size() != mesh.nodes().size())
            throw std::runtime_error("HEX20 result node count differs");
        const std::array<std::string, 3> displacement_names = {"displacement_x", "displacement_y", "displacement_z"};
        for (std::size_t region = 0; region < spatial.region_count(); ++region) {
            const auto& local_mesh = spatial.hex20_region_mesh(region);
            for (std::size_t local = 0; local < local_mesh.nodes().size(); ++local) {
                const auto source = local_mesh.source_node_ids()[local];
                const auto global = spatial.global_node(region, local);
                for (std::size_t component = 0; component < 3; ++component)
                    state.at(spatial.field_layout()[component + 1].begin + global) =
                        output.nodal(displacement_names[component]).at(source);
            }
            for (std::size_t local = 0; local < local_mesh.nodes().size(); ++local) {
                if (!local_mesh.temperature_nodes()[local])
                    continue;
                const auto source = local_mesh.source_node_ids()[local];
                state.at(spatial.dof(fuelsim::Field::temperature, spatial.global_temperature_node(region, local))) =
                    output.nodal("temperature").at(source);
            }
        }
        problem.validate_state(state);
        if (name == "quadratic_friction"
            && fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0).size() != 37)
            throw std::runtime_error("H20.35 requires one history slot for each of its 37 contact constraints");

        const auto geometry = inspect_contact_geometry(mesh);
        const auto numerical = inspect_contact_numerics(problem, spatial, state);
        const bool faceted_geometry = geometry.minimum_contact_radius < 0.99
                                      && geometry.maximum_face_planarity_error < 1.0e-12
                                      && geometry.minimum_within_face_normal_dot > 1.0 - 1.0e-12;
        const bool quadratic_geometry = geometry.maximum_contact_radius_error < 1.0e-12
                                        && geometry.maximum_face_planarity_error > 1.0e-3
                                        && geometry.minimum_within_face_normal_dot < 0.99;
        std::cout << "h20_30_" << name << "_minimum_contact_node_radius=" << geometry.minimum_contact_radius << '\n'
                  << "maximum_contact_node_radius_error=" << geometry.maximum_contact_radius_error << '\n'
                  << "maximum_face_planarity_error=" << geometry.maximum_face_planarity_error << '\n'
                  << "minimum_between_face_normal_dot=" << geometry.minimum_face_normal_dot << '\n'
                  << "minimum_within_face_normal_dot=" << geometry.minimum_within_face_normal_dot << '\n'
                  << "contact_residual_balance=" << numerical.residual_balance[0] << ','
                  << numerical.residual_balance[1] << ',' << numerical.residual_balance[2] << '\n'
                  << "contact_jacobian_directional_error=" << numerical.maximum_jacobian_directional_error << '\n';
        const bool passed = check(geometry.minimum_face_normal_dot < 0.95, "H20.30 normals vary across the cylinder")
                            && check(name == "faceted_cylinder" ? faceted_geometry : quadratic_geometry,
                                "H20.30 preserves the requested faceted or genuinely quadratic surface")
                            && check(std::abs(numerical.residual_balance[0]) < 1.0e-8
                                         && std::abs(numerical.residual_balance[1]) < 1.0e-8
                                         && std::abs(numerical.residual_balance[2]) < 1.0e-8,
                                "H20.30 contact residual is action-reaction conservative")
                            && check(numerical.maximum_jacobian_directional_error < 1.0e-7,
                                "H20.30 contact Jacobian matches a centered directional difference");
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
