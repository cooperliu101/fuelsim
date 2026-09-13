#include "c3d8t.hpp"
#include "support/c3d8_common_tests.hpp"

namespace {
using namespace fuelsim::test::c3d8;

bool test_element_average_thermal_expansion_temperature() {
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(unit_cube());
    const fuelsim::IsotropicThermoelasticMaterial material(properties());
    const fuelsim::CartesianTestData data{material, 0.0, 0.0};
    fuelsim::Hex8LocalValues state{};
    const std::array<double, 8> nodal_temperatures{{360.0, 410.0, 445.0, 385.0, 470.0, 430.0, 515.0, 455.0}};
    double element_temperature = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = nodal_temperatures[node];
        element_temperature += nodal_temperatures[node] / 8.0;
    }
    const fuelsim::SymmetricTensor3 imposed = material.eigenstrain(element_temperature);
    const std::array<fuelsim::SymmetricTensor3Values, 8> stresses = fuelsim::compute_c3d8_stress(data, geometry, state);
    double recovered_expansion_maximum_difference = 0.0, stress_range = 0.0;
    double minimum_stress = stresses[0].xx, maximum_stress = stresses[0].xx;
    constexpr std::array<std::size_t, 8> gauss_to_material_node = {0, 1, 3, 2, 4, 5, 7, 6};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const double point_temperature = nodal_temperatures[gauss_to_material_node[q]];
        const fuelsim::ActiveThermoelasticProperties active = material.active_properties(point_temperature);
        const double bulk_factor = 3.0 * active.lame_lambda.value() + 2.0 * active.shear_modulus.value();
        const double recovered_expansion = -stresses[q].xx / bulk_factor;
        recovered_expansion_maximum_difference =
            std::max(recovered_expansion_maximum_difference, std::abs(recovered_expansion - imposed.xx.value()));
        minimum_stress = std::min(minimum_stress, stresses[q].xx);
        maximum_stress = std::max(maximum_stress, stresses[q].xx);
        if (!check(near(stresses[q].xx, stresses[q].yy, 1.0e-13) && near(stresses[q].xx, stresses[q].zz, 1.0e-13)
                       && std::abs(stresses[q].xy) < 1.0e-10 && std::abs(stresses[q].yz) < 1.0e-10
                       && std::abs(stresses[q].xz) < 1.0e-10,
                "nonuniform-temperature HEX8 thermal stress remains isotropic at each integration point"))
            return false;
    }
    stress_range = maximum_stress - minimum_stress;
    std::cout << "hex8_element_expansion_temperature=" << element_temperature << '\n'
              << "hex8_recovered_expansion_maximum_difference=" << recovered_expansion_maximum_difference << '\n'
              << "hex8_point_elasticity_stress_range=" << stress_range << '\n';
    return check(recovered_expansion_maximum_difference < 1.0e-15,
               "HEX8 thermal expansion uses one arithmetic-average nodal temperature")
           && check(stress_range > 1.0e6,
               "HEX8 elasticity retains integration-point temperature dependence while expansion is element constant");
}

double equivalent_stress(const fuelsim::SymmetricTensor3& stress) {
    const double mean = (stress.xx.value() + stress.yy.value() + stress.zz.value()) / 3.0;
    const double xx = stress.xx.value() - mean, yy = stress.yy.value() - mean, zz = stress.zz.value() - mean;
    return std::sqrt(1.5
                     * (xx * xx + yy * yy + zz * zz
                         + 2.0
                               * (stress.xy.value() * stress.xy.value() + stress.yz.value() * stress.yz.value()
                                   + stress.xz.value() * stress.xz.value())));
}

bool test_geometry_and_constant_strain() {
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(coordinates);
    double volume = 0.0;
    for (const fuelsim::Hex8QuadraturePoint& point : geometry.points) {
        volume += point.weighted_measure;
        double shape_sum = 0.0;
        std::array<double, 3> gradient_sum{};
        for (std::size_t node = 0; node < 8; ++node) {
            shape_sum += point.shape[node];
            for (std::size_t component = 0; component < 3; ++component)
                gradient_sum[component] += point.gradient[node][component];
        }
        if (!check(near(shape_sum, 1.0, 1.0e-14), "HEX8 shape functions form a partition of unity")
            || !check(near(gradient_sum[0], 0.0, 1.0e-14) && near(gradient_sum[1], 0.0, 1.0e-14)
                          && near(gradient_sum[2], 0.0, 1.0e-14),
                "HEX8 shape gradients sum to zero"))
            return false;
    }
    if (!check(near(volume, 1.0, 1.0e-14), "HEX8 eight-point integration recovers unit volume"))
        return false;
    fuelsim::Hex8LocalValues state{};
    const double exx = 0.01;
    const double eyy = -0.02;
    const double ezz = 0.03;
    const double exy = 0.004;
    const double eyz = -0.006;
    const double exz = 0.008;
    for (std::size_t node = 0; node < 8; ++node) {
        const fuelsim::CartesianPoint3& point = coordinates[node];
        state[node] = 300.0;
        state[8 + node] = exx * point.x + exy * point.y + exz * point.z;
        state[16 + node] = exy * point.x + eyy * point.y + eyz * point.z;
        state[24 + node] = exz * point.x + eyz * point.y + ezz * point.z;
    }
    const fuelsim::CartesianTestData data{fuelsim::IsotropicThermoelasticMaterial(properties()), 0.0, 0.0};
    const auto stresses = fuelsim::compute_c3d8_stress(data, geometry, state);
    const double lambda = 2.0e11 * 0.25 / (1.25 * 0.5);
    const double shear = 2.0e11 / 2.5;
    const double trace = exx + eyy + ezz;
    for (const fuelsim::SymmetricTensor3Values& stress : stresses)
        if (!check(near(stress.xx, lambda * trace + 2.0 * shear * exx, 2.0e-13)
                       && near(stress.yy, lambda * trace + 2.0 * shear * eyy, 2.0e-13)
                       && near(stress.zz, lambda * trace + 2.0 * shear * ezz, 2.0e-13)
                       && near(stress.xy, 2.0 * shear * exy, 2.0e-13) && near(stress.yz, 2.0 * shear * eyz, 2.0e-13)
                       && near(stress.xz, 2.0 * shear * exz, 2.0e-13),
                "HEX8 reproduces all six constant-strain stress components"))
            return false;
    return true;
}

bool test_distorted_selective_volumetric_integration() {
    fuelsim::Hex8Coordinates coordinates = unit_cube();
    coordinates[1] = {1.08, -0.03, 0.02};
    coordinates[2] = {1.16, 1.05, -0.04};
    coordinates[5] = {0.94, 0.06, 1.12};
    coordinates[6] = {1.24, 1.13, 1.28};
    coordinates[7] = {-0.08, 0.91, 1.06};
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(coordinates);
    constexpr double young_modulus = 2.0e11, poisson_ratio = 0.499;
    const fuelsim::CartesianTestData data{
        fuelsim::IsotropicThermoelasticMaterial(
            fuelsim::test::thermoelastic(0.0, 1.0, young_modulus, poisson_ratio, 0.0, 300.0)),
        0.0,
        0.0};
    fuelsim::Hex8LocalValues state{};
    double average_trace = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        const fuelsim::CartesianPoint3& point = coordinates[node];
        state[node] = 300.0;
        state[8 + node] = 1.0e-3 * (point.x + 0.35 * point.x * point.y - 0.12 * point.z);
        state[16 + node] = 1.0e-3 * (-0.20 * point.x + 0.55 * point.y * point.z + 0.08 * point.z);
        state[24 + node] = 1.0e-3 * (0.14 * point.x * point.z - 0.30 * point.y + 0.65 * point.z);
        for (std::size_t component = 0; component < 3; ++component)
            average_trace += geometry.average_shape_gradient[node][component] * state[8 * (component + 1) + node];
    }
    const double expected_stress_trace = young_modulus * average_trace / (1.0 - 2.0 * poisson_ratio);
    double stress_trace_error = 0.0;
    for (const fuelsim::SymmetricTensor3Values& stress : fuelsim::compute_c3d8_stress(data, geometry, state))
        stress_trace_error =
            std::max(stress_trace_error, std::abs(stress.xx + stress.yy + stress.zz - expected_stress_trace));

    fuelsim::Hex8LocalJacobian jacobian{};
    (void)fuelsim::compute_c3d8_thermoelastic(data, geometry, state, nullptr, 0.0, &jacobian);
    fuelsim::Hex8LocalValues direction{};
    for (std::size_t dof = 0; dof < direction.size(); ++dof)
        direction[dof] = dof < 8 ? 0.0 : std::sin(0.29 * static_cast<double>(dof + 1));
    constexpr double step = 1.0e-8;
    fuelsim::Hex8LocalValues plus = state, minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Hex8LocalResidual plus_residual = fuelsim::compute_c3d8_thermoelastic(data, geometry, plus),
                                     minus_residual = fuelsim::compute_c3d8_thermoelastic(data, geometry, minus);
    double jacobian_error = 0.0, jacobian_scale = 0.0;
    for (std::size_t row = 8; row < 32; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 8; column < 32; ++column)
            analytic += jacobian[row * 32 + column] * direction[column];
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        jacobian_error = std::max(jacobian_error, std::abs(analytic - numerical));
        jacobian_scale = std::max({jacobian_scale, std::abs(analytic), std::abs(numerical)});
    }
    std::cout << "hex8_selective_distorted_stress_trace_maximum_error=" << stress_trace_error << '\n'
              << "hex8_selective_near_incompressible_jacobian_relative_error=" << jacobian_error / jacobian_scale
              << '\n';
    return check(stress_trace_error / std::abs(expected_stress_trace) < 1.0e-12,
               "distorted HEX8 uses one volume-average strain trace at all eight integration points")
           && check(jacobian_error / jacobian_scale < 5.0e-7,
               "near-incompressible distorted HEX8 selective-integration Jacobian matches centered differences");
}

bool test_selective_integration_constrained_face_rank() {
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(unit_cube());
    const fuelsim::CartesianTestData data{
        fuelsim::IsotropicThermoelasticMaterial(fuelsim::test::thermoelastic(0.0, 1.0, 1.0e9, 0.0, 0.0, 300.0)),
        0.0,
        0.0};
    fuelsim::Hex8LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node)
        state[node] = 300.0;
    fuelsim::Hex8LocalJacobian jacobian{};
    (void)fuelsim::compute_c3d8_thermoelastic(data, geometry, state, nullptr, 0.0, &jacobian);
    constexpr std::array<std::size_t, 4> free_nodes = {1, 2, 5, 6};
    std::array<std::array<double, 12>, 12> reduced{};
    double maximum_entry = 0.0;
    for (std::size_t row_component = 0; row_component < 3; ++row_component)
        for (std::size_t row_node = 0; row_node < free_nodes.size(); ++row_node) {
            const std::size_t reduced_row = 4 * row_component + row_node;
            const std::size_t row = 8 * (row_component + 1) + free_nodes[row_node];
            for (std::size_t column_component = 0; column_component < 3; ++column_component)
                for (std::size_t column_node = 0; column_node < free_nodes.size(); ++column_node) {
                    const std::size_t reduced_column = 4 * column_component + column_node;
                    const std::size_t column = 8 * (column_component + 1) + free_nodes[column_node];
                    reduced[reduced_row][reduced_column] = jacobian[row * 32 + column];
                    maximum_entry = std::max(maximum_entry, std::abs(reduced[reduced_row][reduced_column]));
                }
        }
    double minimum_pivot = maximum_entry;
    for (std::size_t column = 0; column < reduced.size(); ++column) {
        std::size_t pivot_row = column;
        for (std::size_t row = column + 1; row < reduced.size(); ++row)
            if (std::abs(reduced[row][column]) > std::abs(reduced[pivot_row][column]))
                pivot_row = row;
        std::swap(reduced[column], reduced[pivot_row]);
        const double pivot = std::abs(reduced[column][column]);
        minimum_pivot = std::min(minimum_pivot, pivot);
        if (pivot == 0.0)
            break;
        for (std::size_t row = column + 1; row < reduced.size(); ++row) {
            const double factor = reduced[row][column] / reduced[column][column];
            for (std::size_t entry = column; entry < reduced.size(); ++entry)
                reduced[row][entry] -= factor * reduced[column][entry];
        }
    }
    std::cout << "hex8_selective_constrained_face_minimum_pivot_ratio=" << minimum_pivot / maximum_entry << '\n';
    return check(minimum_pivot / maximum_entry > 1.0e-10,
        "selective-integration HEX8 retains full stiffness rank after one face is fixed");
}

bool test_free_thermal_expansion_and_jacobian() {
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(coordinates);
    const fuelsim::CartesianTestData data{fuelsim::IsotropicThermoelasticMaterial(properties()), 7.0e5, 0.0};
    fuelsim::Hex8LocalValues state{};
    const double temperature = 650.0;
    const double active_alpha = 1.2e-5 + 1.0e-8 * (temperature - 300.0);
    const double strain = active_alpha * (temperature - 300.0);
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = temperature;
        state[8 + node] = strain * coordinates[node].x;
        state[16 + node] = strain * coordinates[node].y;
        state[24 + node] = strain * coordinates[node].z;
    }
    for (const fuelsim::SymmetricTensor3Values& stress : fuelsim::compute_c3d8_stress(data, geometry, state))
        if (!check(std::max({std::abs(stress.xx),
                       std::abs(stress.yy),
                       std::abs(stress.zz),
                       std::abs(stress.xy),
                       std::abs(stress.yz),
                       std::abs(stress.xz)})
                       < 1.0e-4,
                "uniform three-dimensional thermal expansion is stress free"))
            return false;
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        state[dof] += dof < 8 ? 2.0 * static_cast<double>(dof) : 1.0e-5 * static_cast<double>(dof + 1);
    std::array<double, 32> direction{};
    for (std::size_t dof = 0; dof < direction.size(); ++dof)
        direction[dof] = std::sin(0.37 * static_cast<double>(dof + 1));
    fuelsim::Hex8LocalJacobian jacobian{};
    (void)fuelsim::compute_c3d8_thermoelastic(data, geometry, state, nullptr, 0.0, &jacobian);
    const double epsilon = 1.0e-7;
    fuelsim::Hex8LocalValues plus = state;
    fuelsim::Hex8LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += epsilon * direction[dof];
        minus[dof] -= epsilon * direction[dof];
    }
    const fuelsim::Hex8LocalResidual plus_residual = fuelsim::compute_c3d8_thermoelastic(data, geometry, plus);
    const fuelsim::Hex8LocalResidual minus_residual = fuelsim::compute_c3d8_thermoelastic(data, geometry, minus);
    double maximum_error = 0.0;
    double scale = 0.0;
    for (std::size_t row = 0; row < 32; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 32; ++column)
            analytic += jacobian[row * 32 + column] * direction[column];
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * epsilon);
        maximum_error = std::max(maximum_error, std::abs(analytic - numerical));
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    return check(maximum_error / scale < 3.0e-7,
        "full 32-DOF HEX8 automatic-differentiation Jacobian matches a centered directional difference");
}

bool test_transient_capacity_and_faces() {
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(unit_cube());
    const fuelsim::CartesianTestData data{fuelsim::IsotropicThermoelasticMaterial(properties()), 1.2e7, 0.0};
    fuelsim::Hex8LocalValues old_state{};
    fuelsim::Hex8LocalValues state{};
    for (std::size_t node = 0; node < 8; ++node) {
        old_state[node] = 300.0;
        state[node] = 302.0;
    }
    const fuelsim::Hex8LocalResidual residual =
        fuelsim::compute_c3d8_thermoelastic(data, geometry, state, &old_state, 1.0);
    for (std::size_t node = 0; node < 8; ++node)
        if (!check(std::abs(residual[node]) < 1.0e-7,
                "Backward Euler lumped heat capacity balances uniform volumetric heating; residual="
                    + std::to_string(residual[node])))
            return false;
    const fuelsim::CartesianMaterialHistory history(geometry.points.size());
    fuelsim::Hex8LocalValues nonuniform_state = old_state;
    for (std::size_t node = 0; node < 8; ++node)
        nonuniform_state[node] += 0.25 * static_cast<double>(node + 1);
    fuelsim::Hex8LocalJacobian with_capacity{}, without_capacity{};
    const fuelsim::Hex8LocalResidual with_capacity_residual =
        fuelsim::compute_c3d8_transient(data, geometry, nonuniform_state, old_state, history, 2.0, &with_capacity);
    const fuelsim::Hex8LocalResidual without_capacity_residual = fuelsim::compute_c3d8_transient(data,
        geometry,
        nonuniform_state,
        old_state,
        history,
        2.0,
        &without_capacity,
        false);
    constexpr double volumetric_capacity = 2000.0 * 3000.0;
    double capacity_residual_error = 0.0, capacity_jacobian_error = 0.0;
    for (std::size_t row = 0; row < 8; ++row) {
        const double expected_residual = geometry.capacity_points[row].weighted_measure * volumetric_capacity
                                         * (nonuniform_state[row] - old_state[row]) / 2.0;
        capacity_residual_error = std::max(capacity_residual_error,
            std::abs((with_capacity_residual[row] - without_capacity_residual[row]) - expected_residual));
        for (std::size_t column = 0; column < 8; ++column) {
            const double expected =
                row == column ? geometry.capacity_points[row].weighted_measure * volumetric_capacity / 2.0 : 0.0;
            capacity_jacobian_error = std::max(capacity_jacobian_error,
                std::abs((with_capacity[row * 32 + column] - without_capacity[row * 32 + column]) - expected));
        }
    }
    if (!check(capacity_residual_error < 1.0e-7,
            "HEX8 nodal heat-capacity residual uses the Abaqus corner integration weights")
        || !check(capacity_jacobian_error < 1.0e-7,
            "HEX8 nodal heat-capacity Jacobian is diagonal with the Abaqus corner weights"))
        return false;
    const fuelsim::CartesianTestData nonlinear_capacity_data{
        fuelsim::IsotropicThermoelasticMaterial(capacity_properties()),
        0.0,
        0.0};
    fuelsim::Hex8LocalJacobian nonlinear_with{}, nonlinear_without{};
    (void)fuelsim::compute_c3d8_transient(nonlinear_capacity_data,
        geometry,
        nonuniform_state,
        old_state,
        history,
        2.0,
        &nonlinear_with);
    (void)fuelsim::compute_c3d8_transient(nonlinear_capacity_data,
        geometry,
        nonuniform_state,
        old_state,
        history,
        2.0,
        &nonlinear_without,
        false);
    double nonlinear_capacity_jacobian_error = 0.0;
    for (std::size_t row = 0; row < 8; ++row)
        for (std::size_t column = 0; column < 8; ++column) {
            const double temperature = nonuniform_state[row], increment = temperature - old_state[row];
            const double expected = row == column ? geometry.capacity_points[row].weighted_measure * 2000.0
                                                        * (1000.0 + 2.0 * temperature + 2.0 * increment) / 2.0
                                                  : 0.0;
            nonlinear_capacity_jacobian_error = std::max(nonlinear_capacity_jacobian_error,
                std::abs((nonlinear_with[row * 32 + column] - nonlinear_without[row * 32 + column]) - expected));
        }
    if (!check(nonlinear_capacity_jacobian_error < 1.0e-7,
            "temperature-dependent lumped heat capacity retains an exact diagonal Jacobian"))
        return false;
    fuelsim::Hex8Coordinates distorted_coordinates = unit_cube();
    distorted_coordinates[6] = {1.20, 1.10, 1.30};
    const fuelsim::Hex8Geometry distorted = fuelsim::elements::make_c3d8t_geometry(distorted_coordinates);
    double minimum_capacity_weight = distorted.capacity_points[0].weighted_measure,
           maximum_capacity_weight = minimum_capacity_weight, capacity_gauss_mapping_error = 0.0;
    constexpr std::array<std::size_t, 8> node_to_gauss = {0, 1, 3, 2, 4, 5, 7, 6};
    for (std::size_t node = 0; node < distorted.capacity_points.size(); ++node) {
        const fuelsim::Hex8CapacityPoint& point = distorted.capacity_points[node];
        minimum_capacity_weight = std::min(minimum_capacity_weight, point.weighted_measure);
        maximum_capacity_weight = std::max(maximum_capacity_weight, point.weighted_measure);
        capacity_gauss_mapping_error = std::max(capacity_gauss_mapping_error,
            std::abs(point.weighted_measure - distorted.points[node_to_gauss[node]].weighted_measure));
    }
    std::cout << "hex8_distorted_capacity_weight_minimum=" << minimum_capacity_weight << '\n'
              << "hex8_distorted_capacity_weight_maximum=" << maximum_capacity_weight << '\n'
              << "hex8_distorted_capacity_gauss_mapping_maximum_error=" << capacity_gauss_mapping_error << '\n'
              << "hex8_nonlinear_capacity_jacobian_maximum_error=" << nonlinear_capacity_jacobian_error << '\n';
    if (!check(maximum_capacity_weight - minimum_capacity_weight > 1.0e-3,
            "distorted HEX8 nodal heat capacity uses distinct associated standard Gauss-point volume weights")
        || !check(capacity_gauss_mapping_error < 1.0e-15,
            "distorted HEX8 nodal heat capacity maps every node to the Abaqus-associated standard Gauss point"))
        return false;
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Quad4FaceCoordinates face_coordinates = {
        {coordinates[1], coordinates[2], coordinates[6], coordinates[5]}};
    const fuelsim::Quad4FaceGeometry face = fuelsim::make_quad4_face_geometry(face_coordinates);
    fuelsim::Quad4FaceLocalValues face_state{};
    const fuelsim::Quad4FaceBoundaryData pressure = {fuelsim::Quad4FaceBoundaryKind::pressure,
        fuelsim::CartesianTractionComponent::x,
        5.0,
        0.0};
    const fuelsim::Quad4FaceLocalResidual pressure_residual =
        fuelsim::compute_quad4_face_boundary(pressure, face, face_state);
    double force_x = 0.0;
    double force_y = 0.0;
    double force_z = 0.0;
    for (std::size_t node = 0; node < 4; ++node) {
        force_x += pressure_residual[4 + node];
        force_y += pressure_residual[8 + node];
        force_z += pressure_residual[12 + node];
    }
    if (!check(near(force_x, 5.0, 1.0e-14) && near(force_y, 0.0, 1.0e-14) && near(force_z, 0.0, 1.0e-14),
            "reference pressure uses the outward three-dimensional face area vector and exact total force"))
        return false;
    face_state.fill(0.0);
    for (std::size_t node = 0; node < 4; ++node) {
        face_state[4 + node] = 0.2 * face_coordinates[node].y;
        face_state[8 + node] = 0.1 * face_coordinates[node].y;
        face_state[12 + node] = -0.05 * face_coordinates[node].z;
    }
    const fuelsim::Quad4FaceBoundaryData follower = {fuelsim::Quad4FaceBoundaryKind::pressure,
        fuelsim::CartesianTractionComponent::x,
        5.0,
        0.0,
        true};
    fuelsim::Quad4FaceLocalJacobian follower_jacobian{};
    const auto follower_residual = fuelsim::compute_quad4_face_boundary(follower, face, face_state, &follower_jacobian);
    force_x = 0.0;
    force_y = 0.0;
    force_z = 0.0;
    for (std::size_t node = 0; node < 4; ++node) {
        force_x += follower_residual[4 + node];
        force_y += follower_residual[8 + node];
        force_z += follower_residual[12 + node];
    }
    if (!check(near(force_x, 5.225, 1.0e-14) && near(force_y, -0.95, 1.0e-14) && near(force_z, 0.0, 1.0e-14),
            "three-dimensional follower pressure uses the complete current area vector and exact resultant"))
        return false;
    std::array<double, 16> follower_direction{};
    for (std::size_t dof = 4; dof < 16; ++dof)
        follower_direction[dof] = std::cos(0.41 * static_cast<double>(dof + 1));
    constexpr double follower_step = 1.0e-7;
    auto follower_plus = face_state, follower_minus = face_state;
    for (std::size_t dof = 0; dof < 16; ++dof) {
        follower_plus[dof] += follower_step * follower_direction[dof];
        follower_minus[dof] -= follower_step * follower_direction[dof];
    }
    const auto follower_plus_residual = fuelsim::compute_quad4_face_boundary(follower, face, follower_plus);
    const auto follower_minus_residual = fuelsim::compute_quad4_face_boundary(follower, face, follower_minus);
    double follower_error = 0.0;
    for (std::size_t row = 4; row < 16; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 16; ++column)
            analytic += follower_jacobian[row * 16 + column] * follower_direction[column];
        const double numerical = (follower_plus_residual[row] - follower_minus_residual[row]) / (2.0 * follower_step);
        follower_error = std::max(follower_error, std::abs(analytic - numerical));
    }
    if (!check(follower_error < 2.0e-9,
            "three-dimensional follower-pressure geometric Jacobian matches centered difference"))
        return false;
    const fuelsim::Quad4FaceBoundaryData current_traction = {fuelsim::Quad4FaceBoundaryKind::traction,
        fuelsim::CartesianTractionComponent::x,
        5.0,
        0.0,
        true};
    const fuelsim::Quad4FaceBoundaryData reference_traction = {fuelsim::Quad4FaceBoundaryKind::traction,
        fuelsim::CartesianTractionComponent::x,
        5.0,
        0.0,
        false};
    fuelsim::Quad4FaceLocalJacobian current_traction_jacobian{};
    const auto current_traction_residual =
        fuelsim::compute_quad4_face_boundary(current_traction, face, face_state, &current_traction_jacobian);
    const auto reference_traction_residual = fuelsim::compute_quad4_face_boundary(reference_traction, face, face_state);
    double current_force = 0.0, reference_force = 0.0;
    for (std::size_t node = 0; node < 4; ++node) {
        current_force += current_traction_residual[4 + node];
        reference_force += reference_traction_residual[4 + node];
    }
    auto current_traction_plus = face_state, current_traction_minus = face_state;
    for (std::size_t dof = 0; dof < 16; ++dof) {
        current_traction_plus[dof] += follower_step * follower_direction[dof];
        current_traction_minus[dof] -= follower_step * follower_direction[dof];
    }
    const auto current_traction_plus_residual =
        fuelsim::compute_quad4_face_boundary(current_traction, face, current_traction_plus);
    const auto current_traction_minus_residual =
        fuelsim::compute_quad4_face_boundary(current_traction, face, current_traction_minus);
    double current_traction_error = 0.0;
    for (std::size_t row = 4; row < 8; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 16; ++column)
            analytic += current_traction_jacobian[row * 16 + column] * follower_direction[column];
        const double numerical =
            (current_traction_plus_residual[row] - current_traction_minus_residual[row]) / (2.0 * follower_step);
        current_traction_error = std::max(current_traction_error, std::abs(analytic - numerical));
    }
    std::cout << "current_traction_force=" << current_force << '\n'
              << "reference_traction_force=" << reference_force << '\n'
              << "current_traction_jacobian_error=" << current_traction_error << '\n';
    if (!check(std::abs(current_force - reference_force) > 1.0e-6 && current_traction_error < 2.0e-8,
            "three-dimensional current-configuration component traction uses the current face measure and consistent "
            "geometric Jacobian"))
        return false;
    for (std::size_t node = 0; node < 4; ++node)
        face_state[node] = 350.0;
    const fuelsim::Quad4FaceBoundaryData convection = {fuelsim::Quad4FaceBoundaryKind::convection,
        fuelsim::CartesianTractionComponent::x,
        20.0,
        300.0};
    fuelsim::Quad4FaceLocalJacobian convection_jacobian{};
    const fuelsim::Quad4FaceLocalResidual convection_residual =
        fuelsim::compute_quad4_face_boundary(convection, face, face_state, &convection_jacobian);
    double heat = 0.0;
    double tangent_sum = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        heat += convection_residual[row];
        for (std::size_t column = 0; column < 4; ++column)
            tangent_sum += convection_jacobian[row * 16 + column];
    }
    if (!check(near(heat, 1000.0, 1.0e-14) && near(tangent_sum, 20.0, 1.0e-14),
            "three-dimensional convection has the exact face heat rate and consistent temperature tangent"))
        return false;

    const fuelsim::Quad4FaceCoordinates warped_coordinates = {
        {{0.0, 0.0, 0.0}, {1.2, 0.1, 0.0}, {1.0, 1.1, 0.35}, {-0.1, 0.9, -0.05}}};
    const fuelsim::Quad4FaceGeometry warped_face = fuelsim::make_quad4_face_geometry(warped_coordinates);
    fuelsim::Quad4FaceLocalValues warped_state{};
    std::array<double, 16> warped_direction{};
    for (std::size_t node = 0; node < 4; ++node) {
        warped_state[node] = 335.0 + 17.0 * static_cast<double>(node);
        warped_state[4 + node] = 0.03 * static_cast<double>(node + 1);
        warped_state[8 + node] = -0.02 * static_cast<double>(node * node + 1);
        warped_state[12 + node] = 0.015 * static_cast<double>((node + 1) * (node + 2));
    }
    for (std::size_t dof = 0; dof < warped_direction.size(); ++dof)
        warped_direction[dof] = std::cos(0.37 * static_cast<double>(dof + 1));
    const fuelsim::Quad4FaceBoundaryData current_flux = {fuelsim::Quad4FaceBoundaryKind::surface_heat_flux,
        fuelsim::CartesianTractionComponent::x,
        40.0,
        0.0,
        true};
    const fuelsim::Quad4FaceBoundaryData current_convection = {fuelsim::Quad4FaceBoundaryKind::convection,
        fuelsim::CartesianTractionComponent::x,
        20.0,
        300.0,
        true};
    fuelsim::Quad4FaceLocalJacobian current_flux_jacobian{}, current_convection_jacobian{};
    (void)fuelsim::compute_quad4_face_boundary(current_flux, warped_face, warped_state, &current_flux_jacobian);
    (void)fuelsim::compute_quad4_face_boundary(current_convection,
        warped_face,
        warped_state,
        &current_convection_jacobian);
    constexpr double thermal_face_step = 1.0e-7;
    auto warped_plus = warped_state, warped_minus = warped_state;
    for (std::size_t dof = 0; dof < warped_state.size(); ++dof) {
        warped_plus[dof] += thermal_face_step * warped_direction[dof];
        warped_minus[dof] -= thermal_face_step * warped_direction[dof];
    }
    const fuelsim::Quad4FaceLocalResidual current_flux_plus =
        fuelsim::compute_quad4_face_boundary(current_flux, warped_face, warped_plus);
    const fuelsim::Quad4FaceLocalResidual current_flux_minus =
        fuelsim::compute_quad4_face_boundary(current_flux, warped_face, warped_minus);
    const fuelsim::Quad4FaceLocalResidual current_convection_plus =
        fuelsim::compute_quad4_face_boundary(current_convection, warped_face, warped_plus);
    const fuelsim::Quad4FaceLocalResidual current_convection_minus =
        fuelsim::compute_quad4_face_boundary(current_convection, warped_face, warped_minus);
    double current_flux_jacobian_error = 0.0, current_convection_jacobian_error = 0.0,
           current_flux_displacement_coupling = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        double flux_analytic = 0.0, convection_analytic = 0.0;
        for (std::size_t column = 0; column < 16; ++column) {
            flux_analytic += current_flux_jacobian[row * 16 + column] * warped_direction[column];
            convection_analytic += current_convection_jacobian[row * 16 + column] * warped_direction[column];
            if (column >= 4)
                current_flux_displacement_coupling =
                    std::max(current_flux_displacement_coupling, std::abs(current_flux_jacobian[row * 16 + column]));
        }
        const double flux_numerical = (current_flux_plus[row] - current_flux_minus[row]) / (2.0 * thermal_face_step);
        const double convection_numerical =
            (current_convection_plus[row] - current_convection_minus[row]) / (2.0 * thermal_face_step);
        current_flux_jacobian_error = std::max(current_flux_jacobian_error, std::abs(flux_analytic - flux_numerical));
        current_convection_jacobian_error =
            std::max(current_convection_jacobian_error, std::abs(convection_analytic - convection_numerical));
    }
    std::cout << "hex8_current_surface_flux_jacobian_maximum_error=" << current_flux_jacobian_error << '\n'
              << "hex8_current_convection_jacobian_maximum_error=" << current_convection_jacobian_error << '\n'
              << "hex8_current_surface_flux_displacement_coupling_maximum=" << current_flux_displacement_coupling
              << '\n';
    return check(current_flux_jacobian_error < 2.0e-7 && current_convection_jacobian_error < 2.0e-6,
               "current-configuration surface heat flux and convection Jacobians match centered differences")
           && check(current_flux_displacement_coupling > 1.0e-6,
               "current-configuration surface heat flux has a nonzero displacement coupling block");
}

bool test_cartesian_inelastic_material() {
    const fuelsim::SymmetricTensor3 strain{0.20, -0.04, -0.03, 0.02, -0.015, 0.01};
    const fuelsim::CartesianMaterialPointState committed{};
    bool passed = true;
    for (const std::array<bool, 2> branch : {std::array<bool, 2>{false, true}, {true, false}, {true, true}}) {
        const fuelsim::IsotropicThermoelasticMaterial material(inelastic_properties(branch[0], branch[1]));
        const fuelsim::CartesianInelasticStressResponse response = material.response(strain, 300.0, 1.0, committed);
        const fuelsim::CartesianMaterialPointState state = response.trial_state;
        const double plastic_trace = state.plastic_strain[0] + state.plastic_strain[1] + state.plastic_strain[2];
        const double creep_trace = state.creep_strain[0] + state.creep_strain[1] + state.creep_strain[2];
        passed = check((branch[1] ? state.equivalent_plastic_strain > 0.0 : state.equivalent_plastic_strain == 0.0)
                           && (branch[0] ? state.equivalent_creep_strain > 0.0 : state.equivalent_creep_strain == 0.0)
                           && std::abs(plastic_trace) < 2.0e-15 && std::abs(creep_trace) < 2.0e-15
                           && std::isfinite(equivalent_stress(response.stress)),
                     "Cartesian plastic, creep, and coupled updates activate the requested traceless branches")
                 && passed;
        const adlite::Scalar active_xx = adlite::Scalar::independent(strain.xx.value(), 0, 1);
        const fuelsim::CartesianInelasticStressResponse active =
            material.response({active_xx, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz},
                300.0,
                1.0,
                committed);
        constexpr double step = 1.0e-7;
        const double plus =
            material
                .response({strain.xx.value() + step, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz},
                    300.0,
                    1.0,
                    committed)
                .stress.xx.value();
        const double minus =
            material
                .response({strain.xx.value() - step, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz},
                    300.0,
                    1.0,
                    committed)
                .stress.xx.value();
        const double numerical = (plus - minus) / (2.0 * step);
        passed =
            check(std::abs(active.stress.xx.derivative(0) - numerical) / std::max(1.0, std::abs(numerical)) < 2.0e-7,
                "Cartesian inelastic material automatic-differentiation tangent matches centered difference")
            && passed;
    }
    return passed;
}

bool test_nonaffine_finite_capacity_weights() {
    // Reference Gauss weights independently evaluated in the nonaffine capacity analysis.
    // See verification/abaqus/b523_diagnosis/nonaffine_capacity_candidates.tsv (reference_gauss).
    // The Abaqus 2018 current-volume capacity remains a separate external comparison.
    const std::array<fuelsim::Hex8Coordinates, 2> coordinates{{{{{0.0, 0.0, 0.0},
                                                                   {0.5, 0.0, 0.0},
                                                                   {0.5, 1.0, 0.0},
                                                                   {0.0, 1.0, 0.0},
                                                                   {0.0, 0.0, 1.0},
                                                                   {0.5, 0.0, 1.0},
                                                                   {0.5, 1.0, 1.0},
                                                                   {0.0, 1.0, 1.0}}},
        {{{0.0, 0.0, 0.0},
            {1.2, 0.1, -0.05},
            {1.1, 1.0, 0.1},
            {-0.1, 0.9, 0.0},
            {0.05, -0.1, 1.0},
            {1.3, 0.0, 1.1},
            {1.0, 1.2, 0.9},
            {-0.2, 1.0, 1.2}}}}};
    const fuelsim::Hex8Coordinates displacement{{{0.0, 0.0, 0.0},
        {0.12, 0.025, 0.04},
        {0.06, 0.10, -0.035},
        {-0.015, 0.03, 0.01},
        {0.03, -0.02, 0.08},
        {0.14, 0.015, 0.02},
        {0.09, 0.045, 0.15},
        {-0.04, 0.11, 0.025}}};
    const std::array<std::array<double, 8>, 2> reference_weights{{{{0.062499999999999979,
                                                                      0.062499999999999979,
                                                                      0.062500000000000028,
                                                                      0.062499999999999979,
                                                                      0.062500000000000056,
                                                                      0.062499999999999951,
                                                                      0.062499999999999979,
                                                                      0.062500000000000028}},
        {{0.15233523564835924,
            0.15730243156300722,
            0.13523398686603305,
            0.15778250719944623,
            0.17779552702285578,
            0.18759422891166486,
            0.16325677824052978,
            0.18090763788143716}}}};
    const fuelsim::CartesianTestData data{
        fuelsim::IsotropicThermoelasticMaterial(
            fuelsim::test::thermoelastic(0.0, 1.0e-20, 2.0e11, 0.25, 0.0, 300.0, 0.0, 0.0, 0.0, 2000.0, 3000.0)),
        0.0,
        1.0,
        fuelsim::StrainFormulation::finite};
    constexpr double volumetric_capacity = 2000.0 * 3000.0;
    double maximum_reference_relative_error = 0.0, maximum_tangent_relative_error = 0.0;
    double maximum_temperature_block_error = 0.0, maximum_residual_path_error = 0.0;
    double maximum_displacement_coupling = 0.0;
    for (std::size_t shape = 0; shape < coordinates.size(); ++shape) {
        const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(coordinates[shape]);
        const fuelsim::CartesianMaterialHistory history(8);
        fuelsim::Hex8LocalValues state{}, old_state{};
        for (std::size_t node = 0; node < 8; ++node) {
            old_state[node] = 300.0;
            state[node] = 301.0 + 0.25 * static_cast<double>(node);
            state[8 + node] = displacement[node].x;
            state[16 + node] = displacement[node].y;
            state[24 + node] = displacement[node].z;
        }
        fuelsim::Hex8LocalJacobian jacobian{};
        const auto residual =
            fuelsim::compute_c3d8_transient(data, geometry, state, old_state, history, 1.0, &jacobian);
        const auto ordinary = fuelsim::compute_c3d8_transient(data, geometry, state, old_state, history, 1.0);
        for (std::size_t row = 0; row < 8; ++row) {
            const double expected = reference_weights[shape][row] * volumetric_capacity * (state[row] - old_state[row]);
            maximum_reference_relative_error =
                std::max(maximum_reference_relative_error, std::abs(residual[row] - expected) / expected);
            maximum_residual_path_error =
                std::max(maximum_residual_path_error, std::abs(residual[row] - ordinary[row]) / expected);
            for (std::size_t column = 0; column < 8; ++column) {
                const double expected_tangent =
                    row == column ? reference_weights[shape][row] * volumetric_capacity : 0.0;
                maximum_temperature_block_error = std::max(maximum_temperature_block_error,
                    std::abs(jacobian[row * 32 + column] - expected_tangent)
                        / (reference_weights[shape][row] * volumetric_capacity));
            }
        }
        constexpr double difference_step = 2.0e-7;
        for (std::size_t column = 0; column < 32; ++column) {
            fuelsim::Hex8LocalValues plus = state, minus = state;
            plus[column] += difference_step;
            minus[column] -= difference_step;
            const auto plus_residual = fuelsim::compute_c3d8_transient(data, geometry, plus, old_state, history, 1.0);
            const auto minus_residual = fuelsim::compute_c3d8_transient(data, geometry, minus, old_state, history, 1.0);
            double column_error = 0.0, column_scale = 0.0;
            for (std::size_t row = 0; row < 8; ++row) {
                const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * difference_step);
                const double exact = jacobian[row * 32 + column];
                column_error = std::max(column_error, std::abs(exact - numerical));
                column_scale = std::max({column_scale, std::abs(exact), std::abs(numerical)});
                if (column >= 8)
                    maximum_displacement_coupling =
                        std::max({maximum_displacement_coupling, std::abs(exact), std::abs(numerical)});
            }
            if (column < 8)
                maximum_tangent_relative_error = std::max(maximum_tangent_relative_error,
                    column_scale > 0.0 ? column_error / column_scale : column_error);
        }
    }
    std::cout << "hex8_nonaffine_capacity_reference_maximum_relative_error=" << maximum_reference_relative_error << '\n'
              << "hex8_nonaffine_capacity_tangent_maximum_relative_error=" << maximum_tangent_relative_error << '\n'
              << "hex8_nonaffine_capacity_temperature_block_error=" << maximum_temperature_block_error << '\n';
    return check(maximum_reference_relative_error < 2.0e-13,
               "nonaffine finite capacity uses independently evaluated reference Gauss weights")
           && check(maximum_temperature_block_error < 2.0e-13,
               "nonaffine finite capacity retains the initial-mass diagonal temperature block")
           && check(maximum_residual_path_error < 2.0e-14,
               "ordinary and tangent residual paths agree for nonaffine finite capacity")
           && check(maximum_tangent_relative_error < 2.0e-6 && maximum_displacement_coupling < 1.0e-12,
               "all nonaffine capacity temperature and displacement columns agree with centered differences");
}

bool test_finite_strain_kinematics_and_coupled_jacobian() {
    const fuelsim::Hex8Coordinates coordinates = unit_cube();
    const fuelsim::Hex8Geometry geometry = fuelsim::elements::make_c3d8t_geometry(coordinates);
    fuelsim::Hex8LocalValues state{};
    const double stretch_x = 1.12, stretch_y = 0.94, stretch_z = 1.03;
    for (std::size_t node = 0; node < 8; ++node) {
        state[node] = 300.0;
        state[8 + node] = (stretch_x - 1.0) * coordinates[node].x;
        state[16 + node] = (stretch_y - 1.0) * coordinates[node].y;
        state[24 + node] = (stretch_z - 1.0) * coordinates[node].z;
    }
    const auto diagnostics = fuelsim::test::recover_c3d8t(geometry, state, {}, fuelsim::StrainFormulation::finite);
    const auto& kinematics = diagnostics.points[0];
    const auto hughes_winget = [](double stretch) {
        return 2.0 * (stretch - 1.0) / (stretch + 1.0);
    };
    bool passed = check(near(kinematics.strain_increment.xx, hughes_winget(stretch_x), 2.0e-14)
                            && near(kinematics.strain_increment.yy, hughes_winget(stretch_y), 2.0e-14)
                            && near(kinematics.strain_increment.zz, hughes_winget(stretch_z), 2.0e-14)
                            && near(kinematics.current_weighted_measure,
                                geometry.points[0].weighted_measure * stretch_x * stretch_y * stretch_z,
                                2.0e-14),
        "finite-strain HEX8 recovers the Abaqus Hughes-Winget diagonal increment and current volume measure");
    for (std::size_t node = 0; node < 8; ++node) {
        const fuelsim::CartesianPoint3& point = coordinates[node];
        state[node] = 302.0 + 0.75 * static_cast<double>(node);
        state[8 + node] = 0.20 * point.x + 0.08 * point.y - 0.03 * point.z;
        state[16 + node] = -0.02 * point.x - 0.04 * point.y + 0.06 * point.z;
        state[24 + node] = 0.04 * point.x - 0.05 * point.y - 0.03 * point.z;
    }
    const fuelsim::CartesianTestData data{fuelsim::IsotropicThermoelasticMaterial(inelastic_properties(true, true)),
        8.0e4,
        1.0,
        fuelsim::StrainFormulation::finite};
    const fuelsim::Hex8LocalValues committed_state = [] {
        fuelsim::Hex8LocalValues value{};
        for (std::size_t node = 0; node < 8; ++node)
            value[node] = 298.0 + 0.5 * static_cast<double>(node);
        return value;
    }();
    const fuelsim::CartesianMaterialHistory committed_material(8);
    fuelsim::Hex8LocalJacobian jacobian{};
    const fuelsim::Hex8LocalResidual residual_with_jacobian =
        fuelsim::compute_c3d8_transient(data, geometry, state, committed_state, committed_material, 1.0, &jacobian);
    const fuelsim::Hex8LocalResidual residual_without_jacobian =
        fuelsim::compute_c3d8_transient(data, geometry, state, committed_state, committed_material, 1.0);
    double residual_path_maximum_difference = 0.0, residual_path_scale = 0.0;
    for (std::size_t row = 0; row < residual_with_jacobian.size(); ++row) {
        residual_path_maximum_difference = std::max(residual_path_maximum_difference,
            std::abs(residual_with_jacobian[row] - residual_without_jacobian[row]));
        residual_path_scale = std::max(
            {residual_path_scale, std::abs(residual_with_jacobian[row]), std::abs(residual_without_jacobian[row])});
    }
    std::cout << "hex8_c3d8t_finite_residual_path_maximum_absolute_difference=" << residual_path_maximum_difference
              << '\n'
              << "hex8_c3d8t_finite_residual_path_relative_difference="
              << residual_path_maximum_difference / residual_path_scale << '\n';
    passed = check(residual_path_maximum_difference / residual_path_scale < 2.0e-14,
                 "finite-strain C3D8T Jacobian and ordinary-double residual paths agree to roundoff")
             && passed;
    fuelsim::Hex8LocalValues heated_state = state;
    for (std::size_t node = 0; node < 8; ++node)
        heated_state[node] = 305.0;
    const fuelsim::Hex8LocalResidual with_capacity =
        fuelsim::compute_c3d8_transient(data, geometry, heated_state, committed_state, committed_material, 1.0);
    const fuelsim::Hex8LocalResidual without_capacity = fuelsim::compute_c3d8_transient(data,
        geometry,
        heated_state,
        committed_state,
        committed_material,
        1.0,
        nullptr,
        false);
    double maximum_capacity_difference = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        maximum_capacity_difference =
            std::max(maximum_capacity_difference, std::abs(with_capacity[node] - without_capacity[node]));
    passed = check(maximum_capacity_difference > 1.0e-6,
                 "transient HEX8 thermal time term can be disabled independently of heat conduction")
             && passed;
    std::array<double, 32> direction{};
    for (std::size_t dof = 0; dof < direction.size(); ++dof)
        direction[dof] = std::sin(0.29 * static_cast<double>(dof + 1));
    constexpr double step = 2.0e-7;
    fuelsim::Hex8LocalValues plus = state, minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const auto plus_residual =
        fuelsim::compute_c3d8_transient(data, geometry, plus, committed_state, committed_material, 1.0);
    const auto minus_residual =
        fuelsim::compute_c3d8_transient(data, geometry, minus, committed_state, committed_material, 1.0);
    double thermal_error = 0.0, thermal_scale = 0.0, mechanical_error = 0.0, mechanical_scale = 0.0;
    for (std::size_t row = 0; row < 32; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 32; ++column)
            analytic += jacobian[row * 32 + column] * direction[column];
        const double numerical = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        double& error = row < 8 ? thermal_error : mechanical_error;
        double& scale = row < 8 ? thermal_scale : mechanical_scale;
        error = std::max(error, std::abs(analytic - numerical));
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    double thermal_displacement_coupling = 0.0;
    for (std::size_t row = 0; row < 8; ++row)
        for (std::size_t column = 8; column < 32; ++column)
            thermal_displacement_coupling =
                std::max(thermal_displacement_coupling, std::abs(jacobian[row * 32 + column]));
    std::cout << "hex8_finite_heat_ktu_maximum_absolute=" << thermal_displacement_coupling << '\n'
              << "hex8_finite_heat_jacobian_relative_error=" << thermal_error / thermal_scale << '\n';
    passed = check(thermal_error / thermal_scale < 2.0e-6 && mechanical_error / mechanical_scale < 2.0e-6,
                 "finite-strain coupled HEX8 thermal and mechanical Jacobian rows match centered differences")
             && check(thermal_displacement_coupling > 1.0e-6,
                 "current-configuration HEX8 heat conduction has a nonzero displacement coupling block")
             && passed;
    std::array<double, 4> block_relative_errors{}, block_maximum_entries{};
    for (std::size_t block = 0; block < 4; ++block) {
        const bool thermal_rows = block < 2;
        const bool thermal_columns = block == 0 || block == 2;
        const double block_step = thermal_columns ? 1.0e-4 : step;
        fuelsim::Hex8LocalValues block_plus = state, block_minus = state;
        for (std::size_t column = 0; column < 32; ++column) {
            if ((column < 8) != thermal_columns)
                continue;
            block_plus[column] += block_step * direction[column];
            block_minus[column] -= block_step * direction[column];
        }
        const fuelsim::Hex8LocalResidual block_plus_residual =
            fuelsim::compute_c3d8_transient(data, geometry, block_plus, committed_state, committed_material, 1.0);
        const fuelsim::Hex8LocalResidual block_minus_residual =
            fuelsim::compute_c3d8_transient(data, geometry, block_minus, committed_state, committed_material, 1.0);
        double error_squared = 0.0, scale_squared = 0.0;
        for (std::size_t row = 0; row < 32; ++row) {
            if ((row < 8) != thermal_rows)
                continue;
            double analytic = 0.0;
            for (std::size_t column = 0; column < 32; ++column) {
                if ((column < 8) != thermal_columns)
                    continue;
                const double entry = jacobian[row * 32 + column];
                analytic += entry * direction[column];
                block_maximum_entries[block] = std::max(block_maximum_entries[block], std::abs(entry));
            }
            const double numerical = (block_plus_residual[row] - block_minus_residual[row]) / (2.0 * block_step);
            const double error = analytic - numerical;
            error_squared += error * error;
            scale_squared += numerical * numerical;
        }
        block_relative_errors[block] = std::sqrt(error_squared / scale_squared);
    }
    std::cout << "hex8_finite_ktt_directional_relative_error=" << block_relative_errors[0] << '\n'
              << "hex8_finite_ktu_directional_relative_error=" << block_relative_errors[1] << '\n'
              << "hex8_finite_kut_directional_relative_error=" << block_relative_errors[2] << '\n'
              << "hex8_finite_kuu_directional_relative_error=" << block_relative_errors[3] << '\n'
              << "hex8_finite_ktu_maximum_entry=" << block_maximum_entries[1] << '\n'
              << "hex8_finite_kut_maximum_entry=" << block_maximum_entries[2] << '\n';
    passed = check(*std::max_element(block_relative_errors.begin(), block_relative_errors.end()) < 2.0e-6,
                 "all four finite-strain thermo-mechanical Jacobian blocks match centered directional differences")
             && check(block_maximum_entries[1] > 1.0e-6 && block_maximum_entries[2] > 1.0e-6,
                 "both off-diagonal finite-strain thermo-mechanical Jacobian blocks are active")
             && passed;
    fuelsim::Hex8LocalValues inverted = state;
    for (std::size_t node = 0; node < 8; ++node)
        inverted[8 + node] = -2.0 * coordinates[node].x;
    try {
        (void)fuelsim::compute_c3d8_transient(data, geometry, inverted, committed_state, committed_material, 1.0);
        passed = check(false, "finite-strain HEX8 rejects a nonpositive deformation Jacobian") && passed;
    } catch (const std::domain_error&) {
    }
    fuelsim::CartesianMaterialPointState history;
    history.elastic_strain = {0.01, -0.02, 0.03, 0.004, -0.005, 0.006};
    history.plastic_strain = {0.02, -0.01, -0.01, 0.007, 0.008, -0.009};
    history.creep_strain = {-0.03, 0.01, 0.02, -0.011, 0.012, 0.013};
    history.equivalent_plastic_strain = 0.04;
    history.equivalent_creep_strain = 0.05;
    const fuelsim::CartesianRotation quarter_turn{0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    const auto rotated =
        fuelsim::IsotropicThermoelasticMaterial(inelastic_properties(false, false))
            .incremental_response({0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, quarter_turn, 300.0, 300.0, 1.0, history)
            .trial_state;
    const std::array<std::array<double, 6>, 3> expected = {{{-0.02, 0.01, 0.03, -0.004, 0.006, 0.005},
        {-0.01, 0.02, -0.01, -0.007, -0.009, -0.008},
        {0.01, -0.03, 0.02, 0.011, 0.013, -0.012}}};
    const std::array<const std::array<double, 6>*, 3> actual = {&rotated.elastic_strain,
        &rotated.plastic_strain,
        &rotated.creep_strain};
    bool objective = true;
    for (std::size_t tensor_index = 0; tensor_index < actual.size(); ++tensor_index)
        for (std::size_t component = 0; component < 6; ++component)
            objective =
                objective && near((*actual[tensor_index])[component], expected[tensor_index][component], 1.0e-14);
    passed = check(objective && rotated.equivalent_plastic_strain == history.equivalent_plastic_strain
                       && rotated.equivalent_creep_strain == history.equivalent_creep_strain,
                 "finite-strain Cartesian elastic, plastic, and creep tensors rotate objectively while scalars do not")
             && passed;
    return passed;
}

int run_c3d8t_tests() {
    bool passed = test_history_geometry(false);
    passed = test_fixed_initial_mass_capacity(false) && passed;
    passed = test_geometry_and_constant_strain() && passed;
    passed = test_distorted_selective_volumetric_integration() && passed;
    passed = test_selective_integration_constrained_face_rank() && passed;
    passed = test_free_thermal_expansion_and_jacobian() && passed;
    passed = test_element_average_thermal_expansion_temperature() && passed;
    passed = test_transient_capacity_and_faces() && passed;
    passed = test_nonaffine_finite_capacity_weights() && passed;
    passed = test_cartesian_inelastic_material() && passed;
    passed = test_finite_strain_kinematics_and_coupled_jacobian() && passed;
    if (!passed)
        return 1;
    std::cout << "HEX8 thermo-mechanics tests passed\n";
    return 0;
}
} // namespace

int main() {
    return run_c3d8t_tests();
}
