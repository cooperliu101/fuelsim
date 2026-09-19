#include "support/c3d8_common_tests.hpp"

namespace {
using namespace fuelsim::test::c3d8;

bool test_cartesian_surface_contact_kernels() {
    const fuelsim::Quad4FaceCoordinates
        secondary = {{{0.99, 0.0, 0.0}, {0.99, 1.0, 0.0}, {0.99, 1.0, 1.0}, {0.99, 0.0, 1.0}}},
        primary = {{{1.0, -0.5, -0.5}, {1.0, 1.5, -0.5}, {1.0, 1.5, 1.5}, {1.0, -0.5, 1.5}}};
    const fuelsim::Quad4FaceGeometry secondary_face = fuelsim::make_quad4_face_geometry(secondary);
    const fuelsim::Quad4FaceQuadraturePoint& point = secondary_face.points[0];
    const fuelsim::Quad4ToQuad4HeatGeometry heat_geometry{secondary,
        primary,
        point.shape,
        point.derivative_xi,
        point.derivative_eta,
        1.0};
    fuelsim::Quad4SurfaceContactLocalValues state{}, committed{};
    for (std::size_t node = 0; node < 4; ++node) {
        state[node] = 400.0;
        state[4 + node] = 300.0;
        committed[node] = 400.0;
        committed[4 + node] = 300.0;
        state[8 + node] = 0.02;
        state[16 + node] = 0.001;
    }
    fuelsim::Quad4SurfaceContactLocalJacobian heat_jacobian{};
    const fuelsim::Quad4SurfaceContactLocalResidual heat =
        fuelsim::compute_quad4_to_quad4_gap_heat({0.2, 0.001}, heat_geometry, state, &heat_jacobian);
    auto translated_geometry = heat_geometry;
    for (auto* coordinates : {&translated_geometry.primary_coordinates, &translated_geometry.secondary_coordinates})
        for (auto& coordinate : *coordinates) {
            coordinate.x += 128.0;
            coordinate.y -= 64.0;
            coordinate.z += 32.0;
        }
    fuelsim::Quad4SurfaceContactLocalJacobian translated_jacobian{};
    const auto translated_heat =
        fuelsim::compute_quad4_to_quad4_gap_heat({0.2, 0.001}, translated_geometry, state, &translated_jacobian);
    for (std::size_t row = 0; row < heat.size(); ++row) {
        if (!check(std::abs(translated_heat[row] - heat[row]) < 1e-9 * std::max(1.0, std::abs(heat[row])),
                "Global translation preserves Q4 contact heat residuals"))
            return false;
        for (std::size_t column = 0; column < heat.size(); ++column)
            if (!check(std::abs(
                           translated_jacobian[row * heat.size() + column] - heat_jacobian[row * heat.size() + column])
                           < 1e-9 * std::max(1.0, std::abs(heat_jacobian[row * heat.size() + column])),
                    "Global translation preserves Q4 contact heat Jacobians"))
                return false;
    }
    double secondary_heat = 0.0, primary_heat = 0.0;
    for (std::size_t node = 0; node < 4; ++node) {
        secondary_heat += heat[node];
        primary_heat += heat[4 + node];
    }
    bool passed = check(near(secondary_heat, 5000.0, 1.0e-12) && near(primary_heat, -5000.0, 1.0e-12),
        "three-dimensional gap heat transfer is exactly conservative on the secondary quadrature point");
    fuelsim::GapHeatProperties affine_heat_properties{};
    affine_heat_properties.law = fuelsim::GapHeatConductanceLaw::affine;
    affine_heat_properties.conductance = 50.0;
    affine_heat_properties.clearance_derivative = -1000.0;
    affine_heat_properties.pressure_derivative = 2.0e-6;
    affine_heat_properties.temperature_derivative = 0.1;
    affine_heat_properties.reference_temperature = 350.0;
    affine_heat_properties.contact_penalty = 1.0e5;
    fuelsim::Quad4SurfaceContactLocalJacobian affine_heat_jacobian{};
    const fuelsim::Quad4SurfaceContactLocalResidual affine_heat =
        fuelsim::compute_quad4_to_quad4_gap_heat(affine_heat_properties, heat_geometry, state, &affine_heat_jacobian);
    secondary_heat = 0.0;
    primary_heat = 0.0;
    for (std::size_t node = 0; node < 4; ++node) {
        secondary_heat += affine_heat[node];
        primary_heat += affine_heat[4 + node];
    }
    passed = check(near(secondary_heat, 1500.05, 1.0e-12) && near(primary_heat, -1500.05, 1.0e-12),
                 "clearance-, pressure-, and temperature-dependent gap conductance is exactly conservative")
             && passed;
    std::array<double, 32> affine_heat_direction{};
    for (std::size_t dof = 0; dof < affine_heat_direction.size(); ++dof)
        affine_heat_direction[dof] = std::sin(0.29 * static_cast<double>(dof + 1));
    constexpr double affine_heat_step = 1.0e-6;
    fuelsim::Quad4SurfaceContactLocalValues affine_heat_plus = state, affine_heat_minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        affine_heat_plus[dof] += affine_heat_step * affine_heat_direction[dof];
        affine_heat_minus[dof] -= affine_heat_step * affine_heat_direction[dof];
    }
    const fuelsim::Quad4SurfaceContactLocalResidual affine_heat_plus_residual =
                                                        fuelsim::compute_quad4_to_quad4_gap_heat(affine_heat_properties,
                                                            heat_geometry,
                                                            affine_heat_plus),
                                                    affine_heat_minus_residual =
                                                        fuelsim::compute_quad4_to_quad4_gap_heat(affine_heat_properties,
                                                            heat_geometry,
                                                            affine_heat_minus);
    double affine_heat_maximum_error = 0.0, affine_heat_scale = 0.0;
    for (std::size_t row = 0; row < 8; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 32; ++column)
            analytic += affine_heat_jacobian[row * 32 + column] * affine_heat_direction[column];
        const double numerical =
            (affine_heat_plus_residual[row] - affine_heat_minus_residual[row]) / (2.0 * affine_heat_step);
        affine_heat_maximum_error = std::max(affine_heat_maximum_error, std::abs(analytic - numerical));
        affine_heat_scale = std::max({affine_heat_scale, std::abs(analytic), std::abs(numerical)});
    }
    std::cout << "hex8_affine_gap_heat_jacobian_relative_error=" << affine_heat_maximum_error / affine_heat_scale
              << '\n';
    passed = check(affine_heat_maximum_error / affine_heat_scale < 1.0e-7,
                 "clearance-, pressure-, and temperature-dependent gap heat automatic-differentiation Jacobian "
                 "matches centered difference")
             && passed;
    std::array<std::array<double, 4>, 4> shapes{}, derivatives_xi{}, derivatives_eta{};
    for (std::size_t q = 0; q < 4; ++q) {
        shapes[q] = secondary_face.points[q].shape;
        derivatives_xi[q] = secondary_face.points[q].derivative_xi;
        derivatives_eta[q] = secondary_face.points[q].derivative_eta;
    }
    const fuelsim::NodeToQuad4ContactGeometry
        contact_geometry{secondary, primary, shapes, derivatives_xi, derivatives_eta, 0, 1.0};
    fuelsim::Quad4SurfaceContactLocalJacobian contact_jacobian{};
    const fuelsim::NormalContactProperties stick_properties{1000.0, 1.0, false};
    const fuelsim::ContactPointHistory history{};
    const fuelsim::Quad4SurfaceContactLocalResidual contact = fuelsim::compute_node_to_quad4_contact(stick_properties,
        contact_geometry,
        state,
        committed,
        history,
        &contact_jacobian);
    const fuelsim::CartesianContactPointValue stick =
        fuelsim::compute_node_to_quad4_contact_value(stick_properties, contact_geometry, state, committed, history);
    std::array<double, 3> resultant{};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t node = 0; node < 8; ++node)
            resultant[component] += contact[8 * (component + 1) + node];
    passed =
        check(stick.projected && near(stick.gap, -0.01, 1.0e-12) && near(stick.pressure, 10.0, 1.0e-12)
                  && near(stick.tributary_area, 0.25, 1.0e-12) && near(stick.contact_force, 2.5, 1.0e-12)
                  && near(stick.tangential_traction, 1.0, 1.0e-12) && !stick.sliding && near(resultant[0], 0.0, 1.0e-12)
                  && near(resultant[1], 0.0, 1.0e-12) && near(resultant[2], 0.0, 1.0e-12),
            "three-dimensional node-to-face contact recovers pressure, tributary area, stick traction, and "
            "equal-and-opposite force")
        && passed;
    const fuelsim::CartesianContactPointValue sliding =
        fuelsim::compute_node_to_quad4_contact_value({1000.0, 0.05, false},
            contact_geometry,
            state,
            committed,
            history);
    passed = check(sliding.sliding && near(sliding.tangential_traction, 0.5, 1.0e-12)
                       && near(sliding.elastic_tangential_slip[1], 0.0005, 1.0e-12)
                       && near(sliding.friction_dissipation, 6.25e-5, 1.0e-12),
                 "three-dimensional Coulomb contact caps sliding traction, stores the vector elastic slip, and "
                 "integrates the dissipated sliding work")
             && passed;
    std::array<double, 32> direction{};
    for (std::size_t dof = 8; dof < direction.size(); ++dof)
        direction[dof] = std::sin(0.37 * static_cast<double>(dof + 1));
    constexpr double step = 1.0e-7;
    fuelsim::Quad4SurfaceContactLocalValues plus = state, minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Quad4SurfaceContactLocalResidual
        plus_contact =
            fuelsim::compute_node_to_quad4_contact(stick_properties, contact_geometry, plus, committed, history),
        minus_contact =
            fuelsim::compute_node_to_quad4_contact(stick_properties, contact_geometry, minus, committed, history);
    double maximum_error = 0.0, scale = 0.0;
    std::size_t maximum_row = 0;
    double maximum_analytic = 0.0, maximum_numerical = 0.0;
    for (std::size_t row = 8; row < 32; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 32; ++column)
            analytic += contact_jacobian[row * 32 + column] * direction[column];
        const double numerical = (plus_contact[row] - minus_contact[row]) / (2.0 * step);
        if (std::abs(analytic - numerical) > maximum_error) {
            maximum_error = std::abs(analytic - numerical);
            maximum_row = row;
            maximum_analytic = analytic;
            maximum_numerical = numerical;
        }
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    std::cout << "cartesian_contact_jacobian_relative_error=" << maximum_error / scale << '\n'
              << "cartesian_contact_jacobian_maximum_row=" << maximum_row << '\n'
              << "cartesian_contact_jacobian_maximum_analytic=" << maximum_analytic << '\n'
              << "cartesian_contact_jacobian_maximum_numerical=" << maximum_numerical << '\n';
    passed = check(maximum_error / scale < 2.0e-6,
                 "three-dimensional sticking contact automatic-differentiation Jacobian matches centered difference")
             && passed;

    const fuelsim::Quad4FaceQuadraturePoint finite_point =
        fuelsim::make_quad4_face_quadrature_point(secondary, -0.5, -0.5, 1.0);
    const fuelsim::Quad4FaceQuadraturePoint finite_normal_point =
        fuelsim::make_quad4_face_quadrature_point(secondary, -2.0 / 3.0, -2.0 / 3.0, 1.0);
    const fuelsim::Quad4ToQuad4MechanicalGeometry finite_geometry{secondary,
        primary,
        finite_point.shape,
        finite_point.derivative_xi,
        finite_point.derivative_eta,
        finite_normal_point.derivative_xi,
        finite_normal_point.derivative_eta,
        1.0,
        1.0,
        1.0};
    fuelsim::Quad4SurfaceContactLocalJacobian finite_jacobian{};
    const fuelsim::Quad4SurfaceContactLocalResidual finite_contact =
        fuelsim::compute_quad4_to_quad4_contact(stick_properties,
            finite_geometry,
            state,
            committed,
            history,
            &finite_jacobian);
    const fuelsim::CartesianContactPointValue finite_stick =
        fuelsim::compute_quad4_to_quad4_contact_value(stick_properties, finite_geometry, state, committed, history);
    const fuelsim::CartesianContactPointValue finite_sliding =
        fuelsim::compute_quad4_to_quad4_contact_value({1000.0, 0.05, false},
            finite_geometry,
            state,
            committed,
            history);
    resultant = {};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t node = 0; node < 8; ++node)
            resultant[component] += finite_contact[8 * (component + 1) + node];
    passed =
        check(finite_stick.projected && near(finite_stick.gap, -0.01, 1.0e-12)
                  && near(finite_stick.pressure, 10.0, 1.0e-12) && near(finite_stick.tributary_area, 0.25, 1.0e-12)
                  && near(finite_stick.contact_force, 2.5, 1.0e-12)
                  && near(finite_stick.tangential_traction, 1.0, 1.0e-12) && !finite_stick.sliding
                  && near(finite_stick.friction_dissipation, 0.0, 1.0e-12) && finite_sliding.sliding
                  && near(finite_sliding.friction_dissipation, 6.25e-5, 1.0e-12) && near(resultant[0], 0.0, 1.0e-12)
                  && near(resultant[1], 0.0, 1.0e-12) && near(resultant[2], 0.0, 1.0e-12),
            "HEX8 finite-sliding surface contact uses the node-centered area, current projection, and "
            "equal-and-opposite three-component force")
        && passed;
    fuelsim::ContactPointHistory accumulated_slip_history;
    accumulated_slip_history.cartesian_total_tangential_slip = {0.0, 2.0e-4, -3.0e-4};
    fuelsim::Quad4SurfaceContactLocalValues open_state = committed;
    const fuelsim::CartesianContactPointValue open_surface =
                                                  fuelsim::compute_quad4_to_quad4_contact_value(stick_properties,
                                                      finite_geometry,
                                                      open_state,
                                                      open_state,
                                                      accumulated_slip_history),
                                              open_node = fuelsim::compute_node_to_quad4_contact_value(stick_properties,
                                                  contact_geometry,
                                                  open_state,
                                                  open_state,
                                                  accumulated_slip_history);
    bool open_slip_is_retained =
        open_surface.projected && open_node.projected && open_surface.pressure == 0.0 && open_node.pressure == 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        open_slip_is_retained = open_slip_is_retained
                                && near(open_surface.tangential_slip[component],
                                    accumulated_slip_history.cartesian_total_tangential_slip[component],
                                    1.0e-12)
                                && near(open_node.tangential_slip[component],
                                    accumulated_slip_history.cartesian_total_tangential_slip[component],
                                    1.0e-12);
    passed = check(open_slip_is_retained,
                 "HEX8 accumulated total tangential slip remains constant while projected contact is open")
             && passed;
    plus = state;
    minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Quad4SurfaceContactLocalResidual
        finite_plus =
            fuelsim::compute_quad4_to_quad4_contact(stick_properties, finite_geometry, plus, committed, history),
        finite_minus =
            fuelsim::compute_quad4_to_quad4_contact(stick_properties, finite_geometry, minus, committed, history);
    maximum_error = 0.0;
    scale = 0.0;
    for (std::size_t row = 8; row < 32; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < 32; ++column)
            analytic += finite_jacobian[row * 32 + column] * direction[column];
        const double numerical = (finite_plus[row] - finite_minus[row]) / (2.0 * step);
        maximum_error = std::max(maximum_error, std::abs(analytic - numerical));
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    std::cout << "hex8_finite_sliding_contact_jacobian_relative_error=" << maximum_error / scale << '\n';
    passed = check(maximum_error / scale < 2.0e-6,
                 "HEX8 finite-sliding surface-contact automatic-differentiation Jacobian matches centered difference")
             && passed;

    fuelsim::Quad4FiniteRegionNormalGeometryJacobian finite_region_jacobian{};
    const fuelsim::Quad4FiniteRegionNormalGeometryValue finite_region =
        fuelsim::compute_quad4_finite_region_normal_geometry(finite_geometry, state, &finite_region_jacobian);
    const fuelsim::Quad4FiniteRegionNormalGeometryValue
        finite_region_plus = fuelsim::compute_quad4_finite_region_normal_geometry(finite_geometry, plus),
        finite_region_minus = fuelsim::compute_quad4_finite_region_normal_geometry(finite_geometry, minus);
    maximum_error = 0.0;
    scale = 0.0;
    for (std::size_t row = 0; row < fuelsim::quad4_finite_region_normal_geometry_output_count; ++row) {
        double analytic = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            analytic += finite_region_jacobian[row * direction.size() + column] * direction[column];
        const auto output = [row](const fuelsim::Quad4FiniteRegionNormalGeometryValue& value) {
            if (row == 0)
                return value.area;
            if (row == 1)
                return value.gap_integral;
            return value.unit_pressure_residual[row - 2];
        };
        const double numerical = (output(finite_region_plus) - output(finite_region_minus)) / (2.0 * step);
        maximum_error = std::max(maximum_error, std::abs(analytic - numerical));
        scale = std::max({scale, std::abs(analytic), std::abs(numerical)});
    }
    resultant = {};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t node = 0; node < 8; ++node)
            resultant[component] += finite_region.unit_pressure_residual[8 * (component + 1) + node];
    std::cout << "hex8_finite_region_normal_geometry_jacobian_relative_error=" << maximum_error / scale << '\n';
    passed = check(finite_region.projected && near(finite_region.area, 0.25, 1.0e-12)
                       && near(finite_region.gap_integral, -0.0025, 1.0e-12) && near(resultant[0], 0.0, 1.0e-12)
                       && near(resultant[1], 0.0, 1.0e-12) && near(resultant[2], 0.0, 1.0e-12)
                       && maximum_error / scale < 2.0e-6,
                 "HEX8 finite-region normal geometry integrates area and gap, conserves force, and matches a "
                 "centered directional difference")
             && passed;

    fuelsim::Quad4SurfaceContactLocalValues objective_committed{}, objective_current{};
    for (std::size_t node = 0; node < 4; ++node)
        objective_committed[8 + node] = 0.02;
    constexpr double angle = 0.55, slip_first = 2.0e-4, slip_second = -3.0e-4;
    const double cosine = std::cos(angle), sine = std::sin(angle);
    for (std::size_t node = 0; node < 8; ++node) {
        const fuelsim::CartesianPoint3& reference = node < 4 ? secondary[node] : primary[node - 4];
        const double committed_x = reference.x + objective_committed[8 + node],
                     committed_y = reference.y + objective_committed[16 + node];
        objective_current[8 + node] = cosine * committed_x - sine * committed_y - reference.x;
        objective_current[16 + node] = sine * committed_x + cosine * committed_y - reference.y;
        objective_current[24 + node] = objective_committed[24 + node];
    }
    fuelsim::ContactPointHistory objective_history;
    objective_history.cartesian_elastic_tangential_slip = {0.0, slip_first, slip_second};
    objective_history.cartesian_tangent_basis_initialized = true;
    objective_history.cartesian_contact_normal = {1.0, 0.0, 0.0};
    objective_history.cartesian_contact_tangent_first = {0.0, 1.0, 0.0};
    const fuelsim::CartesianContactPointValue objective =
        fuelsim::compute_quad4_to_quad4_contact_value({1000.0, 1.0, false, 1.0e-2},
            finite_geometry,
            objective_current,
            objective_committed,
            objective_history);
    const std::array<double, 3> expected_normal = {cosine, sine, 0.0}, expected_tangent = {-sine, cosine, 0.0},
                                expected_slip = {-slip_first * sine, slip_first * cosine, slip_second};
    bool objective_rotation = objective.projected && !objective.sliding;
    for (std::size_t component = 0; component < 3; ++component)
        objective_rotation = objective_rotation
                             && near(objective.normal[component], expected_normal[component], 1.0e-12)
                             && near(objective.tangent_first[component], expected_tangent[component], 1.0e-12)
                             && near(objective.elastic_tangential_slip[component], expected_slip[component], 1.0e-12);
    passed =
        check(objective_rotation
                  && near(std::hypot(objective.elastic_tangential_slip[0],
                              std::hypot(objective.elastic_tangential_slip[1], objective.elastic_tangential_slip[2])),
                      std::hypot(slip_first, slip_second),
                      1.0e-12),
            "HEX8 finite-sliding nonzero two-component elastic-slip history rotates objectively through a "
            "31.5-degree tangent-plane rotation")
        && passed;
    auto reversed_history = objective_history;
    for (double& component : reversed_history.cartesian_contact_tangent_first)
        component = -component;
    const auto reversed = fuelsim::compute_quad4_to_quad4_contact_value({1000.0, 1.0, false, 1.0e-2},
        finite_geometry,
        objective_current,
        objective_committed,
        reversed_history);
    const auto original_residual = fuelsim::compute_quad4_to_quad4_contact({1000.0, 1.0, false, 1.0e-2},
        finite_geometry,
        objective_current,
        objective_committed,
        objective_history);
    const auto reversed_residual = fuelsim::compute_quad4_to_quad4_contact({1000.0, 1.0, false, 1.0e-2},
        finite_geometry,
        objective_current,
        objective_committed,
        reversed_history);
    double residual_difference = 0.0, traction_difference = 0.0, history_difference = 0.0;
    for (std::size_t i = 0; i < original_residual.size(); ++i)
        residual_difference += std::abs(original_residual[i] - reversed_residual[i]);
    for (std::size_t i = 0; i < 3; ++i) {
        traction_difference +=
            std::abs(objective.tangential_traction_vector[i] - reversed.tangential_traction_vector[i]);
        history_difference += std::abs(objective.elastic_tangential_slip[i] - reversed.elastic_tangential_slip[i]);
    }
    passed =
        check(residual_difference > 1.0e-6 && traction_difference > 1.0e-6 && history_difference > 1.0e-6,
            "Reversing only the stored directed tangent changes residual, friction traction and updated nonzero slip")
        && passed;
    return passed;
}
} // namespace

int main() {
    return (test_cartesian_surface_contact_kernels()) ? 0 : 1;
}
