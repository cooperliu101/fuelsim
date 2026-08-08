#include "fuelsim/boundary.hpp"
#include "fuelsim/dof_map.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/material.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/quad4_rz.hpp"
#include "fuelsim/time_table.hpp"
#include "support/steady_fuel_cladding_problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double scaled_error(double actual, double expected) {
    return std::abs(actual - expected) /
           (1.0 + std::max(std::abs(actual), std::abs(expected)));
}

double relative_difference(double actual, double expected) {
    return std::abs(actual - expected) / std::abs(expected);
}

fuelsim::ThermoelasticProperties properties() {
    return {
        3824.0, 0.61, 2.0e11, 0.316, 1.0e-5, 600.0,
    };
}

bool test_mesh_and_geometry() {
    const double inner = 0.0;
    const double outer = 0.004;
    const double length = 0.01;
    const fuelsim::StructuredRzMesh mesh =
        fuelsim::StructuredRzMesh::make_annulus(inner, outer, length, 4, 3);

    bool passed = true;
    passed = check(mesh.nodes().size() == 20, "structured mesh node count") &&
             passed;
    passed =
        check(mesh.elements().size() == 12, "structured mesh element count") &&
        passed;

    double integrated_volume = 0.0;
    for (const fuelsim::Quad4Element& element : mesh.elements()) {
        fuelsim::Quad4Coordinates coordinates{};
        for (std::size_t node = 0; node < element.nodes.size(); ++node)
            coordinates[node] = mesh.nodes()[element.nodes[node]];
        const fuelsim::Quad4RzGeometry geometry =
            fuelsim::make_quad4_rz_geometry(coordinates);

        for (const fuelsim::RzQuadraturePoint& point : geometry.points) {
            double shape_sum = 0.0;
            double gradient_r_sum = 0.0;
            double gradient_z_sum = 0.0;
            for (std::size_t node = 0; node < 4; ++node) {
                shape_sum += point.shape[node];
                gradient_r_sum += point.gradient_r[node];
                gradient_z_sum += point.gradient_z[node];
            }
            passed = check(std::abs(shape_sum - 1.0) < 1.0e-14,
                           "Quad4 shape functions form a partition of unity") &&
                     passed;
            passed = check(std::abs(gradient_r_sum) < 1.0e-11 &&
                               std::abs(gradient_z_sum) < 1.0e-11,
                           "Quad4 physical gradients sum to zero") &&
                     passed;
            passed = check(point.radius > 0.0 && point.weighted_measure > 0.0,
                           "RZ quadrature measure is positive") &&
                     passed;
            integrated_volume += point.weighted_measure;
        }
    }

    const double exact_volume = pi * (outer * outer - inner * inner) * length;
    const double volume_error =
        std::abs(integrated_volume - exact_volume) / exact_volume;
    std::cout << "geometry_volume_relative_error=" << volume_error << '\n';
    passed = check(volume_error < 1.0e-13,
                   "RZ quadrature integrates annular volume") &&
             passed;
    return passed;
}

bool test_element_jacobian() {
    const fuelsim::Quad4Coordinates coordinates = {{
        {0.001, 0.0},
        {0.003, 0.0},
        {0.003, 0.004},
        {0.001, 0.004},
    }};
    const fuelsim::Quad4RzGeometry geometry =
        fuelsim::make_quad4_rz_geometry(coordinates);
    const fuelsim::Quad4RzThermoelasticKernel kernel(
        fuelsim::IsotropicThermoelasticMaterial(properties()), 2.0e8,
        fuelsim::StrainFormulation::small);

    const fuelsim::LocalValues state = {
        710.0,  680.0,  650.0, 690.0,   0.0,    2.0e-6,
        2.5e-6, 0.4e-6, 0.0,   -0.2e-6, 3.0e-6, 2.5e-6,
    };
    const fuelsim::LocalValues direction = {
        0.7,    -0.4,    0.3,     -0.6,   0.2e-6,  -0.4e-6,
        0.5e-6, -0.1e-6, -0.3e-6, 0.6e-6, -0.2e-6, 0.4e-6,
    };

    const fuelsim::LocalSystem system = kernel.linearize(geometry, state);

    constexpr double step = 1.0e-4;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        kernel.residual(geometry, plus);
    const fuelsim::LocalResidual minus_residual =
        kernel.residual(geometry, minus);

    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column) {
            ad_direction +=
                system.jacobian[row * fuelsim::quad4_local_dof_count + column] *
                direction[column];
        }
        const double finite_difference =
            (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                       "AD element Jacobian row " + std::to_string(row) +
                           " matches centered finite difference") &&
                 passed;
    }
    std::cout << "element_jacobian_maximum_scaled_error="
              << maximum_jacobian_error << '\n';

    double heat_to_displacement = 0.0;
    double temperature_to_mechanics = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        for (std::size_t column = 4; column < 12; ++column) {
            heat_to_displacement +=
                std::abs(system.jacobian[row * 12 + column]);
        }
    }
    for (std::size_t row = 4; row < 12; ++row) {
        for (std::size_t column = 0; column < 4; ++column) {
            temperature_to_mechanics +=
                std::abs(system.jacobian[row * 12 + column]);
        }
    }

    passed =
        check(heat_to_displacement == 0.0, "reference-mesh thermal residual "
                                           "does not depend on displacement") &&
        passed;
    passed =
        check(temperature_to_mechanics > 0.0,
              "thermal expansion produces temperature-mechanics coupling") &&
        passed;
    return passed;
}

bool test_finite_strain_kinematics_and_jacobian() {
    const fuelsim::Quad4Coordinates coordinates = {{
        {1.0, 0.0},
        {2.0, 0.0},
        {2.0, 1.0},
        {1.0, 1.0},
    }};
    const fuelsim::Quad4RzGeometry geometry =
        fuelsim::make_quad4_rz_geometry(coordinates);
    const fuelsim::Quad4RzThermoelasticKernel kernel(
        fuelsim::IsotropicThermoelasticMaterial(properties()), 0.0,
        fuelsim::StrainFormulation::finite);

    constexpr double radial_stretch = 1.08;
    constexpr double axial_stretch = 0.96;
    fuelsim::LocalValues uniform_state{};
    for (std::size_t node = 0; node < 4; ++node) {
        uniform_state[node] = 600.0;
        uniform_state[4 + node] =
            (radial_stretch - 1.0) * coordinates[node].r;
        uniform_state[8 + node] =
            (axial_stretch - 1.0) * coordinates[node].z;
    }
    fuelsim::LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < uniform_state.size(); ++dof)
        passive_state[dof] = uniform_state[dof];

    bool passed = true;
    double maximum_strain_error = 0.0;
    double maximum_measure_error = 0.0;
    for (const fuelsim::RzQuadraturePoint& point : geometry.points) {
        const fuelsim::AxisymmetricKinematics kinematics =
            fuelsim::evaluate_axisymmetric_kinematics(
                point, passive_state, fuelsim::StrainFormulation::finite);
        const auto taylor_increment = [](double stretch) {
            const double cinv = 1.0 / (stretch * stretch) - 1.0;
            return -0.5 * cinv + 0.25 * cinv * cinv;
        };
        maximum_strain_error =
            std::max({maximum_strain_error,
                      std::abs(kinematics.strain_rr.value() -
                               taylor_increment(radial_stretch)),
                      std::abs(kinematics.strain_zz.value() -
                               taylor_increment(axial_stretch)),
                      std::abs(kinematics.strain_hoop.value() -
                               taylor_increment(radial_stretch)),
                      std::abs(kinematics.strain_rz.value())});
        const double expected_measure =
            point.weighted_measure * radial_stretch * radial_stretch *
            axial_stretch;
        maximum_measure_error =
            std::max(maximum_measure_error,
                     std::abs(kinematics.weighted_measure.value() -
                              expected_measure) /
                         expected_measure);
    }
    passed = check(maximum_strain_error < 1.0e-14,
                   "finite RZ uniform stretches give MOOSE Taylor strain "
                   "increments") &&
             passed;
    passed = check(maximum_measure_error < 1.0e-14,
                   "finite RZ current measure follows the deformation "
                   "Jacobian") &&
             passed;

    constexpr double old_radial_stretch = 1.03;
    constexpr double old_axial_stretch = 0.98;
    fuelsim::LocalValues committed_state{};
    for (std::size_t node = 0; node < 4; ++node) {
        committed_state[node] = 600.0;
        committed_state[4 + node] =
            (old_radial_stretch - 1.0) * coordinates[node].r;
        committed_state[8 + node] =
            (old_axial_stretch - 1.0) * coordinates[node].z;
    }
    double maximum_incremental_error = 0.0;
    for (const fuelsim::RzQuadraturePoint& point : geometry.points) {
        const fuelsim::AxisymmetricKinematics kinematics =
            fuelsim::evaluate_axisymmetric_incremental_kinematics(
                point, passive_state, committed_state,
                fuelsim::StrainFormulation::finite);
        const auto taylor_increment = [](double stretch) {
            const double cinv = 1.0 / (stretch * stretch) - 1.0;
            return -0.5 * cinv + 0.25 * cinv * cinv;
        };
        maximum_incremental_error =
            std::max({maximum_incremental_error,
                      std::abs(kinematics.strain_rr.value() -
                               taylor_increment(radial_stretch /
                                                old_radial_stretch)),
                      std::abs(kinematics.strain_zz.value() -
                               taylor_increment(axial_stretch /
                                                old_axial_stretch)),
                      std::abs(kinematics.strain_hoop.value() -
                               taylor_increment(radial_stretch /
                                                old_radial_stretch)),
                      std::abs(kinematics.strain_rz.value())});
    }
    passed = check(maximum_incremental_error < 1.0e-14,
                   "finite RZ strain uses current-to-committed MOOSE "
                   "incremental deformation") &&
             passed;

    fuelsim::LocalValues state = uniform_state;
    state[5] += 0.025;
    state[6] += 0.015;
    state[10] -= 0.020;
    state[11] += 0.010;
    const fuelsim::LocalValues direction = {
        0.2,  -0.3, 0.4,  -0.1, 0.3,  -0.5,
        0.2,  0.4,  -0.2, 0.35, -0.45, 0.25,
    };
    const fuelsim::LocalSystem system = kernel.linearize(geometry, state);
    constexpr double step = 1.0e-5;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual = kernel.residual(geometry, plus);
    const fuelsim::LocalResidual minus_residual =
        kernel.residual(geometry, minus);
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction +=
                system.jacobian[row * fuelsim::quad4_local_dof_count + column] *
                direction[column];
        const double finite_difference =
            (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        maximum_jacobian_error =
            std::max(maximum_jacobian_error,
                     scaled_error(ad_direction, finite_difference));
    }
    passed = check(maximum_jacobian_error < 2.0e-7,
                   "finite RZ AD Jacobian matches centered finite "
                   "difference") &&
             passed;

    fuelsim::LocalValues inverted = uniform_state;
    for (std::size_t node = 0; node < 4; ++node)
        inverted[4 + node] = -1.1 * coordinates[node].r + 2.0;
    bool inversion_rejected = false;
    try {
        (void)kernel.residual(geometry, inverted);
    } catch (const std::domain_error&) {
        inversion_rejected = true;
    }
    passed = check(inversion_rejected,
                   "finite RZ rejects a nonpositive in-plane Jacobian while "
                   "the current radius remains positive") &&
             passed;
    fuelsim::LocalValues collapsed_radius = uniform_state;
    for (std::size_t node = 0; node < 4; ++node)
        collapsed_radius[4 + node] = -2.1;
    bool radius_rejected = false;
    try {
        (void)kernel.residual(geometry, collapsed_radius);
    } catch (const std::domain_error&) {
        radius_rejected = true;
    }
    passed = check(radius_rejected,
                   "finite RZ rejects nonpositive hoop stretch and current "
                   "radius while the in-plane Jacobian remains positive") &&
             passed;
    std::cout << "finite_strain_uniform_taylor_increment_error="
              << maximum_strain_error << '\n';
    std::cout << "finite_strain_current_measure_relative_error="
              << maximum_measure_error << '\n';
    std::cout << "finite_strain_incremental_taylor_error="
              << maximum_incremental_error << '\n';
    std::cout << "finite_strain_directional_jacobian_error="
              << maximum_jacobian_error << '\n';
    return passed;
}

bool test_heat_interface_case(const std::string& name,
                              const fuelsim::Line2RzGapHeatKernel& kernel,
                              const fuelsim::Line2RzHeatGeometry& geometry,
                              const fuelsim::LocalValues& state) {
    const fuelsim::LocalValues direction = {
        0.7,    -0.4,    0.3,     -0.6,   0.2e-6,  -0.4e-6,
        0.5e-6, -0.1e-6, -0.3e-6, 0.6e-6, -0.2e-6, 0.4e-6,
    };
    const fuelsim::LocalSystem system = kernel.linearize(geometry, state);

    constexpr double step = 1.0e-4;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }

    const fuelsim::LocalResidual plus_residual =
        kernel.residual(geometry, plus);
    const fuelsim::LocalResidual minus_residual =
        kernel.residual(geometry, minus);

    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column) {
            ad_direction +=
                system.jacobian[row * fuelsim::local_dof_count + column] *
                direction[column];
        }
        const double finite_difference =
            (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed =
            check(error < 1.0e-7, name + " interface AD Jacobian row " +
                                      std::to_string(row) +
                                      " matches centered finite difference") &&
            passed;
    }

    const fuelsim::LocalResidual residual = kernel.residual(geometry, state);
    double thermal_sum = 0.0;
    double thermal_scale = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        thermal_sum += residual[row];
        thermal_scale += std::abs(residual[row]);
    }
    passed = check(std::abs(thermal_sum) < 1.0e-13 * (1.0 + thermal_scale),
                   name + " interface conserves heat") &&
             passed;

    for (std::size_t row = 4; row < fuelsim::local_dof_count; ++row) {
        passed = check(residual[row] == 0.0,
                       name + " gap heat kernel has no mechanical residual") &&
                 passed;
    }

    std::cout << name << "_heat_jacobian_maximum_scaled_error="
              << maximum_jacobian_error << '\n';
    return passed;
}

bool test_contact_interface_case(
    const std::string& name, const fuelsim::NodeToLineRzContactKernel& kernel,
    const fuelsim::NodeToLineRzContactGeometry& geometry,
    const fuelsim::LocalValues& state,
    const fuelsim::ContactPointHistory& history = {}) {
    const fuelsim::LocalValues direction = {
        0.7,    -0.4,    0.3,     -0.6,   0.2e-6,  -0.4e-6,
        0.5e-6, -0.1e-6, -0.3e-6, 0.6e-6, -0.2e-6, 0.4e-6,
    };
    const fuelsim::LocalSystem system =
        kernel.linearize(geometry, state, state, history);

    constexpr double step = 1.0e-4;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        kernel.residual(geometry, plus, state, history);
    const fuelsim::LocalResidual minus_residual =
        kernel.residual(geometry, minus, state, history);

    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column) {
            ad_direction +=
                system.jacobian[row * fuelsim::local_dof_count + column] *
                direction[column];
        }
        const double finite_difference =
            (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed =
            check(error < 1.0e-7, name + " contact AD Jacobian row " +
                                      std::to_string(row) +
                                      " matches centered finite difference") &&
            passed;
    }

    const fuelsim::LocalResidual residual =
        kernel.residual(geometry, state, state, history);
    double radial_sum = 0.0;
    double radial_scale = 0.0;
    for (std::size_t row = 4; row < 8; ++row) {
        radial_sum += residual[row];
        radial_scale += std::abs(residual[row]);
    }
    passed = check(std::abs(radial_sum) < 1.0e-13 * (1.0 + radial_scale),
                   name + " contact conserves radial force") &&
             passed;
    double axial_sum = 0.0;
    double axial_scale = 0.0;
    for (std::size_t row = 8; row < 12; ++row) {
        axial_sum += residual[row];
        axial_scale += std::abs(residual[row]);
    }
    passed = check(std::abs(axial_sum) < 1.0e-13 * (1.0 + axial_scale),
                   name + " contact conserves axial force") &&
             passed;
    for (std::size_t row = 0; row < 4; ++row)
        passed = check(residual[row] == 0.0,
                       name + " contact has no thermal residual") &&
                 passed;
    std::cout << name << "_contact_jacobian_maximum_scaled_error="
              << maximum_jacobian_error << '\n';
    return passed;
}

bool test_friction_contact_case(
    const std::string& name,
    const fuelsim::NodeToLineRzContactKernel& kernel,
    const fuelsim::NodeToLineRzContactGeometry& geometry,
    const fuelsim::LocalValues& state,
    const fuelsim::LocalValues& committed_state,
    const fuelsim::ContactPointHistory& history) {
    const fuelsim::LocalValues direction = {
        0.0, 0.0, 0.0, 0.0, 0.2e-6,  -0.4e-6,
        0.5e-6, -0.1e-6, 0.0, 0.6e-6, -0.2e-6, 0.4e-6,
    };
    const fuelsim::LocalSystem system =
        kernel.linearize(geometry, state, committed_state, history);
    constexpr double step = 1.0e-4;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        kernel.residual(geometry, plus, committed_state, history);
    const fuelsim::LocalResidual minus_residual =
        kernel.residual(geometry, minus, committed_state, history);

    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction +=
                system.jacobian[row * fuelsim::local_dof_count + column] *
                direction[column];
        const double finite_difference =
            (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                       name + " friction AD Jacobian row " +
                           std::to_string(row) +
                           " matches centered finite difference") &&
                 passed;
    }
    const fuelsim::LocalResidual residual =
        kernel.residual(geometry, state, committed_state, history);
    double radial_sum = 0.0;
    double axial_sum = 0.0;
    double force_scale = 0.0;
    for (std::size_t row = 4; row < 8; ++row) {
        radial_sum += residual[row];
        force_scale += std::abs(residual[row]);
    }
    for (std::size_t row = 8; row < 12; ++row) {
        axial_sum += residual[row];
        force_scale += std::abs(residual[row]);
    }
    passed = check(std::abs(radial_sum) < 1.0e-13 * (1.0 + force_scale) &&
                       std::abs(axial_sum) < 1.0e-13 * (1.0 + force_scale),
                   name + " friction reaction is discretely conservative") &&
             passed;
    std::cout << name << "_friction_jacobian_maximum_scaled_error="
              << maximum_jacobian_error << '\n';
    return passed;
}

bool test_gap_heat_and_normal_contact() {
    const fuelsim::Line2InterfaceSideCoordinates fuel = {{
        {0.004120, 0.0},
        {0.004120, 0.001},
    }};
    const fuelsim::Line2InterfaceSideCoordinates cladding = {{
        {0.004122, 0.0},
        {0.004122, 0.001002},
    }};
    const fuelsim::Line2RzHeatGeometry heat_geometry =
        fuelsim::make_line2_rz_heat_geometry(fuel, cladding, 0.0);
    const fuelsim::Line2RzGapHeatKernel heat_kernel({0.4, 1.0e-6});
    const fuelsim::NodeToLineRzContactGeometry contact_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding, 1,
                                                       true, false, 0.0);
    const fuelsim::NodeToLineRzContactKernel contact_kernel({1.0e14});

    const fuelsim::LocalValues open_state = {
        750.0,  740.0, 610.0, 620.0, 0.2e-6, 0.3e-6,
        0.1e-6, 0.0,   0.0,   0.0,   0.0,    0.0,
    };
    const fuelsim::LocalValues minimum_gap_state = {
        750.0, 740.0, 610.0, 620.0, 1.4e-6, 1.5e-6,
        0.0,   0.0,   0.0,   0.0,   0.0,    0.0,
    };
    const fuelsim::LocalValues closed_state = {
        750.0, 740.0, 610.0, 620.0, 3.0e-6, 3.2e-6,
        0.0,   0.0,   0.0,   0.0,   0.0,    0.0,
    };

    bool passed = test_heat_interface_case("open", heat_kernel, heat_geometry,
                                           open_state);
    passed = test_heat_interface_case("minimum_gap", heat_kernel, heat_geometry,
                                      minimum_gap_state) &&
             passed;
    passed = test_contact_interface_case("closed", contact_kernel,
                                         contact_geometry, closed_state) &&
             passed;

    const fuelsim::HeatQuadratureValues open_values =
        heat_kernel.quadrature_values(heat_geometry, open_state);
    for (const fuelsim::HeatQuadratureValue& value : open_values) {
        passed = check(value.gap > 1.0e-6 && value.weighted_measure > 0.0,
                       "open gap heat point has positive gap and measure") &&
                 passed;
    }

    const fuelsim::ContactPointValue open_contact =
        contact_kernel.value(contact_geometry, open_state, open_state, {});
    const fuelsim::ContactPointValue closed_contact =
        contact_kernel.value(contact_geometry, closed_state, closed_state, {});
    passed = check(open_contact.projected && open_contact.gap > 0.0 &&
                       open_contact.pressure == 0.0,
                   "open NTS node projects with zero pressure") &&
             passed;
    passed = check(closed_contact.projected && closed_contact.gap < 0.0 &&
                       closed_contact.pressure > 0.0 &&
                       closed_contact.contact_force > 0.0,
                   "closed NTS node develops pressure and nodal force") &&
             passed;

    fuelsim::LocalValues outside_state = closed_state;
    outside_state[9] = 5.0e-6;
    const fuelsim::ContactPointValue outside_contact =
        contact_kernel.value(contact_geometry, outside_state, outside_state,
                             {});
    passed = check(!outside_contact.projected &&
                       outside_contact.contact_force == 0.0,
                   "out-of-segment NTS projection is inactive") &&
             passed;

    const fuelsim::NodeToLineRzContactGeometry radial_endpoint_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding, 0,
                                                       true, true, 0.0);
    fuelsim::LocalValues radial_endpoint_state = closed_state;
    radial_endpoint_state[8] = -5.0e-16;
    const fuelsim::ContactPointValue radial_endpoint =
        contact_kernel.value(radial_endpoint_geometry,
                             radial_endpoint_state, radial_endpoint_state, {});
    passed = check(radial_endpoint.projected,
                   "radial NTS reference endpoint retains projection after "
                   "roundoff-scale axial motion") &&
             passed;

    // A side that is vertical in the reference mesh need not remain vertical.
    // Unequal radial displacement of its primary endpoints must therefore use
    // the current inclined normal instead of a reference-geometry shortcut.
    fuelsim::LocalValues current_sloped_state = closed_state;
    current_sloped_state[4] = 20.0e-6;
    current_sloped_state[5] = 20.0e-6;
    current_sloped_state[6] = 0.0;
    current_sloped_state[7] = 10.0e-6;
    passed = test_heat_interface_case("current_sloped_reference_vertical",
                                      heat_kernel, heat_geometry,
                                      current_sloped_state) &&
             test_contact_interface_case(
                 "current_sloped_reference_vertical", contact_kernel,
                 contact_geometry, current_sloped_state) &&
             passed;
    const fuelsim::HeatQuadratureValues current_sloped_heat_values =
        heat_kernel.quadrature_values(heat_geometry, current_sloped_state);
    const fuelsim::LocalResidual current_sloped_heat_residual =
        heat_kernel.residual(heat_geometry, current_sloped_state);
    double expected_current_primary_node_1 = 0.0;
    double fixed_reference_primary_node_1 = 0.0;
    for (std::size_t q = 0; q < heat_geometry.points.size(); ++q) {
        const fuelsim::Line2RzHeatQuadraturePoint& point =
            heat_geometry.points[q];
        const double secondary_r =
            point.secondary_shape[0] *
                (heat_geometry.secondary_coordinates[0].r +
                 current_sloped_state[4]) +
            point.secondary_shape[1] *
                (heat_geometry.secondary_coordinates[1].r +
                 current_sloped_state[5]);
        const double secondary_z =
            point.secondary_shape[0] *
                (heat_geometry.secondary_coordinates[0].z +
                 current_sloped_state[8]) +
            point.secondary_shape[1] *
                (heat_geometry.secondary_coordinates[1].z +
                 current_sloped_state[9]);
        const double primary_r_0 =
            heat_geometry.primary_coordinates[0].r + current_sloped_state[6];
        const double primary_r_1 =
            heat_geometry.primary_coordinates[1].r + current_sloped_state[7];
        const double primary_z_0 =
            heat_geometry.primary_coordinates[0].z + current_sloped_state[10];
        const double primary_z_1 =
            heat_geometry.primary_coordinates[1].z + current_sloped_state[11];
        const double tangent_r = primary_r_1 - primary_r_0;
        const double tangent_z = primary_z_1 - primary_z_0;
        const double fraction =
            ((secondary_r - primary_r_0) * tangent_r +
             (secondary_z - primary_z_0) * tangent_z) /
            (tangent_r * tangent_r + tangent_z * tangent_z);
        const double heat_rate = current_sloped_heat_values[q].heat_flux *
                                 current_sloped_heat_values[q].weighted_measure;
        expected_current_primary_node_1 -= fraction * heat_rate;
        fixed_reference_primary_node_1 -= point.primary_shape[1] * heat_rate;
    }
    passed =
        check(scaled_error(current_sloped_heat_residual[3],
                           expected_current_primary_node_1) < 1.0e-13 &&
                  scaled_error(current_sloped_heat_residual[3],
                               fixed_reference_primary_node_1) > 1.0e-6,
              "thermal contact uses current primary projection rather than "
              "the reference interpolation fraction") &&
        passed;

    fuelsim::LocalValues thermal_projection_lost_state = open_state;
    thermal_projection_lost_state[10] = 0.45e-3;
    thermal_projection_lost_state[11] = -0.45e-3;
    passed = test_heat_interface_case("current_projection_endpoint_clamp",
                                      heat_kernel, heat_geometry,
                                      thermal_projection_lost_state) &&
             passed;
    const fuelsim::LocalResidual clamped_heat_residual = heat_kernel.residual(
        heat_geometry, thermal_projection_lost_state);
    passed = check(clamped_heat_residual[2] < 0.0 &&
                       clamped_heat_residual[3] < 0.0,
                   "thermal quadrature points that cross opposite ends remain "
                   "uniquely assigned to the owned primary endpoints") &&
             passed;
    const fuelsim::ContactPointValue current_sloped_contact =
        contact_kernel.value(contact_geometry, current_sloped_state,
                             current_sloped_state, {});
    const fuelsim::LocalResidual current_sloped_residual =
        contact_kernel.residual(contact_geometry, current_sloped_state,
                                current_sloped_state, {});
    passed = check(current_sloped_contact.projected &&
                       current_sloped_contact.pressure > 0.0 &&
                       std::abs(current_sloped_residual[9]) > 0.0,
                   "reference-vertical contact follows the inclined current "
                   "primary normal in both RZ equations") &&
             passed;

    const fuelsim::Line2InterfaceSideCoordinates vertex_secondary = {{
        {0.004000, 0.000500},
        {0.004000, 0.001000},
    }};
    const fuelsim::Line2InterfaceSideCoordinates vertex_primary_lower = {{
        {0.004002, 0.000000},
        {0.004002, 0.001000},
    }};
    const fuelsim::Line2InterfaceSideCoordinates vertex_primary_upper = {{
        {0.004002, 0.001000},
        {0.004002, 0.002000},
    }};
    const fuelsim::NodeToLineRzContactGeometry vertex_lower_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            vertex_secondary, vertex_primary_lower, 1, true, false, 0.0);
    const fuelsim::NodeToLineRzContactGeometry vertex_upper_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            vertex_secondary, vertex_primary_upper, 1, false, false, 0.0);
    fuelsim::LocalValues vertex_state{};
    vertex_state[5] = 3.0e-6;
    const fuelsim::ContactPointValue vertex_lower_reference =
        contact_kernel.value(vertex_lower_geometry, vertex_state,
                             vertex_state, {});
    const fuelsim::ContactPointValue vertex_upper_reference =
        contact_kernel.value(vertex_upper_geometry, vertex_state,
                             vertex_state, {});
    passed = check(!vertex_lower_reference.projected &&
                       vertex_upper_reference.projected,
                   "internal primary vertex has one reference owner") &&
             passed;
    vertex_state[9] = -1.0e-8;
    const fuelsim::ContactPointValue vertex_lower_slid =
        contact_kernel.value(vertex_lower_geometry, vertex_state,
                             vertex_state, {});
    const fuelsim::ContactPointValue vertex_upper_slid =
        contact_kernel.value(vertex_upper_geometry, vertex_state,
                             vertex_state, {});
    passed = check(vertex_lower_slid.projected &&
                       !vertex_upper_slid.projected &&
                       vertex_lower_slid.contact_force > 0.0,
                   "internal primary vertex slide transfers unique NTS "
                   "ownership without double force") &&
             passed;
    vertex_state[9] = -2.0e-3;
    const fuelsim::ContactPointValue vertex_upper_far =
        contact_kernel.value(vertex_upper_geometry, vertex_state,
                             vertex_state, {});
    passed = check(!vertex_upper_far.projected,
                   "reference endpoint does not mask a stale far projection") &&
             passed;

    const fuelsim::Line2InterfaceSideCoordinates general_secondary = {{
        {1.000000, 5.000000},
        {2.000000, 5.500000},
    }};
    const fuelsim::Line2InterfaceSideCoordinates general_primary_lower = {{
        {1.000000, 0.000000},
        {4.000000, 4.000000},
    }};
    const fuelsim::Line2InterfaceSideCoordinates general_primary_upper = {{
        {4.000000, 4.000000},
        {7.000000, 8.000000},
    }};
    const fuelsim::NodeToLineRzContactGeometry general_lower_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            general_secondary, general_primary_lower, 1, true, false, 0.0);
    const fuelsim::NodeToLineRzContactGeometry general_upper_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            general_secondary, general_primary_upper, 1, false, true, 0.0);
    fuelsim::LocalValues general_vertex_state{};
    const fuelsim::ContactPointValue general_lower_reference =
        contact_kernel.value(general_lower_geometry, general_vertex_state,
                             general_vertex_state, {});
    const fuelsim::ContactPointValue general_upper_reference =
        contact_kernel.value(general_upper_geometry, general_vertex_state,
                             general_vertex_state, {});
    passed = check(!general_lower_reference.projected &&
                       general_upper_reference.projected,
                   "sloped internal primary vertex has one reference owner") &&
             passed;
    general_vertex_state[5] = -1.0e-8;
    general_vertex_state[9] = -1.0e-8;
    const fuelsim::ContactPointValue general_lower_slid =
        contact_kernel.value(general_lower_geometry, general_vertex_state,
                             general_vertex_state, {});
    const fuelsim::ContactPointValue general_upper_slid =
        contact_kernel.value(general_upper_geometry, general_vertex_state,
                             general_vertex_state, {});
    passed = check(general_lower_slid.projected &&
                       !general_upper_slid.projected,
                   "sloped internal vertex slide keeps unique NTS ownership") &&
             passed;

    const fuelsim::Line2InterfaceSideCoordinates lower_pellet = {{
        {0.0, 0.001000},
        {0.004, 0.001000},
    }};
    const fuelsim::Line2InterfaceSideCoordinates upper_pellet = {{
        {0.0, 0.001002},
        {0.0041, 0.001002},
    }};
    const fuelsim::Line2RzHeatGeometry axial_heat_geometry =
        fuelsim::make_line2_rz_heat_geometry(lower_pellet, upper_pellet, 0.0);
    const fuelsim::NodeToLineRzContactGeometry axial_contact_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            lower_pellet, upper_pellet, 1, true, true, 0.0);
    const fuelsim::LocalValues axial_open_state = {
        750.0, 740.0, 610.0, 620.0, 0.0, 0.0,
        0.0,   0.0,   0.2e-6, 0.3e-6, 0.0, 0.0,
    };
    const fuelsim::LocalValues axial_closed_state = {
        750.0, 740.0, 610.0, 620.0, 0.0, 0.0,
        0.0,   0.0,   3.0e-6, 3.2e-6, 0.0, 0.0,
    };
    passed = test_heat_interface_case("axial_open", heat_kernel,
                                      axial_heat_geometry,
                                      axial_open_state) &&
             passed;
    passed = test_contact_interface_case(
                 "axial_closed", contact_kernel, axial_contact_geometry,
                 axial_closed_state) &&
             passed;
    const fuelsim::ContactPointValue axial_contact = contact_kernel.value(
        axial_contact_geometry, axial_closed_state, axial_closed_state, {});
    passed = check(axial_contact.projected && axial_contact.gap < 0.0 &&
                       axial_contact.pressure > 0.0,
                   "horizontal pellet faces develop axial contact") &&
             passed;

    const fuelsim::Line2InterfaceSideCoordinates sloped_secondary = {{
        {0.004200, 0.000200},
        {0.004800, 0.000800},
    }};
    const fuelsim::Line2InterfaceSideCoordinates sloped_primary = {{
        {0.004002, -0.000002},
        {0.005002, 0.000998},
    }};
    const fuelsim::Line2RzHeatGeometry sloped_heat_geometry =
        fuelsim::make_line2_rz_heat_geometry(sloped_secondary,
                                             sloped_primary, 0.0);
    const fuelsim::NodeToLineRzContactGeometry sloped_contact_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            sloped_secondary, sloped_primary, 1, true, true, 0.0);
    const fuelsim::LocalValues sloped_closed_state = {
        750.0, 740.0, 610.0, 620.0, 3.0e-6, 3.0e-6,
        0.0,   0.0,   -3.0e-6, -3.0e-6, 0.0, 0.0,
    };
    passed = test_heat_interface_case("sloped_open", heat_kernel,
                                      sloped_heat_geometry, open_state) &&
             passed;
    passed = test_contact_interface_case(
                 "sloped_closed", contact_kernel, sloped_contact_geometry,
                 sloped_closed_state) &&
             passed;
    const fuelsim::LocalResidual sloped_residual = contact_kernel.residual(
        sloped_contact_geometry, sloped_closed_state, sloped_closed_state, {});
    const fuelsim::ContactPointValue sloped_contact = contact_kernel.value(
        sloped_contact_geometry, sloped_closed_state, sloped_closed_state, {});
    passed = check(sloped_contact.projected && sloped_contact.gap < 0.0 &&
                       sloped_contact.pressure > 0.0 &&
                       std::abs(sloped_residual[5]) > 0.0 &&
                       std::abs(sloped_residual[9]) > 0.0,
                   "45-degree contact activates both normal components") &&
             passed;

    const fuelsim::NodeToLineRzContactKernel augmented_kernel(
        {1.0e14, 0.0, true});
    fuelsim::ContactPointHistory augmented_history;
    augmented_history.normal_multiplier = 2.0e6;
    const fuelsim::ContactPointValue augmented = augmented_kernel.value(
        contact_geometry, closed_state, closed_state, augmented_history);
    const double expected_augmented_pressure =
        augmented_history.normal_multiplier - 1.0e14 * augmented.gap;
    passed = check(
                 relative_difference(augmented.pressure,
                                     expected_augmented_pressure) < 1.0e-13,
                 "augmented contact adds the committed normal multiplier to "
                 "the penalty traction") &&
             test_contact_interface_case(
                 "augmented_normal", augmented_kernel, contact_geometry,
                 closed_state, augmented_history) &&
             passed;

    const fuelsim::NodeToLineRzContactKernel friction_kernel({1.0e14, 0.3});
    fuelsim::LocalValues sticking_state = closed_state;
    sticking_state[9] = 1.0e-7;
    const fuelsim::ContactPointValue sticking = friction_kernel.value(
        contact_geometry, sticking_state, closed_state, {});
    passed = check(!sticking.sliding &&
                       relative_difference(sticking.tangential_traction,
                                           1.0e7) <
                           1.0e-13,
                   "Coulomb contact matches the closed-form sticking "
                   "traction") &&
             passed;
    passed = test_friction_contact_case(
                 "sticking", friction_kernel, contact_geometry,
                 sticking_state, closed_state, {}) &&
             passed;

    fuelsim::LocalValues sliding_state = closed_state;
    sliding_state[9] = 6.0e-7;
    const fuelsim::ContactPointValue sliding = friction_kernel.value(
        contact_geometry, sliding_state, closed_state, {});
    const double sliding_limit = 0.3 * sliding.pressure;
    passed = check(sliding.sliding &&
                       relative_difference(sliding.tangential_traction,
                                           sliding_limit) < 1.0e-13 &&
                       relative_difference(sliding.elastic_tangential_slip,
                                           sliding_limit / 1.0e14) < 1.0e-13,
                   "Coulomb contact caps sliding traction and returns the "
                   "elastic slip") &&
             passed;
    passed = test_friction_contact_case(
                 "sliding", friction_kernel, contact_geometry, sliding_state,
                 closed_state, {}) &&
             passed;

    fuelsim::LocalValues resticking_state = closed_state;
    resticking_state[9] = -1.0e-7;
    const fuelsim::ContactPointHistory sliding_history =
        friction_kernel.trial_history(contact_geometry, sliding_state,
                                      closed_state, {});
    const fuelsim::ContactPointValue resticking = friction_kernel.value(
        contact_geometry, resticking_state, closed_state, sliding_history);
    passed = check(!resticking.sliding &&
                       resticking.tangential_traction > 0.0 &&
                       resticking.tangential_traction <
                           0.3 * resticking.pressure,
                   "Coulomb contact returns from sliding to sticking under "
                   "reverse tangential motion") &&
             passed;
    passed = test_friction_contact_case(
                 "sliding_to_sticking", friction_kernel, contact_geometry,
                 resticking_state, closed_state, sliding_history) &&
             passed;

    const fuelsim::NodeToLineRzContactKernel sloped_friction_kernel(
        {1.0e14, 0.3});
    const double sloped_tangent_component = 1.0 / std::sqrt(2.0);
    fuelsim::LocalValues sloped_sticking_state = sloped_closed_state;
    sloped_sticking_state[5] += 1.0e-7 * sloped_tangent_component;
    sloped_sticking_state[9] += 1.0e-7 * sloped_tangent_component;
    const fuelsim::ContactPointValue sloped_sticking =
        sloped_friction_kernel.value(sloped_contact_geometry,
                                     sloped_sticking_state,
                                     sloped_closed_state, {});
    passed = check(!sloped_sticking.sliding &&
                       relative_difference(sloped_sticking.tangential_traction,
                                           1.0e7) < 1.0e-13,
                   "sloped Coulomb contact matches the closed-form sticking "
                   "traction") &&
             passed;
    passed = test_friction_contact_case(
                 "sloped_sticking", sloped_friction_kernel,
                 sloped_contact_geometry, sloped_sticking_state,
                 sloped_closed_state, {}) &&
             passed;

    fuelsim::LocalValues sloped_sliding_state = sloped_closed_state;
    sloped_sliding_state[5] += 6.0e-7 * sloped_tangent_component;
    sloped_sliding_state[9] += 6.0e-7 * sloped_tangent_component;
    const fuelsim::ContactPointValue sloped_sliding =
        sloped_friction_kernel.value(sloped_contact_geometry,
                                     sloped_sliding_state,
                                     sloped_closed_state, {});
    const double sloped_sliding_limit = 0.3 * sloped_sliding.pressure;
    passed = check(sloped_sliding.sliding &&
                       relative_difference(sloped_sliding.tangential_traction,
                                           sloped_sliding_limit) < 1.0e-13 &&
                       relative_difference(
                           sloped_sliding.elastic_tangential_slip,
                           sloped_sliding_limit / 1.0e14) < 1.0e-13,
                   "sloped Coulomb contact caps sliding traction and returns "
                   "the elastic slip") &&
             passed;
    passed = test_friction_contact_case(
                 "sloped_sliding", sloped_friction_kernel,
                 sloped_contact_geometry, sloped_sliding_state,
                 sloped_closed_state, {}) &&
             passed;

    fuelsim::LocalValues sloped_resticking_state = sloped_closed_state;
    sloped_resticking_state[5] -= 1.0e-7 * sloped_tangent_component;
    sloped_resticking_state[9] -= 1.0e-7 * sloped_tangent_component;
    const fuelsim::ContactPointHistory sloped_sliding_history =
        sloped_friction_kernel.trial_history(sloped_contact_geometry,
                                             sloped_sliding_state,
                                             sloped_closed_state, {});
    const fuelsim::ContactPointValue sloped_resticking =
        sloped_friction_kernel.value(sloped_contact_geometry,
                                     sloped_resticking_state,
                                     sloped_closed_state,
                                     sloped_sliding_history);
    passed = check(!sloped_resticking.sliding &&
                       sloped_resticking.tangential_traction > 0.0 &&
                       sloped_resticking.tangential_traction <
                           0.3 * sloped_resticking.pressure,
                   "sloped Coulomb contact returns from sliding to sticking "
                   "under reverse tangential motion") &&
             passed;
    passed = test_friction_contact_case(
                 "sloped_sliding_to_sticking", sloped_friction_kernel,
                 sloped_contact_geometry, sloped_resticking_state,
                 sloped_closed_state, sloped_sliding_history) &&
             passed;

    const fuelsim::NodeToLineRzContactKernel explicit_zero_friction(
        {1.0e14, 0.0});
    const fuelsim::LocalResidual legacy_normal = contact_kernel.residual(
        contact_geometry, closed_state, closed_state, {});
    const fuelsim::LocalResidual zero_friction =
        explicit_zero_friction.residual(contact_geometry, closed_state,
                                        closed_state, {});
    passed = check(legacy_normal == zero_friction,
                   "mu equal to zero preserves every normal-contact residual "
                   "entry exactly") &&
             passed;
    return passed;
}

bool test_zero_gap_contact_orientation() {
    // Coincident fuel and cladding surfaces: the secondary node rides exactly
    // on the primary segment, so the raw reference normal gap is exactly
    // zero and the orientation must come from the material-side hint. The
    // fuel parent element sits at smaller radii, so its centroid is on the
    // negative side of the primary base normal (tangent_z, -tangent_r)/length
    // and the hint is negative.
    const fuelsim::Line2InterfaceSideCoordinates fuel = {{
        {0.004120, 0.0},
        {0.004120, 0.001},
    }};
    const fuelsim::Line2InterfaceSideCoordinates cladding_coincident = {{
        {0.004120, 0.0},
        {0.004120, 0.001002},
    }};
    constexpr double epsilon_gap = 1.0e-9;
    const fuelsim::Line2InterfaceSideCoordinates cladding_open = {{
        {0.004120 + epsilon_gap, 0.0},
        {0.004120 + epsilon_gap, 0.001002},
    }};
    const fuelsim::Line2InterfaceSideCoordinates cladding_inward = {{
        {0.004120 - epsilon_gap, 0.0},
        {0.004120 - epsilon_gap, 0.001002},
    }};

    bool passed = true;
    const fuelsim::NodeToLineRzContactGeometry zero_gap_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            fuel, cladding_coincident, 1, true, false, -1.0e-4);
    const fuelsim::NodeToLineRzContactGeometry opened_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            fuel, cladding_open, 1, true, false, 0.0);
    passed = check(zero_gap_geometry.normal_orientation == 1.0 &&
                       zero_gap_geometry.normal_orientation ==
                           opened_geometry.normal_orientation,
                   "zero-gap hint orientation matches the same geometry with "
                   "the gap opened by 1e-9 m") &&
             passed;

    const fuelsim::NodeToLineRzContactGeometry flipped_hint_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            fuel, cladding_coincident, 1, true, false, 1.0e-4);
    const fuelsim::NodeToLineRzContactGeometry inward_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            fuel, cladding_inward, 1, true, false, 0.0);
    passed = check(flipped_hint_geometry.normal_orientation == -1.0 &&
                       flipped_hint_geometry.normal_orientation ==
                           inward_geometry.normal_orientation,
                   "flipped material side flips the orientation, matching a "
                   "1e-9 m inward offset of the primary surface") &&
             passed;

    bool missing_hint_rejected = false;
    try {
        (void)fuelsim::make_node_to_line_rz_contact_geometry(
            fuel, cladding_coincident, 1, true, false, 0.0);
    } catch (const std::invalid_argument&) {
        missing_hint_rejected = true;
    }
    passed = check(missing_hint_rejected,
                   "on-segment zero gap without a hint remains an explicit "
                   "error") &&
             passed;

    const fuelsim::Line2RzHeatGeometry zero_gap_heat =
        fuelsim::make_line2_rz_heat_geometry(fuel, cladding_coincident,
                                             -1.0e-4);
    const fuelsim::Line2RzHeatGeometry opened_heat =
        fuelsim::make_line2_rz_heat_geometry(fuel, cladding_open, 0.0);
    for (std::size_t point = 0; point < zero_gap_heat.points.size(); ++point) {
        passed = check(zero_gap_heat.points[point].normal_orientation == 1.0 &&
                           zero_gap_heat.points[point].normal_orientation ==
                               opened_heat.points[point].normal_orientation,
                       "zero-gap heat quadrature orientation matches the "
                       "1e-9 m opened geometry") &&
                 passed;
    }
    bool heat_missing_hint_rejected = false;
    try {
        (void)fuelsim::make_line2_rz_heat_geometry(fuel, cladding_coincident,
                                                   0.0);
    } catch (const std::invalid_argument&) {
        heat_missing_hint_rejected = true;
    }
    passed = check(heat_missing_hint_rejected,
                   "zero-gap heat geometry without a hint remains an "
                   "explicit error") &&
             passed;

    // Branch checks avoid the nondifferentiable kink at gap exactly zero:
    // the open state leaves a +1e-6 m gap and the closed state penetrates
    // 3e-6 m, while the 1e-4 scaled finite-difference direction moves the
    // gap by at most 6e-11 m and cannot cross the kink.
    const fuelsim::NodeToLineRzContactKernel contact_kernel({1.0e14});
    const fuelsim::LocalValues zero_gap_open_state = {
        750.0, 740.0, 610.0, 620.0, 0.2e-6, -1.0e-6,
        0.0,   0.0,   0.0,   0.0,   0.0,    0.0,
    };
    const fuelsim::LocalValues zero_gap_closed_state = {
        750.0, 740.0, 610.0, 620.0, 3.0e-6, 3.2e-6,
        0.0,   0.0,   0.0,   0.0,   0.0,    0.0,
    };
    passed = test_contact_interface_case("zero_gap_open", contact_kernel,
                                         zero_gap_geometry,
                                         zero_gap_open_state) &&
             passed;
    passed = test_contact_interface_case("zero_gap_closed", contact_kernel,
                                         zero_gap_geometry,
                                         zero_gap_closed_state) &&
             passed;

    const fuelsim::ContactPointValue zero_gap_open =
        contact_kernel.value(zero_gap_geometry, zero_gap_open_state,
                             zero_gap_open_state, {});
    const fuelsim::ContactPointValue zero_gap_closed =
        contact_kernel.value(zero_gap_geometry, zero_gap_closed_state,
                             zero_gap_closed_state, {});
    passed =
        check(zero_gap_open.projected && zero_gap_open.gap > 0.0 &&
                  zero_gap_open.pressure == 0.0,
              "zero-gap geometry opens with zero pressure at +1e-6 m gap") &&
        check(zero_gap_closed.projected && zero_gap_closed.gap < 0.0 &&
                  scaled_error(zero_gap_closed.pressure, 3.2e8) < 1.0e-13,
              "zero-gap geometry develops penalty pressure over a 3.2e-6 m "
              "penetration") &&
        passed;

    // Gap continuity oracle for the general current-normal formula. Opening
    // the reference surfaces by 1e-9 m must increase the measured gap by
    // exactly 1e-9 m at the same state, which pins both the orientation and
    // the gap sign convention.
    const fuelsim::Line2InterfaceSideCoordinates sloped_fuel = {{
        {0.004120, 0.0},
        {0.0041205, 0.001},
    }};
    const fuelsim::NodeToLineRzContactGeometry sloped_zero_gap_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            sloped_fuel, cladding_coincident, 0, true, false, -1.0e-4);
    const fuelsim::NodeToLineRzContactGeometry sloped_opened_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(
            sloped_fuel, cladding_open, 0, true, false, 0.0);
    passed = check(sloped_zero_gap_geometry.normal_orientation == 1.0 &&
                       sloped_zero_gap_geometry.normal_orientation ==
                           sloped_opened_geometry.normal_orientation,
                   "sloped zero-gap hint orientation matches the 1e-9 m "
                   "opened geometry") &&
             passed;

    const fuelsim::LocalValues oracle_closed_state = {
        750.0, 740.0, 610.0, 620.0, 3.0e-6, 3.0e-6,
        0.0,   0.0,   1.0e-6, 0.0,   0.0,    0.0,
    };
    passed = test_contact_interface_case("zero_gap_sloped_closed",
                                         contact_kernel,
                                         sloped_zero_gap_geometry,
                                         oracle_closed_state) &&
             passed;
    const fuelsim::ContactPointValue sloped_zero_value =
        contact_kernel.value(sloped_zero_gap_geometry, oracle_closed_state,
                             oracle_closed_state, {});
    const fuelsim::ContactPointValue sloped_opened_value =
        contact_kernel.value(sloped_opened_geometry, oracle_closed_state,
                             oracle_closed_state, {});
    const double gap_oracle_error =
        std::abs((sloped_opened_value.gap - sloped_zero_value.gap) -
                 epsilon_gap);
    std::cout << "zero_gap_oracle_coincident_gap=" << sloped_zero_value.gap
              << '\n'
              << "zero_gap_oracle_opened_gap=" << sloped_opened_value.gap
              << '\n'
              << "zero_gap_oracle_gap_continuity_error=" << gap_oracle_error
              << '\n';
    passed =
        check(sloped_zero_value.projected && sloped_opened_value.projected &&
                  gap_oracle_error < 1.0e-15,
              "opening the reference surfaces by 1e-9 m increases the "
              "measured gap by exactly 1e-9 m at the same state") &&
        check(scaled_error(sloped_zero_value.pressure,
                           -1.0e14 * sloped_zero_value.gap) < 1.0e-13 &&
                  scaled_error(sloped_opened_value.pressure,
                               -1.0e14 * sloped_opened_value.gap) < 1.0e-13 &&
                  sloped_zero_value.pressure > 2.9e8 &&
                  sloped_opened_value.pressure > 2.9e8,
              "both geometries follow the penalty law on their own closed "
              "branch gaps near 3e8 Pa") &&
        passed;

    // The gap heat kernel has no kink at zero gap; run its two smooth
    // branches on the coincident geometry: +2e-6 m gap above the 1e-6 m
    // minimum gap and a -3e-6 m penetration on the minimum-gap branch.
    const fuelsim::Line2RzGapHeatKernel heat_kernel({0.4, 1.0e-6});
    const fuelsim::LocalValues zero_gap_heat_open_state = {
        750.0, 740.0, 610.0, 620.0, -2.0e-6, -2.0e-6,
        0.0,   0.0,   0.0,    0.0,   0.0,     0.0,
    };
    passed = test_heat_interface_case("zero_gap_heat_open", heat_kernel,
                                      zero_gap_heat,
                                      zero_gap_heat_open_state) &&
             passed;
    passed = test_heat_interface_case("zero_gap_heat_minimum", heat_kernel,
                                      zero_gap_heat, zero_gap_closed_state) &&
             passed;
    const fuelsim::HeatQuadratureValues zero_gap_heat_values =
        heat_kernel.quadrature_values(zero_gap_heat,
                                      zero_gap_heat_open_state);
    for (const fuelsim::HeatQuadratureValue& value : zero_gap_heat_values) {
        passed = check(value.gap > 1.0e-6 && value.weighted_measure > 0.0,
                       "zero-gap heat point reports the +2e-6 m open gap") &&
                 passed;
    }
    return passed;
}

bool test_m1_dof_layout() {
    const fuelsim::SteadyFuelCladdingParameters parameters = {
        0.004,
        0.0041,
        0.0046,
        0.010,
        0.01002,
        2,
        1,
        2,
        properties(),
        {
            0.0,
            16.0,
            75.0e9,
            0.3,
            5.0e-6,
            600.0,
        },
        1.0e8,
        600.0,
        600.0,
        0.4,
        1.0e-6,
        1.0e14,
    };
    fuelsim::SteadyFuelCladdingProblem problem(parameters);

    bool passed = true;
    passed =
        check(problem.thermal_interface_count() == parameters.axial_elements,
              "M1 has one STS heat contribution per fuel axial element") &&
        passed;
    passed = check(problem.contact_contribution_count() >
                       problem.thermal_interface_count(),
                   "M1 has local NTS candidate contributions") &&
             passed;
    passed = check(problem.dof_count() ==
                       3 * (problem.fuel_mesh().nodes().size() +
                            problem.cladding_mesh().nodes().size()),
                   "M1 uses one field-major map for both independent meshes") &&
             passed;

    const fuelsim::LocalDofs interface = problem.thermal_interface_dofs(0);
    const std::size_t fuel_outer =
        problem.fuel_mesh().node_id(problem.fuel_mesh().radial_elements(), 0);
    const std::size_t cladding_inner = problem.cladding_mesh().node_id(0, 0);
    passed = check(interface[0] == problem.dof_map().temperature(
                                       problem.fuel_global_node(fuel_outer)) &&
                       interface[2] ==
                           problem.dof_map().temperature(
                               problem.cladding_global_node(cladding_inner)),
                   "M1 interface DOFs preserve fuel/cladding node ownership") &&
             passed;
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        problem.summarize_contact_nodes(problem.initial_state());
    passed = check(std::all_of(contact_nodes.begin(), contact_nodes.end(),
                               [](const fuelsim::ContactNodeSummary& node) {
                                   return node.projected;
                               }),
                   "taller cladding contains every initial NTS projection") &&
             passed;

    const fuelsim::LocalValues initial_element_state =
        problem.contribution_state(0, problem.initial_state());
    const fuelsim::LocalResidual source_residual =
        problem.contribution_residual(0, initial_element_state);
    problem.set_volumetric_heat_source(2.0 * parameters.volumetric_heat_source);
    const fuelsim::LocalResidual doubled_source_residual =
        problem.contribution_residual(0, initial_element_state);
    for (std::size_t node = 0; node < fuelsim::quad4_node_count; ++node) {
        passed = check(std::abs(doubled_source_residual[node] -
                                2.0 * source_residual[node]) <
                           1.0e-12 * (1.0 + std::abs(source_residual[node])),
                       "M1 updates heat loading without rebuilding geometry") &&
                 passed;
    }
    passed = check(problem.parameters().volumetric_heat_source ==
                           2.0 * parameters.volumetric_heat_source &&
                       problem.fuel_kernel().volumetric_heat_source() ==
                           problem.parameters().volumetric_heat_source,
                   "M1 heat-source parameter and kernel remain synchronized") &&
             passed;
    return passed;
}

bool test_time_table_and_convection() {
    const fuelsim::PiecewiseLinearTimeTable table("power", {0.0, 2.0, 5.0},
                                                  {0.0, 1.0, 0.4});
    bool passed =
        check(table.value(0.0) == 0.0 && table.value(1.0) == 0.5 &&
                  table.value(3.0) == 0.8 && table.value(8.0) == 0.4,
              "piecewise-linear table interpolates and holds endpoints");

    const fuelsim::Line2RzConvectionGeometry geometry =
        fuelsim::make_line2_rz_convection_geometry(
            {{{0.005, 0.0}, {0.005, 0.01}}}, {{1, 2}});
    const fuelsim::Line2RzConvectionKernel kernel({1000.0, 500.0});
    const fuelsim::LocalValues state = {
        590.0, 600.0, 600.0, 610.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    };
    const fuelsim::LocalValues direction = {
        0.2, -0.7, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    };
    const fuelsim::LocalSystem system = kernel.linearize(geometry, state);
    const double expected_heat = 1000.0 * 100.0 * 2.0 * pi * 0.005 * 0.01;
    passed = check(scaled_error(system.residual[1] + system.residual[2],
                                expected_heat) < 1.0e-13,
                   "convection integrates the RZ surface heat loss") &&
             passed;

    constexpr double step = 1.0e-4;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        kernel.residual(geometry, plus);
    const fuelsim::LocalResidual minus_residual =
        kernel.residual(geometry, minus);
    double maximum_error = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row) {
        double tangent = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column)
            tangent += system.jacobian[row * state.size() + column] *
                       direction[column];
        const double finite_difference =
            (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        maximum_error =
            std::max(maximum_error, scaled_error(tangent, finite_difference));
    }
    std::cout << "convection_directional_jacobian_error=" << maximum_error
              << '\n';
    passed = check(maximum_error < 1.0e-10,
                   "convection AD Jacobian matches centered differences") &&
             passed;
    return passed;
}

bool test_follower_pressure() {
    const fuelsim::Line2RzPressureGeometry geometry =
        fuelsim::make_line2_rz_pressure_geometry(
            {{{0.005, 0.0}, {0.005, 0.01}}}, {{1, 2}});
    constexpr double pressure = 3.0e6;
    const fuelsim::Line2RzPressureKernel follower({pressure, true});
    const fuelsim::LocalValues state = {
        600.0, 600.0, 600.0, 600.0,
        0.0,   0.001, 0.002, 0.0,
        0.0,   0.0004, -0.0002, 0.0,
    };
    const fuelsim::LocalValues direction = {
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.3, -0.2, 0.0,
        0.0, -0.4, 0.5, 0.0,
    };
    const fuelsim::LocalSystem system = follower.linearize(geometry, state);
    const double first_radius = 0.006;
    const double second_radius = 0.007;
    const double delta_radius = second_radius - first_radius;
    const double delta_axial = 0.0098 - 0.0004;
    const double expected_radial =
        pi * pressure * delta_axial * (first_radius + second_radius);
    const double expected_axial =
        -pi * pressure * delta_radius * (first_radius + second_radius);
    bool passed =
        check(scaled_error(system.residual[5] + system.residual[6],
                           expected_radial) < 1.0e-13 &&
                  scaled_error(system.residual[9] + system.residual[10],
                               expected_axial) < 1.0e-13,
              "follower pressure uses current RZ radius and outward normal");

    constexpr double step = 1.0e-7;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual = follower.residual(geometry, plus);
    const fuelsim::LocalResidual minus_residual =
        follower.residual(geometry, minus);
    double maximum_error = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row) {
        double tangent = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column)
            tangent += system.jacobian[row * state.size() + column] *
                       direction[column];
        const double finite_difference =
            (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        maximum_error =
            std::max(maximum_error, scaled_error(tangent, finite_difference));
    }
    std::cout << "follower_pressure_directional_jacobian_error="
              << maximum_error << '\n';
    passed = check(maximum_error < 1.0e-8,
                   "follower-pressure AD Jacobian matches centered "
                   "differences") &&
             passed;

    const fuelsim::Line2RzPressureKernel dead({pressure, false});
    const fuelsim::LocalSystem dead_system = dead.linearize(geometry, state);
    const double maximum_dead_tangent = *std::max_element(
        dead_system.jacobian.begin(), dead_system.jacobian.end(),
        [](double left, double right) {
            return std::abs(left) < std::abs(right);
        });
    passed = check(maximum_dead_tangent == 0.0,
                   "reference pressure has an exactly zero geometric "
                   "tangent") &&
             passed;
    return passed;
}

bool test_current_configuration_traction() {
    const fuelsim::Line2RzTractionGeometry geometry =
        fuelsim::make_line2_rz_traction_geometry(
            {{{0.005, 0.0}, {0.005, 0.01}}}, {{1, 2}});
    constexpr double traction = 2.0e6;
    const fuelsim::Line2RzTractionKernel current(
        {fuelsim::TractionComponent::axial, traction, true});
    const fuelsim::LocalValues state = {
        600.0, 600.0, 600.0, 600.0,
        0.0,   0.001, 0.002, 0.0,
        0.0,   0.0004, -0.0002, 0.0,
    };
    const fuelsim::LocalValues direction = {
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.3, -0.2, 0.0,
        0.0, -0.4, 0.5, 0.0,
    };
    const fuelsim::LocalSystem system = current.linearize(geometry, state);
    const double current_length = std::hypot(0.001, 0.0094);
    const double expected_axial =
        -2.0 * pi * 0.0065 * current_length * traction;
    bool passed =
        check(scaled_error(system.residual[9] + system.residual[10],
                           expected_axial) < 1.0e-13,
              "current-configuration component traction uses current RZ "
              "surface measure");

    constexpr double step = 1.0e-7;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        current.residual(geometry, plus);
    const fuelsim::LocalResidual minus_residual =
        current.residual(geometry, minus);
    double maximum_error = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row) {
        double tangent = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column)
            tangent += system.jacobian[row * state.size() + column] *
                       direction[column];
        const double finite_difference =
            (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        maximum_error =
            std::max(maximum_error, scaled_error(tangent, finite_difference));
    }
    std::cout << "current_traction_directional_jacobian_error="
              << maximum_error << '\n';
    passed = check(maximum_error < 1.0e-8,
                   "current-configuration traction AD Jacobian matches "
                   "centered differences") &&
             passed;

    const fuelsim::Line2RzTractionKernel reference(
        {fuelsim::TractionComponent::axial, traction, false});
    const fuelsim::LocalSystem reference_system =
        reference.linearize(geometry, state);
    const double maximum_reference_tangent = *std::max_element(
        reference_system.jacobian.begin(), reference_system.jacobian.end(),
        [](double left, double right) {
            return std::abs(left) < std::abs(right);
        });
    return check(maximum_reference_tangent == 0.0,
                 "reference-configuration traction has zero geometric "
                 "tangent") &&
           passed;
}

bool test_temperature_active_thermoelastic_properties() {
    fuelsim::ThermoelasticProperties active_properties = properties();
    active_properties.young_modulus_temperature_coefficient = -8.0e7;
    active_properties.poisson_ratio_temperature_coefficient = 2.0e-5;
    active_properties.thermal_expansion_temperature_coefficient = 3.0e-9;
    const fuelsim::IsotropicThermoelasticMaterial material(active_properties);
    constexpr double temperature = 725.0;
    const adlite::Scalar active_temperature =
        adlite::Scalar::independent(temperature, 0, 1);
    const fuelsim::AxisymmetricStress active = material.stress(
        1.1e-3, -0.4e-3, 0.2e-3, 0.3e-3, active_temperature);
    constexpr double step = 1.0e-3;
    const fuelsim::AxisymmetricStress plus = material.stress(
        1.1e-3, -0.4e-3, 0.2e-3, 0.3e-3, temperature + step);
    const fuelsim::AxisymmetricStress minus = material.stress(
        1.1e-3, -0.4e-3, 0.2e-3, 0.3e-3, temperature - step);
    const std::array<double, 4> analytic = {
        active.rr.derivative(0), active.zz.derivative(0),
        active.hoop.derivative(0), active.rz.derivative(0)};
    const std::array<double, 4> finite_difference = {
        (plus.rr.value() - minus.rr.value()) / (2.0 * step),
        (plus.zz.value() - minus.zz.value()) / (2.0 * step),
        (plus.hoop.value() - minus.hoop.value()) / (2.0 * step),
        (plus.rz.value() - minus.rz.value()) / (2.0 * step)};
    double maximum_error = 0.0;
    for (std::size_t component = 0; component < analytic.size(); ++component)
        maximum_error =
            std::max(maximum_error,
                     scaled_error(analytic[component],
                                  finite_difference[component]));
    std::cout << "active_thermoelastic_temperature_tangent_error="
              << maximum_error << '\n';
    return check(maximum_error < 1.0e-8,
                 "temperature-dependent thermoelastic AD tangent matches "
                 "centered differences");
}

} // namespace

int main() {
    std::cout << std::scientific << std::setprecision(12);
    bool passed = true;
    passed = test_mesh_and_geometry() && passed;
    passed = test_element_jacobian() && passed;
    passed = test_finite_strain_kinematics_and_jacobian() && passed;
    passed = test_gap_heat_and_normal_contact() && passed;
    passed = test_zero_gap_contact_orientation() && passed;
    passed = test_m1_dof_layout() && passed;
    passed = test_time_table_and_convection() && passed;
    passed = test_follower_pressure() && passed;
    passed = test_current_configuration_traction() && passed;
    passed = test_temperature_active_thermoelastic_properties() && passed;

    if (!passed)
        return 1;
    std::cout << "[PASS] fuelsim core geometry, DOF, and AD Jacobian tests\n";
    return 0;
}
#include "fuelsim/boundary.hpp"
