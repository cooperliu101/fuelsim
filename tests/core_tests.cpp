#include "fuelsim/dof_map.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/material.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/quad4_rz.hpp"
#include "support/steady_fuel_cladding_problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
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
        fuelsim::IsotropicThermoelasticMaterial(properties()), 2.0e8);

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
            check(error < 1.0e-7, name + " contact AD Jacobian row " +
                                      std::to_string(row) +
                                      " matches centered finite difference") &&
            passed;
    }

    const fuelsim::LocalResidual residual = kernel.residual(geometry, state);
    double radial_sum = 0.0;
    double radial_scale = 0.0;
    for (std::size_t row = 4; row < 8; ++row) {
        radial_sum += residual[row];
        radial_scale += std::abs(residual[row]);
    }
    passed = check(std::abs(radial_sum) < 1.0e-13 * (1.0 + radial_scale),
                   name + " contact conserves radial force") &&
             passed;
    for (std::size_t row = 0; row < 4; ++row)
        passed = check(residual[row] == 0.0,
                       name + " contact has no thermal residual") &&
                 passed;
    for (std::size_t row = 8; row < fuelsim::local_dof_count; ++row)
        passed = check(residual[row] == 0.0,
                       name + " frictionless contact has no axial residual") &&
                 passed;

    std::cout << name << "_contact_jacobian_maximum_scaled_error="
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
        fuelsim::make_line2_rz_heat_geometry(fuel, cladding);
    const fuelsim::Line2RzGapHeatKernel heat_kernel({0.4, 1.0e-6});
    const fuelsim::NodeToLineRzContactGeometry contact_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding, 1,
                                                       false);
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
        contact_kernel.value(contact_geometry, open_state);
    const fuelsim::ContactPointValue closed_contact =
        contact_kernel.value(contact_geometry, closed_state);
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
        contact_kernel.value(contact_geometry, outside_state);
    passed = check(!outside_contact.projected &&
                       outside_contact.contact_force == 0.0,
                   "out-of-segment NTS projection is inactive") &&
             passed;
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

} // namespace

int main() {
    std::cout << std::scientific << std::setprecision(12);
    bool passed = true;
    passed = test_mesh_and_geometry() && passed;
    passed = test_element_jacobian() && passed;
    passed = test_gap_heat_and_normal_contact() && passed;
    passed = test_m1_dof_layout() && passed;

    if (!passed)
        return 1;
    std::cout << "[PASS] fuelsim core geometry, DOF, and AD Jacobian tests\n";
    return 0;
}
