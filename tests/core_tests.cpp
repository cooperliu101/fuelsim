#include "fuelsim/dof_map.hpp"
#include "fuelsim/interface.hpp"
#include "fuelsim/m1_problem.hpp"
#include "fuelsim/material.hpp"
#include "fuelsim/mesh.hpp"
#include "fuelsim/quad4_rz.hpp"

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

bool test_interface_case(const std::string& name,
                         const fuelsim::Line2RzGapContactKernel& kernel,
                         const fuelsim::Line2RzInterfaceGeometry& geometry,
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
    double radial_sum = 0.0;
    double radial_scale = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        thermal_sum += residual[row];
        thermal_scale += std::abs(residual[row]);
    }
    for (std::size_t row = 4; row < 8; ++row) {
        radial_sum += residual[row];
        radial_scale += std::abs(residual[row]);
    }
    passed = check(std::abs(thermal_sum) < 1.0e-13 * (1.0 + thermal_scale),
                   name + " interface conserves heat") &&
             passed;
    passed = check(std::abs(radial_sum) < 1.0e-13 * (1.0 + radial_scale),
                   name + " interface conserves radial force") &&
             passed;

    for (std::size_t row = 8; row < fuelsim::local_dof_count; ++row) {
        passed = check(residual[row] == 0.0,
                       name + " frictionless interface has zero axial "
                              "residual") &&
                 passed;
    }

    std::cout << name << "_interface_jacobian_maximum_scaled_error="
              << maximum_jacobian_error << '\n';
    return passed;
}

bool test_gap_contact_interface() {
    const fuelsim::Line2InterfaceSideCoordinates fuel = {{
        {0.004120, 0.0},
        {0.004120, 0.010},
    }};
    const fuelsim::Line2InterfaceSideCoordinates cladding = {{
        {0.004122, 0.0},
        {0.004122, 0.010},
    }};
    const fuelsim::Line2RzInterfaceGeometry geometry =
        fuelsim::make_line2_rz_interface_geometry(fuel, cladding);
    const fuelsim::Line2RzGapContactKernel kernel({0.4, 1.0e-6, 1.0e14});

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

    bool passed = test_interface_case("open", kernel, geometry, open_state);
    passed = test_interface_case("minimum_gap", kernel, geometry,
                                 minimum_gap_state) &&
             passed;
    passed =
        test_interface_case("closed", kernel, geometry, closed_state) && passed;

    const fuelsim::InterfaceQuadratureValues open_values =
        kernel.quadrature_values(geometry, open_state);
    const fuelsim::InterfaceQuadratureValues closed_values =
        kernel.quadrature_values(geometry, closed_state);
    for (const fuelsim::InterfaceQuadratureValue& value : open_values) {
        passed = check(value.gap > 1.0e-6 && value.pressure == 0.0,
                       "open interface has positive gap and zero pressure") &&
                 passed;
    }
    for (const fuelsim::InterfaceQuadratureValue& value : closed_values) {
        passed = check(value.gap < 0.0 && value.pressure > 0.0,
                       "closed interface has penetration pressure") &&
                 passed;
    }
    return passed;
}

bool test_m1_dof_layout() {
    const fuelsim::M1Parameters parameters = {
        0.004,
        0.0041,
        0.0046,
        0.010,
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
    const fuelsim::M1Problem problem(parameters);

    bool passed = true;
    passed = check(problem.interface_count() == parameters.axial_elements,
                   "M1 has one interface contribution per axial element") &&
             passed;
    passed = check(problem.dof_count() ==
                       3 * (problem.fuel_mesh().nodes().size() +
                            problem.cladding_mesh().nodes().size()),
                   "M1 uses one field-major map for both independent meshes") &&
             passed;

    const fuelsim::LocalDofs interface = problem.interface_dofs(0);
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
    return passed;
}

} // namespace

int main() {
    std::cout << std::scientific << std::setprecision(12);
    bool passed = true;
    passed = test_mesh_and_geometry() && passed;
    passed = test_element_jacobian() && passed;
    passed = test_gap_contact_interface() && passed;
    passed = test_m1_dof_layout() && passed;

    if (!passed)
        return 1;
    std::cout << "[PASS] fuelsim core geometry, DOF, and AD Jacobian tests\n";
    return 0;
}
