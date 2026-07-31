#include "fuelsim/m2_problem.hpp"
#include "fuelsim/m2_solver.hpp"
#include "fuelsim/petsc_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double relative_error(double actual, double expected) {
    if (!std::isfinite(actual) || !std::isfinite(expected) || expected == 0.0)
        return std::numeric_limits<double>::infinity();
    return std::abs(actual - expected) / std::abs(expected);
}

fuelsim::M2Parameters pcmi_parameters() {
    const fuelsim::M1Parameters base = {
        0.004120,
        0.004121,
        0.004692,
        0.010,
        0.010020,
        6,
        2,
        4,
        {
            3824.0,
            0.61,
            2.0e11,
            0.316,
            10.0e-6,
            600.0,
        },
        {
            0.0,
            16.0,
            75.0e9,
            0.3,
            0.0,
            600.0,
        },
        2.0e8,
        600.0,
        600.0,
        0.4,
        1.0e-6,
        1.0e14,
    };
    const fuelsim::TransientInelasticProperties fuel = {
        10970.0,         300.0,      fuelsim::InelasticBehavior::elastic,
        {0.0, 1.0, 1.0}, {1.0, 0.0},
    };
    const fuelsim::TransientInelasticProperties cladding = {
        6500.0,
        330.0,
        fuelsim::InelasticBehavior::norton_creep_j2_plasticity,
        {1.0e-5, 5.0e6, 3.0},
        {5.0e6, 2.0e9},
    };
    return {base, fuel, cladding};
}

struct WeightedHistory final {
    double average_plastic;
    double average_creep;
    double average_equivalent_stress;
    double maximum_plastic;
    double maximum_creep;
};

WeightedHistory cladding_history(const fuelsim::M2Problem& problem,
                                 const std::vector<double>& state) {
    double measure = 0.0;
    double weighted_plastic = 0.0;
    double weighted_creep = 0.0;
    double weighted_equivalent_stress = 0.0;
    double maximum_plastic = 0.0;
    double maximum_creep = 0.0;
    const fuelsim::ThermoelasticProperties& properties =
        problem.parameters().base.cladding;
    const double shear_modulus =
        properties.young_modulus / (2.0 * (1.0 + properties.poisson_ratio));
    const double lame_lambda = properties.young_modulus *
                               properties.poisson_ratio /
                               ((1.0 + properties.poisson_ratio) *
                                (1.0 - 2.0 * properties.poisson_ratio));
    for (std::size_t element = 0;
         element < problem.base_problem().cladding_element_count(); ++element) {
        const fuelsim::Quad4RzGeometry& geometry =
            problem.base_problem().cladding_element_geometry(element);
        const fuelsim::LocalDofs dofs =
            problem.base_problem().cladding_element_dofs(element);
        fuelsim::LocalValues local_state{};
        for (std::size_t dof = 0; dof < dofs.size(); ++dof)
            local_state[dof] = state[dofs[dof]];
        const fuelsim::Quad4MaterialHistory& history =
            problem.cladding_material_history(element);
        for (std::size_t q = 0; q < history.size(); ++q) {
            const fuelsim::RzQuadraturePoint& point = geometry.points[q];
            double temperature = 0.0;
            double radial_displacement = 0.0;
            double strain_rr = 0.0;
            double strain_zz = 0.0;
            double strain_rz = 0.0;
            for (std::size_t node = 0; node < 4; ++node) {
                temperature += point.shape[node] * local_state[node];
                radial_displacement +=
                    point.shape[node] * local_state[4 + node];
                strain_rr += point.gradient_r[node] * local_state[4 + node];
                strain_zz += point.gradient_z[node] * local_state[8 + node];
                strain_rz +=
                    0.5 * (point.gradient_z[node] * local_state[4 + node] +
                           point.gradient_r[node] * local_state[8 + node]);
            }
            const double thermal_strain =
                properties.thermal_expansion *
                (temperature - properties.reference_temperature);
            const double elastic_rr = strain_rr - history[q].plastic_strain[0] -
                                      history[q].creep_strain[0] -
                                      thermal_strain;
            const double elastic_zz = strain_zz - history[q].plastic_strain[1] -
                                      history[q].creep_strain[1] -
                                      thermal_strain;
            const double elastic_hoop = radial_displacement / point.radius -
                                        history[q].plastic_strain[2] -
                                        history[q].creep_strain[2] -
                                        thermal_strain;
            const double elastic_rz = strain_rz - history[q].plastic_strain[3] -
                                      history[q].creep_strain[3];
            const double trace = elastic_rr + elastic_zz + elastic_hoop;
            const double stress_rr =
                2.0 * shear_modulus * elastic_rr + lame_lambda * trace;
            const double stress_zz =
                2.0 * shear_modulus * elastic_zz + lame_lambda * trace;
            const double stress_hoop =
                2.0 * shear_modulus * elastic_hoop + lame_lambda * trace;
            const double stress_rz = 2.0 * shear_modulus * elastic_rz;
            const double mean_stress =
                (stress_rr + stress_zz + stress_hoop) / 3.0;
            const double deviator_rr = stress_rr - mean_stress;
            const double deviator_zz = stress_zz - mean_stress;
            const double deviator_hoop = stress_hoop - mean_stress;
            const double equivalent_stress = std::sqrt(
                1.5 *
                (deviator_rr * deviator_rr + deviator_zz * deviator_zz +
                 deviator_hoop * deviator_hoop + 2.0 * stress_rz * stress_rz));
            const double weight = point.weighted_measure;
            measure += weight;
            weighted_plastic += weight * history[q].equivalent_plastic_strain;
            weighted_creep += weight * history[q].equivalent_creep_strain;
            weighted_equivalent_stress += weight * equivalent_stress;
            maximum_plastic =
                std::max(maximum_plastic, history[q].equivalent_plastic_strain);
            maximum_creep =
                std::max(maximum_creep, history[q].equivalent_creep_strain);
        }
    }
    if (!(measure > 0.0))
        throw std::runtime_error(
            "PCMI cladding history requires positive integration measure");
    return {
        weighted_plastic / measure,
        weighted_creep / measure,
        weighted_equivalent_stress / measure,
        maximum_plastic,
        maximum_creep,
    };
}

bool fuel_history_is_elastic(const fuelsim::M2Problem& problem) {
    for (std::size_t element = 0;
         element < problem.base_problem().fuel_element_count(); ++element) {
        for (const fuelsim::MaterialPointState& point :
             problem.fuel_material_history(element)) {
            if (point.equivalent_plastic_strain != 0.0 ||
                point.equivalent_creep_strain != 0.0)
                return false;
        }
    }
    return true;
}

bool test_pcmi_coupled_cladding() {
    fuelsim::M2Problem problem(pcmi_parameters());
    const fuelsim::M2TimeOptions time_options = {
        20.0, 1.0, 0.125, 1.0, 1.0, 0.5, 3, 20.0,
    };
    fuelsim::SolverOptions solver_options;
    solver_options.maximum_iterations = 80;

    const fuelsim::M2TimeStepper time_stepper;
    const fuelsim::M2TransientResult result =
        time_stepper.solve(problem, time_options, solver_options);
    if (!result.completed)
        return check(false, "PCMI transient completes all twenty time steps");

    const fuelsim::M1Problem& base = problem.base_problem();
    const fuelsim::DofMap& dofs = base.dof_map();
    const std::size_t axial_mid = base.fuel_mesh().axial_elements() / 2;
    const std::size_t fuel_center = base.fuel_mesh().node_id(0, axial_mid);
    const std::size_t fuel_surface =
        base.fuel_mesh().node_id(base.fuel_mesh().radial_elements(), axial_mid);
    const std::size_t cladding_inner =
        base.cladding_mesh().node_id(0, axial_mid);
    const std::size_t cladding_outer = base.cladding_mesh().node_id(
        base.cladding_mesh().radial_elements(), axial_mid);
    const std::size_t fuel_axis_top =
        base.fuel_mesh().node_id(0, base.fuel_mesh().axial_elements());

    const std::vector<double>& state = result.committed_state;
    const double temperature_center =
        state[dofs.temperature(base.fuel_global_node(fuel_center))];
    const double temperature_fuel_surface =
        state[dofs.temperature(base.fuel_global_node(fuel_surface))];
    const double temperature_cladding_inner =
        state[dofs.temperature(base.cladding_global_node(cladding_inner))];
    const double radial_fuel_surface =
        state[dofs.radial_displacement(base.fuel_global_node(fuel_surface))];
    const double radial_cladding_inner = state[dofs.radial_displacement(
        base.cladding_global_node(cladding_inner))];
    const double radial_cladding_outer = state[dofs.radial_displacement(
        base.cladding_global_node(cladding_outer))];
    const double axial_fuel_top =
        state[dofs.axial_displacement(base.fuel_global_node(fuel_axis_top))];
    const fuelsim::InterfaceSummary interface =
        problem.summarize_interface(state);
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        base.summarize_contact_nodes(state);
    const WeightedHistory history = cladding_history(problem, state);

    // Generated by verification/moose/m23_pcmi_coupled_cladding_rz.i. The
    // fuel is elastic, while the cladding uses coupled Norton creep and J2
    // plasticity in the same material update.
    constexpr double expected_temperature_center = 737.32437038895;
    constexpr double expected_temperature_fuel_surface = 613.74977861649;
    constexpr double expected_temperature_cladding_inner = 612.79216505292;
    constexpr double expected_radial_fuel_surface = 2.9177714961082e-6;
    constexpr double expected_radial_cladding_inner = 1.9098157198910e-6;
    constexpr double expected_radial_cladding_outer = 1.7860488488808e-6;
    constexpr double expected_axial_fuel_top = 9.2609839250304e-6;
    constexpr double expected_average_plastic = 2.4156288403990e-4;
    constexpr double expected_average_creep = 1.2769576281297e-4;
    constexpr double expected_average_equivalent_stress = 5483125.7680798;
    constexpr double expected_total_contact_force = 191.06119157558877;
    constexpr std::array<double, 5> expected_contact_pressure = {
        730945.89357970, 735854.19529622, 777525.10021447,
        640537.88144925, 858282.65485111,
    };

    const double temperature_center_error =
        relative_error(temperature_center, expected_temperature_center);
    const double temperature_fuel_surface_error = relative_error(
        temperature_fuel_surface, expected_temperature_fuel_surface);
    const double temperature_cladding_inner_error = relative_error(
        temperature_cladding_inner, expected_temperature_cladding_inner);
    const double radial_fuel_surface_error =
        relative_error(radial_fuel_surface, expected_radial_fuel_surface);
    const double radial_cladding_inner_error =
        relative_error(radial_cladding_inner, expected_radial_cladding_inner);
    const double radial_cladding_outer_error =
        relative_error(radial_cladding_outer, expected_radial_cladding_outer);
    const double axial_fuel_top_error =
        relative_error(axial_fuel_top, expected_axial_fuel_top);
    const double average_plastic_error =
        relative_error(history.average_plastic, expected_average_plastic);
    const double average_creep_error =
        relative_error(history.average_creep, expected_average_creep);
    const double average_equivalent_stress_error = relative_error(
        history.average_equivalent_stress, expected_average_equivalent_stress);
    const double total_contact_force_error = relative_error(
        interface.total_contact_force, expected_total_contact_force);
    double pressure_difference_squared = 0.0;
    double pressure_reference_squared = 0.0;
    double pressure_maximum_point_error = 0.0;
    if (contact_nodes.size() == expected_contact_pressure.size()) {
        for (std::size_t node = 0; node < expected_contact_pressure.size();
             ++node) {
            const double difference =
                contact_nodes[node].pressure - expected_contact_pressure[node];
            pressure_difference_squared += difference * difference;
            pressure_reference_squared += expected_contact_pressure[node] *
                                          expected_contact_pressure[node];
            pressure_maximum_point_error =
                std::max(pressure_maximum_point_error,
                         relative_error(contact_nodes[node].pressure,
                                        expected_contact_pressure[node]));
        }
    }
    const double contact_pressure_relative_l2 =
        pressure_reference_squared > 0.0
            ? std::sqrt(pressure_difference_squared /
                        pressure_reference_squared)
            : std::numeric_limits<double>::infinity();

    bool passed = check(result.accepted_steps.size() == 20,
                        "PCMI transient commits twenty accepted steps");
    passed = check(result.aggregate_timing.workspace_setups == 1,
                   "PCMI transient reuses one PETSc workspace") &&
             passed;
    passed = check(result.total_cutbacks == 0,
                   "PCMI reference path converges without cutbacks") &&
             passed;
    passed = check(fuel_history_is_elastic(problem),
                   "PCMI fuel remains on the elastic material path") &&
             passed;
    passed =
        check(interface.projected_contact_nodes ==
                  problem.parameters().base.axial_elements + 1,
              "PCMI taller cladding contains every fuel-node projection") &&
        passed;
    passed = check(interface.active_contact_nodes > 0 &&
                       interface.minimum_contact_gap < 0.0 &&
                       interface.maximum_contact_pressure > 0.0,
                   "fuel thermal expansion closes the gap and develops "
                   "contact pressure") &&
             passed;
    passed =
        check(history.average_plastic > 0.0 && history.average_creep > 0.0 &&
                  history.maximum_plastic > 0.0 && history.maximum_creep > 0.0,
              "PCMI cladding accumulates plastic and creep history") &&
        passed;
    passed = check(temperature_center_error < 1.0e-3 &&
                       temperature_fuel_surface_error < 1.0e-3 &&
                       temperature_cladding_inner_error < 1.0e-3,
                   "PCMI temperatures match MOOSE within 0.1%") &&
             passed;
    passed = check(radial_fuel_surface_error < 1.0e-3 &&
                       radial_cladding_inner_error < 1.0e-3 &&
                       radial_cladding_outer_error < 1.0e-3 &&
                       axial_fuel_top_error < 1.0e-3,
                   "PCMI displacements match MOOSE within 0.1%") &&
             passed;
    passed = check(average_plastic_error < 1.0e-3,
                   "PCMI average plastic strain matches MOOSE within 0.1%") &&
             passed;
    passed = check(average_creep_error < 1.0e-3,
                   "PCMI average creep strain matches MOOSE within 0.1%") &&
             passed;
    passed =
        check(average_equivalent_stress_error < 1.0e-3,
              "PCMI average equivalent stress matches MOOSE within 0.1%") &&
        passed;
    passed =
        check(contact_nodes.size() == expected_contact_pressure.size(),
              "PCMI and MOOSE pressure vectors have the same node count") &&
        passed;
    passed = check(contact_pressure_relative_l2 < 1.0e-2,
                   "PCMI contact-pressure L2 error is below 1%") &&
             passed;
    passed = check(pressure_maximum_point_error < 1.0e-2,
                   "PCMI maximum nodal pressure error is below 1%") &&
             passed;
    passed = check(total_contact_force_error < 1.0e-2,
                   "PCMI total contact-force error is below 1%") &&
             passed;

    std::cout << "pcmi_temperature_center=" << temperature_center << '\n';
    std::cout << "pcmi_temperature_fuel_surface=" << temperature_fuel_surface
              << '\n';
    std::cout << "pcmi_temperature_cladding_inner="
              << temperature_cladding_inner << '\n';
    std::cout << "pcmi_radial_fuel_surface=" << radial_fuel_surface << '\n';
    std::cout << "pcmi_radial_cladding_inner=" << radial_cladding_inner << '\n';
    std::cout << "pcmi_radial_cladding_outer=" << radial_cladding_outer << '\n';
    std::cout << "pcmi_axial_fuel_top=" << axial_fuel_top << '\n';
    std::cout << "pcmi_minimum_contact_gap=" << interface.minimum_contact_gap
              << '\n';
    std::cout << "pcmi_maximum_contact_pressure="
              << interface.maximum_contact_pressure << '\n';
    std::cout << "pcmi_total_contact_force=" << interface.total_contact_force
              << '\n';
    std::cout << "pcmi_projected_contact_nodes="
              << interface.projected_contact_nodes << '\n';
    std::cout << "pcmi_active_contact_nodes=" << interface.active_contact_nodes
              << '\n';
    std::cout << "pcmi_cladding_average_equivalent_plastic_strain="
              << history.average_plastic << '\n';
    std::cout << "pcmi_cladding_average_equivalent_creep_strain="
              << history.average_creep << '\n';
    std::cout << "pcmi_cladding_average_equivalent_stress="
              << history.average_equivalent_stress << '\n';
    std::cout << "pcmi_cladding_maximum_equivalent_plastic_strain="
              << history.maximum_plastic << '\n';
    std::cout << "pcmi_cladding_maximum_equivalent_creep_strain="
              << history.maximum_creep << '\n';
    for (std::size_t node = 0; node < contact_nodes.size(); ++node) {
        std::cout << "pcmi_contact_pressure_" << node << '='
                  << contact_nodes[node].pressure << '\n';
    }
    std::cout << "pcmi_moose_temperature_center_relative_error="
              << temperature_center_error << '\n';
    std::cout << "pcmi_moose_temperature_fuel_surface_relative_error="
              << temperature_fuel_surface_error << '\n';
    std::cout << "pcmi_moose_temperature_cladding_inner_relative_error="
              << temperature_cladding_inner_error << '\n';
    std::cout << "pcmi_moose_radial_fuel_surface_relative_error="
              << radial_fuel_surface_error << '\n';
    std::cout << "pcmi_moose_radial_cladding_inner_relative_error="
              << radial_cladding_inner_error << '\n';
    std::cout << "pcmi_moose_radial_cladding_outer_relative_error="
              << radial_cladding_outer_error << '\n';
    std::cout << "pcmi_moose_axial_fuel_top_relative_error="
              << axial_fuel_top_error << '\n';
    std::cout << "pcmi_moose_average_plastic_relative_error="
              << average_plastic_error << '\n';
    std::cout << "pcmi_moose_average_creep_relative_error="
              << average_creep_error << '\n';
    std::cout << "pcmi_moose_average_equivalent_stress_relative_error="
              << average_equivalent_stress_error << '\n';
    std::cout << "pcmi_moose_contact_pressure_relative_l2="
              << contact_pressure_relative_l2 << '\n';
    std::cout << "pcmi_moose_contact_pressure_maximum_point_relative_error="
              << pressure_maximum_point_error << '\n';
    std::cout << "pcmi_moose_total_contact_force_relative_error="
              << total_contact_force_error << '\n';
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv,
            "fuelsim M2 PCMI coupled cladding MOOSE comparison tests\n");
        if (!test_pcmi_coupled_cladding())
            return 1;
        std::cout << "[PASS] fuelsim M2 PCMI coupled cladding tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] PCMI tests raised: " << error.what() << '\n';
        return 1;
    }
}
