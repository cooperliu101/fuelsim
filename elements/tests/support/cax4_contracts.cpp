#include "boundary_types.hpp"
#include "cax4rt.hpp"
#include "contact_types.hpp"
#include "line2_rz_boundary.hpp"
#include "line2_rz_contact.hpp"
#include "support/cax4_evaluation.hpp"
#include "support/element_test_data.hpp"

#include "axisymmetric_types.hpp"
#include "support/material_factory.hpp"
#include "support/rz_linearization.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

bool check(bool condition, const std::string& message) {
    if (!condition)
        std::cerr << "[FAIL] " << message << '\n';
    return condition;
}

double scaled_error(double actual, double expected) {
    return std::abs(actual - expected) / (1.0 + std::max(std::abs(actual), std::abs(expected)));
}

double relative_difference(double actual, double expected) {
    return std::abs(actual - expected) / std::abs(expected);
}

fuelsim::ThermoelasticProperties properties() {
    return fuelsim::test::thermoelastic(3824.0, 0.61, 2.0e11, 0.316, 1.0e-5, 600.0);
}

bool test_cax_kinematics_and_jacobian(bool reduced) {
    const std::string name = reduced ? "CAX4RT" : "CAX4T";
    const std::size_t points = reduced ? 1 : 4;
    const fuelsim::Quad4Coordinates coordinates = {{{1.0, 0.0}, {2.1, 0.1}, {2.0, 1.2}, {0.9, 1.0}}};
    const auto geometry = fuelsim::make_quad4_rz_geometry(coordinates);
    fuelsim::AxisymmetricTestData data{fuelsim::IsotropicThermoelasticMaterial(properties())};
    data.element_formulation = reduced ? fuelsim::RzElementFormulation::cax4rt : fuelsim::RzElementFormulation::cax4t;
    const fuelsim::Cax4LocalValues direction = {0.2, -0.3, 0.4, -0.1, 0.3, -0.5, 0.2, 0.4, -0.2, 0.35, -0.45, 0.25};
    bool passed = true;
    if (!reduced) {
        // Fully restrained heating has the same hydrostatic stress at every
        // point, even on a distorted element with nonuniform corner temperatures.
        for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
            data.strain_formulation = formulation;
            fuelsim::Cax4LocalValues initial{}, old{}, heated{};
            for (std::size_t n = 0; n < 4; ++n)
                initial[n] = 600.0;
            old[0] = 610.0;
            old[1] = 630.0;
            old[2] = 650.0;
            old[3] = 710.0;
            heated[0] = 650.0;
            heated[1] = 690.0;
            heated[2] = 710.0;
            heated[3] = 750.0;
            const auto history = fuelsim::compute_cax4_transient_update(data, geometry, old, initial, {}, 0.1);
            const auto next = fuelsim::compute_cax4_transient_update(data, geometry, heated, old, history, 0.1);
            for (const auto& point : next) {
                const std::array<double, 3> stress = {point.stress.rr, point.stress.zz, point.stress.hoop};
                for (std::size_t component = 0; component < 3; ++component) {
                    passed = check(std::abs(point.elastic_strain[component] + 1e-3) < 1e-14,
                                 "CAX4T expansion uses arithmetic corner temperature across thermal history")
                             && passed;
                    const double expected = -2e11 * 1e-3 / (1.0 - 2.0 * 0.316);
                    passed = check(relative_difference(stress[component], expected) < 1e-12,
                                 "CAX4T restrained nonuniform heating agrees with analytical hydrostatic stress")
                             && passed;
                }
            }
        }
    }
    for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        data.strain_formulation = formulation;
        data.volumetric_heat_source = formulation == fuelsim::StrainFormulation::finite ? 2e6 : 0.0;
        fuelsim::Cax4LocalValues initial{}, old{}, state{};
        for (std::size_t n = 0; n < 4; ++n) {
            initial[n] = old[n] = state[n] = 600.0;
            old[4 + n] = 0.03 * coordinates[n].r;
            old[8 + n] = -0.02 * coordinates[n].z;
            state[4 + n] = 0.08 * coordinates[n].r;
            state[8 + n] = -0.04 * coordinates[n].z;
        }
        const auto old_history = fuelsim::compute_cax4_transient_update(data, geometry, old, initial, {}, 0.1);
        const auto affine = fuelsim::compute_cax4_transient_update(data, geometry, state, old, old_history, 0.1);
        const bool finite = formulation == fuelsim::StrainFormulation::finite;
        const double radial = finite ? 2.0 * 0.03 / 2.03 + 2.0 * 0.05 / 2.11 : 0.08;
        const double axial = finite ? -2.0 * 0.02 / 1.98 - 2.0 * 0.02 / 1.94 : -0.04;
        double affine_error = 0.0;
        for (std::size_t q = 0; q < points; ++q) {
            const auto& h = affine[q];
            affine_error = std::max({affine_error,
                std::abs(h.elastic_strain[0] - radial),
                std::abs(h.elastic_strain[1] - axial),
                std::abs(h.elastic_strain[2] - radial),
                std::abs(h.elastic_strain[3])});
        }
        passed = check(affine_error < 2e-14, name + " homogeneous stretch follows the analytical incremental strain")
                 && passed;
        if (reduced)
            passed = check(fuelsim::elements::cax4rt_hourglass_energy({data.material,
                               geometry,
                               state,
                               state,
                               nullptr,
                               0.0,
                               data.time,
                               data.volumetric_heat_source,
                               data.strain_formulation,
                               false,
                               data.initial_temperature})
                               < 1e-20,
                         "CAX4RT affine deformation has no artificial hourglass energy")
                     && passed;
        state[5] += 0.025;
        state[6] += 0.015;
        state[10] -= 0.020;
        state[11] += 0.010;
        if (reduced) {
            auto heated = state;
            auto unheated = state;
            for (std::size_t n = 0; n < 4; ++n) {
                unheated[n] = 600.0;
                heated[n] = 601.0;
            }
            const double saved_source = data.volumetric_heat_source;
            data.volumetric_heat_source = data.material.heat_capacity(601.0, {}).value() / 0.1;
            const auto residual = fuelsim::compute_cax4_transient(data, geometry, heated, unheated, {}, 0.1);
            for (std::size_t n = 0; n < 4; ++n)
                passed = check(std::abs(residual[n]) < 1e-10 * data.volumetric_heat_source,
                             "CAX4RT uniform heating balances uniform source on a distorted element")
                         && passed;
            data.volumetric_heat_source = saved_source;
        }
        state[0] += 3.0;
        state[2] -= 2.0;
        fuelsim::Cax4LocalJacobian jacobian{};
        const auto active = fuelsim::compute_cax4_transient(data, geometry, state, old, old_history, 0.1, &jacobian);
        const auto passive = fuelsim::compute_cax4_transient(data, geometry, state, old, old_history, 0.1);
        if (!reduced) {
            const auto no_capacity =
                fuelsim::compute_cax4_transient(data, geometry, state, old, old_history, 0.1, nullptr, false);
            auto thermal_coordinates = coordinates;
            if (finite)
                for (std::size_t n = 0; n < 4; ++n) {
                    thermal_coordinates[n].r += state[4 + n];
                    thermal_coordinates[n].z += state[8 + n];
                }
            const auto thermal_geometry = fuelsim::make_quad4_rz_geometry(thermal_coordinates);
            for (std::size_t n = 0; n < 4; ++n) {
                double weight = 0.0;
                for (const auto& p : thermal_geometry.points)
                    weight += p.weighted_measure * p.shape[n];
                const auto& x = coordinates[n];
                const double expected = weight * data.material.heat_capacity(state[n], {data.time, x.r, 0, x.z}).value()
                                        * (state[n] - old[n]) / 0.1;
                passed = check(scaled_error(active[n] - no_capacity[n], expected) < 1e-11,
                             "CAX4T nodal capacity uses row-sum weight and each node's temperature rate")
                         && passed;
            }
        }
        passed =
            check(active == passive, name + " residual and Jacobian evaluations return identical residuals") && passed;
        constexpr double step = 1e-5;
        auto plus = state, minus = state;
        for (std::size_t j = 0; j < 12; ++j) {
            plus[j] += step * direction[j];
            minus[j] -= step * direction[j];
        }
        const auto rp = fuelsim::compute_cax4_transient(data, geometry, plus, old, old_history, 0.1);
        const auto rm = fuelsim::compute_cax4_transient(data, geometry, minus, old, old_history, 0.1);
        double derivative_error = 0.0;
        for (std::size_t i = 0; i < 12; ++i) {
            double ad = 0.0;
            for (std::size_t j = 0; j < 12; ++j)
                ad += jacobian[12 * i + j] * direction[j];
            derivative_error = std::max(derivative_error, scaled_error(ad, (rp[i] - rm[i]) / (2.0 * step)));
        }
        passed = check(derivative_error < 2e-7, name + " coupled volume and hoop Jacobian matches centered differences")
                 && passed;
        std::cout << name + "_" << (finite ? "finite" : "small") << "_directional_jacobian_error=" << derivative_error
                  << '\n';
        // Axial rigid translation is an exact axisymmetric rigid motion, including with committed stress.
        auto translated = state;
        for (std::size_t n = 0; n < 4; ++n)
            translated[8 + n] += 0.5;
        const auto shifted = fuelsim::compute_cax4_transient(data, geometry, translated, old, old_history, 0.1);
        double translation_error = 0.0;
        for (std::size_t i = 0; i < 12; ++i)
            translation_error = std::max(translation_error, scaled_error(shifted[i], passive[i]));
        passed = check(translation_error < 1e-12, name + " internal forces are invariant under axial rigid translation")
                 && passed;
        if (reduced) {
            const double energy = fuelsim::elements::cax4rt_hourglass_energy({data.material,
                geometry,
                state,
                state,
                nullptr,
                0.0,
                data.time,
                data.volumetric_heat_source,
                data.strain_formulation,
                false,
                data.initial_temperature});
            const double shifted_energy = fuelsim::elements::cax4rt_hourglass_energy({data.material,
                geometry,
                translated,
                translated,
                nullptr,
                0.0,
                data.time,
                data.volumetric_heat_source,
                data.strain_formulation,
                false,
                data.initial_temperature});
            passed = check(energy > 0.0 && std::abs(energy - shifted_energy) < 1e-12 * energy,
                         "CAX4RT hourglass energy is positive and invariant under axial translation")
                     && passed;
        }
        if (finite) {
            auto invalid_midpoint = initial;
            for (std::size_t n = 0; n < 4; ++n) {
                invalid_midpoint[4 + n] = 5.0 - 3.0 * coordinates[n].r;
                invalid_midpoint[8 + n] = -1.5 * coordinates[n].z;
            }
            // det(F_new)=1 and all current radii are positive, but det(F_mid)<0.
            bool rejected = false;
            try {
                (void)fuelsim::compute_cax4_transient(data, geometry, invalid_midpoint, initial, {}, 0.1);
            } catch (const std::domain_error& error) {
                rejected = std::string(error.what()).find("midpoint") != std::string::npos;
            }
            passed = check(rejected, name + " rejects an invalid midpoint even when its current volume is positive")
                     && passed;
        }
    }
    for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
        data.strain_formulation = formulation;
        data.volumetric_heat_source = reduced && formulation == fuelsim::StrainFormulation::finite ? 2e6 : 0.0;
        for (int mechanism = 0; mechanism < 3; ++mechanism) {
            auto material = properties();
            if (mechanism != 0)
                material = fuelsim::test::with_norton(material, 1e-4, 1e8, 3.0);
            if (mechanism != 1)
                material = fuelsim::test::with_plasticity(material, 1e8, 2e10);
            data.material = fuelsim::IsotropicThermoelasticMaterial(material);
            fuelsim::Cax4LocalValues initial{}, old{}, state{};
            for (std::size_t n = 0; n < 4; ++n) {
                initial[n] = old[n] = state[n] = 600.0;
                old[4 + n] = 0.001 * coordinates[n].r;
                old[8 + n] = -0.001 * coordinates[n].z;
                state[4 + n] = 0.004 * coordinates[n].r;
                state[8 + n] = -0.005 * coordinates[n].z;
            }
            state[5] += 0.0003;
            state[10] -= 0.0002;
            const auto history = fuelsim::compute_cax4_transient_update(data, geometry, old, initial, {}, 0.1);
            fuelsim::Cax4LocalJacobian jacobian{};
            const auto residual = fuelsim::compute_cax4_transient(data, geometry, state, old, history, 0.1, &jacobian);
            passed = check(residual == fuelsim::compute_cax4_transient(data, geometry, state, old, history, 0.1),
                         name + " inelastic material residual agrees exactly between passive and Jacobian paths")
                     && passed;
            auto plus = state, minus = state;
            constexpr double step = 1e-5;
            for (std::size_t j = 0; j < 12; ++j) {
                plus[j] += step * direction[j];
                minus[j] -= step * direction[j];
            }
            const auto rp = fuelsim::compute_cax4_transient(data, geometry, plus, old, history, 0.1);
            const auto rm = fuelsim::compute_cax4_transient(data, geometry, minus, old, history, 0.1);
            double error = 0.0;
            for (std::size_t i = 0; i < 12; ++i) {
                double ad = 0.0;
                for (std::size_t j = 0; j < 12; ++j)
                    ad += jacobian[12 * i + j] * direction[j];
                error = std::max(error, scaled_error(ad, (rp[i] - rm[i]) / (2.0 * step)));
            }
            const auto next = fuelsim::compute_cax4_transient_update(data, geometry, state, old, history, 0.1);
            for (std::size_t q = 0; q < points; ++q) {
                const auto& point = next[q];
                if (mechanism != 0)
                    passed = check(point.equivalent_creep_strain > 1e-6, name + " creep branch is active") && passed;
                if (mechanism != 1)
                    passed =
                        check(point.equivalent_plastic_strain > 1e-6, name + " plastic branch is active") && passed;
            }
            passed = check(error < 2e-6, name + " active inelastic Jacobian matches centered differences") && passed;
            std::cout << name + "_inelastic_" << static_cast<int>(formulation) << '_' << mechanism
                      << "_jacobian_error=" << error << '\n';
        }
    }
    return passed;
}

using HeatPointGeometries =
    std::array<fuelsim::Line2RzHeatPointGeometry, fuelsim::line2_interface_quadrature_point_count>;

HeatPointGeometries make_heat_point_geometries(const fuelsim::Line2InterfaceSideCoordinates& secondary,
    const fuelsim::Line2InterfaceSideCoordinates& primary,
    double zero_gap_orientation_hint = 0.0) {
    const auto integration =
        fuelsim::make_line2_rz_heat_quadrature(secondary, primary, -1.0, 1.0, zero_gap_orientation_hint);
    HeatPointGeometries result{};
    for (std::size_t q = 0; q < result.size(); ++q) {
        result[q] = fuelsim::make_line2_rz_heat_point_geometry(secondary,
            primary,
            integration[q].secondary_shape,
            integration[q].integration_weight,
            true,
            zero_gap_orientation_hint);
    }
    return result;
}

fuelsim::Cax4LocalResidual summed_heat_residual(const fuelsim::GapHeatProperties& properties,
    const HeatPointGeometries& geometries,
    const fuelsim::Cax4LocalValues& state) {
    fuelsim::Cax4LocalResidual result{};
    for (const fuelsim::Line2RzHeatPointGeometry& geometry : geometries) {
        const fuelsim::Cax4LocalResidual point = fuelsim::compute_line2_rz_gap_heat(properties, geometry, state);
        for (std::size_t row = 0; row < result.size(); ++row)
            result[row] += point[row];
    }
    return result;
}

bool test_heat_interface_case(const std::string& name,
    const fuelsim::GapHeatProperties& properties,
    const HeatPointGeometries& geometries,
    const fuelsim::Cax4LocalValues& state) {
    const fuelsim::Cax4LocalValues direction = {
        0.7,
        -0.4,
        0.3,
        -0.6,
        0.2e-6,
        -0.4e-6,
        0.5e-6,
        -0.1e-6,
        -0.3e-6,
        0.6e-6,
        -0.2e-6,
        0.4e-6,
    };
    fuelsim::rz::LocalLinearization system{};
    for (const fuelsim::Line2RzHeatPointGeometry& geometry : geometries) {
        const fuelsim::rz::LocalLinearization point =
            fuelsim::rz::linearize_line2_rz_gap_heat(properties, geometry, state);
        for (std::size_t row = 0; row < system.residual.size(); ++row)
            system.residual[row] += point.residual[row];
        for (std::size_t entry = 0; entry < system.jacobian.size(); ++entry)
            system.jacobian[entry] += point.jacobian[entry];
    }
    constexpr double step = 1.0e-4;
    fuelsim::Cax4LocalValues plus = state;
    fuelsim::Cax4LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Cax4LocalResidual plus_residual = summed_heat_residual(properties, geometries, plus);
    const fuelsim::Cax4LocalResidual minus_residual = summed_heat_residual(properties, geometries, minus);
    bool passed = true;
    for (const fuelsim::Line2RzHeatPointGeometry& geometry : geometries) {
        passed = check(fuelsim::compute_line2_rz_gap_heat_value(properties, geometry, state).projected,
                     name + " current integration point projects onto its candidate segment")
                 && passed;
    }
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::cax4_local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                     name + " interface AD Jacobian row " + std::to_string(row) + " matches centered finite difference")
                 && passed;
    }
    const fuelsim::Cax4LocalResidual residual = summed_heat_residual(properties, geometries, state);
    double thermal_sum = 0.0;
    double thermal_scale = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        thermal_sum += residual[row];
        thermal_scale += std::abs(residual[row]);
    }
    passed =
        check(std::abs(thermal_sum) < 1.0e-13 * (1.0 + thermal_scale), name + " interface conserves heat") && passed;
    for (std::size_t row = 4; row < fuelsim::cax4_local_dof_count; ++row)
        passed = check(residual[row] == 0.0, name + " gap heat kernel has no mechanical residual") && passed;
    std::cout << name << "_heat_jacobian_maximum_scaled_error=" << maximum_jacobian_error << '\n';
    return passed;
}

bool test_heat_point_interface_case(const std::string& name,
    const fuelsim::GapHeatProperties& properties,
    const fuelsim::Line2RzHeatPointGeometry& geometry,
    const fuelsim::Cax4LocalValues& state) {
    const fuelsim::Cax4LocalValues direction = {
        0.7,
        -0.4,
        0.3,
        -0.6,
        0.2e-6,
        -0.4e-6,
        0.5e-6,
        -0.1e-6,
        -0.3e-6,
        0.6e-6,
        -0.2e-6,
        0.4e-6,
    };
    const fuelsim::rz::LocalLinearization system =
        fuelsim::rz::linearize_line2_rz_gap_heat(properties, geometry, state);
    constexpr double step = 1.0e-4;
    fuelsim::Cax4LocalValues plus = state;
    fuelsim::Cax4LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Cax4LocalResidual plus_residual = fuelsim::compute_line2_rz_gap_heat(properties, geometry, plus);
    const fuelsim::Cax4LocalResidual minus_residual = fuelsim::compute_line2_rz_gap_heat(properties, geometry, minus);
    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::cax4_local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                     name + " point AD Jacobian row " + std::to_string(row) + " matches centered finite difference")
                 && passed;
    }
    const fuelsim::Cax4LocalResidual residual = fuelsim::compute_line2_rz_gap_heat(properties, geometry, state);
    double thermal_sum = 0.0;
    double thermal_scale = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        thermal_sum += residual[row];
        thermal_scale += std::abs(residual[row]);
    }
    passed = check(std::abs(thermal_sum) < 1.0e-13 * (1.0 + thermal_scale), name + " point contribution conserves heat")
             && passed;
    std::cout << name << "_heat_point_jacobian_maximum_scaled_error=" << maximum_jacobian_error << '\n';
    return passed;
}

bool test_contact_interface_case(const std::string& name,
    const fuelsim::NormalContactProperties& properties,
    const fuelsim::NodeToLineRzContactGeometry& geometry,
    const fuelsim::Cax4LocalValues& state,
    const fuelsim::ContactPointHistory& history = {}) {
    const fuelsim::Cax4LocalValues direction = {
        0.7,
        -0.4,
        0.3,
        -0.6,
        0.2e-6,
        -0.4e-6,
        0.5e-6,
        -0.1e-6,
        -0.3e-6,
        0.6e-6,
        -0.2e-6,
        0.4e-6,
    };
    const fuelsim::rz::LocalLinearization system =
        fuelsim::rz::linearize_node_to_line_rz_contact(properties, geometry, state, state, history);
    constexpr double step = 1.0e-4;
    fuelsim::Cax4LocalValues plus = state;
    fuelsim::Cax4LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Cax4LocalResidual plus_residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, plus, state, history);
    const fuelsim::Cax4LocalResidual minus_residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, minus, state, history);
    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::cax4_local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                     name + " contact AD Jacobian row " + std::to_string(row) + " matches centered finite difference")
                 && passed;
    }
    const fuelsim::Cax4LocalResidual residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, state, state, history);
    double radial_sum = 0.0;
    double radial_scale = 0.0;
    for (std::size_t row = 4; row < 8; ++row) {
        radial_sum += residual[row];
        radial_scale += std::abs(residual[row]);
    }
    passed = check(std::abs(radial_sum) < 1.0e-13 * (1.0 + radial_scale), name + " contact conserves radial force")
             && passed;
    double axial_sum = 0.0;
    double axial_scale = 0.0;
    for (std::size_t row = 8; row < 12; ++row) {
        axial_sum += residual[row];
        axial_scale += std::abs(residual[row]);
    }
    passed =
        check(std::abs(axial_sum) < 1.0e-13 * (1.0 + axial_scale), name + " contact conserves axial force") && passed;
    for (std::size_t row = 0; row < 4; ++row)
        passed = check(residual[row] == 0.0, name + " contact has no thermal residual") && passed;
    std::cout << name << "_contact_jacobian_maximum_scaled_error=" << maximum_jacobian_error << '\n';
    return passed;
}

bool test_friction_contact_case(const std::string& name,
    const fuelsim::NormalContactProperties& properties,
    const fuelsim::NodeToLineRzContactGeometry& geometry,
    const fuelsim::Cax4LocalValues& state,
    const fuelsim::Cax4LocalValues& committed_state,
    const fuelsim::ContactPointHistory& history) {
    const fuelsim::Cax4LocalValues direction = {
        0.0,
        0.0,
        0.0,
        0.0,
        0.2e-6,
        -0.4e-6,
        0.5e-6,
        -0.1e-6,
        0.0,
        0.6e-6,
        -0.2e-6,
        0.4e-6,
    };
    const fuelsim::rz::LocalLinearization system =
        fuelsim::rz::linearize_node_to_line_rz_contact(properties, geometry, state, committed_state, history);
    constexpr double step = 1.0e-4;
    fuelsim::Cax4LocalValues plus = state;
    fuelsim::Cax4LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Cax4LocalResidual plus_residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, plus, committed_state, history);
    const fuelsim::Cax4LocalResidual minus_residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, minus, committed_state, history);
    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::cax4_local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                     name + " friction AD Jacobian row " + std::to_string(row) + " matches centered finite difference")
                 && passed;
    }
    const fuelsim::Cax4LocalResidual residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, state, committed_state, history);
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
    passed = check(std::abs(radial_sum) < 1.0e-13 * (1.0 + force_scale)
                       && std::abs(axial_sum) < 1.0e-13 * (1.0 + force_scale),
                 name + " friction reaction is discretely conservative")
             && passed;
    std::cout << name << "_friction_jacobian_maximum_scaled_error=" << maximum_jacobian_error << '\n';
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
    const HeatPointGeometries heat_geometries = make_heat_point_geometries(fuel, cladding);
    const fuelsim::GapHeatProperties heat_properties{0.4, 1.0e-6};
    const fuelsim::NodeToLineRzContactGeometry contact_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding, 1, true, false, 0.0);
    const fuelsim::NormalContactProperties contact_properties{1.0e14};
    const fuelsim::Cax4LocalValues open_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        0.2e-6,
        0.3e-6,
        0.1e-6,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    const fuelsim::Cax4LocalValues minimum_gap_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        1.4e-6,
        1.5e-6,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    const fuelsim::Cax4LocalValues closed_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        3.0e-6,
        3.2e-6,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    bool passed = test_heat_interface_case("open", heat_properties, heat_geometries, open_state);
    passed = test_heat_interface_case("minimum_gap", heat_properties, heat_geometries, minimum_gap_state) && passed;
    passed = test_contact_interface_case("closed", contact_properties, contact_geometry, closed_state) && passed;
    for (const fuelsim::Line2RzHeatPointGeometry& geometry : heat_geometries) {
        const fuelsim::HeatQuadratureValue value =
            fuelsim::compute_line2_rz_gap_heat_value(heat_properties, geometry, open_state);
        const fuelsim::ContactProjectionValue projection =
            fuelsim::compute_line2_rz_heat_projection(geometry, open_state);
        passed = check(value.gap > 1.0e-6 && value.weighted_measure > 0.0 && projection.projected
                           && scaled_error(projection.gap, value.gap) < 1.0e-14,
                     "double-only thermal search projection matches the full interface value")
                 && passed;
    }
    const fuelsim::ContactPointValue open_contact = fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
        contact_geometry,
        open_state,
        open_state,
        {});
    const fuelsim::ContactPointValue closed_contact = fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
        contact_geometry,
        closed_state,
        closed_state,
        {});
    const fuelsim::ContactProjectionValue closed_projection =
        fuelsim::compute_node_to_line_rz_contact_projection(contact_geometry, closed_state);
    passed = check(open_contact.projected && open_contact.gap > 0.0 && open_contact.pressure == 0.0,
                 "open NTS node projects with zero pressure")
             && passed;
    passed = check(closed_contact.projected && closed_contact.gap < 0.0 && closed_contact.pressure > 0.0
                       && closed_contact.contact_force > 0.0 && closed_projection.projected
                       && scaled_error(closed_projection.gap, closed_contact.gap) < 1.0e-14,
                 "double-only mechanical search projection matches the full contact value")
             && passed;
    fuelsim::Cax4LocalValues outside_state = closed_state;
    outside_state[9] = 5.0e-6;
    const fuelsim::ContactPointValue outside_contact =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            contact_geometry,
            outside_state,
            outside_state,
            {});
    passed = check(!outside_contact.projected && outside_contact.contact_force == 0.0,
                 "out-of-segment NTS projection is inactive")
             && passed;
    const fuelsim::NodeToLineRzContactGeometry radial_endpoint_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding, 0, true, true, 0.0);
    fuelsim::Cax4LocalValues radial_endpoint_state = closed_state;
    radial_endpoint_state[8] = -5.0e-16;
    const fuelsim::ContactPointValue radial_endpoint =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            radial_endpoint_geometry,
            radial_endpoint_state,
            radial_endpoint_state,
            {});
    passed = check(radial_endpoint.projected,
                 "radial NTS reference endpoint retains projection after "
                 "roundoff-scale axial motion")
             && passed;
    // A side that is vertical in the reference mesh need not remain vertical.
    // Unequal radial displacement of its primary endpoints must therefore use
    // the current inclined normal instead of a reference-geometry shortcut.
    fuelsim::Cax4LocalValues current_sloped_state = closed_state;
    current_sloped_state[4] = 20.0e-6;
    current_sloped_state[5] = 20.0e-6;
    current_sloped_state[6] = 0.0;
    current_sloped_state[7] = 10.0e-6;
    passed = test_heat_interface_case("current_sloped_reference_vertical",
                 heat_properties,
                 heat_geometries,
                 current_sloped_state)
             && test_contact_interface_case("current_sloped_reference_vertical",
                 contact_properties,
                 contact_geometry,
                 current_sloped_state)
             && passed;
    const fuelsim::Cax4LocalResidual current_sloped_heat_residual =
        summed_heat_residual(heat_properties, heat_geometries, current_sloped_state);
    double expected_current_primary_node_1 = 0.0;
    double fixed_reference_primary_node_1 = 0.0;
    for (const fuelsim::Line2RzHeatPointGeometry& geometry : heat_geometries) {
        const fuelsim::Line2RzHeatQuadraturePoint& point = geometry.point;
        const double secondary_r =
            point.secondary_shape[0] * (geometry.secondary_coordinates[0].r + current_sloped_state[4])
            + point.secondary_shape[1] * (geometry.secondary_coordinates[1].r + current_sloped_state[5]);
        const double secondary_z =
            point.secondary_shape[0] * (geometry.secondary_coordinates[0].z + current_sloped_state[8])
            + point.secondary_shape[1] * (geometry.secondary_coordinates[1].z + current_sloped_state[9]);
        const double primary_r_0 = geometry.primary_coordinates[0].r + current_sloped_state[6];
        const double primary_r_1 = geometry.primary_coordinates[1].r + current_sloped_state[7];
        const double primary_z_0 = geometry.primary_coordinates[0].z + current_sloped_state[10];
        const double primary_z_1 = geometry.primary_coordinates[1].z + current_sloped_state[11];
        const double tangent_r = primary_r_1 - primary_r_0;
        const double tangent_z = primary_z_1 - primary_z_0;
        const double fraction = ((secondary_r - primary_r_0) * tangent_r + (secondary_z - primary_z_0) * tangent_z)
                                / (tangent_r * tangent_r + tangent_z * tangent_z);
        const fuelsim::HeatQuadratureValue value =
            fuelsim::compute_line2_rz_gap_heat_value(heat_properties, geometry, current_sloped_state);
        const double heat_rate = value.heat_flux * value.weighted_measure;
        expected_current_primary_node_1 -= fraction * heat_rate;
        fixed_reference_primary_node_1 -= point.primary_shape[1] * heat_rate;
    }
    passed = check(scaled_error(current_sloped_heat_residual[3], expected_current_primary_node_1) < 1.0e-13
                       && scaled_error(current_sloped_heat_residual[3], fixed_reference_primary_node_1) > 1.0e-6,
                 "thermal contact uses current primary projection rather than "
                 "the reference interpolation fraction")
             && passed;
    fuelsim::Cax4LocalValues thermal_projection_lost_state = open_state;
    thermal_projection_lost_state[10] = 0.45e-3;
    thermal_projection_lost_state[11] = -0.45e-3;
    for (const fuelsim::Line2RzHeatPointGeometry& geometry : heat_geometries) {
        const fuelsim::HeatQuadratureValue value =
            fuelsim::compute_line2_rz_gap_heat_value(heat_properties, geometry, thermal_projection_lost_state);
        const fuelsim::Cax4LocalResidual residual =
            fuelsim::compute_line2_rz_gap_heat(heat_properties, geometry, thermal_projection_lost_state);
        const fuelsim::rz::LocalLinearization system =
            fuelsim::rz::linearize_line2_rz_gap_heat(heat_properties, geometry, thermal_projection_lost_state);
        passed = check(!value.projected
                           && std::all_of(residual.begin(), residual.end(), [](double entry) { return entry == 0.0; })
                           && std::all_of(system.jacobian.begin(),
                               system.jacobian.end(),
                               [](double entry) { return entry == 0.0; }),
                     "thermal quadrature point outside the complete primary segment is not clamped")
                 && passed;
    }
    const fuelsim::ContactPointValue current_sloped_contact =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            contact_geometry,
            current_sloped_state,
            current_sloped_state,
            {});
    const fuelsim::Cax4LocalResidual current_sloped_residual =
        fuelsim::compute_node_to_line_rz_contact(contact_properties,
            contact_geometry,
            current_sloped_state,
            current_sloped_state,
            {});
    passed = check(current_sloped_contact.projected && current_sloped_contact.pressure > 0.0
                       && std::abs(current_sloped_residual[9]) > 0.0,
                 "reference-vertical contact follows the inclined current "
                 "primary normal in both RZ equations")
             && passed;
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
        fuelsim::make_node_to_line_rz_contact_geometry(vertex_secondary, vertex_primary_lower, 1, true, false, 0.0);
    const fuelsim::NodeToLineRzContactGeometry vertex_upper_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(vertex_secondary, vertex_primary_upper, 1, false, false, 0.0);
    fuelsim::Cax4LocalValues vertex_state{};
    vertex_state[5] = 3.0e-6;
    const fuelsim::ContactPointValue vertex_lower_reference =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            vertex_lower_geometry,
            vertex_state,
            vertex_state,
            {});
    const fuelsim::ContactPointValue vertex_upper_reference =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            vertex_upper_geometry,
            vertex_state,
            vertex_state,
            {});
    passed = check(!vertex_lower_reference.projected && vertex_upper_reference.projected,
                 "internal primary vertex has one reference owner")
             && passed;
    vertex_state[9] = -1.0e-8;
    const fuelsim::ContactPointValue vertex_lower_slid =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            vertex_lower_geometry,
            vertex_state,
            vertex_state,
            {});
    const fuelsim::ContactPointValue vertex_upper_slid =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            vertex_upper_geometry,
            vertex_state,
            vertex_state,
            {});
    passed = check(vertex_lower_slid.projected && !vertex_upper_slid.projected && vertex_lower_slid.contact_force > 0.0,
                 "internal primary vertex slide transfers unique NTS "
                 "ownership without double force")
             && passed;
    vertex_state[9] = -2.0e-3;
    const fuelsim::ContactPointValue vertex_upper_far =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            vertex_upper_geometry,
            vertex_state,
            vertex_state,
            {});
    passed = check(!vertex_upper_far.projected, "reference endpoint does not mask a stale far projection") && passed;
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
        fuelsim::make_node_to_line_rz_contact_geometry(general_secondary, general_primary_lower, 1, true, false, 0.0);
    const fuelsim::NodeToLineRzContactGeometry general_upper_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(general_secondary, general_primary_upper, 1, false, true, 0.0);
    fuelsim::Cax4LocalValues general_vertex_state{};
    const fuelsim::ContactPointValue general_lower_reference =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            general_lower_geometry,
            general_vertex_state,
            general_vertex_state,
            {});
    const fuelsim::ContactPointValue general_upper_reference =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            general_upper_geometry,
            general_vertex_state,
            general_vertex_state,
            {});
    passed = check(!general_lower_reference.projected && general_upper_reference.projected,
                 "sloped internal primary vertex has one reference owner")
             && passed;
    general_vertex_state[5] = -1.0e-8;
    general_vertex_state[9] = -1.0e-8;
    const fuelsim::ContactPointValue general_lower_slid =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            general_lower_geometry,
            general_vertex_state,
            general_vertex_state,
            {});
    const fuelsim::ContactPointValue general_upper_slid =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            general_upper_geometry,
            general_vertex_state,
            general_vertex_state,
            {});
    passed = check(general_lower_slid.projected && !general_upper_slid.projected,
                 "sloped internal vertex slide keeps unique NTS ownership")
             && passed;
    const fuelsim::Line2InterfaceSideCoordinates lower_pellet = {{
        {0.0, 0.001000},
        {0.004, 0.001000},
    }};
    const fuelsim::Line2InterfaceSideCoordinates upper_pellet = {{
        {0.0, 0.001002},
        {0.0041, 0.001002},
    }};
    const HeatPointGeometries axial_heat_geometries = make_heat_point_geometries(lower_pellet, upper_pellet);
    const fuelsim::NodeToLineRzContactGeometry axial_contact_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(lower_pellet, upper_pellet, 1, true, true, 0.0);
    const fuelsim::Cax4LocalValues axial_open_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.2e-6,
        0.3e-6,
        0.0,
        0.0,
    };
    const fuelsim::Cax4LocalValues axial_closed_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        0.0,
        0.0,
        0.0,
        0.0,
        3.0e-6,
        3.2e-6,
        0.0,
        0.0,
    };
    passed = test_heat_interface_case("axial_open", heat_properties, axial_heat_geometries, axial_open_state) && passed;
    passed = test_contact_interface_case("axial_closed", contact_properties, axial_contact_geometry, axial_closed_state)
             && passed;
    const fuelsim::ContactPointValue axial_contact = fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
        axial_contact_geometry,
        axial_closed_state,
        axial_closed_state,
        {});
    passed = check(axial_contact.projected && axial_contact.gap < 0.0 && axial_contact.pressure > 0.0,
                 "horizontal pellet faces develop axial contact")
             && passed;
    const fuelsim::Line2InterfaceSideCoordinates sloped_secondary = {{
        {0.004200, 0.000200},
        {0.004800, 0.000800},
    }};
    const fuelsim::Line2InterfaceSideCoordinates sloped_primary = {{
        {0.004002, -0.000002},
        {0.005002, 0.000998},
    }};
    const HeatPointGeometries sloped_heat_geometries = make_heat_point_geometries(sloped_secondary, sloped_primary);
    const fuelsim::NodeToLineRzContactGeometry sloped_contact_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(sloped_secondary, sloped_primary, 1, true, true, 0.0);
    const fuelsim::Cax4LocalValues sloped_closed_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        3.0e-6,
        3.0e-6,
        0.0,
        0.0,
        -3.0e-6,
        -3.0e-6,
        0.0,
        0.0,
    };
    passed = test_heat_interface_case("sloped_open", heat_properties, sloped_heat_geometries, open_state) && passed;
    passed =
        test_contact_interface_case("sloped_closed", contact_properties, sloped_contact_geometry, sloped_closed_state)
        && passed;
    const fuelsim::Cax4LocalResidual sloped_residual = fuelsim::compute_node_to_line_rz_contact(contact_properties,
        sloped_contact_geometry,
        sloped_closed_state,
        sloped_closed_state,
        {});
    const fuelsim::ContactPointValue sloped_contact = fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
        sloped_contact_geometry,
        sloped_closed_state,
        sloped_closed_state,
        {});
    passed = check(sloped_contact.projected && sloped_contact.gap < 0.0 && sloped_contact.pressure > 0.0
                       && std::abs(sloped_residual[5]) > 0.0 && std::abs(sloped_residual[9]) > 0.0,
                 "45-degree contact activates both normal components")
             && passed;
    const fuelsim::NormalContactProperties augmented_properties{1.0e14, 0.0, true};
    fuelsim::ContactPointHistory augmented_history;
    augmented_history.normal_multiplier = 2.0e6;
    const fuelsim::ContactPointValue augmented = fuelsim::compute_node_to_line_rz_contact_value(augmented_properties,
        contact_geometry,
        closed_state,
        closed_state,
        augmented_history);
    const double expected_augmented_pressure = augmented_history.normal_multiplier - 1.0e14 * augmented.gap;
    passed = check(relative_difference(augmented.pressure, expected_augmented_pressure) < 1.0e-13,
                 "augmented contact adds the committed normal multiplier to "
                 "the penalty traction")
             && test_contact_interface_case("augmented_normal",
                 augmented_properties,
                 contact_geometry,
                 closed_state,
                 augmented_history)
             && passed;
    const fuelsim::NormalContactProperties friction_properties{1.0e14, 0.3};
    fuelsim::Cax4LocalValues sticking_state = closed_state;
    sticking_state[9] = 1.0e-7;
    const fuelsim::ContactPointValue sticking = fuelsim::compute_node_to_line_rz_contact_value(friction_properties,
        contact_geometry,
        sticking_state,
        closed_state,
        {});
    passed = check(!sticking.sliding && relative_difference(sticking.tangential_traction, 1.0e7) < 1.0e-13,
                 "Coulomb contact matches the closed-form sticking "
                 "traction")
             && passed;
    passed =
        test_friction_contact_case("sticking", friction_properties, contact_geometry, sticking_state, closed_state, {})
        && passed;
    fuelsim::Cax4LocalValues sliding_state = closed_state;
    sliding_state[9] = 6.0e-7;
    const fuelsim::ContactPointValue sliding = fuelsim::compute_node_to_line_rz_contact_value(friction_properties,
        contact_geometry,
        sliding_state,
        closed_state,
        {});
    const double sliding_limit = 0.3 * sliding.pressure;
    passed = check(sliding.sliding && relative_difference(sliding.tangential_traction, sliding_limit) < 1.0e-13
                       && relative_difference(sliding.elastic_tangential_slip, sliding_limit / 1.0e14) < 1.0e-13,
                 "Coulomb contact caps sliding traction and returns the "
                 "elastic slip")
             && passed;
    passed =
        test_friction_contact_case("sliding", friction_properties, contact_geometry, sliding_state, closed_state, {})
        && passed;
    fuelsim::Cax4LocalValues resticking_state = closed_state;
    resticking_state[9] = -1.0e-7;
    const fuelsim::ContactPointHistory sliding_history = {sliding.elastic_tangential_slip, sliding.sliding, 0.0};
    const fuelsim::ContactPointValue resticking = fuelsim::compute_node_to_line_rz_contact_value(friction_properties,
        contact_geometry,
        resticking_state,
        closed_state,
        sliding_history);
    passed = check(!resticking.sliding && resticking.tangential_traction > 0.0
                       && resticking.tangential_traction < 0.3 * resticking.pressure,
                 "Coulomb contact returns from sliding to sticking under "
                 "reverse tangential motion")
             && passed;
    passed = test_friction_contact_case("sliding_to_sticking",
                 friction_properties,
                 contact_geometry,
                 resticking_state,
                 closed_state,
                 sliding_history)
             && passed;
    const fuelsim::NormalContactProperties sloped_friction_properties{1.0e14, 0.3};
    const fuelsim::NormalContactProperties elastic_slip_properties{1.0e14, 0.3, false, 5.0e-7};
    const auto abaqus_stick = fuelsim::compute_node_to_line_rz_contact_value(elastic_slip_properties,
        contact_geometry,
        sticking_state,
        closed_state,
        {});
    passed = check(!abaqus_stick.sliding
                       && relative_difference(abaqus_stick.tangential_traction,
                              0.3 * abaqus_stick.pressure * abaqus_stick.elastic_tangential_slip / 5e-7)
                              < 1e-13,
                 "RZ elastic-slip stiffness depends on the current normal pressure")
             && passed;
    passed = test_friction_contact_case("pressure_dependent_stick",
                 elastic_slip_properties,
                 contact_geometry,
                 sticking_state,
                 closed_state,
                 {})
             && passed;
    passed = test_friction_contact_case("pressure_dependent_slide",
                 elastic_slip_properties,
                 contact_geometry,
                 sliding_state,
                 closed_state,
                 {})
             && passed;
    fuelsim::ContactPointHistory accumulated_history;
    accumulated_history.total_tangential_slip = 4e-6;
    passed = check(fuelsim::compute_node_to_line_rz_contact(elastic_slip_properties,
                       contact_geometry,
                       sticking_state,
                       closed_state,
                       accumulated_history)
                       == fuelsim::compute_node_to_line_rz_contact(elastic_slip_properties,
                           contact_geometry,
                           sticking_state,
                           closed_state,
                           {}),
                 "RZ accumulated-slip output history does not change the contact residual")
             && passed;
    const auto accumulated = fuelsim::compute_node_to_line_rz_contact_value(elastic_slip_properties,
        contact_geometry,
        sticking_state,
        closed_state,
        accumulated_history);
    passed = check(std::abs(accumulated.total_tangential_slip - 4.1e-6) < 1e-18,
                 "RZ total slip adds active relative motion independently of elastic slip")
             && passed;
    auto open_slip_state = sticking_state;
    open_slip_state[5] -= 1e-3;
    const auto open_slip = fuelsim::compute_node_to_line_rz_contact_value(elastic_slip_properties,
        contact_geometry,
        open_slip_state,
        closed_state,
        accumulated_history);
    passed =
        check(open_slip.pressure == 0.0 && open_slip.total_tangential_slip == accumulated_history.total_tangential_slip,
            "RZ accumulated slip is frozen while contact is open")
        && passed;
    const double sloped_tangent_component = 1.0 / std::sqrt(2.0);
    fuelsim::Cax4LocalValues sloped_sticking_state = sloped_closed_state;
    sloped_sticking_state[5] += 1.0e-7 * sloped_tangent_component;
    sloped_sticking_state[9] += 1.0e-7 * sloped_tangent_component;
    const fuelsim::ContactPointValue sloped_sticking =
        fuelsim::compute_node_to_line_rz_contact_value(sloped_friction_properties,
            sloped_contact_geometry,
            sloped_sticking_state,
            sloped_closed_state,
            {});
    passed =
        check(!sloped_sticking.sliding && relative_difference(sloped_sticking.tangential_traction, 1.0e7) < 1.0e-13,
            "sloped Coulomb contact matches the closed-form sticking "
            "traction")
        && passed;
    passed = test_friction_contact_case("sloped_sticking",
                 sloped_friction_properties,
                 sloped_contact_geometry,
                 sloped_sticking_state,
                 sloped_closed_state,
                 {})
             && passed;
    fuelsim::Cax4LocalValues sloped_sliding_state = sloped_closed_state;
    sloped_sliding_state[5] += 6.0e-7 * sloped_tangent_component;
    sloped_sliding_state[9] += 6.0e-7 * sloped_tangent_component;
    const fuelsim::ContactPointValue sloped_sliding =
        fuelsim::compute_node_to_line_rz_contact_value(sloped_friction_properties,
            sloped_contact_geometry,
            sloped_sliding_state,
            sloped_closed_state,
            {});
    const double sloped_sliding_limit = 0.3 * sloped_sliding.pressure;
    passed = check(sloped_sliding.sliding
                       && relative_difference(sloped_sliding.tangential_traction, sloped_sliding_limit) < 1.0e-13
                       && relative_difference(sloped_sliding.elastic_tangential_slip, sloped_sliding_limit / 1.0e14)
                              < 1.0e-13,
                 "sloped Coulomb contact caps sliding traction and returns "
                 "the elastic slip")
             && passed;
    passed = test_friction_contact_case("sloped_sliding",
                 sloped_friction_properties,
                 sloped_contact_geometry,
                 sloped_sliding_state,
                 sloped_closed_state,
                 {})
             && passed;
    fuelsim::Cax4LocalValues sloped_resticking_state = sloped_closed_state;
    sloped_resticking_state[5] -= 1.0e-7 * sloped_tangent_component;
    sloped_resticking_state[9] -= 1.0e-7 * sloped_tangent_component;
    const fuelsim::ContactPointHistory sloped_sliding_history = {sloped_sliding.elastic_tangential_slip,
        sloped_sliding.sliding,
        0.0};
    const fuelsim::ContactPointValue sloped_resticking =
        fuelsim::compute_node_to_line_rz_contact_value(sloped_friction_properties,
            sloped_contact_geometry,
            sloped_resticking_state,
            sloped_closed_state,
            sloped_sliding_history);
    passed = check(!sloped_resticking.sliding && sloped_resticking.tangential_traction > 0.0
                       && sloped_resticking.tangential_traction < 0.3 * sloped_resticking.pressure,
                 "sloped Coulomb contact returns from sliding to sticking "
                 "under reverse tangential motion")
             && passed;
    passed = test_friction_contact_case("sloped_sliding_to_sticking",
                 sloped_friction_properties,
                 sloped_contact_geometry,
                 sloped_resticking_state,
                 sloped_closed_state,
                 sloped_sliding_history)
             && passed;
    const fuelsim::NormalContactProperties explicit_zero_friction{1.0e14, 0.0};
    const fuelsim::Cax4LocalResidual legacy_normal =
        fuelsim::compute_node_to_line_rz_contact(contact_properties, contact_geometry, closed_state, closed_state, {});
    const fuelsim::Cax4LocalResidual zero_friction = fuelsim::compute_node_to_line_rz_contact(explicit_zero_friction,
        contact_geometry,
        closed_state,
        closed_state,
        {});
    passed = check(legacy_normal == zero_friction,
                 "mu equal to zero preserves every normal-contact residual "
                 "entry exactly")
             && passed;
    return passed;
}

bool test_heat_point_primary_owner() {
    const fuelsim::Line2InterfaceSideCoordinates sloped_secondary = {{{1.0, .2}, {1.2, 1.2}}};
    const fuelsim::Line2InterfaceSideCoordinates sloped_primary = {{{1.4, 0.0}, {1.8, 2.0}}};
    const fuelsim::Cax4LocalValues nodal_state = {500, 650, 300, 350, .01, .02, .03, -.01, .01, .02, -.02, .01};
    bool nodal_passed = true;
    for (std::size_t node = 0; node < 2; ++node) {
        const std::array<double, 2> shape = {node == 0 ? 1.0 : 0.0, node == 1 ? 1.0 : 0.0};
        auto geometry =
            fuelsim::make_line2_rz_heat_point_geometry(sloped_secondary, sloped_primary, shape, 1.0, true, 0.0);
        geometry.point.nodal = true;
        nodal_passed = test_heat_point_interface_case("nts_sloped_" + std::to_string(node),
                           fuelsim::GapHeatProperties{.2, 1e-5},
                           geometry,
                           nodal_state)
                       && nodal_passed;
        const auto value = fuelsim::compute_line2_rz_gap_heat_value({.2, 1e-5}, geometry, nodal_state);
        const double r0 = 1.01, r1 = 1.22, length = std::hypot(.21, 1.01);
        const double expected = pi * length * (node == 0 ? 2 * r0 + r1 : r0 + 2 * r1) / 3;
        nodal_passed = check(relative_difference(value.weighted_measure, expected) < 1e-13,
                           "NTS heat area equals analytical linear-shape integral on current sloped edge")
                       && nodal_passed;
    }
    const fuelsim::Line2InterfaceSideCoordinates secondary = {{{1.0, 0.5}, {1.0, 1.5}}};
    const fuelsim::Line2InterfaceSideCoordinates primary_lower = {{{1.1, 0.0}, {1.1, 1.0}}};
    const fuelsim::Line2InterfaceSideCoordinates primary_upper = {{{1.1, 1.0}, {1.1, 2.0}}};
    const std::array<double, 2> secondary_shape = {0.5, 0.5};
    const fuelsim::Line2RzHeatPointGeometry lower =
        fuelsim::make_line2_rz_heat_point_geometry(secondary, primary_lower, secondary_shape, 1.0, false, 0.0);
    const fuelsim::Line2RzHeatPointGeometry upper =
        fuelsim::make_line2_rz_heat_point_geometry(secondary, primary_upper, secondary_shape, 1.0, true, 0.0);
    const fuelsim::GapHeatProperties properties{0.2, 1.0e-5};
    fuelsim::Cax4LocalValues state = {
        500.0,
        500.0,
        300.0,
        300.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    bool passed = nodal_passed;
    const fuelsim::HeatQuadratureValue lower_vertex =
        fuelsim::compute_line2_rz_gap_heat_value(properties, lower, state);
    const fuelsim::HeatQuadratureValue upper_vertex =
        fuelsim::compute_line2_rz_gap_heat_value(properties, upper, state);
    passed = check(!lower_vertex.projected && upper_vertex.projected,
                 "thermal point on an internal primary vertex belongs only to the following half-open segment")
             && passed;
    state[8] = 0.1;
    state[9] = 0.1;
    const fuelsim::HeatQuadratureValue lower_after = fuelsim::compute_line2_rz_gap_heat_value(properties, lower, state);
    const fuelsim::HeatQuadratureValue upper_after = fuelsim::compute_line2_rz_gap_heat_value(properties, upper, state);
    passed = check(!lower_after.projected && upper_after.projected,
                 "thermal point transfers uniquely to the following primary segment after forward sliding")
             && passed;
    passed = test_heat_point_interface_case("current_primary_owner", properties, upper, state) && passed;
    state[8] = -0.1;
    state[9] = -0.1;
    const fuelsim::HeatQuadratureValue lower_before =
        fuelsim::compute_line2_rz_gap_heat_value(properties, lower, state);
    const fuelsim::HeatQuadratureValue upper_before =
        fuelsim::compute_line2_rz_gap_heat_value(properties, upper, state);
    passed = check(lower_before.projected && !upper_before.projected,
                 "thermal point transfers uniquely to the preceding primary segment after reverse sliding")
             && passed;
    state[8] = 1.1;
    state[9] = 1.1;
    passed = check(!fuelsim::compute_line2_rz_gap_heat_value(properties, lower, state).projected
                       && !fuelsim::compute_line2_rz_gap_heat_value(properties, upper, state).projected,
                 "thermal point outside the complete primary chain is not clamped to a distant endpoint")
             && passed;
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
        fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding_coincident, 1, true, false, -1.0e-4);
    const fuelsim::NodeToLineRzContactGeometry opened_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding_open, 1, true, false, 0.0);
    passed = check(zero_gap_geometry.normal_orientation == 1.0
                       && zero_gap_geometry.normal_orientation == opened_geometry.normal_orientation,
                 "zero-gap hint orientation matches the same geometry with "
                 "the gap opened by 1e-9 m")
             && passed;
    const fuelsim::NodeToLineRzContactGeometry flipped_hint_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding_coincident, 1, true, false, 1.0e-4);
    const fuelsim::NodeToLineRzContactGeometry inward_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding_inward, 1, true, false, 0.0);
    passed = check(flipped_hint_geometry.normal_orientation == -1.0
                       && flipped_hint_geometry.normal_orientation == inward_geometry.normal_orientation,
                 "flipped material side flips the orientation, matching a "
                 "1e-9 m inward offset of the primary surface")
             && passed;
    bool missing_hint_rejected = false;
    try {
        (void)fuelsim::make_node_to_line_rz_contact_geometry(fuel, cladding_coincident, 1, true, false, 0.0);
    } catch (const std::invalid_argument&) {
        missing_hint_rejected = true;
    }
    passed = check(missing_hint_rejected,
                 "on-segment zero gap without a hint remains an explicit "
                 "error")
             && passed;
    const HeatPointGeometries zero_gap_heat = make_heat_point_geometries(fuel, cladding_coincident, -1.0e-4);
    const HeatPointGeometries opened_heat = make_heat_point_geometries(fuel, cladding_open);
    for (std::size_t point = 0; point < zero_gap_heat.size(); ++point) {
        passed =
            check(zero_gap_heat[point].point.normal_orientation == 1.0
                      && zero_gap_heat[point].point.normal_orientation == opened_heat[point].point.normal_orientation,
                "zero-gap heat quadrature orientation matches the "
                "1e-9 m opened geometry")
            && passed;
    }
    bool heat_missing_hint_rejected = false;
    try {
        (void)fuelsim::make_line2_rz_heat_quadrature(fuel, cladding_coincident, -1.0, 1.0, 0.0);
    } catch (const std::invalid_argument&) {
        heat_missing_hint_rejected = true;
    }
    passed = check(heat_missing_hint_rejected,
                 "zero-gap heat geometry without a hint remains an "
                 "explicit error")
             && passed;
    // Branch checks avoid the nondifferentiable kink at gap exactly zero:
    // the open state leaves a +1e-6 m gap and the closed state penetrates
    // 3e-6 m, while the 1e-4 scaled finite-difference direction moves the
    // gap by at most 6e-11 m and cannot cross the kink.
    const fuelsim::NormalContactProperties contact_properties{1.0e14};
    const fuelsim::Cax4LocalValues zero_gap_open_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        0.2e-6,
        -1.0e-6,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    const fuelsim::Cax4LocalValues zero_gap_closed_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        3.0e-6,
        3.2e-6,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    passed = test_contact_interface_case("zero_gap_open", contact_properties, zero_gap_geometry, zero_gap_open_state)
             && passed;
    passed =
        test_contact_interface_case("zero_gap_closed", contact_properties, zero_gap_geometry, zero_gap_closed_state)
        && passed;
    const fuelsim::ContactPointValue zero_gap_open = fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
        zero_gap_geometry,
        zero_gap_open_state,
        zero_gap_open_state,
        {});
    const fuelsim::ContactPointValue zero_gap_closed =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            zero_gap_geometry,
            zero_gap_closed_state,
            zero_gap_closed_state,
            {});
    passed = check(zero_gap_open.projected && zero_gap_open.gap > 0.0 && zero_gap_open.pressure == 0.0,
                 "zero-gap geometry opens with zero pressure at +1e-6 m gap")
             && check(zero_gap_closed.projected && zero_gap_closed.gap < 0.0
                          && scaled_error(zero_gap_closed.pressure, 3.2e8) < 1.0e-13,
                 "zero-gap geometry develops penalty pressure over a 3.2e-6 m "
                 "penetration")
             && passed;
    // Gap continuity oracle for the general current-normal formula. Opening
    // the reference surfaces by 1e-9 m must increase the measured gap by
    // exactly 1e-9 m at the same state, which pins both the orientation and
    // the gap sign convention.
    const fuelsim::Line2InterfaceSideCoordinates sloped_fuel = {{
        {0.004120, 0.0},
        {0.0041205, 0.001},
    }};
    const fuelsim::NodeToLineRzContactGeometry sloped_zero_gap_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(sloped_fuel, cladding_coincident, 0, true, false, -1.0e-4);
    const fuelsim::NodeToLineRzContactGeometry sloped_opened_geometry =
        fuelsim::make_node_to_line_rz_contact_geometry(sloped_fuel, cladding_open, 0, true, false, 0.0);
    passed = check(sloped_zero_gap_geometry.normal_orientation == 1.0
                       && sloped_zero_gap_geometry.normal_orientation == sloped_opened_geometry.normal_orientation,
                 "sloped zero-gap hint orientation matches the 1e-9 m "
                 "opened geometry")
             && passed;
    const fuelsim::Cax4LocalValues oracle_closed_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        3.0e-6,
        3.0e-6,
        0.0,
        0.0,
        1.0e-6,
        0.0,
        0.0,
        0.0,
    };
    passed = test_contact_interface_case("zero_gap_sloped_closed",
                 contact_properties,
                 sloped_zero_gap_geometry,
                 oracle_closed_state)
             && passed;
    const fuelsim::ContactPointValue sloped_zero_value =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            sloped_zero_gap_geometry,
            oracle_closed_state,
            oracle_closed_state,
            {});
    const fuelsim::ContactPointValue sloped_opened_value =
        fuelsim::compute_node_to_line_rz_contact_value(contact_properties,
            sloped_opened_geometry,
            oracle_closed_state,
            oracle_closed_state,
            {});
    const double gap_oracle_error = std::abs((sloped_opened_value.gap - sloped_zero_value.gap) - epsilon_gap);
    std::cout << "zero_gap_oracle_coincident_gap=" << sloped_zero_value.gap << '\n'
              << "zero_gap_oracle_opened_gap=" << sloped_opened_value.gap << '\n'
              << "zero_gap_oracle_gap_continuity_error=" << gap_oracle_error << '\n';
    passed = check(sloped_zero_value.projected && sloped_opened_value.projected && gap_oracle_error < 1.0e-15,
                 "opening the reference surfaces by 1e-9 m increases the "
                 "measured gap by exactly 1e-9 m at the same state")
             && check(scaled_error(sloped_zero_value.pressure, -1.0e14 * sloped_zero_value.gap) < 1.0e-13
                          && scaled_error(sloped_opened_value.pressure, -1.0e14 * sloped_opened_value.gap) < 1.0e-13
                          && sloped_zero_value.pressure > 2.9e8 && sloped_opened_value.pressure > 2.9e8,
                 "both geometries follow the penalty law on their own closed "
                 "branch gaps near 3e8 Pa")
             && passed;
    // The gap heat kernel has no kink at zero gap; run its two smooth
    // branches on the coincident geometry: +2e-6 m gap above the 1e-6 m
    // minimum gap and a -3e-6 m penetration on the minimum-gap branch.
    const fuelsim::GapHeatProperties heat_properties{0.4, 1.0e-6};
    const fuelsim::Cax4LocalValues zero_gap_heat_open_state = {
        750.0,
        740.0,
        610.0,
        620.0,
        -2.0e-6,
        -2.0e-6,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    passed = test_heat_interface_case("zero_gap_heat_open", heat_properties, zero_gap_heat, zero_gap_heat_open_state)
             && passed;
    passed = test_heat_interface_case("zero_gap_heat_minimum", heat_properties, zero_gap_heat, zero_gap_closed_state)
             && passed;
    for (const fuelsim::Line2RzHeatPointGeometry& geometry : zero_gap_heat) {
        const fuelsim::HeatQuadratureValue value =
            fuelsim::compute_line2_rz_gap_heat_value(heat_properties, geometry, zero_gap_heat_open_state);
        passed = check(value.gap > 1.0e-6 && value.weighted_measure > 0.0,
                     "zero-gap heat point reports the +2e-6 m open gap")
                 && passed;
    }
    return passed;
}

bool test_follower_pressure() {
    const fuelsim::Line2RzBoundaryGeometry geometry =
        fuelsim::make_line2_rz_boundary_geometry({{{0.005, 0.0}, {0.005, 0.01}}}, {{1, 2}});
    constexpr double pressure = 3.0e6;
    const fuelsim::Line2RzBoundaryData follower = {fuelsim::Line2RzBoundaryKind::pressure,
        fuelsim::TractionComponent::radial,
        pressure,
        0.0,
        true};
    const fuelsim::Cax4LocalValues state = {
        600.0,
        600.0,
        600.0,
        600.0,
        0.0,
        0.001,
        0.002,
        0.0,
        0.0,
        0.0004,
        -0.0002,
        0.0,
    };
    const fuelsim::Cax4LocalValues direction = {
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.3,
        -0.2,
        0.0,
        0.0,
        -0.4,
        0.5,
        0.0,
    };
    const fuelsim::rz::LocalLinearization system = fuelsim::rz::linearize_line2_rz_boundary(follower, geometry, state);
    const double first_radius = 0.006;
    const double second_radius = 0.007;
    const double delta_radius = second_radius - first_radius;
    const double delta_axial = 0.0098 - 0.0004;
    const double expected_radial = pi * pressure * delta_axial * (first_radius + second_radius);
    const double expected_axial = -pi * pressure * delta_radius * (first_radius + second_radius);
    bool passed = check(scaled_error(system.residual[5] + system.residual[6], expected_radial) < 1.0e-13
                            && scaled_error(system.residual[9] + system.residual[10], expected_axial) < 1.0e-13,
        "follower pressure uses current RZ radius and outward normal");
    constexpr double step = 1.0e-7;
    fuelsim::Cax4LocalValues plus = state;
    fuelsim::Cax4LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Cax4LocalResidual plus_residual = fuelsim::compute_line2_rz_boundary(follower, geometry, plus);
    const fuelsim::Cax4LocalResidual minus_residual = fuelsim::compute_line2_rz_boundary(follower, geometry, minus);
    double maximum_error = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row) {
        double tangent = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column)
            tangent += system.jacobian[row * state.size() + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        maximum_error = std::max(maximum_error, scaled_error(tangent, finite_difference));
    }
    std::cout << "follower_pressure_directional_jacobian_error=" << maximum_error << '\n';
    passed = check(maximum_error < 1.0e-8,
                 "follower-pressure AD Jacobian matches centered "
                 "differences")
             && passed;
    const fuelsim::Line2RzBoundaryData dead = {fuelsim::Line2RzBoundaryKind::pressure,
        fuelsim::TractionComponent::radial,
        pressure,
        0.0,
        false};
    const fuelsim::rz::LocalLinearization dead_system = fuelsim::rz::linearize_line2_rz_boundary(dead, geometry, state);
    const double maximum_dead_tangent = *std::max_element(dead_system.jacobian.begin(),
        dead_system.jacobian.end(),
        [](double left, double right) { return std::abs(left) < std::abs(right); });
    passed = check(maximum_dead_tangent == 0.0,
                 "reference pressure has an exactly zero geometric "
                 "tangent")
             && passed;
    return passed;
}

bool test_current_configuration_traction() {
    const fuelsim::Line2RzBoundaryGeometry geometry =
        fuelsim::make_line2_rz_boundary_geometry({{{0.005, 0.0}, {0.005, 0.01}}}, {{1, 2}});
    constexpr double traction = 2.0e6;
    const fuelsim::Line2RzBoundaryData current = {fuelsim::Line2RzBoundaryKind::traction,
        fuelsim::TractionComponent::axial,
        traction,
        0.0,
        true};
    const fuelsim::Cax4LocalValues state = {
        600.0,
        600.0,
        600.0,
        600.0,
        0.0,
        0.001,
        0.002,
        0.0,
        0.0,
        0.0004,
        -0.0002,
        0.0,
    };
    const fuelsim::Cax4LocalValues direction = {
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.3,
        -0.2,
        0.0,
        0.0,
        -0.4,
        0.5,
        0.0,
    };
    const fuelsim::rz::LocalLinearization system = fuelsim::rz::linearize_line2_rz_boundary(current, geometry, state);
    const double current_length = std::hypot(0.001, 0.0094);
    const double expected_axial = -2.0 * pi * 0.0065 * current_length * traction;
    bool passed = check(scaled_error(system.residual[9] + system.residual[10], expected_axial) < 1.0e-13,
        "current-configuration component traction uses current RZ "
        "surface measure");
    constexpr double step = 1.0e-7;
    fuelsim::Cax4LocalValues plus = state;
    fuelsim::Cax4LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Cax4LocalResidual plus_residual = fuelsim::compute_line2_rz_boundary(current, geometry, plus);
    const fuelsim::Cax4LocalResidual minus_residual = fuelsim::compute_line2_rz_boundary(current, geometry, minus);
    double maximum_error = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row) {
        double tangent = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column)
            tangent += system.jacobian[row * state.size() + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        maximum_error = std::max(maximum_error, scaled_error(tangent, finite_difference));
    }
    std::cout << "current_traction_directional_jacobian_error=" << maximum_error << '\n';
    passed = check(maximum_error < 1.0e-8,
                 "current-configuration traction AD Jacobian matches "
                 "centered differences")
             && passed;
    const fuelsim::Line2RzBoundaryData reference = {fuelsim::Line2RzBoundaryKind::traction,
        fuelsim::TractionComponent::axial,
        traction,
        0.0,
        false};
    const fuelsim::rz::LocalLinearization reference_system =
        fuelsim::rz::linearize_line2_rz_boundary(reference, geometry, state);
    const double maximum_reference_tangent = *std::max_element(reference_system.jacobian.begin(),
        reference_system.jacobian.end(),
        [](double left, double right) { return std::abs(left) < std::abs(right); });
    return check(maximum_reference_tangent == 0.0,
               "reference-configuration traction has zero geometric "
               "tangent")
           && passed;
}

bool test_temperature_active_thermoelastic_properties() {
    const fuelsim::ThermoelasticProperties active_properties =
        fuelsim::test::thermoelastic(3824.0, 0.61, 2.0e11, 0.316, 1.0e-5, 600.0, -8.0e7, 2.0e-5, 3.0e-9);
    const fuelsim::IsotropicThermoelasticMaterial material(active_properties);
    constexpr double temperature = 725.0;
    const adlite::Scalar active_temperature = adlite::Scalar::independent(temperature, 0, 1);
    const fuelsim::AxisymmetricStress active = material.stress(1.1e-3, -0.4e-3, 0.2e-3, 0.3e-3, active_temperature);
    constexpr double step = 1.0e-3;
    const fuelsim::AxisymmetricStress plus = material.stress(1.1e-3, -0.4e-3, 0.2e-3, 0.3e-3, temperature + step);
    const fuelsim::AxisymmetricStress minus = material.stress(1.1e-3, -0.4e-3, 0.2e-3, 0.3e-3, temperature - step);
    const std::array<double, 4> analytic = {active.rr.derivative(0),
        active.zz.derivative(0),
        active.hoop.derivative(0),
        active.rz.derivative(0)};
    const std::array<double, 4> finite_difference = {(plus.rr.value() - minus.rr.value()) / (2.0 * step),
        (plus.zz.value() - minus.zz.value()) / (2.0 * step),
        (plus.hoop.value() - minus.hoop.value()) / (2.0 * step),
        (plus.rz.value() - minus.rz.value()) / (2.0 * step)};
    double maximum_error = 0.0;
    for (std::size_t component = 0; component < analytic.size(); ++component)
        maximum_error = std::max(maximum_error, scaled_error(analytic[component], finite_difference[component]));
    std::cout << "active_thermoelastic_temperature_tangent_error=" << maximum_error << '\n';
    return check(maximum_error < 1.0e-8,
        "temperature-dependent thermoelastic AD tangent matches "
        "centered differences");
}

bool test_line2_convection() {
    bool passed = true;
    const fuelsim::Line2RzBoundaryGeometry geometry =
        fuelsim::make_line2_rz_boundary_geometry({{{0.005, 0.0}, {0.005, 0.01}}}, {{1, 2}});
    const fuelsim::Line2RzBoundaryData data = {fuelsim::Line2RzBoundaryKind::convection,
        fuelsim::TractionComponent::radial,
        1000.0,
        500.0,
        false};
    const fuelsim::Cax4LocalValues state = {
        590.0,
        600.0,
        600.0,
        610.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    const fuelsim::Cax4LocalValues direction = {
        0.2,
        -0.7,
        0.4,
        0.3,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    };
    const fuelsim::rz::LocalLinearization system = fuelsim::rz::linearize_line2_rz_boundary(data, geometry, state);
    const double expected_heat = 1000.0 * 100.0 * 2.0 * pi * 0.005 * 0.01;
    passed = check(scaled_error(system.residual[1] + system.residual[2], expected_heat) < 1.0e-13,
                 "convection integrates the RZ surface heat loss")
             && passed;
    constexpr double step = 1.0e-4;
    fuelsim::Cax4LocalValues plus = state;
    fuelsim::Cax4LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::Cax4LocalResidual plus_residual = fuelsim::compute_line2_rz_boundary(data, geometry, plus);
    const fuelsim::Cax4LocalResidual minus_residual = fuelsim::compute_line2_rz_boundary(data, geometry, minus);
    double maximum_error = 0.0;
    for (std::size_t row = 0; row < state.size(); ++row) {
        double tangent = 0.0;
        for (std::size_t column = 0; column < state.size(); ++column)
            tangent += system.jacobian[row * state.size() + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        maximum_error = std::max(maximum_error, scaled_error(tangent, finite_difference));
    }
    std::cout << "convection_directional_jacobian_error=" << maximum_error << '\n';
    passed = check(maximum_error < 1.0e-10, "convection AD Jacobian matches centered differences") && passed;
    return passed;
}
} // namespace

int run_cax4t_contract_tests() {
    std::cout << std::scientific << std::setprecision(12);
    return test_cax_kinematics_and_jacobian(false) && test_temperature_active_thermoelastic_properties() ? 0 : 1;
}

int run_cax4rt_tests() {
    std::cout << std::scientific << std::setprecision(12);
    return test_cax_kinematics_and_jacobian(true) ? 0 : 1;
}

int run_axisymmetric_interface_tests() {
    std::cout << std::scientific << std::setprecision(12);
    bool passed = true;
    passed = test_gap_heat_and_normal_contact() && passed;
    passed = test_heat_point_primary_owner() && passed;
    passed = test_zero_gap_contact_orientation() && passed;
    passed = test_follower_pressure() && passed;
    passed = test_current_configuration_traction() && passed;
    passed = test_line2_convection() && passed;
    if (passed)
        std::cout << "[PASS] independent RZ element and interface tests\n";
    return passed ? 0 : 1;
}
