#include "core/detail/cax4rt.hpp"
#include "core/spatial_layout.hpp"
#include "fuelsim/core/contact.hpp"
#include "fuelsim/core/mesh.hpp"
#include "fuelsim/core/rz_quad4.hpp"
#include "fuelsim/core/spatial_definition.hpp"
#include "fuelsim/core/steady_problem.hpp"
#include "fuelsim/elements/material.hpp"
#include "support/material_factory.hpp"
#include "support/mesh_fixture.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double scaled_error(double actual, double expected) {
    return std::abs(actual - expected) / (1.0 + std::max(std::abs(actual), std::abs(expected)));
}

double relative_difference(double actual, double expected) {
    return std::abs(actual - expected) / std::abs(expected);
}

bool test_contact_search_tree() {
    std::vector<fuelsim::spatial_detail::ContactSearchBox> boxes;
    for (std::size_t item = 0; item < 97; ++item) {
        const double coordinate = 0.001 * static_cast<double>((37 * item) % 97);
        boxes.push_back({{coordinate, -0.01, 0.0}, {coordinate, 0.01, 0.0}, item});
    }
    fuelsim::spatial_detail::ContactSearchTree tree;
    tree.build(std::move(boxes));
    fuelsim::spatial_detail::ContactSearchQuery query;
    bool passed = true;
    std::size_t maximum_visited_candidates = 0;
    for (std::size_t sample = 0; sample < 151; ++sample) {
        const std::array<double, 3> point = {0.00064 * static_cast<double>(sample),
            0.007 * std::sin(static_cast<double>(sample)),
            0.0};
        double exhaustive_distance = std::numeric_limits<double>::infinity();
        std::size_t exhaustive_item = std::numeric_limits<std::size_t>::max();
        for (std::size_t item = 0; item < 97; ++item) {
            const double coordinate = 0.001 * static_cast<double>((37 * item) % 97);
            const double distance = std::abs(point[0] - coordinate);
            if (distance < exhaustive_distance || (distance == exhaustive_distance && item < exhaustive_item)) {
                exhaustive_distance = distance;
                exhaustive_item = item;
            }
        }
        tree.begin_query(point, query);
        double tree_distance = std::numeric_limits<double>::infinity();
        std::size_t tree_item = std::numeric_limits<std::size_t>::max(), candidate = 0;
        std::size_t visited_candidates = 0;
        while (tree.next_candidate(query, tree_distance, candidate)) {
            ++visited_candidates;
            const double coordinate = 0.001 * static_cast<double>((37 * candidate) % 97);
            const double distance = std::abs(point[0] - coordinate);
            if (distance < tree_distance || (distance == tree_distance && candidate < tree_item)) {
                tree_distance = distance;
                tree_item = candidate;
            }
        }
        passed = check(tree_item == exhaustive_item && tree_distance == exhaustive_distance,
                     "contact search tree matches exhaustive nearest-candidate selection")
                 && passed;
        maximum_visited_candidates = std::max(maximum_visited_candidates, visited_candidates);
    }
    std::cout << "contact_search_maximum_visited_candidates=" << maximum_visited_candidates << '\n';
    passed = check(maximum_visited_candidates <= 2,
                 "contact search tree prunes the 97-segment exact search to at most two candidates")
             && passed;
    tree.build({{{0.03, -0.01, 0.0}, {0.03, 0.01, 0.0}, 8}, {{0.03, -0.01, 0.0}, {0.03, 0.01, 0.0}, 2}});
    tree.begin_query({0.04, 0.0, 0.0}, query);
    double tied_distance = std::numeric_limits<double>::infinity();
    std::size_t tied_item = std::numeric_limits<std::size_t>::max(), candidate = 0;
    while (tree.next_candidate(query, tied_distance, candidate)) {
        constexpr double distance = 0.01;
        if (distance < tied_distance || (distance == tied_distance && candidate < tied_item)) {
            tied_distance = distance;
            tied_item = candidate;
        }
    }
    passed = check(tied_item == 2, "contact search tree retains equal-distance candidates for deterministic ownership")
             && passed;
    return passed;
}

fuelsim::ThermoelasticProperties properties() {
    return fuelsim::test::thermoelastic(3824.0, 0.61, 2.0e11, 0.316, 1.0e-5, 600.0);
}

bool test_mesh_and_geometry() {
    const double inner = 0.0;
    const double outer = 0.004;
    const double length = 0.01;
    const fuelsim::UnstructuredQuad4Mesh source =
        fuelsim::test::make_disconnected_annular_mesh({{1, "solid", inner, outer, length, 4, 3}});
    const fuelsim::RegionMesh mesh =
        fuelsim::RegionMesh::from_unstructured_block(source, source.element_block("solid").id);
    bool passed = true;
    passed = check(mesh.nodes().size() == 20, "structured mesh node count") && passed;
    passed = check(mesh.elements().size() == 12, "structured mesh element count") && passed;
    double integrated_volume = 0.0;
    for (const fuelsim::Quad4Element& element : mesh.elements()) {
        fuelsim::Quad4Coordinates coordinates{};
        for (std::size_t node = 0; node < element.nodes.size(); ++node)
            coordinates[node] = mesh.nodes()[element.nodes[node]];
        const fuelsim::Quad4RzGeometry geometry = fuelsim::make_quad4_rz_geometry(coordinates);
        for (const fuelsim::RzQuadraturePoint& point : geometry.points) {
            double shape_sum = 0.0;
            double gradient_r_sum = 0.0;
            double gradient_z_sum = 0.0;
            for (std::size_t node = 0; node < 4; ++node) {
                shape_sum += point.shape[node];
                gradient_r_sum += point.gradient_r[node];
                gradient_z_sum += point.gradient_z[node];
            }
            passed =
                check(std::abs(shape_sum - 1.0) < 1.0e-14, "Quad4 shape functions form a partition of unity") && passed;
            passed = check(std::abs(gradient_r_sum) < 1.0e-11 && std::abs(gradient_z_sum) < 1.0e-11,
                         "Quad4 physical gradients sum to zero")
                     && passed;
            passed = check(point.radius > 0.0 && point.weighted_measure > 0.0, "RZ quadrature measure is positive")
                     && passed;
            integrated_volume += point.weighted_measure;
        }
    }
    const double exact_volume = pi * (outer * outer - inner * inner) * length;
    const double volume_error = std::abs(integrated_volume - exact_volume) / exact_volume;
    std::cout << "geometry_volume_relative_error=" << volume_error << '\n';
    passed = check(volume_error < 1.0e-13, "RZ quadrature integrates annular volume") && passed;
    return passed;
}

bool test_element_jacobian() {
    const fuelsim::Quad4Coordinates coordinates = {{
        {0.001, 0.0},
        {0.003, 0.0},
        {0.003, 0.004},
        {0.001, 0.004},
    }};
    const fuelsim::Quad4RzGeometry geometry = fuelsim::make_quad4_rz_geometry(coordinates);
    const fuelsim::Quad4RzData data{fuelsim::IsotropicThermoelasticMaterial(properties()),
        2.0e8,
        1.25,
        fuelsim::StrainFormulation::small};
    const fuelsim::LocalValues state = {
        710.0,
        680.0,
        650.0,
        690.0,
        0.0,
        2.0e-6,
        2.5e-6,
        0.4e-6,
        0.0,
        -0.2e-6,
        3.0e-6,
        2.5e-6,
    };
    const fuelsim::LocalValues direction = {
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
    const fuelsim::rz::LocalLinearization system = fuelsim::rz::linearize_quad4_rz_thermoelastic(data, geometry, state);
    constexpr double step = 1.0e-4;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual = fuelsim::compute_quad4_rz_thermoelastic(data, geometry, plus);
    const fuelsim::LocalResidual minus_residual = fuelsim::compute_quad4_rz_thermoelastic(data, geometry, minus);
    bool passed = true;
    const std::array<fuelsim::AxisymmetricStressValues, 4> stresses =
        fuelsim::compute_quad4_rz_thermoelastic_stress(data, geometry, state);
    passed = check(std::all_of(stresses.begin(),
                       stresses.end(),
                       [](const fuelsim::AxisymmetricStressValues& stress) {
                           return std::isfinite(stress.rr) && std::isfinite(stress.zz) && std::isfinite(stress.hoop)
                                  && std::isfinite(stress.rz);
                       }),
                 "procedural Quad4 RZ kernel returns finite stress at every quadrature point")
             && passed;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::quad4_local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                     "AD element Jacobian row " + std::to_string(row) + " matches centered finite difference")
                 && passed;
    }
    std::cout << "element_jacobian_maximum_scaled_error=" << maximum_jacobian_error << '\n';
    double heat_to_displacement = 0.0;
    double temperature_to_mechanics = 0.0;
    for (std::size_t row = 0; row < 4; ++row)
        for (std::size_t column = 4; column < 12; ++column)
            heat_to_displacement += std::abs(system.jacobian[row * 12 + column]);
    for (std::size_t row = 4; row < 12; ++row)
        for (std::size_t column = 0; column < 4; ++column)
            temperature_to_mechanics += std::abs(system.jacobian[row * 12 + column]);
    passed = check(heat_to_displacement == 0.0,
                 "reference-mesh thermal residual "
                 "does not depend on displacement")
             && passed;
    passed =
        check(temperature_to_mechanics > 0.0, "thermal expansion produces temperature-mechanics coupling") && passed;
    return passed;
}

bool test_cax_kinematics_and_jacobian(bool reduced) {
    const std::string name = reduced ? "CAX4RT" : "CAX4T";
    const std::size_t points = reduced ? 1 : 4;
    const fuelsim::Quad4Coordinates coordinates = {{{1.0, 0.0}, {2.1, 0.1}, {2.0, 1.2}, {0.9, 1.0}}};
    const auto geometry = fuelsim::make_quad4_rz_geometry(coordinates);
    fuelsim::Quad4RzData data{fuelsim::IsotropicThermoelasticMaterial(properties())};
    data.element_formulation = reduced ? fuelsim::RzElementFormulation::cax4rt : fuelsim::RzElementFormulation::cax4t;
    const fuelsim::LocalValues direction = {0.2, -0.3, 0.4, -0.1, 0.3, -0.5, 0.2, 0.4, -0.2, 0.35, -0.45, 0.25};
    bool passed = true;
    if (!reduced) {
        // Fully restrained heating has the same hydrostatic stress at every
        // point, even on a distorted element with nonuniform corner temperatures.
        for (const auto formulation : {fuelsim::StrainFormulation::small, fuelsim::StrainFormulation::finite}) {
            data.strain_formulation = formulation;
            fuelsim::LocalValues initial{}, old{}, heated{};
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
            const auto history = fuelsim::compute_quad4_rz_transient_update(data, geometry, old, initial, {}, 0.1);
            const auto next = fuelsim::compute_quad4_rz_transient_update(data, geometry, heated, old, history, 0.1);
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
        fuelsim::LocalValues initial{}, old{}, state{};
        for (std::size_t n = 0; n < 4; ++n) {
            initial[n] = old[n] = state[n] = 600.0;
            old[4 + n] = 0.03 * coordinates[n].r;
            old[8 + n] = -0.02 * coordinates[n].z;
            state[4 + n] = 0.08 * coordinates[n].r;
            state[8 + n] = -0.04 * coordinates[n].z;
        }
        const auto old_history = fuelsim::compute_quad4_rz_transient_update(data, geometry, old, initial, {}, 0.1);
        const auto affine = fuelsim::compute_quad4_rz_transient_update(data, geometry, state, old, old_history, 0.1);
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
            passed = check(fuelsim::rz::cax4rt_hourglass_energy(data, geometry, state) < 1e-20,
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
            const auto residual = fuelsim::compute_quad4_rz_transient(data, geometry, heated, unheated, {}, 0.1);
            for (std::size_t n = 0; n < 4; ++n)
                passed = check(std::abs(residual[n]) < 1e-10 * data.volumetric_heat_source,
                             "CAX4RT uniform heating balances uniform source on a distorted element")
                         && passed;
            data.volumetric_heat_source = saved_source;
        }
        state[0] += 3.0;
        state[2] -= 2.0;
        fuelsim::LocalJacobian jacobian{};
        const auto active =
            fuelsim::compute_quad4_rz_transient(data, geometry, state, old, old_history, 0.1, &jacobian);
        const auto passive = fuelsim::compute_quad4_rz_transient(data, geometry, state, old, old_history, 0.1);
        if (!reduced) {
            const auto no_capacity =
                fuelsim::compute_quad4_rz_transient(data, geometry, state, old, old_history, 0.1, nullptr, false);
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
        const auto rp = fuelsim::compute_quad4_rz_transient(data, geometry, plus, old, old_history, 0.1);
        const auto rm = fuelsim::compute_quad4_rz_transient(data, geometry, minus, old, old_history, 0.1);
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
        const auto shifted = fuelsim::compute_quad4_rz_transient(data, geometry, translated, old, old_history, 0.1);
        double translation_error = 0.0;
        for (std::size_t i = 0; i < 12; ++i)
            translation_error = std::max(translation_error, scaled_error(shifted[i], passive[i]));
        passed = check(translation_error < 1e-12, name + " internal forces are invariant under axial rigid translation")
                 && passed;
        if (reduced) {
            const double energy = fuelsim::rz::cax4rt_hourglass_energy(data, geometry, state);
            const double shifted_energy = fuelsim::rz::cax4rt_hourglass_energy(data, geometry, translated);
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
                (void)fuelsim::compute_quad4_rz_transient(data, geometry, invalid_midpoint, initial, {}, 0.1);
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
            fuelsim::LocalValues initial{}, old{}, state{};
            for (std::size_t n = 0; n < 4; ++n) {
                initial[n] = old[n] = state[n] = 600.0;
                old[4 + n] = 0.001 * coordinates[n].r;
                old[8 + n] = -0.001 * coordinates[n].z;
                state[4 + n] = 0.004 * coordinates[n].r;
                state[8 + n] = -0.005 * coordinates[n].z;
            }
            state[5] += 0.0003;
            state[10] -= 0.0002;
            const auto history = fuelsim::compute_quad4_rz_transient_update(data, geometry, old, initial, {}, 0.1);
            fuelsim::LocalJacobian jacobian{};
            const auto residual =
                fuelsim::compute_quad4_rz_transient(data, geometry, state, old, history, 0.1, &jacobian);
            passed = check(residual == fuelsim::compute_quad4_rz_transient(data, geometry, state, old, history, 0.1),
                         name + " inelastic material residual agrees exactly between passive and Jacobian paths")
                     && passed;
            auto plus = state, minus = state;
            constexpr double step = 1e-5;
            for (std::size_t j = 0; j < 12; ++j) {
                plus[j] += step * direction[j];
                minus[j] -= step * direction[j];
            }
            const auto rp = fuelsim::compute_quad4_rz_transient(data, geometry, plus, old, history, 0.1);
            const auto rm = fuelsim::compute_quad4_rz_transient(data, geometry, minus, old, history, 0.1);
            double error = 0.0;
            for (std::size_t i = 0; i < 12; ++i) {
                double ad = 0.0;
                for (std::size_t j = 0; j < 12; ++j)
                    ad += jacobian[12 * i + j] * direction[j];
                error = std::max(error, scaled_error(ad, (rp[i] - rm[i]) / (2.0 * step)));
            }
            const auto next = fuelsim::compute_quad4_rz_transient_update(data, geometry, state, old, history, 0.1);
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

bool test_finite_strain_kinematics_and_jacobian() {
    const fuelsim::Quad4Coordinates coordinates = {{
        {1.0, 0.0},
        {2.0, 0.0},
        {2.0, 1.0},
        {1.0, 1.0},
    }};
    const fuelsim::Quad4RzGeometry geometry = fuelsim::make_quad4_rz_geometry(coordinates);
    const fuelsim::Quad4RzData data{fuelsim::IsotropicThermoelasticMaterial(properties()),
        0.0,
        0.0,
        fuelsim::StrainFormulation::finite};
    constexpr double radial_stretch = 1.08;
    constexpr double axial_stretch = 0.96;
    fuelsim::LocalValues uniform_state{};
    for (std::size_t node = 0; node < 4; ++node) {
        uniform_state[node] = 600.0;
        uniform_state[4 + node] = (radial_stretch - 1.0) * coordinates[node].r;
        uniform_state[8 + node] = (axial_stretch - 1.0) * coordinates[node].z;
    }
    fuelsim::LocalAdValues passive_state{};
    for (std::size_t dof = 0; dof < uniform_state.size(); ++dof)
        passive_state[dof] = uniform_state[dof];
    bool passed = true;
    double maximum_strain_error = 0.0;
    double maximum_measure_error = 0.0;
    for (const fuelsim::RzQuadraturePoint& point : geometry.points) {
        const fuelsim::AxisymmetricKinematics kinematics =
            fuelsim::evaluate_axisymmetric_kinematics(point, passive_state, fuelsim::StrainFormulation::finite);
        const auto taylor_increment = [](double stretch) {
            const double cinv = 1.0 / (stretch * stretch) - 1.0;
            return -0.5 * cinv + 0.25 * cinv * cinv;
        };
        maximum_strain_error = std::max({maximum_strain_error,
            std::abs(kinematics.strain_rr.value() - taylor_increment(radial_stretch)),
            std::abs(kinematics.strain_zz.value() - taylor_increment(axial_stretch)),
            std::abs(kinematics.strain_hoop.value() - taylor_increment(radial_stretch)),
            std::abs(kinematics.strain_rz.value())});
        const double expected_measure = point.weighted_measure * radial_stretch * radial_stretch * axial_stretch;
        maximum_measure_error = std::max(maximum_measure_error,
            std::abs(kinematics.weighted_measure.value() - expected_measure) / expected_measure);
    }
    passed = check(maximum_strain_error < 1.0e-14,
                 "finite RZ uniform stretches give MOOSE Taylor strain "
                 "increments")
             && passed;
    passed = check(maximum_measure_error < 1.0e-14,
                 "finite RZ current measure follows the deformation "
                 "Jacobian")
             && passed;
    constexpr double old_radial_stretch = 1.03;
    constexpr double old_axial_stretch = 0.98;
    fuelsim::LocalValues committed_state{};
    for (std::size_t node = 0; node < 4; ++node) {
        committed_state[node] = 600.0;
        committed_state[4 + node] = (old_radial_stretch - 1.0) * coordinates[node].r;
        committed_state[8 + node] = (old_axial_stretch - 1.0) * coordinates[node].z;
    }
    double maximum_incremental_error = 0.0;
    for (const fuelsim::RzQuadraturePoint& point : geometry.points) {
        const fuelsim::AxisymmetricKinematics kinematics = fuelsim::evaluate_axisymmetric_incremental_kinematics(point,
            passive_state,
            committed_state,
            fuelsim::StrainFormulation::finite);
        const auto taylor_increment = [](double stretch) {
            const double cinv = 1.0 / (stretch * stretch) - 1.0;
            return -0.5 * cinv + 0.25 * cinv * cinv;
        };
        maximum_incremental_error = std::max({maximum_incremental_error,
            std::abs(kinematics.strain_rr.value() - taylor_increment(radial_stretch / old_radial_stretch)),
            std::abs(kinematics.strain_zz.value() - taylor_increment(axial_stretch / old_axial_stretch)),
            std::abs(kinematics.strain_hoop.value() - taylor_increment(radial_stretch / old_radial_stretch)),
            std::abs(kinematics.strain_rz.value())});
    }
    passed = check(maximum_incremental_error < 1.0e-14,
                 "finite RZ strain uses current-to-committed MOOSE "
                 "incremental deformation")
             && passed;
    fuelsim::LocalValues state = uniform_state;
    state[5] += 0.025;
    state[6] += 0.015;
    state[10] -= 0.020;
    state[11] += 0.010;
    const fuelsim::LocalValues direction = {
        0.2,
        -0.3,
        0.4,
        -0.1,
        0.3,
        -0.5,
        0.2,
        0.4,
        -0.2,
        0.35,
        -0.45,
        0.25,
    };
    const fuelsim::rz::LocalLinearization system = fuelsim::rz::linearize_quad4_rz_thermoelastic(data, geometry, state);
    constexpr double step = 1.0e-5;
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual = fuelsim::compute_quad4_rz_thermoelastic(data, geometry, plus);
    const fuelsim::LocalResidual minus_residual = fuelsim::compute_quad4_rz_thermoelastic(data, geometry, minus);
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::quad4_local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        maximum_jacobian_error = std::max(maximum_jacobian_error, scaled_error(ad_direction, finite_difference));
    }
    passed = check(maximum_jacobian_error < 2.0e-7,
                 "finite RZ AD Jacobian matches centered finite "
                 "difference")
             && passed;
    fuelsim::LocalValues inverted = uniform_state;
    for (std::size_t node = 0; node < 4; ++node)
        inverted[4 + node] = -1.1 * coordinates[node].r + 2.0;
    bool inversion_rejected = false;
    try {
        (void)fuelsim::compute_quad4_rz_thermoelastic(data, geometry, inverted);
    } catch (const std::domain_error&) {
        inversion_rejected = true;
    }
    passed = check(inversion_rejected,
                 "finite RZ rejects a nonpositive in-plane Jacobian while "
                 "the current radius remains positive")
             && passed;
    fuelsim::LocalValues collapsed_radius = uniform_state;
    for (std::size_t node = 0; node < 4; ++node)
        collapsed_radius[4 + node] = -2.1;
    bool radius_rejected = false;
    try {
        (void)fuelsim::compute_quad4_rz_thermoelastic(data, geometry, collapsed_radius);
    } catch (const std::domain_error&) {
        radius_rejected = true;
    }
    passed = check(radius_rejected,
                 "finite RZ rejects nonpositive hoop stretch and current "
                 "radius while the in-plane Jacobian remains positive")
             && passed;
    std::cout << "finite_strain_uniform_taylor_increment_error=" << maximum_strain_error << '\n';
    std::cout << "finite_strain_current_measure_relative_error=" << maximum_measure_error << '\n';
    std::cout << "finite_strain_incremental_taylor_error=" << maximum_incremental_error << '\n';
    std::cout << "finite_strain_directional_jacobian_error=" << maximum_jacobian_error << '\n';
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

fuelsim::LocalResidual summed_heat_residual(const fuelsim::GapHeatProperties& properties,
    const HeatPointGeometries& geometries,
    const fuelsim::LocalValues& state) {
    fuelsim::LocalResidual result{};
    for (const fuelsim::Line2RzHeatPointGeometry& geometry : geometries) {
        const fuelsim::LocalResidual point = fuelsim::compute_line2_rz_gap_heat(properties, geometry, state);
        for (std::size_t row = 0; row < result.size(); ++row)
            result[row] += point[row];
    }
    return result;
}

bool test_heat_interface_case(const std::string& name,
    const fuelsim::GapHeatProperties& properties,
    const HeatPointGeometries& geometries,
    const fuelsim::LocalValues& state) {
    const fuelsim::LocalValues direction = {
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
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual = summed_heat_residual(properties, geometries, plus);
    const fuelsim::LocalResidual minus_residual = summed_heat_residual(properties, geometries, minus);
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
            ad_direction += system.jacobian[row * fuelsim::local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                     name + " interface AD Jacobian row " + std::to_string(row) + " matches centered finite difference")
                 && passed;
    }
    const fuelsim::LocalResidual residual = summed_heat_residual(properties, geometries, state);
    double thermal_sum = 0.0;
    double thermal_scale = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
        thermal_sum += residual[row];
        thermal_scale += std::abs(residual[row]);
    }
    passed =
        check(std::abs(thermal_sum) < 1.0e-13 * (1.0 + thermal_scale), name + " interface conserves heat") && passed;
    for (std::size_t row = 4; row < fuelsim::local_dof_count; ++row)
        passed = check(residual[row] == 0.0, name + " gap heat kernel has no mechanical residual") && passed;
    std::cout << name << "_heat_jacobian_maximum_scaled_error=" << maximum_jacobian_error << '\n';
    return passed;
}

bool test_heat_point_interface_case(const std::string& name,
    const fuelsim::GapHeatProperties& properties,
    const fuelsim::Line2RzHeatPointGeometry& geometry,
    const fuelsim::LocalValues& state) {
    const fuelsim::LocalValues direction = {
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
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual = fuelsim::compute_line2_rz_gap_heat(properties, geometry, plus);
    const fuelsim::LocalResidual minus_residual = fuelsim::compute_line2_rz_gap_heat(properties, geometry, minus);
    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                     name + " point AD Jacobian row " + std::to_string(row) + " matches centered finite difference")
                 && passed;
    }
    const fuelsim::LocalResidual residual = fuelsim::compute_line2_rz_gap_heat(properties, geometry, state);
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
    const fuelsim::LocalValues& state,
    const fuelsim::ContactPointHistory& history = {}) {
    const fuelsim::LocalValues direction = {
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
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, plus, state, history);
    const fuelsim::LocalResidual minus_residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, minus, state, history);
    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                     name + " contact AD Jacobian row " + std::to_string(row) + " matches centered finite difference")
                 && passed;
    }
    const fuelsim::LocalResidual residual =
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
    const fuelsim::LocalValues& state,
    const fuelsim::LocalValues& committed_state,
    const fuelsim::ContactPointHistory& history) {
    const fuelsim::LocalValues direction = {
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
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, plus, committed_state, history);
    const fuelsim::LocalResidual minus_residual =
        fuelsim::compute_node_to_line_rz_contact(properties, geometry, minus, committed_state, history);
    bool passed = true;
    double maximum_jacobian_error = 0.0;
    for (std::size_t row = 0; row < system.residual.size(); ++row) {
        double ad_direction = 0.0;
        for (std::size_t column = 0; column < direction.size(); ++column)
            ad_direction += system.jacobian[row * fuelsim::local_dof_count + column] * direction[column];
        const double finite_difference = (plus_residual[row] - minus_residual[row]) / (2.0 * step);
        const double error = scaled_error(ad_direction, finite_difference);
        maximum_jacobian_error = std::max(maximum_jacobian_error, error);
        passed = check(error < 1.0e-7,
                     name + " friction AD Jacobian row " + std::to_string(row) + " matches centered finite difference")
                 && passed;
    }
    const fuelsim::LocalResidual residual =
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
    const fuelsim::LocalValues open_state = {
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
    const fuelsim::LocalValues minimum_gap_state = {
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
    const fuelsim::LocalValues closed_state = {
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
    fuelsim::LocalValues outside_state = closed_state;
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
    fuelsim::LocalValues radial_endpoint_state = closed_state;
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
    fuelsim::LocalValues current_sloped_state = closed_state;
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
    const fuelsim::LocalResidual current_sloped_heat_residual =
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
    fuelsim::LocalValues thermal_projection_lost_state = open_state;
    thermal_projection_lost_state[10] = 0.45e-3;
    thermal_projection_lost_state[11] = -0.45e-3;
    for (const fuelsim::Line2RzHeatPointGeometry& geometry : heat_geometries) {
        const fuelsim::HeatQuadratureValue value =
            fuelsim::compute_line2_rz_gap_heat_value(heat_properties, geometry, thermal_projection_lost_state);
        const fuelsim::LocalResidual residual =
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
    const fuelsim::LocalResidual current_sloped_residual = fuelsim::compute_node_to_line_rz_contact(contact_properties,
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
    fuelsim::LocalValues vertex_state{};
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
    fuelsim::LocalValues general_vertex_state{};
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
    const fuelsim::LocalValues axial_open_state = {
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
    const fuelsim::LocalValues axial_closed_state = {
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
    const fuelsim::LocalValues sloped_closed_state = {
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
    const fuelsim::LocalResidual sloped_residual = fuelsim::compute_node_to_line_rz_contact(contact_properties,
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
    fuelsim::LocalValues sticking_state = closed_state;
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
    fuelsim::LocalValues sliding_state = closed_state;
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
    fuelsim::LocalValues resticking_state = closed_state;
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
    fuelsim::LocalValues sloped_sticking_state = sloped_closed_state;
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
    fuelsim::LocalValues sloped_sliding_state = sloped_closed_state;
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
    fuelsim::LocalValues sloped_resticking_state = sloped_closed_state;
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
    const fuelsim::LocalResidual legacy_normal =
        fuelsim::compute_node_to_line_rz_contact(contact_properties, contact_geometry, closed_state, closed_state, {});
    const fuelsim::LocalResidual zero_friction = fuelsim::compute_node_to_line_rz_contact(explicit_zero_friction,
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
    const fuelsim::LocalValues nodal_state = {500, 650, 300, 350, .01, .02, .03, -.01, .01, .02, -.02, .01};
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
    fuelsim::LocalValues state = {
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

fuelsim::UnstructuredQuad4Mesh thermal_owner_transfer_mesh() {
    return fuelsim::UnstructuredQuad4Mesh(
        {
            {0.0, 0.5},
            {1.0, 0.5},
            {1.0, 1.5},
            {0.0, 1.5},
            {1.1, 0.0},
            {1.2, 0.0},
            {1.2, 1.0},
            {1.1, 1.0},
            {1.2, 2.0},
            {1.1, 2.0},
        },
        {{{{0, 1, 2, 3}}}, {{{4, 5, 6, 7}}}, {{{7, 6, 8, 9}}}},
        {1, 2, 2},
        {{1, "secondary"}, {2, "primary"}},
        {},
        {{11, "secondary_face", {{{0, 1}}}}, {21, "primary_face", {{{1, 3}, {2, 3}}}}});
}

bool test_thermal_owner_transfer_assembly() {
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back({"secondary", "secondary", properties(), 0.0, 500.0});
    definition.regions.push_back({"primary", "primary", properties(), 0.0, 300.0});
    definition.contacts.push_back(
        {"thermal_slide", "primary_face", "secondary_face", true, false, 0.2, 1.0e-5, 1.0e14});
    fuelsim::SteadyProblem problem(std::move(definition), thermal_owner_transfer_mesh());
    std::vector<double> state = problem.initial_state();
    const auto active_primary_centers = [&](const std::vector<double>& current) {
        problem.validate_state(current);
        std::vector<double> centers;
        for (std::size_t contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
            contribution < problem.contribution_count();
            ++contribution) {
            if (fuelsim::rz::ProblemAccess::contribution_type(problem, contribution)
                != fuelsim::SpatialContributionType::thermal_contact)
                continue;
            const fuelsim::LocalResidual residual = fuelsim::rz::ProblemAccess::contribution_residual(problem,
                contribution,
                fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, current));
            double magnitude = 0.0;
            for (std::size_t row = 0; row < 4; ++row)
                magnitude += std::abs(residual[row]);
            if (magnitude == 0.0)
                continue;
            const fuelsim::LocalDofs dofs = fuelsim::rz::ProblemAccess::contribution_dofs(problem, contribution);
            const std::size_t primary_offset = fuelsim::rz::ProblemAccess::region_node_offset(problem, 1);
            const std::size_t first = dofs[2] - primary_offset;
            const std::size_t second = dofs[3] - primary_offset;
            centers.push_back(0.5
                              * (fuelsim::rz::ProblemAccess::region_mesh(problem, 1).nodes()[first].z
                                  + fuelsim::rz::ProblemAccess::region_mesh(problem, 1).nodes()[second].z));
        }
        return centers;
    };
    const std::vector<double> initial_centers = active_primary_centers(state);
    bool passed =
        check(initial_centers.size() == 4 && std::count(initial_centers.begin(), initial_centers.end(), 0.5) == 2
                  && std::count(initial_centers.begin(), initial_centers.end(), 1.5) == 2,
            "each initial thermal integration point has exactly one primary-segment owner");
    const std::size_t secondary_offset = fuelsim::rz::ProblemAccess::region_node_offset(problem, 0);
    for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes().size(); ++node)
        state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement,
            secondary_offset + node)] = 0.25;
    const std::vector<double> split_centers = active_primary_centers(state);
    passed = check(split_centers.size() == 4 && std::count(split_centers.begin(), split_centers.end(), 0.5) == 1
                       && std::count(split_centers.begin(), split_centers.end(), 1.5) == 3,
                 "two integration points from one reference fragment may select different current primary segments")
             && passed;
    double assembled_heat_sum = 0.0;
    double assembled_heat_scale = 0.0;
    for (std::size_t contribution = fuelsim::rz::ProblemAccess::volume_contribution_count(problem);
        contribution < problem.contribution_count();
        ++contribution) {
        if (fuelsim::rz::ProblemAccess::contribution_type(problem, contribution)
            != fuelsim::SpatialContributionType::thermal_contact)
            continue;
        const fuelsim::LocalResidual residual = fuelsim::rz::ProblemAccess::contribution_residual(problem,
            contribution,
            fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, state));
        for (std::size_t row = 0; row < 4; ++row) {
            assembled_heat_sum += residual[row];
            assembled_heat_scale += std::abs(residual[row]);
        }
    }
    passed = check(std::abs(assembled_heat_sum) < 1.0e-13 * (1.0 + assembled_heat_scale),
                 "split-owner thermal candidates conserve the assembled interface heat rate")
             && passed;
    for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes().size(); ++node)
        state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement,
            secondary_offset + node)] = 0.6;
    const std::vector<double> slid_centers = active_primary_centers(state);
    passed =
        check(slid_centers.size() == 4
                  && std::all_of(slid_centers.begin(), slid_centers.end(), [](double center) { return center == 1.5; }),
            "all thermal integration points transfer to the current upper primary segment after large sliding")
        && passed;
    const fuelsim::InterfaceSummary summary = fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, state);
    passed = check(std::isfinite(summary.total_heat_rate) && summary.total_heat_rate > 0.0,
                 "large-sliding thermal owner transfer retains a finite positive conservative heat rate")
             && passed;
    for (std::size_t node = 0; node < fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes().size(); ++node)
        state[fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::axial_displacement,
            secondary_offset + node)] = 2.0;
    bool rejected = false;
    try {
        problem.validate_state(state);
    } catch (const std::domain_error&) {
        rejected = true;
    }
    bool summary_rejected = false;
    try {
        (void)fuelsim::rz::ProblemAccess::summarize_interface(problem, 0, state);
    } catch (const std::domain_error&) {
        summary_rejected = true;
    }
    passed = check(rejected && summary_rejected,
                 "validation and interface diagnostics reject a thermal integration point that leaves the complete "
                 "chain")
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
    const fuelsim::LocalValues zero_gap_open_state = {
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
    const fuelsim::LocalValues zero_gap_closed_state = {
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
    const fuelsim::LocalValues oracle_closed_state = {
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
    const fuelsim::LocalValues zero_gap_heat_open_state = {
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

bool test_m1_dof_layout() {
    constexpr std::size_t fuel_radial_elements = 2;
    constexpr std::size_t cladding_radial_elements = 1;
    constexpr std::size_t axial_elements = 2;
    constexpr double heat_source = 1.0e8;
    const fuelsim::UnstructuredQuad4Mesh source = fuelsim::test::make_disconnected_annular_mesh(
        {{1, "fuel", 0.0, 0.004, 0.010, fuel_radial_elements, axial_elements},
            {2, "clad", 0.0041, 0.0046, 0.01002, cladding_radial_elements, axial_elements}});
    fuelsim::SpatialDefinition definition;
    definition.regions.push_back({"fuel", "fuel", properties(), heat_source, 600.0});
    definition.regions.push_back(
        {"clad", "clad", fuelsim::test::thermoelastic(0.0, 16.0, 75.0e9, 0.3, 5.0e-6, 600.0), 0.0, 600.0});
    definition.contacts.push_back({"fuel_clad", "clad_inner", "fuel_outer", true, true, 0.4, 1.0e-6, 1.0e14});
    definition.boundary_conditions.push_back({"fuel_axis",
        fuelsim::BoundaryConditionType::dirichlet,
        "fuel_inner",
        fuelsim::Field::radial_displacement,
        0.0});
    definition.boundary_conditions.push_back({"fuel_bottom",
        fuelsim::BoundaryConditionType::dirichlet,
        "fuel_bottom",
        fuelsim::Field::axial_displacement,
        0.0});
    definition.boundary_conditions.push_back({"clad_bottom",
        fuelsim::BoundaryConditionType::dirichlet,
        "clad_bottom",
        fuelsim::Field::axial_displacement,
        0.0});
    definition.boundary_conditions.push_back({"clad_temperature",
        fuelsim::BoundaryConditionType::dirichlet,
        "clad_outer",
        fuelsim::Field::temperature,
        600.0});
    definition.boundary_conditions.push_back({"clad_pressure",
        fuelsim::BoundaryConditionType::pressure,
        "clad_outer",
        fuelsim::Field::radial_displacement,
        1.0e5});
    definition.boundary_conditions.push_back({"clad_traction",
        fuelsim::BoundaryConditionType::traction,
        "clad_outer",
        fuelsim::Field::axial_displacement,
        1.0e3});
    fuelsim::BoundaryConditionDefinition convection{"clad_convection",
        fuelsim::BoundaryConditionType::convection,
        "clad_outer",
        fuelsim::Field::temperature,
        0.0};
    convection.heat_transfer_coefficient = 100.0;
    convection.ambient_temperature = 600.0;
    definition.boundary_conditions.push_back(std::move(convection));
    fuelsim::SteadyProblem problem(std::move(definition), source);
    bool passed = true;
    constexpr std::size_t contribution_type_count =
        static_cast<std::size_t>(fuelsim::SpatialContributionType::convection) + 1;
    std::array<std::size_t, contribution_type_count> contribution_counts{};
    std::size_t thermal_contributions = 0;
    std::size_t mechanical_contributions = 0;
    std::size_t first_thermal = problem.contribution_count();
    std::size_t previous_type = 0;
    const std::vector<double> initial_state = problem.initial_state();
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        const fuelsim::SpatialContributionType type =
            fuelsim::rz::ProblemAccess::contribution_type(problem, contribution);
        const std::size_t type_index = static_cast<std::size_t>(type);
        ++contribution_counts.at(type_index);
        passed =
            check(contribution == 0 || type_index >= previous_type, "spatial contribution categories are contiguous")
            && passed;
        previous_type = type_index;
        const fuelsim::LocalValues local =
            fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, initial_state);
        const fuelsim::rz::LocalLinearization system =
            fuelsim::rz::ProblemAccess::linearize_contribution(problem, contribution, local);
        const bool finite_residual = std::all_of(system.residual.begin(), system.residual.end(), [](double value) {
            return std::isfinite(value);
        });
        const bool finite_jacobian = std::all_of(system.jacobian.begin(), system.jacobian.end(), [](double value) {
            return std::isfinite(value);
        });
        passed = check(finite_residual && finite_jacobian,
                     "every spatial contribution routes to finite residual "
                     "and Jacobian values")
                 && passed;
        if (type == fuelsim::SpatialContributionType::thermal_contact) {
            first_thermal = std::min(first_thermal, contribution);
            ++thermal_contributions;
        } else if (type == fuelsim::SpatialContributionType::mechanical_contact) {
            ++mechanical_contributions;
        }
    }
    constexpr std::array expected_types = {fuelsim::SpatialContributionType::volume,
        fuelsim::SpatialContributionType::thermal_contact,
        fuelsim::SpatialContributionType::mechanical_contact,
        fuelsim::SpatialContributionType::pressure,
        fuelsim::SpatialContributionType::traction,
        fuelsim::SpatialContributionType::convection};
    passed = check(std::all_of(expected_types.begin(),
                       expected_types.end(),
                       [&](fuelsim::SpatialContributionType type) {
                           return contribution_counts.at(static_cast<std::size_t>(type)) > 0;
                       }),
                 "volume, thermal contact, mechanical contact, pressure, traction, and convection contributions are "
                 "all routed")
             && passed;
    bool rejected_out_of_range = false;
    try {
        static_cast<void>(fuelsim::rz::ProblemAccess::contribution_dofs(problem, problem.contribution_count()));
    } catch (const std::out_of_range&) {
        rejected_out_of_range = true;
    }
    passed = check(rejected_out_of_range, "spatial contribution routing rejects the end index") && passed;
    passed = check(thermal_contributions == 2 * (axial_elements + 1),
                 "every M1 STS integration point forms one active thermal contribution")
             && passed;
    passed = check(mechanical_contributions == 2 * axial_elements,
                 "M1 has one active NTS contribution per secondary half-node")
             && passed;
    passed = check(problem.sparsity_contribution_count() > problem.contribution_count(),
                 "M1 reserves every primary candidate block outside the active contribution loop")
             && passed;
    passed = check(problem.dof_count()
                       == 3
                              * (fuelsim::rz::ProblemAccess::region_mesh(problem, 0).nodes().size()
                                  + fuelsim::rz::ProblemAccess::region_mesh(problem, 1).nodes().size()),
                 "M1 uses one field-major map for both independent meshes")
             && passed;
    const fuelsim::LocalDofs interface = fuelsim::rz::ProblemAccess::contribution_dofs(problem, first_thermal);
    const std::size_t fuel_outer = fuelsim::test::annular_node_id(fuel_radial_elements, fuel_radial_elements, 0);
    const std::size_t cladding_inner = fuelsim::test::annular_node_id(cladding_radial_elements, 0, 0);
    passed = check(interface[0]
                           == fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature,
                               fuelsim::rz::ProblemAccess::region_node_offset(problem, 0) + fuel_outer)
                       && interface[2]
                              == fuelsim::rz::ProblemAccess::dof_map(problem).dof(fuelsim::Field::temperature,
                                  fuelsim::rz::ProblemAccess::region_node_offset(problem, 1) + cladding_inner),
                 "M1 interface DOFs preserve fuel/cladding node ownership")
             && passed;
    const std::vector<fuelsim::ContactNodeSummary> contact_nodes =
        fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, initial_state);
    passed = check(std::all_of(contact_nodes.begin(),
                       contact_nodes.end(),
                       [](const fuelsim::ContactNodeSummary& node) { return node.projected; }),
                 "taller cladding contains every initial NTS projection")
             && passed;
    const fuelsim::LocalValues initial_element_state =
        fuelsim::rz::ProblemAccess::contribution_state(problem, 0, initial_state);
    const fuelsim::LocalResidual source_residual =
        fuelsim::rz::ProblemAccess::contribution_residual(problem, 0, initial_element_state);
    problem.set_load_factor(2.0);
    const fuelsim::LocalResidual doubled_source_residual =
        fuelsim::rz::ProblemAccess::contribution_residual(problem, 0, initial_element_state);
    for (std::size_t node = 0; node < fuelsim::quad4_node_count; ++node) {
        passed = check(std::abs(doubled_source_residual[node] - 2.0 * source_residual[node])
                           < 1.0e-12 * (1.0 + std::abs(source_residual[node])),
                     "M1 updates heat loading without rebuilding geometry")
                 && passed;
    }
    passed = check(fuelsim::rz::ProblemAccess::definition(problem).regions[0].volumetric_heat_source == heat_source
                       && fuelsim::rz::ProblemAccess::region_kernel_data(problem, 0).volumetric_heat_source
                              == 2.0 * heat_source,
                 "M1 load factor updates the production heat-source kernel")
             && passed;
    return passed;
}

bool test_time_table_and_convection() {
    const fuelsim::PiecewiseLinearTimeTable table("power", {0.0, 2.0, 5.0}, {0.0, 1.0, 0.4});
    bool passed =
        check(table.value(0.0) == 0.0 && table.value(1.0) == 0.5 && table.value(3.0) == 0.8 && table.value(8.0) == 0.4,
            "piecewise-linear table interpolates and holds endpoints");
    passed = check(table.average_value(0.0, 1.0) == 0.25 && std::abs(table.average_value(1.0, 3.0) - 0.825) < 1.0e-15
                       && std::abs(table.average_value(4.0, 8.0) - 0.425) < 1.0e-15,
                 "piecewise-linear table averages exactly across knots and held endpoints")
             && passed;
    const fuelsim::Line2RzBoundaryGeometry geometry =
        fuelsim::make_line2_rz_boundary_geometry({{{0.005, 0.0}, {0.005, 0.01}}}, {{1, 2}});
    const fuelsim::Line2RzBoundaryData data = {fuelsim::Line2RzBoundaryKind::convection,
        fuelsim::TractionComponent::radial,
        1000.0,
        500.0,
        false};
    const fuelsim::LocalValues state = {
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
    const fuelsim::LocalValues direction = {
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
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual = fuelsim::compute_line2_rz_boundary(data, geometry, plus);
    const fuelsim::LocalResidual minus_residual = fuelsim::compute_line2_rz_boundary(data, geometry, minus);
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

bool test_follower_pressure() {
    const fuelsim::Line2RzBoundaryGeometry geometry =
        fuelsim::make_line2_rz_boundary_geometry({{{0.005, 0.0}, {0.005, 0.01}}}, {{1, 2}});
    constexpr double pressure = 3.0e6;
    const fuelsim::Line2RzBoundaryData follower = {fuelsim::Line2RzBoundaryKind::pressure,
        fuelsim::TractionComponent::radial,
        pressure,
        0.0,
        true};
    const fuelsim::LocalValues state = {
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
    const fuelsim::LocalValues direction = {
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
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual = fuelsim::compute_line2_rz_boundary(follower, geometry, plus);
    const fuelsim::LocalResidual minus_residual = fuelsim::compute_line2_rz_boundary(follower, geometry, minus);
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
    const fuelsim::LocalValues state = {
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
    const fuelsim::LocalValues direction = {
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
    fuelsim::LocalValues plus = state;
    fuelsim::LocalValues minus = state;
    for (std::size_t dof = 0; dof < state.size(); ++dof) {
        plus[dof] += step * direction[dof];
        minus[dof] -= step * direction[dof];
    }
    const fuelsim::LocalResidual plus_residual = fuelsim::compute_line2_rz_boundary(current, geometry, plus);
    const fuelsim::LocalResidual minus_residual = fuelsim::compute_line2_rz_boundary(current, geometry, minus);
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
} // namespace

int main() {
    std::cout << std::scientific << std::setprecision(12);
    bool passed = true;
    passed = test_contact_search_tree() && passed;
    passed = test_mesh_and_geometry() && passed;
    passed = test_element_jacobian() && passed;
    passed = test_finite_strain_kinematics_and_jacobian() && passed;
    passed = test_cax_kinematics_and_jacobian(false) && passed;
    passed = test_cax_kinematics_and_jacobian(true) && passed;
    passed = test_gap_heat_and_normal_contact() && passed;
    passed = test_heat_point_primary_owner() && passed;
    passed = test_thermal_owner_transfer_assembly() && passed;
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

#include "fuelsim/core/rz_quad4.hpp"
