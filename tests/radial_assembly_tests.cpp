#include "core/radial_assembly.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
using namespace fuelsim;
constexpr double pi = 3.141592653589793238462643383279502884;

void require(bool value, const char* message) {
    if (!value)
        throw std::runtime_error(message);
}

void close(double actual, double expected, const char* message, double tolerance = 2.0e-12) {
    require(std::abs(actual - expected) <= tolerance * std::max({1.0, std::abs(actual), std::abs(expected)}), message);
}

void reject(const std::function<void()>& operation, const char* message) {
    try {
        operation();
    } catch (const std::invalid_argument&) {
        return;
    } catch (const std::domain_error&) {
        return;
    }
    throw std::runtime_error(message);
}

ThermoelasticProperties properties() {
    const auto registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>();
    functions->name = "radial_assembly_test";
    functions->thermal = registry.bind_thermal("constant_thermophysical",
        {{"conductivity", 5.0}, {"density", 1000.0}, {"specific_heat", 500.0}});
    functions->elasticity = registry.bind_elasticity("linear_temperature_isotropic",
        {{"young_modulus", 2.0e11},
            {"poisson_ratio", 0.3},
            {"reference_temperature", 600.0},
            {"young_modulus_temperature_coefficient", 0.0},
            {"poisson_ratio_temperature_coefficient", 0.0}});
    functions->eigenstrains.push_back(registry.bind_eigenstrain("thermal",
        "linear_temperature_isotropic_thermal_expansion",
        {{"thermal_expansion", 1.0e-5},
            {"reference_temperature", 600.0},
            {"thermal_expansion_temperature_coefficient", 0.0}}));
    return {functions, 2.0e11};
}

UnstructuredBar2Mesh layered_mesh(bool matching = true) {
    const double clad_mid = matching ? 1.0 : 1.1;
    return {{{0.0, 0.5},
                {0.004, 0.5},
                {0.0045, 0.5 * clad_mid},
                {0.005, 0.5 * clad_mid},
                {0.0, 1.5},
                {0.004, 1.5},
                {0.0045, 0.5 * (clad_mid + 2.0)},
                {0.005, 0.5 * (clad_mid + 2.0)},
                {0.0, 0.0},
                {0.0, 1.0},
                {0.0, 2.0},
                {0.0, 0.0},
                {0.0, clad_mid},
                {0.0, 2.0},
                {0.002, 0.5}},
        {{{0, 14}, {8, 9}}, {{14, 1}, {8, 9}}, {{2, 3}, {11, 12}}, {{4, 5}, {9, 10}}, {{6, 7}, {12, 13}}},
        {1, 1, 2, 1, 2},
        {{1, "fuel"}, {2, "clad"}},
        {{1, "fuel_bottom", {8}}, {2, "clad_bottom", {11}}, {3, "fuel_top", {10}}},
        {{1, "fuel_outer", {{1, 1}, {3, 1}}},
            {2, "clad_inner", {{2, 0}, {4, 0}}},
            {3, "clad_outer", {{2, 1}, {4, 1}}}}};
}

SpatialDefinition layered_definition() {
    SpatialDefinition definition;
    definition.regions.push_back({"fuel", "fuel", properties(), 1.0e7, 600.0});
    definition.regions.push_back({"clad", "clad", properties(), 0.0, 600.0});
    for (auto& region : definition.regions)
        region.radial_gps = true;
    definition.contacts.push_back({"gap", "clad_inner", "fuel_outer", true, true, 0.1, 1.0e-6, 1.0e12, 0.3, 1.0e-5});
    definition.boundary_conditions.push_back(
        {"fuel_anchor", BoundaryConditionType::dirichlet, "fuel_bottom", Field::axial_displacement, 0.0});
    definition.boundary_conditions.push_back(
        {"clad_anchor", BoundaryConditionType::dirichlet, "clad_bottom", Field::axial_displacement, 0.0});
    definition.boundary_conditions.push_back(
        {"end_force", BoundaryConditionType::axial_force, "fuel_top", Field::axial_displacement, 17.5});
    return definition;
}

std::vector<double>
local_state(const radial::SpatialAssembly& spatial, std::size_t index, const std::vector<double>& global) {
    std::vector<std::size_t> dofs;
    spatial.contribution_dofs(index, dofs);
    std::vector<double> state;
    for (const auto dof : dofs)
        state.push_back(global[dof]);
    return state;
}

bool same_histories(const std::vector<std::vector<ContactPointHistory>>& a,
    const std::vector<std::vector<ContactPointHistory>>& b) {
    if (a.size() != b.size())
        return false;
    for (std::size_t c = 0; c < a.size(); ++c) {
        if (a[c].size() != b[c].size())
            return false;
        for (std::size_t q = 0; q < a[c].size(); ++q) {
            const auto& x = a[c][q];
            const auto& y = b[c][q];
            if (x.elastic_tangential_slip != y.elastic_tangential_slip
                || x.total_tangential_slip != y.total_tangential_slip || x.sliding != y.sliding
                || x.normal_multiplier != y.normal_multiplier
                || x.cartesian_elastic_tangential_slip != y.cartesian_elastic_tangential_slip
                || x.cartesian_total_tangential_slip != y.cartesian_total_tangential_slip
                || x.cartesian_tangent_basis_initialized != y.cartesian_tangent_basis_initialized
                || x.cartesian_contact_normal != y.cartesian_contact_normal
                || x.cartesian_contact_tangent_first != y.cartesian_contact_tangent_first)
                return false;
        }
    }
    return true;
}

void test_layout_and_controls() {
    auto definition = layered_definition();
    definition.time_tables.emplace_back("source", std::vector<double>{0.0, 1.0}, std::vector<double>{0.0, 2.0});
    definition.regions[0].heat_source_function = "source";
    definition.regions[0].heat_source_time_evaluation = HeatSourceTimeEvaluation::interval_average;
    radial::SpatialAssembly spatial(definition, layered_mesh());
    require(spatial.node_count() == 9 && spatial.temperature_node_count() == 9 && spatial.dof_count() == 24,
        "Radial and axial field cardinalities must be independent");
    require(spatial.field_layout()[2].begin == 18 && spatial.field_layout()[2].end == 24,
        "Axial controls must follow the two radial nodal fields");
    require(spatial.region_axial_dofs(0, 0) == spatial.region_axial_dofs(0, 1),
        "Radial elements in one layer must share both axial controls");
    require(spatial.region_axial_dofs(0, 1)[1] == spatial.region_axial_dofs(0, 2)[0],
        "Adjacent fuel layers must share one axial control");
    require(spatial.axial_dof(8) != spatial.axial_dof(11),
        "Coincident fuel and cladding axial controls must remain independent");
    require(spatial.dirichlet_conditions().size() == 4,
        "Two axis radial constraints must accompany the two prescribed axial anchors");
    const auto first = spatial.contribution_partition(0, 2), second = spatial.contribution_partition(1, 2);
    require(first.first == 0 && first.second == second.first && second.second == spatial.contribution_count(),
        "Two partitions must cover every contribution exactly once");
    require(spatial.contact_source_elements(0) == std::vector<std::pair<std::size_t, std::size_t>>{{2, 1}, {4, 3}},
        "Contact pair mapping must preserve source element numbers and both axial layers");
    spatial.set_time(1.0);
    spatial.set_heat_source_interval(0.0, 1.0);
    const auto state = spatial.initial_state();
    double generated = 0.0;
    for (std::size_t e = 0; e < spatial.region_element_count(0); ++e)
        generated += spatial.evaluate_volume(0, e, state, state, nullptr, 0.0, 1.0, false).generated_heat_rate;
    close(generated,
        pi * 0.004 * 0.004 * 2.0 * 1.0e7,
        "One-time body assembly must preserve the full cylinder heat source and interval average");
    spatial.set_time(0.0);
    close(spatial.evaluate_volume(0, 0, state, state, nullptr, 0.0, 0.0, false).generated_heat_rate,
        0.0,
        "Restoring time must restore the instantaneous source control");
    std::vector<double> residual;
    spatial.compute_boundary(spatial.contribution_count() - 1,
        local_state(spatial, spatial.contribution_count() - 1, state),
        residual,
        nullptr);
    require(residual == std::vector<double>{-17.5},
        "Axial force must be applied once in newtons without an area factor");
}

void test_contact_history_transaction() {
    radial::SpatialAssembly spatial(layered_definition(), layered_mesh());
    const auto old = spatial.initial_state();
    const auto initial_histories = spatial.committed_contact_histories();
    auto state = old;
    for (const auto source : spatial.region_mesh(0).source_node_ids()) {
        const auto node = spatial.radial_node(source);
        state[spatial.dof(Field::temperature, node)] = 700.0;
        state[spatial.dof(Field::radial_displacement, node)] = 0.15 * spatial.source_mesh().nodes()[source].r;
    }
    state[spatial.axial_dof(9)] = 1.0e-4;
    state[spatial.axial_dof(10)] = 2.0e-4;
    spatial.validate_state(state);
    for (std::size_t pair = 0; pair < 2; ++pair) {
        const auto value = spatial.evaluate_contact(0, pair, state);
        close(value.normal_force,
            1.0e8 * 2.0 * pi * 0.004,
            "Contact normal force must include exactly one cylindrical layer area");
        require(value.points[0].sliding && value.points[1].sliding,
            "The two axial contact points must independently reach sliding");
        require(value.points[0].total_tangential_slip != value.points[1].total_tangential_slip,
            "Contact output must preserve the two axial slip values");
        const auto index = spatial.volume_contribution_count() + pair;
        const auto local = local_state(spatial, index, state);
        radial::LocalContribution first, repeated;
        spatial.compute_contribution(index, local, old, nullptr, 0.0, 0.0, false, true, first);
        spatial.compute_contribution(index, local, old, nullptr, 0.0, 0.0, false, false, repeated);
        require(first.residual == repeated.residual,
            "Contact residual calls must repeat exactly without committing history");
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(index, dofs);
        double heat = 0.0, radial_force = 0.0, axial_force = 0.0;
        for (std::size_t i = 0; i < dofs.size(); ++i) {
            if (dofs[i] < spatial.temperature_node_count())
                heat += first.residual[i];
            else if (dofs[i] < 2 * spatial.node_count())
                radial_force += first.residual[i];
            else
                axial_force += first.residual[i];
        }
        close(heat, 0.0, "Thermal contact must conserve heat across its two sides");
        close(radial_force, 0.0, "Mechanical contact must conserve radial force across its two sides");
        close(axial_force, 0.0, "Friction must conserve axial force across both axial chains", 1.0e-8);
    }
    require(same_histories(initial_histories, spatial.committed_contact_histories()),
        "Evaluating contact residuals and output must not change committed history");
    require(spatial.commit_contact_state(state) > 0.0, "Accepting sliding must return positive friction dissipation");
    const auto accepted = spatial.committed_contact_histories();
    require(!same_histories(initial_histories, accepted), "Accepted sliding must commit both quadrature histories");
    auto corrupt = accepted;
    corrupt[0][0].elastic_tangential_slip = std::numeric_limits<double>::quiet_NaN();
    reject([&]() { spatial.restore_contact_state(old, corrupt); }, "Invalid restored contact history must be rejected");
    require(same_histories(accepted, spatial.committed_contact_histories()),
        "Failed history restoration must preserve the previous accepted state");
    spatial.restore_contact_state(old, initial_histories);
    require(same_histories(initial_histories, spatial.committed_contact_histories()),
        "Restoring a complete state must recover contact history exactly");
    close(spatial.summarize_interface(0, old).total_contact_force,
        0.0,
        "Restoring the initial open contact must recover its zero force");
}

void test_current_boundary_jacobians() {
    const UnstructuredBar2Mesh mesh({{1.0, 1.0}, {2.0, 1.0}, {0.0, 0.0}, {0.0, 2.0}},
        {{{0, 1}, {2, 3}}},
        {1},
        {{1, "body"}},
        {{1, "lower", {2}}, {2, "upper", {3}}},
        {{1, "inner", {{0, 0}}}, {2, "outer", {{0, 1}}}});
    SpatialDefinition definition;
    definition.regions.push_back({"body", "body", properties(), 0.0, 600.0});
    definition.regions[0].radial_gps = true;
    definition.regions[0].strain_formulation = StrainFormulation::finite;
    definition.boundary_conditions = {
        {"pressure", BoundaryConditionType::pressure, "outer", Field::radial_displacement, 1.0e6},
        {"radial_traction", BoundaryConditionType::traction, "inner", Field::radial_displacement, 2.0e6},
        {"axial_traction", BoundaryConditionType::traction, "outer", Field::axial_displacement, 3.0e6},
        {"flux", BoundaryConditionType::heat_flux, "outer", Field::temperature, 1234.0},
        {"cooling", BoundaryConditionType::convection, "outer", Field::temperature, 0.0}};
    definition.boundary_conditions.back().heat_transfer_coefficient = 25.0;
    definition.boundary_conditions.back().ambient_temperature = 300.0;
    radial::SpatialAssembly spatial(definition, mesh);
    auto state = spatial.initial_state();
    state[0] = 620.0;
    state[1] = 630.0;
    state[2] = 0.05;
    state[3] = 0.1;
    state[4] = 0.2;
    state[5] = 0.5;
    spatial.validate_state(state);
    const double outer_area = 2.0 * pi * 2.1 * 2.3, inner_area = 2.0 * pi * 1.05 * 2.3;
    for (std::size_t b = 0; b < 5; ++b) {
        const auto index = spatial.volume_contribution_count() + b;
        const auto local = local_state(spatial, index, state);
        std::vector<double> residual, jacobian;
        spatial.compute_boundary(index, local, residual, &jacobian);
        if (b == 0)
            close(residual[1], 1.0e6 * outer_area, "Finite-strain pressure must use current cylindrical area");
        else if (b == 1)
            close(residual[1], -2.0e6 * inner_area, "Radial traction must retain its global component direction");
        else if (b == 2)
            close(residual[2] + residual[3],
                -3.0e6 * outer_area,
                "Axial cylindrical traction must split its complete force between axial controls");
        else if (b == 3)
            close(residual[0], -1234.0 * outer_area, "Finite-strain heat flux must use current cylindrical area");
        else
            close(residual[0],
                25.0 * (630.0 - 300.0) * outer_area,
                "Finite-strain convection must use current cylindrical area");
        for (std::size_t column = 0; column < local.size(); ++column) {
            const double step = column == 0 ? 1.0e-3 : 1.0e-6;
            auto plus = local, minus = local;
            plus[column] += step;
            minus[column] -= step;
            std::vector<double> upper, lower;
            spatial.compute_boundary(index, plus, upper, nullptr);
            spatial.compute_boundary(index, minus, lower, nullptr);
            for (std::size_t row = 0; row < local.size(); ++row)
                close(jacobian[row * local.size() + column],
                    (upper[row] - lower[row]) / (2.0 * step),
                    "Current cylindrical boundary tangent must include every geometric derivative",
                    5.0e-8);
        }
    }
    auto bad = state;
    bad[5] = -2.0;
    reject([&]() { spatial.validate_state(bad); }, "A nonpositive current layer height must reject the trial state");

    definition.boundary_conditions[0].configuration_explicit = true;
    definition.boundary_conditions[0].use_displaced_geometry = false;
    radial::SpatialAssembly reference(definition, mesh);
    std::vector<double> residual, jacobian;
    reference.compute_boundary(1, local_state(reference, 1, state), residual, &jacobian);
    close(residual[1], 1.0e6 * 2.0 * pi * 2.0 * 2.0, "Explicit reference pressure must honor the reference area");
    for (const auto value : jacobian)
        require(value == 0.0, "Reference cylindrical pressure has no geometric tangent");
}

void test_dimensionless_friction_slip_tolerance() {
    const auto original = layered_mesh();
    auto nodes = original.nodes();
    for (const auto source : {0U, 1U, 2U, 3U, 14U})
        nodes[source].z = 0.005;
    for (const auto source : {4U, 5U, 6U, 7U})
        nodes[source].z = 0.02;
    nodes[9].z = nodes[12].z = 0.01;
    nodes[10].z = nodes[13].z = 0.03;
    const UnstructuredBar2Mesh mesh(nodes,
        original.elements(),
        original.element_block_ids(),
        original.element_blocks(),
        original.node_sets(),
        original.side_sets());
    auto definition = layered_definition();
    definition.contacts[0].friction_slip_tolerance = 0.1;
    radial::SpatialAssembly spatial(definition, mesh);
    auto state = spatial.initial_state();
    for (const auto source : spatial.region_mesh(0).source_node_ids())
        state[spatial.dof(Field::radial_displacement, spatial.radial_node(source))] = 0.15 * nodes[source].r;
    for (const auto source : {8U, 9U, 10U})
        state[spatial.axial_dof(source)] = 0.005;
    for (std::size_t pair = 0; pair < 2; ++pair)
        for (const auto& point : spatial.evaluate_contact(0, pair, state).points) {
            require(point.sliding, "The prescribed relative motion must reach the Coulomb sliding branch");
            close(std::abs(point.elastic_tangential_slip),
                0.1 * (0.01 + 0.02) * 0.5,
                "Dimensionless slip tolerance must multiply mean primary reference segment length");
        }
}

void test_shared_axial_contact_mapping() {
    const auto original = layered_mesh();
    auto elements = original.elements();
    elements[2].axial_nodes = {{8, 9}};
    elements[4].axial_nodes = {{9, 10}};
    auto sets = original.node_sets();
    sets[1].nodes = {8};
    const UnstructuredBar2Mesh mesh(original.nodes(),
        elements,
        original.element_block_ids(),
        original.element_blocks(),
        sets,
        original.side_sets());
    auto definition = layered_definition();
    definition.boundary_conditions.erase(definition.boundary_conditions.begin() + 1);
    radial::SpatialAssembly spatial(definition, mesh);
    auto state = spatial.initial_state();
    for (const auto source : spatial.region_mesh(0).source_node_ids())
        state[spatial.dof(Field::radial_displacement, spatial.radial_node(source))] = 0.15 * mesh.nodes()[source].r;
    const auto index = spatial.volume_contribution_count();
    std::vector<std::size_t> dofs;
    spatial.contribution_dofs(index, dofs);
    require(dofs.size() == 6 && std::adjacent_find(dofs.begin(), dofs.end()) == dofs.end(),
        "Contact must compress repeated shared axial controls before global assembly");
    radial::LocalContribution contribution;
    spatial.compute_contribution(index,
        local_state(spatial, index, state),
        spatial.initial_state(),
        nullptr,
        0.0,
        0.0,
        false,
        true,
        contribution);
    for (std::size_t i = 0; i < dofs.size(); ++i)
        if (dofs[i] >= 2 * spatial.node_count()) {
            close(contribution.residual[i], 0.0, "Shared axial motion must cancel both sides of friction residual");
            for (std::size_t j = 0; j < dofs.size(); ++j)
                close(contribution.jacobian[i * dofs.size() + j],
                    0.0,
                    "Shared axial motion must cancel both sides of friction tangent");
        }
}

void test_unsupported_inputs() {
    auto definition = layered_definition();
    definition.contacts[0].automatic_penalty = true;
    reject([&]() { radial::SpatialAssembly spatial(definition, layered_mesh()); },
        "Automatic penalty must be rejected");
    definition = layered_definition();
    definition.contacts[0].mechanical_sliding = MechanicalContactSliding::finite;
    reject([&]() { radial::SpatialAssembly spatial(definition, layered_mesh()); },
        "Cross-segment finite sliding must be rejected");
    definition = layered_definition();
    definition.contacts[0].mechanical_formulation = MechanicalContactFormulation::augmented_lagrangian;
    reject([&]() { radial::SpatialAssembly spatial(definition, layered_mesh()); },
        "Augmented contact must be rejected");
    definition = layered_definition();
    definition.contacts[0].mechanical_discretization = MechanicalContactDiscretization::node_to_surface;
    reject([&]() { radial::SpatialAssembly spatial(definition, layered_mesh()); },
        "Node-to-surface contact must be rejected");
    definition = layered_definition();
    reject([&]() { radial::SpatialAssembly spatial(definition, layered_mesh(false)); },
        "Unmatched axial segments must be rejected");
    definition.regions[0].strain_formulation = StrainFormulation::finite;
    reject([&]() { radial::SpatialAssembly spatial(definition, layered_mesh()); },
        "Mixed contact strain formulations must be rejected");
    for (const bool thermal : {false, true}) {
        definition = layered_definition();
        definition.contacts[0].thermal = thermal;
        definition.contacts[0].mechanical = !thermal;
        radial::SpatialAssembly spatial(definition, layered_mesh());
        require(spatial.contribution_count() == 8,
            "Thermal-only and mechanical-only contact must both retain two axial pairs");
    }
}

void test_contact_boundary_aliases() {
    const auto original = layered_mesh();
    auto sides = original.side_sets();
    sides.push_back({4, "fuel_outer_alias", sides[0].sides});
    sides.push_back({5, "clad_inner_alias", sides[1].sides});
    const UnstructuredBar2Mesh mesh(original.nodes(),
        original.elements(),
        original.element_block_ids(),
        original.element_blocks(),
        original.node_sets(),
        sides);
    for (const bool thermal : {false, true}) {
        auto definition = layered_definition();
        definition.contacts[0].thermal = thermal;
        definition.contacts[0].mechanical = !thermal;
        auto alias = definition.contacts[0];
        alias.name = "gap_alias";
        alias.primary = "clad_inner_alias";
        alias.secondary = "fuel_outer_alias";
        definition.contacts.push_back(alias);
        reject([&]() { radial::SpatialAssembly spatial(definition, mesh); },
            "Different boundary names must not duplicate a physical secondary contact constraint");
    }
    auto definition = layered_definition();
    definition.contacts[0].mechanical = false;
    auto mechanical = definition.contacts[0];
    mechanical.name = "mechanical_gap";
    mechanical.primary = "clad_inner_alias";
    mechanical.secondary = "fuel_outer_alias";
    mechanical.thermal = false;
    mechanical.mechanical = true;
    definition.contacts.push_back(mechanical);
    radial::SpatialAssembly spatial(definition, mesh);
    require(spatial.contribution_count() == 10,
        "A thermal-only and a mechanical-only contact may intentionally share the same physical interface");
}
} // namespace

int main() {
    try {
        test_layout_and_controls();
        test_contact_history_transaction();
        test_current_boundary_jacobians();
        test_dimensionless_friction_slip_tolerance();
        test_shared_axial_contact_mapping();
        test_unsupported_inputs();
        test_contact_boundary_aliases();
        std::cout << "Radial assembly checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Radial assembly checks failed: " << error.what() << '\n';
        return 1;
    }
}
