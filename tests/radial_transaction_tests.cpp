#include "core/problem_backend_access.hpp"
#include "io/checkpoint.hpp"
#include "solver/solve_workflows.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
using namespace fuelsim;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

bool near(double first, double second) {
    return std::abs(first - second) < 1e-12 * std::max({1.0, std::abs(first), std::abs(second)});
}

// This deliberately small internal fixture supplies admissible states to
// inspect transactions. It does not solve a physical verification case.
UnstructuredBar2Mesh make_mesh() {
    return UnstructuredBar2Mesh({{0.5, 0.5},
                                    {1.0, 0.5},
                                    {1.05, 0.5},
                                    {1.3, 0.5},
                                    {0.5, 1.5},
                                    {1.0, 1.5},
                                    {1.05, 1.5},
                                    {1.3, 1.5},
                                    {0.0, 0.0},
                                    {0.0, 1.0},
                                    {0.0, 2.0},
                                    {0.0, 0.0},
                                    {0.0, 1.0},
                                    {0.0, 2.0}},
        {{{0, 1}, {8, 9}}, {{2, 3}, {11, 12}}, {{4, 5}, {9, 10}}, {{6, 7}, {12, 13}}},
        {1, 2, 1, 2},
        {{1, "fuel"}, {2, "clad"}},
        {},
        {{1, "fuel_outer", {{0, 1}, {2, 1}}}, {2, "clad_inner", {{1, 0}, {3, 0}}}});
}

SpatialDefinition make_definition() {
    const auto registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>();
    functions->name = "radial_transaction";
    functions->thermal = registry.bind_thermal("constant_thermophysical",
        {{"conductivity", 10.0}, {"density", 1000.0}, {"specific_heat", 1000.0}});
    functions->elasticity =
        registry.bind_elasticity("constant_isotropic", {{"young_modulus", 1.0e9}, {"poisson_ratio", 0.3}});
    functions->eigenstrains.push_back(registry.bind_eigenstrain("thermal",
        "isotropic_thermal_expansion",
        {{"thermal_expansion", 1.0e-5}, {"reference_temperature", 600.0}}));
    SpatialDefinition definition;
    for (const char* name : {"fuel", "clad"}) {
        RegionDefinition region{};
        region.name = name;
        region.block = name;
        region.material = {functions, 1.0e9};
        region.initial_temperature = 600.0;
        region.volumetric_heat_source = 1000.0;
        region.heat_source_function = "source";
        region.heat_source_time_evaluation = HeatSourceTimeEvaluation::interval_average;
        region.strain_formulation = StrainFormulation::finite;
        region.radial_gps = true;
        definition.regions.push_back(region);
    }
    definition.time_tables.emplace_back("source",
        std::vector<double>{0.0, 0.1, 0.2, 0.4},
        std::vector<double>{0.0, 10.0, 0.0, 20.0});
    ContactDefinition contact{};
    contact.name = "interface";
    contact.primary = "clad_inner";
    contact.secondary = "fuel_outer";
    contact.thermal = true;
    contact.mechanical = true;
    contact.gap_conductivity = 1.0;
    contact.minimum_gap = 0.001;
    contact.penalty = 1.0e6;
    contact.friction_coefficient = 0.3;
    definition.contacts.push_back(contact);
    return definition;
}

std::string checkpoint_bytes(const std::string& path, const TransientProblem& problem) {
    write_transient_checkpoint(path, problem, 0.05);
    std::ifstream file(path, std::ios::binary);
    if (!file)
        throw std::runtime_error("Cannot read the radial transaction checkpoint");
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

double current_volume(const radial::SpatialAssembly& spatial,
    std::size_t region,
    std::size_t element,
    const std::vector<double>& state) {
    const auto& geometry = spatial.region_element_geometry(region, element);
    const auto values = spatial.volume_state(spatial.region_element_offset(region) + element, state);
    const double inner = geometry.radii[0] + values[2], outer = geometry.radii[1] + values[3];
    const double height = geometry.z_upper - geometry.z_lower + values[5] - values[4];
    return std::acos(-1.0) * (outer * outer - inner * inner) * height;
}

void run(const std::string& path) {
    const auto mesh = make_mesh();
    const auto definition = make_definition();
    TransientProblem problem(definition, mesh);
    problem.track_previous_committed_solution(true);
    const auto& spatial = BackendAccess::radial_spatial(problem);
    const auto& fields = problem.field_layout();
    require(fields.size() == 3 && fields[0].end - fields[0].begin == 8 && fields[1].end - fields[1].begin == 8
                && fields[2].end - fields[2].begin == 6,
        "Radial temperature/displacement nodes and shared axial section fields must retain different sizes");
    require(spatial.region_axial_dofs(0, 0)[1] == spatial.region_axial_dofs(0, 1)[0],
        "Adjacent axial segments of one body must share their common axial end section");
    auto trial = problem.initial_solution();
    for (const auto source : mesh.radial_source_node_ids()) {
        const auto node = spatial.radial_node(source);
        trial[spatial.dof(Field::temperature, node)] = 610.0 + static_cast<double>(source);
        if (source == 0 || source == 1 || source == 4 || source == 5)
            trial[spatial.dof(Field::radial_displacement, node)] = 0.06;
    }
    for (std::size_t i = 0; i < 3; ++i) {
        trial[spatial.axial_dof(8 + i)] = 0.01 + 0.001 * static_cast<double>(i);
        trial[spatial.axial_dof(11 + i)] = 0.0005 * static_cast<double>(i);
    }
    problem.begin_time_step({0.1, 1.0, true});
    problem.commit_time_step(trial);
    const auto baseline = BackendAccess::committed_state(problem);
    const auto snapshot = problem.capture_state();
    const auto bytes = checkpoint_bytes(path, problem);
    require(baseline.radial_material_histories.size() == 2 && baseline.radial_material_histories[0].size() == 2
                && baseline.radial_material_histories[0][0].size() == 2,
        "Every radial element must commit exactly two material points");
    require(baseline.material_histories.empty() && baseline.cartesian_material_histories.empty()
                && baseline.quad8_material_histories.empty(),
        "Radial history must not be stored in an unrelated topology's material container");
    require(baseline.contact_histories.size() == 1 && baseline.contact_histories[0].size() == 4,
        "Two axial contact segments must retain their four independent quadrature histories");
    for (const auto& history : baseline.contact_histories[0])
        require(history.sliding && history.elastic_tangential_slip > 0.0 && history.total_tangential_slip > 0.0,
            "The committed contact fixture must exercise nonzero sliding and elastic tangential history");
    require(problem.last_conservation_summary().friction_dissipation_increment > 0.0,
        "Accepted radial contact must report its frictional dissipation");
    require(problem.has_previous_committed_solution() && problem.previous_committed_time() == 0.0,
        "The previous complete solution must remain available for history-dependent result output");

    TransientTimeOptions options{};
    options.time_error_relative_tolerance = 0.0;
    options.displacement_time_absolute_tolerance = 1e-8;
    options.strain_history_time_absolute_tolerance = 1e-8;
    options.stress_history_time_absolute_tolerance = 1.0;
    require(problem.step_doubling_error(snapshot, snapshot, options).maximum == 0.0,
        "Identical complete snapshots must have zero time error");
    for (std::size_t field = 0; field < fields.size(); ++field) {
        auto changed = baseline;
        changed.solution[fields[field].end - 1] += field == 0 ? 1.0 : 1e-4;
        BackendAccess::restore_committed_state(problem, std::move(changed));
        const auto estimate = problem.step_doubling_error(snapshot, problem.capture_state(), options);
        require(estimate.nodal_fields.size() == 3 && estimate.nodal_fields[field].value > 1.0,
            "Time error must inspect the last active degree of freedom in every independently sized field");
        problem.restore_state(snapshot);
    }
    for (std::size_t component = 0; component < 6; ++component) {
        auto changed = baseline;
        auto& point = changed.radial_material_histories[1][1][1];
        if (component == 0)
            point.elastic_strain[0] += 1e-4;
        if (component == 1)
            point.plastic_strain[0] += 1e-4;
        if (component == 2)
            point.creep_strain[0] += 1e-4;
        if (component == 3)
            point.equivalent_plastic_strain += 1e-4;
        if (component == 4)
            point.equivalent_creep_strain += 1e-4;
        if (component == 5)
            point.stress.zz += 100.0;
        BackendAccess::restore_committed_state(problem, std::move(changed));
        const auto estimate = problem.step_doubling_error(snapshot, problem.capture_state(), options);
        require(estimate.maximum > 1.0,
            "Time error must inspect all material histories and stress at the second point of the last radial element");
        problem.restore_state(snapshot);
    }
    auto changed_contact = baseline;
    auto& contact_history = changed_contact.contact_histories[0].back();
    contact_history.elastic_tangential_slip += 1e-4;
    contact_history.cartesian_elastic_tangential_slip[2] = contact_history.elastic_tangential_slip;
    BackendAccess::restore_committed_state(problem, std::move(changed_contact));
    require(problem.step_doubling_error(snapshot, problem.capture_state(), options).contact_friction > 1.0,
        "Time error must include the last contact quadrature point's force-carrying elastic slip");
    problem.restore_state(snapshot);
    require(checkpoint_bytes(path, problem) == bytes,
        "Snapshot restoration must preserve the complete serialized state");

    problem.begin_time_step({0.3, 1.0, true});
    const auto averaged =
        spatial
            .evaluate_volume(0, 0, trial, baseline.solution, &baseline.radial_material_histories[0][0], 0.2, 0.3, true);
    const double volume = current_volume(spatial, 0, 0, trial);
    require(near(averaged.generated_heat_rate / volume, 5000.0),
        "A step spanning a source-table corner must use its exact interval average");
    auto invalid = trial;
    invalid[spatial.axial_dof(9)] = -2.0;
    bool rejected = false;
    try {
        problem.commit_time_step(invalid);
    } catch (const std::domain_error&) {
        rejected = true;
    }
    require(rejected, "A finite trial with reversed axial geometry must fail before committing state");
    problem.rollback_time_step();
    require(!problem.time_step_active() && problem.committed_time() == 0.1 && checkpoint_bytes(path, problem) == bytes,
        "A failed finite step must restore nodes, both material points, contact history, time and conservation");
    const auto restored_source = spatial.evaluate_volume(0,
        0,
        trial,
        baseline.solution,
        &baseline.radial_material_histories[0][0],
        0.2,
        0.1,
        false);
    require(near(restored_source.generated_heat_rate / volume, 10000.0),
        "Rollback must restore the heat source at the committed physical time");
    problem.begin_time_step({0.3, 1.0, true});
    const auto repeated =
        spatial
            .evaluate_volume(0, 0, trial, baseline.solution, &baseline.radial_material_histories[0][0], 0.2, 0.3, true);
    require(repeated.residual == averaged.residual && repeated.history[1].stress.zz == averaged.history[1].stress.zz,
        "Repeating the same step after failure must preserve exact local residual and material repeatability");
    problem.rollback_time_step();

    auto invalid_restore = baseline;
    invalid_restore.radial_material_histories[1].pop_back();
    rejected = false;
    try {
        BackendAccess::restore_committed_state(problem, std::move(invalid_restore));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected && checkpoint_bytes(path, problem) == bytes,
        "Malformed material layout restoration must fail without partially replacing committed state");
    invalid_restore = baseline;
    invalid_restore.contact_histories[0].back().total_tangential_slip = std::numeric_limits<double>::quiet_NaN();
    rejected = false;
    try {
        BackendAccess::restore_committed_state(problem, std::move(invalid_restore));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected && checkpoint_bytes(path, problem) == bytes,
        "Malformed contact restoration must fail without partially replacing valid earlier contact points");

    TransientProblem restored(definition, mesh);
    require(restore_transient_checkpoint(path, restored) == 0.05
                && checkpoint_bytes(path + ".roundtrip", restored) == bytes,
        "Checkpoint roundtrip must preserve unequal fields, both material points, contact and previous states exactly");
    std::remove(path.c_str());
    std::remove((path + ".roundtrip").c_str());
    std::cout << "radial_gps_transaction_contracts=passed\n";
}
} // namespace

int main(int argc, char** argv) {
    try {
        run(argc > 1 ? argv[1] : "radial_transaction.checkpoint");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
