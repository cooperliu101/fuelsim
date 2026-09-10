#include "cax4rt.hpp"

#include "axisymmetric_types.hpp"
#include "core/mesh.hpp"
#include "core/spatial_definition.hpp"
#include "core/spatial_layout.hpp"
#include "core/steady_problem.hpp"
#include "material.hpp"
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
            const fuelsim::Cax4LocalResidual residual = fuelsim::rz::ProblemAccess::contribution_residual(problem,
                contribution,
                fuelsim::rz::ProblemAccess::contribution_state(problem, contribution, current));
            double magnitude = 0.0;
            for (std::size_t row = 0; row < 4; ++row)
                magnitude += std::abs(residual[row]);
            if (magnitude == 0.0)
                continue;
            const fuelsim::Cax4LocalDofs dofs = fuelsim::rz::ProblemAccess::contribution_dofs(problem, contribution);
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
        const fuelsim::Cax4LocalResidual residual = fuelsim::rz::ProblemAccess::contribution_residual(problem,
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
        const fuelsim::Cax4LocalValues local =
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
    const fuelsim::Cax4LocalDofs interface = fuelsim::rz::ProblemAccess::contribution_dofs(problem, first_thermal);
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
    const fuelsim::Cax4LocalValues initial_element_state =
        fuelsim::rz::ProblemAccess::contribution_state(problem, 0, initial_state);
    const fuelsim::Cax4LocalResidual source_residual =
        fuelsim::rz::ProblemAccess::contribution_residual(problem, 0, initial_element_state);
    problem.set_load_factor(2.0);
    const fuelsim::Cax4LocalResidual doubled_source_residual =
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

bool test_time_table() {
    const fuelsim::PiecewiseLinearTimeTable table("power", {0.0, 2.0, 5.0}, {0.0, 1.0, 0.4});
    bool passed =
        check(table.value(0.0) == 0.0 && table.value(1.0) == 0.5 && table.value(3.0) == 0.8 && table.value(8.0) == 0.4,
            "piecewise-linear table interpolates and holds endpoints");
    passed = check(table.average_value(0.0, 1.0) == 0.25 && std::abs(table.average_value(1.0, 3.0) - 0.825) < 1.0e-15
                       && std::abs(table.average_value(4.0, 8.0) - 0.425) < 1.0e-15,
                 "piecewise-linear table averages exactly across knots and held endpoints")
             && passed;
    return passed;
}

} // namespace

int main() {
    std::cout << std::scientific << std::setprecision(12);
    bool passed = true;
    passed = test_contact_search_tree() && passed;
    passed = test_mesh_and_geometry() && passed;
    passed = test_thermal_owner_transfer_assembly() && passed;
    passed = test_m1_dof_layout() && passed;
    passed = test_time_table() && passed;
    if (!passed)
        return 1;
    std::cout << "[PASS] fuelsim core geometry, DOF, and AD Jacobian tests\n";
    return 0;
}

#include "axisymmetric_types.hpp"
