#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct DisplacementReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> displacement{};
};

struct PressureReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    double pressure;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> values;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) values.push_back(value);
    return values;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete H20.24 reference row in " + path);
    return std::stod(values[index]);
}

std::vector<DisplacementReference> read_displacement(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.24 Abaqus displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.24 displacement header in " + path);
    std::vector<DisplacementReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 3, path)),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
            {number(values, 0, path), number(values, 1, path), number(values, 2, path)}});
    }
    return result;
}

std::vector<PressureReference> read_pressure(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.24 Abaqus pressure reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "pressure,id,x,y,z") throw std::invalid_argument("Unexpected H20.24 pressure header in " + path);
    std::vector<PressureReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)}, number(values, 0, path)});
    }
    return result;
}

double read_normal_reaction(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read H20.24 Abaqus reaction reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "time,reaction_x,reaction_y,reaction_z")
        throw std::invalid_argument("Unexpected H20.24 reaction header in " + path);
    std::vector<std::string> final_values;
    while (std::getline(input, line))
        if (!line.empty()) final_values = split_csv(line);
    if (final_values.empty() || std::abs(number(final_values, 0, path) - 1.0) > 1.0e-12)
        throw std::invalid_argument("H20.24 Abaqus reaction does not end at unit load");
    return std::abs(number(final_values, 1, path));
}

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = definition.solver.absolute_tolerance;
    options.relative_tolerance = definition.solver.relative_tolerance;
    options.step_tolerance = definition.solver.step_tolerance;
    options.maximum_iterations = definition.solver.maximum_iterations;
    options.linear_solver = definition.solver.linear_solver;
    options.direct_factorization = definition.solver.direct_factorization;
    options.field_residual_scaling = definition.solver.field_residual_scaling;
    options.linear_relative_tolerance = definition.solver.linear_relative_tolerance;
    options.maximum_linear_iterations = definition.solver.maximum_linear_iterations;
    return options;
}

double coordinate_difference(const fuelsim::CartesianPoint3& first, const fuelsim::CartesianPoint3& second) {
    return std::max({std::abs(first.x - second.x), std::abs(first.y - second.y), std::abs(first.z - second.z)});
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_h20_24_hex20_nonmatching_abaqus_tests <case.fsi> <displacement.csv> "
                     "<reaction.csv> <pressure.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim H20.24 Abaqus comparison\n");
        const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(argv[1]);
        if (definition.problem != fuelsim::CaseProblem::steady ||
            definition.geometry != fuelsim::CaseGeometry::cartesian_3d || definition.spatial.contacts.size() != 1 ||
            definition.spatial.contacts[0].thermal || !definition.spatial.contacts[0].mechanical ||
            definition.spatial.contacts[0].mechanical_discretization !=
                fuelsim::MechanicalContactDiscretization::surface_to_surface ||
            definition.spatial.contacts[0].friction_coefficient != 0.0)
            throw std::invalid_argument("H20.24 requires isolated frictionless HEX20 surface-to-surface contact");
        const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
        fuelsim::SteadyProblem problem(definition.spatial, mesh);
        const fuelsim::SteadyResult solve =
            fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
        bool passed =
            check(solve.completed && solve.solve.converged && solve.completed_steps == 1,
                "H20.24 Fuelsim solve completes one compression step") &&
            check(solve.aggregate_timing.workspace_setups == 1, "H20.24 Fuelsim solve constructs one PETSc workspace");

        const auto displacement = read_displacement(argv[2]);
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        const auto& fields = spatial.field_layout();
        fuelsim::test::FieldErrorMetrics displacement_x;
        std::vector<bool> present(mesh.nodes().size(), false);
        double maximum_coordinate_difference = 0.0;
        for (std::size_t region = 0; region < spatial.region_count(); ++region) {
            const auto& region_mesh = spatial.hex20_region_mesh(region);
            for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
                const std::size_t source = region_mesh.source_node_ids()[local];
                const auto found = std::find_if(displacement.begin(), displacement.end(),
                    [&](const DisplacementReference& value) { return value.id == source; });
                if (found == displacement.end() || present[source])
                    throw std::invalid_argument("H20.24 displacement source-node mapping is incomplete or repeated");
                present[source] = true;
                maximum_coordinate_difference =
                    std::max(maximum_coordinate_difference, coordinate_difference(mesh.nodes()[source], found->point));
                const std::size_t global = spatial.global_node(region, local);
                displacement_x.add(solve.solve.state[fields[1].begin + global], found->displacement[0]);
            }
        }

        const auto pressure = read_pressure(argv[4]);
        const auto contact = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
        const auto secondary_sources = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
        fuelsim::test::FieldErrorMetrics contact_pressure;
        double maximum_pressure_coordinate_difference = 0.0;
        for (std::size_t node = 0; node < contact.size(); ++node) {
            const auto found = std::find_if(pressure.begin(), pressure.end(),
                [&](const PressureReference& value) { return value.id == secondary_sources.at(node); });
            if (found == pressure.end())
                throw std::invalid_argument("H20.24 pressure source-node mapping is incomplete");
            maximum_pressure_coordinate_difference = std::max(
                maximum_pressure_coordinate_difference, coordinate_difference(mesh.nodes()[found->id], found->point));
            contact_pressure.add(contact[node].pressure, found->pressure);
        }
        const fuelsim::InterfaceSummary interface =
            fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
        const double reference_force = read_normal_reaction(argv[3]);
        const double force_error = std::abs(interface.total_contact_force - reference_force) / reference_force;
        fuelsim::test::print_relative_metrics("h20_24_displacement_x", displacement_x);
        fuelsim::test::print_relative_metrics("h20_24_contact_pressure", contact_pressure);
        std::cout << "h20_24_normal_resultant=" << interface.total_contact_force << '\n'
                  << "h20_24_abaqus_normal_resultant=" << reference_force << '\n'
                  << "h20_24_normal_resultant_relative_error=" << force_error << '\n';

        constexpr double field_tolerance = 6.0e-2;
        constexpr double resultant_tolerance = 5.0e-3;
        passed =
            check(displacement.size() == mesh.nodes().size() && displacement_x.value_count == mesh.nodes().size(),
                "H20.24 compares all eighty-eight normal displacements") &&
            check(pressure.size() == 13 && contact_pressure.value_count == 13,
                "H20.24 compares all thirteen secondary-face nodal pressures") &&
            check(maximum_coordinate_difference < 1.0e-6 && maximum_pressure_coordinate_difference < 1.0e-6,
                "H20.24 Abaqus fields use the tracked nonmatching Exodus coordinates") &&
            check(interface.active_contact_nodes == 13 && interface.unprojected_contact_nodes == 0,
                "H20.24 keeps all thirteen secondary contact nodes active and projected") &&
            check(fuelsim::test::relative_metrics_below(displacement_x, field_tolerance),
                "H20.24 normal-displacement three Abaqus errors are below 6 percent") &&
            check(fuelsim::test::relative_metrics_below(contact_pressure, field_tolerance),
                "H20.24 contact-pressure three Abaqus errors are below 6 percent") &&
            check(force_error < resultant_tolerance, "H20.24 normal resultant agrees with Abaqus below 0.5 percent") &&
            passed;
        if (passed && session.rank() == 0) std::cout << "[PASS] H20.24 nonmatching HEX20 Abaqus comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.24 Abaqus comparison raised: " << error.what() << '\n';
        return 1;
    }
}
