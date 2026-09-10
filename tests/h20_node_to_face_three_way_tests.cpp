#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/field_error_metrics.hpp"
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

struct MooseContactReference final {
    std::size_t id;
    fuelsim::CartesianPoint3 point;
    double pressure;
    double nodal_area;
};

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> values;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        values.push_back(value);
    return values;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size())
        throw std::invalid_argument("Incomplete H20.25 reference row in " + path);
    return std::stod(values[index]);
}

std::vector<DisplacementReference> read_displacement(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.25 displacement reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "disp_x,disp_y,disp_z,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.25 displacement header in " + path);
    std::vector<DisplacementReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 3, path)),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
            {number(values, 0, path), number(values, 1, path), number(values, 2, path)}});
    }
    return result;
}

std::vector<PressureReference> read_abaqus_pressure(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.25 Abaqus pressure reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "pressure,id,x,y,z")
        throw std::invalid_argument("Unexpected H20.25 pressure header in " + path);
    std::vector<PressureReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            number(values, 0, path)});
    }
    return result;
}

std::vector<MooseContactReference> read_moose_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.25 MOOSE contact reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "contact_pressure,id,nodal_area,penetration,x,y,z")
        throw std::invalid_argument("Unexpected H20.25 MOOSE contact header in " + path);
    std::vector<MooseContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const auto values = split_csv(line);
        result.push_back({static_cast<std::size_t>(number(values, 1, path)),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)},
            number(values, 0, path),
            number(values, 2, path)});
    }
    return result;
}

double read_normal_reaction(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read H20.25 reaction reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "time,reaction_x,reaction_y,reaction_z")
        throw std::invalid_argument("Unexpected H20.25 reaction header in " + path);
    std::vector<std::string> final_values;
    while (std::getline(input, line))
        if (!line.empty())
            final_values = split_csv(line);
    if (final_values.empty() || std::abs(number(final_values, 0, path) - 1.0) > 1.0e-12)
        throw std::invalid_argument("H20.25 reaction does not end at unit load");
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

double relative_error(double first, double second) {
    if (!(second > 0.0))
        throw std::invalid_argument("H20.25 relative-error reference must be positive");
    return std::abs(first - second) / second;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "Usage: fuelsim_h20_25_hex20_node_to_face_three_way_tests <case.fsi> "
                     "<moose_displacement.csv> <moose_contact.csv> <moose_reaction.csv> "
                     "<abaqus_displacement.csv> <abaqus_pressure.csv> <abaqus_reaction.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim H20.25 three-way node-to-face comparison\n");
        fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(argv[1]);
        if (definition.problem != fuelsim::CaseProblem::steady
            || definition.geometry != fuelsim::CaseGeometry::cartesian_3d || definition.spatial.contacts.size() != 1
            || definition.spatial.contacts[0].thermal || !definition.spatial.contacts[0].mechanical
            || definition.spatial.contacts[0].friction_coefficient != 0.0)
            throw std::invalid_argument("H20.25 requires isolated frictionless HEX20 mechanical contact");
        definition.spatial.contacts[0].mechanical_discretization =
            fuelsim::MechanicalContactDiscretization::node_to_surface;
        definition.spatial.contacts[0].quad8_nodal_area_rule = fuelsim::Quad8NodalAreaRule::consistent_shape;

        const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
        fuelsim::SteadyProblem problem(definition.spatial, mesh);
        const fuelsim::SteadyResult solve =
            fuelsim::solve_steady(problem, definition.steady_execution, solver_options(definition));
        bool passed = check(solve.completed && solve.solve.converged && solve.completed_steps == 1,
                          "H20.25 Fuelsim node-to-face solve completes one compression step")
                      && check(solve.aggregate_timing.workspace_setups == 1,
                          "H20.25 Fuelsim solve constructs one PETSc workspace");

        const auto moose_displacement = read_displacement(argv[2]);
        const auto moose_contact = read_moose_contact(argv[3]);
        const double moose_force = read_normal_reaction(argv[4]);
        const auto abaqus_displacement = read_displacement(argv[5]);
        const auto abaqus_pressure = read_abaqus_pressure(argv[6]);
        const double abaqus_force = read_normal_reaction(argv[7]);
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        const auto& fields = spatial.field_layout();
        std::array<fuelsim::test::FieldErrorMetrics, 3> fuelsim_moose_displacement;
        std::array<fuelsim::test::FieldErrorMetrics, 3> fuelsim_abaqus_displacement;
        std::array<fuelsim::test::FieldErrorMetrics, 3> moose_abaqus_displacement;
        std::vector<bool> present(mesh.nodes().size(), false);
        double maximum_moose_coordinate_difference = 0.0;
        double maximum_abaqus_coordinate_difference = 0.0;
        for (std::size_t region = 0; region < spatial.region_count(); ++region) {
            const auto& region_mesh = spatial.hex20_region_mesh(region);
            for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
                const std::size_t source = region_mesh.source_node_ids()[local];
                const auto moose = std::find_if(moose_displacement.begin(),
                    moose_displacement.end(),
                    [&](const DisplacementReference& value) { return value.id == source; });
                const auto abaqus = std::find_if(abaqus_displacement.begin(),
                    abaqus_displacement.end(),
                    [&](const DisplacementReference& value) { return value.id == source; });
                if (moose == moose_displacement.end() || abaqus == abaqus_displacement.end() || present[source])
                    throw std::invalid_argument("H20.25 displacement source-node mapping is incomplete or repeated");
                present[source] = true;
                maximum_moose_coordinate_difference = std::max(maximum_moose_coordinate_difference,
                    coordinate_difference(mesh.nodes().at(source), moose->point));
                maximum_abaqus_coordinate_difference = std::max(maximum_abaqus_coordinate_difference,
                    coordinate_difference(mesh.nodes().at(source), abaqus->point));
                const std::size_t global = spatial.global_node(region, local);
                for (std::size_t component = 0; component < 3; ++component) {
                    const double actual = solve.solve.state.at(fields[component + 1].begin + global);
                    fuelsim_moose_displacement[component].add(actual, moose->displacement[component]);
                    fuelsim_abaqus_displacement[component].add(actual, abaqus->displacement[component]);
                    moose_abaqus_displacement[component].add(moose->displacement[component],
                        abaqus->displacement[component]);
                }
            }
        }

        const auto contact = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, solve.solve.state);
        const auto source_nodes = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
        fuelsim::test::FieldErrorMetrics fuelsim_moose_pressure;
        fuelsim::test::FieldErrorMetrics fuelsim_abaqus_pressure;
        fuelsim::test::FieldErrorMetrics moose_abaqus_pressure;
        std::size_t active = 0;
        for (std::size_t index = 0; index < contact.size(); ++index) {
            const std::size_t source = source_nodes.at(index);
            const auto moose = std::find_if(moose_contact.begin(),
                moose_contact.end(),
                [&](const MooseContactReference& value) { return value.id == source; });
            const auto abaqus = std::find_if(abaqus_pressure.begin(),
                abaqus_pressure.end(),
                [&](const PressureReference& value) { return value.id == source; });
            if (moose == moose_contact.end() || abaqus == abaqus_pressure.end())
                throw std::invalid_argument("H20.25 contact source-node mapping is incomplete");
            maximum_moose_coordinate_difference = std::max(maximum_moose_coordinate_difference,
                coordinate_difference({contact[index].x, contact[index].y, contact[index].z}, moose->point));
            maximum_abaqus_coordinate_difference = std::max(maximum_abaqus_coordinate_difference,
                coordinate_difference({contact[index].x, contact[index].y, contact[index].z}, abaqus->point));
            fuelsim_moose_pressure.add(contact[index].pressure, moose->pressure);
            fuelsim_abaqus_pressure.add(contact[index].pressure, abaqus->pressure);
            moose_abaqus_pressure.add(moose->pressure, abaqus->pressure);
            if (contact[index].pressure > 0.0)
                ++active;
            std::cout << "h20_25_contact_node=" << source << ",fuelsim_pressure=" << contact[index].pressure
                      << ",moose_pressure=" << moose->pressure << ",abaqus_pressure=" << abaqus->pressure
                      << ",fuelsim_area=" << contact[index].tributary_area << '\n';
        }
        const auto extra_abaqus =
            std::find_if(abaqus_pressure.begin(), abaqus_pressure.end(), [&](const PressureReference& value) {
                return std::find(source_nodes.begin(), source_nodes.end(), value.id) == source_nodes.end();
            });
        const std::size_t extra_abaqus_count = static_cast<std::size_t>(
            std::count_if(abaqus_pressure.begin(), abaqus_pressure.end(), [&](const PressureReference& value) {
                return std::find(source_nodes.begin(), source_nodes.end(), value.id) == source_nodes.end();
            }));

        const fuelsim::InterfaceSummary interface =
            fuelsim::cartesian::ProblemAccess::summarize_interface(problem, 0, solve.solve.state);
        const double fuelsim_moose_force_error = relative_error(interface.total_contact_force, moose_force);
        const double fuelsim_abaqus_force_error = relative_error(interface.total_contact_force, abaqus_force);
        const double moose_abaqus_force_error = relative_error(moose_force, abaqus_force);
        const std::array<std::string, 3> names{"x", "y", "z"};
        for (std::size_t component = 0; component < 3; ++component) {
            fuelsim::test::print_relative_metrics("h20_25_fuelsim_moose_displacement_" + names[component],
                fuelsim_moose_displacement[component]);
            fuelsim::test::print_relative_metrics("h20_25_fuelsim_abaqus_displacement_" + names[component],
                fuelsim_abaqus_displacement[component]);
            fuelsim::test::print_relative_metrics("h20_25_moose_abaqus_displacement_" + names[component],
                moose_abaqus_displacement[component]);
        }
        fuelsim::test::print_relative_metrics("h20_25_fuelsim_moose_pressure", fuelsim_moose_pressure);
        fuelsim::test::print_relative_metrics("h20_25_fuelsim_abaqus_pressure", fuelsim_abaqus_pressure);
        fuelsim::test::print_relative_metrics("h20_25_moose_abaqus_pressure", moose_abaqus_pressure);
        std::cout << "h20_25_fuelsim_normal_resultant=" << interface.total_contact_force << '\n'
                  << "h20_25_moose_normal_resultant=" << moose_force << '\n'
                  << "h20_25_abaqus_normal_resultant=" << abaqus_force << '\n'
                  << "h20_25_fuelsim_moose_normal_resultant_relative_error=" << fuelsim_moose_force_error << '\n'
                  << "h20_25_fuelsim_abaqus_normal_resultant_relative_error=" << fuelsim_abaqus_force_error << '\n'
                  << "h20_25_moose_abaqus_normal_resultant_relative_error=" << moose_abaqus_force_error << '\n'
                  << "h20_25_abaqus_generated_center_pressure="
                  << (extra_abaqus == abaqus_pressure.end() ? 0.0 : extra_abaqus->pressure) << '\n';

        constexpr double same_discretization_tolerance = 1.0e-4;
        constexpr double abaqus_displacement_tolerance = 5.0e-3;
        constexpr double force_tolerance = 5.0e-3;
        passed =
            check(moose_displacement.size() == 40 && abaqus_displacement.size() == 40,
                "H20.25 compares all forty original displacement nodes")
            && check(maximum_moose_coordinate_difference < 1.0e-12 && maximum_abaqus_coordinate_difference < 1.0e-9,
                "H20.25 uses the same tracked Hex20 reference coordinates")
            && check(contact.size() == 8 && moose_contact.size() == 8 && active == 8,
                "H20.25 keeps all eight original secondary contact nodes active")
            && check(abaqus_pressure.size() == 9 && extra_abaqus_count == 1 && extra_abaqus != abaqus_pressure.end(),
                "H20.25 records the one Abaqus-generated midface contact node separately")
            && check(fuelsim::test::relative_metrics_below(fuelsim_moose_displacement[0], same_discretization_tolerance)
                         && fuelsim::test::relative_metrics_below(fuelsim_moose_pressure, same_discretization_tolerance)
                         && fuelsim_moose_force_error < same_discretization_tolerance,
                "H20.25 Fuelsim and MOOSE node-to-face normal fields and resultant agree below 0.01 percent")
            && check(
                fuelsim::test::relative_metrics_below(fuelsim_abaqus_displacement[0], abaqus_displacement_tolerance)
                    && fuelsim_abaqus_force_error < force_tolerance,
                "H20.25 Fuelsim and converted-element Abaqus normal displacement and resultant agree below 0.5 "
                "percent")
            && check(fuelsim::test::relative_metrics_below(moose_abaqus_displacement[0], abaqus_displacement_tolerance)
                         && moose_abaqus_force_error < force_tolerance,
                "H20.25 MOOSE and converted-element Abaqus normal displacement and resultant agree below 0.5 "
                "percent")
            && check(interface.total_tangential_force < 1.0e-10,
                "H20.25 pure-normal loading transfers no tangential resultant")
            && passed;
        if (passed && session.rank() == 0)
            std::cout << "[PASS] H20.25 Hex20 node-to-face three-way comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] H20.25 Hex20 node-to-face comparison raised: " << error.what() << '\n';
        return 1;
    }
}
