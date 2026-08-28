#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct ContactReference final {
    std::size_t step, id;
    fuelsim::CartesianPoint3 point;
    std::array<double, 3> normal_force, tangential_force;
    double slip_1, slip_2, opening, pressure;
};

struct ReactionReference final {
    std::size_t step;
    std::array<double, 3> reaction;
};

struct StepState final {
    double time;
    std::vector<double> solution;
    std::vector<fuelsim::CartesianContactNodeSummary> contact;
    std::vector<fuelsim::ContactPointHistory> histories;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

fuelsim::UnstructuredHex8Mesh mesh() {
    const std::vector<fuelsim::CartesianPoint3> nodes = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {2.0, 1.0, 1.0},
        {1.0, 1.0, 1.0}};
    const std::vector<fuelsim::Hex8Element> elements = {
        {{{0, 1, 2, 3, 4, 5, 6, 7}}}, {{{8, 9, 10, 11, 12, 13, 14, 15}}}};
    const std::vector<fuelsim::ElementSide> primary_all = {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}},
                                            secondary_all = {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}};
    return fuelsim::UnstructuredHex8Mesh(nodes, elements, {1, 2}, {{1, "primary"}, {2, "secondary"}},
        {{10, "primary_all_nodes", {0, 1, 2, 3, 4, 5, 6, 7}},
            {11, "secondary_all_nodes", {8, 9, 10, 11, 12, 13, 14, 15}},
            {12, "secondary_contact_nodes", {8, 11, 12, 15}}},
        {{20, "primary_all", primary_all}, {21, "primary_contact", {{0, 1}}}, {30, "secondary_all", secondary_all},
            {31, "secondary_contact", {{1, 3}}}});
}

void write_label_set(std::ofstream& output, const std::vector<std::size_t>& labels) {
    for (std::size_t index = 0; index < labels.size(); ++index) {
        output << labels[index] + 1;
        if ((index + 1) % 16 == 0 || index + 1 == labels.size())
            output << '\n';
        else
            output << ", ";
    }
}

void write_abaqus_input(const std::string& path, const fuelsim::UnstructuredHex8Mesh& generated) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write B4.0 Abaqus input: " + path);
    output << std::setprecision(16) << "*Heading\n"
           << "** B4.0: C3D8 small-sliding STS biaxial Coulomb friction history.\n"
           << "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
           << "*Node\n";
    for (std::size_t node = 0; node < generated.nodes().size(); ++node) {
        const fuelsim::CartesianPoint3& point = generated.nodes()[node];
        output << node + 1 << ", " << point.x << ", " << point.y << ", " << point.z << '\n';
    }
    output << "*Element, type=C3D8, elset=PRIMARY\n1";
    for (const std::size_t node : generated.elements()[0].nodes) output << ", " << node + 1;
    output << "\n*Element, type=C3D8, elset=SECONDARY\n2";
    for (const std::size_t node : generated.elements()[1].nodes) output << ", " << node + 1;
    output << "\n*Nset, nset=PRIMARY_ALL\n";
    write_label_set(output, generated.node_sets().at(0).nodes);
    output << "*Nset, nset=SECONDARY_ALL\n";
    write_label_set(output, generated.node_sets().at(1).nodes);
    output << "*Nset, nset=SECONDARY_CONTACT_NODES\n";
    write_label_set(output, generated.node_sets().at(2).nodes);
    output << "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\nPRIMARY, S4\n"
           << "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\nSECONDARY, S6\n"
           << "*Material, name=ELASTIC\n*Elastic\n1.e9, 0.\n"
           << "*Solid Section, elset=PRIMARY, material=ELASTIC\n,\n"
           << "*Solid Section, elset=SECONDARY, material=ELASTIC\n,\n"
           << "*Surface Interaction, name=PENALTY_CONTACT\n"
           << "*Surface Behavior, penalty=LINEAR\n1.e8,\n"
           << "*Friction, slip tolerance=1.e-5\n0.3,\n"
           << "*Contact Pair, interaction=PENALTY_CONTACT, type=SURFACE TO SURFACE, small sliding, adjust=0.\n"
           << "SECONDARY_CONTACT, PRIMARY_CONTACT\n";
    const std::array<std::string, 4> names = {"STICK", "SLIDE_A", "SLIDE_B", "SLIDE_C"};
    const std::array<double, 4> displacement_y = {3.0e-6, 12.0e-6, 20.0e-6, 25.0e-6},
                                displacement_z = {4.0e-6, 16.0e-6, 15.0e-6, 30.0e-6};
    for (std::size_t step = 0; step < names.size(); ++step) {
        output << "*Step, name=" << names[step] << ", nlgeom=NO, inc=100\n"
               << "*Static\n1., 1., 1.e-8, 1.\n"
               << "*Boundary\nPRIMARY_ALL, 1, 3, 0.\n"
               << "SECONDARY_ALL, 1, 1, -1.e-4\n"
               << "SECONDARY_ALL, 2, 2, " << displacement_y[step] << '\n'
               << "SECONDARY_ALL, 3, 3, " << displacement_z[step] << '\n'
               << "*Output, field, frequency=1\n*Node Output\nCOORD, RF, U\n"
               << "*Contact Output\nCSTRESS, CDISP, CFORCE\n*End Step\n";
    }
}

void write_fsi(const std::string& path) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Could not write B4.0 Fuelsim input: " + path);
    output << "[Case]\n  version = 3\n  problem = transient\n  geometry = cartesian_3d\n[]\n\n"
           << "[Mesh]\n  type = exodus\n  file = ../abaqus/b40_hex8_sts_friction_mesh.e\n[]\n\n"
           << "[Materials]\n  [elastic]\n"
           << "    [thermal]\n      function = constant_thermophysical\n      conductivity = 1\n"
           << "      density = 1\n      specific_heat = 1\n    []\n"
           << "    [elasticity]\n      function = constant_isotropic\n      young_modulus = 1e9\n"
           << "      poisson_ratio = 0\n    []\n  []\n[]\n\n"
           << "[Regions]\n"
           << "  [primary]\n    block = primary\n    material = elastic\n    strain = small\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n"
           << "  [secondary]\n    block = secondary\n    material = elastic\n    strain = small\n"
           << "    initial_temperature = 300\n    volumetric_heat_source = 0\n  []\n[]\n\n"
           << "[TimeFunctions]\n"
           << "  [normal_path]\n    type = piecewise_linear\n    times = 0 1 2 3 4\n"
           << "    values = 0 -1e-4 -1e-4 -1e-4 -1e-4\n  []\n"
           << "  [tangent_y_path]\n    type = piecewise_linear\n    times = 0 1 2 3 4\n"
           << "    values = 0 3e-6 12e-6 20e-6 25e-6\n  []\n"
           << "  [tangent_z_path]\n    type = piecewise_linear\n    times = 0 1 2 3 4\n"
           << "    values = 0 4e-6 16e-6 15e-6 30e-6\n  []\n[]\n\n"
           << "[Contact]\n  [interface]\n    primary = primary_contact\n    secondary = secondary_contact\n"
           << "    [mechanical]\n      formulation = penalty\n      discretization = surface_to_surface\n"
           << "      sliding = small\n      penalty = 1e8\n      mu = 0.3\n      slip_tolerance = 1e-5\n"
           << "    []\n  []\n[]\n\n"
           << "[BoundaryConditions]\n"
           << "  [primary_temperature]\n    type = dirichlet\n    boundary = primary_all\n"
           << "    field = temperature\n    value = 300\n  []\n"
           << "  [secondary_temperature]\n    type = dirichlet\n    boundary = secondary_all\n"
           << "    field = temperature\n    value = 300\n  []\n";
    for (const std::string& component : {std::string("x"), std::string("y"), std::string("z")})
        output << "  [primary_" << component << "]\n    type = dirichlet\n    boundary = primary_all\n"
               << "    field = displacement_" << component << "\n    value = 0\n  []\n";
    output << "  [secondary_x]\n    type = dirichlet\n    boundary = secondary_all\n"
           << "    field = displacement_x\n    value = 1\n    function = normal_path\n  []\n"
           << "  [secondary_y]\n    type = dirichlet\n    boundary = secondary_all\n"
           << "    field = displacement_y\n    value = 1\n    function = tangent_y_path\n  []\n"
           << "  [secondary_z]\n    type = dirichlet\n    boundary = secondary_all\n"
           << "    field = displacement_z\n    value = 1\n    function = tangent_z_path\n  []\n[]\n\n"
           << "[Executioner]\n  type = transient\n  end_time = 4\n  initial_time_step = 1\n"
           << "  minimum_time_step = 0.125\n  maximum_time_step = 1\n  growth_factor = 1\n"
           << "  cutback_factor = 0.5\n  maximum_cutbacks = 3\n  load_ramp_time = 0\n[]\n\n"
           << "[Solver]\n  absolute_tolerance = 1e-8\n  relative_tolerance = 1e-11\n"
           << "  step_tolerance = 1e-12\n  maximum_iterations = 40\n  linear_solver = direct\n"
           << "  direct_factorization = mumps\n  field_residual_scaling = true\n"
           << "  linear_relative_tolerance = 1e-11\n  maximum_linear_iterations = 400\n[]\n\n"
           << "[Outputs]\n  console = true\n[]\n";
}

void generate(const std::string& abaqus_directory, const std::string& fuelsim_directory) {
    const fuelsim::UnstructuredHex8Mesh generated = mesh();
    fuelsim::write_exodus_hex8(abaqus_directory + "/b40_hex8_sts_friction_mesh.e", generated);
    write_abaqus_input(abaqus_directory + "/b40_hex8_sts_friction.inp", generated);
    write_fsi(fuelsim_directory + "/transient_hex8_sts_friction_abaqus.fsi");
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    if (column >= values.size()) throw std::invalid_argument("B4.0 CSV row is too short: " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[column], &parsed);
    if (parsed != values[column].size() || !std::isfinite(result))
        throw std::invalid_argument("B4.0 CSV contains an invalid number: " + path);
    return result;
}

std::size_t index_value(const std::vector<std::string>& values, std::size_t column, const std::string& path) {
    const double result = number(values, column, path);
    if (result < 0.0 || result > static_cast<double>(std::numeric_limits<std::size_t>::max()) ||
        std::floor(result) != result)
        throw std::invalid_argument("B4.0 CSV contains an invalid integer: " + path);
    return static_cast<std::size_t>(result);
}

std::vector<ContactReference> read_contact(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.0 contact reference: " + path);
    std::string line;
    if (!std::getline(input, line) ||
        line != "step,id,x,y,z,cnormf_x,cnormf_y,cnormf_z,cshearf_x,cshearf_y,cshearf_z,cslip1,cslip2,copen,cpress")
        throw std::invalid_argument("Unexpected B4.0 contact header: " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        result.push_back({index_value(values, 0, path), index_value(values, 1, path),
            {number(values, 2, path), number(values, 3, path), number(values, 4, path)},
            {number(values, 5, path), number(values, 6, path), number(values, 7, path)},
            {number(values, 8, path), number(values, 9, path), number(values, 10, path)}, number(values, 11, path),
            number(values, 12, path), number(values, 13, path), number(values, 14, path)});
    }
    if (result.size() != 16) throw std::invalid_argument("B4.0 contact reference must contain sixteen rows");
    return result;
}

std::vector<ReactionReference> read_reaction(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read B4.0 reaction reference: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "step,reaction_x,reaction_y,reaction_z")
        throw std::invalid_argument("Unexpected B4.0 reaction header: " + path);
    std::vector<ReactionReference> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        result.push_back({index_value(values, 0, path),
            {number(values, 1, path), number(values, 2, path), number(values, 3, path)}});
    }
    if (result.size() != 4) throw std::invalid_argument("B4.0 reaction reference must contain four rows");
    return result;
}

class Recorder final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        states.push_back({step.time, problem.committed_solution(),
            fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, problem.committed_solution()),
            fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0)});
    }

    std::vector<StepState> states;
};

fuelsim::SolverOptions solver_options(const fuelsim::FuelSimCaseDefinition& input) {
    fuelsim::SolverOptions result;
    result.absolute_tolerance = input.solver.absolute_tolerance;
    result.relative_tolerance = input.solver.relative_tolerance;
    result.step_tolerance = input.solver.step_tolerance;
    result.maximum_iterations = input.solver.maximum_iterations;
    result.linear_solver = input.solver.linear_solver;
    result.direct_factorization = input.solver.direct_factorization;
    result.field_residual_scaling = input.solver.field_residual_scaling;
    result.linear_relative_tolerance = input.solver.linear_relative_tolerance;
    result.maximum_linear_iterations = input.solver.maximum_linear_iterations;
    return result;
}

fuelsim::TransientTimeOptions time_options(
    const fuelsim::FuelSimCaseDefinition& input, double end_time, double initial_time_step) {
    fuelsim::TransientTimeOptions result = input.transient_execution;
    result.end_time = end_time;
    result.initial_time_step = initial_time_step;
    return result;
}

bool histories_equal(
    const std::vector<fuelsim::ContactPointHistory>& first, const std::vector<fuelsim::ContactPointHistory>& second) {
    if (first.size() != second.size()) return false;
    for (std::size_t point = 0; point < first.size(); ++point)
        if (first[point].elastic_tangential_slip != second[point].elastic_tangential_slip ||
            first[point].sliding != second[point].sliding ||
            first[point].normal_multiplier != second[point].normal_multiplier ||
            first[point].cartesian_elastic_tangential_slip != second[point].cartesian_elastic_tangential_slip ||
            first[point].cartesian_total_tangential_slip != second[point].cartesian_total_tangential_slip ||
            first[point].cartesian_tangent_basis_initialized != second[point].cartesian_tangent_basis_initialized ||
            first[point].cartesian_contact_normal != second[point].cartesian_contact_normal ||
            first[point].cartesian_contact_tangent_first != second[point].cartesian_contact_tangent_first)
            return false;
    return true;
}

double jacobian_error(const fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    double maximum = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index) local[index] = state[dofs[index]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local.size()), plus = local, minus = local;
        constexpr double perturbation = 1.0e-9;
        for (std::size_t index = 0; index < local.size(); ++index) {
            direction[index] = std::sin(static_cast<double>(index + 1));
            plus[index] += perturbation * direction[index];
            minus[index] -= perturbation * direction[index];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t row = 0; row < local.size(); ++row) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local.size(); ++column)
                analytic += jacobian[row * local.size() + column] * direction[column];
            const double reference = (plus_residual[row] - minus_residual[row]) / (2.0 * perturbation);
            difference_squared += std::pow(analytic - reference, 2);
            reference_squared += reference * reference;
        }
        maximum = std::max(maximum, std::sqrt(difference_squared / reference_squared));
    }
    return maximum;
}

bool compare(const fuelsim::TransientProblem& problem, const fuelsim::UnstructuredHex8Mesh& generated,
    const std::vector<StepState>& states, const std::string& contact_path, const std::string& reaction_path) {
    const std::vector<ContactReference> reference = read_contact(contact_path);
    const std::vector<ReactionReference> reactions = read_reaction(reaction_path);
    if (states.size() != 4) throw std::invalid_argument("B4.0 Fuelsim path must contain four states");
    const std::vector<std::size_t> sources =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    fuelsim::test::FieldErrorMetrics normal_force, tangential_y, tangential_z, slip_y, slip_z, opening, pressure;
    double maximum_coordinate_error = 0.0, maximum_zero_force = 0.0, maximum_resultant_error = 0.0;
    bool all_active = true, first_sticks = true, later_path_reaches_sliding = true, biaxial_history = true,
         later_path_is_on_coulomb_circle = true;
    for (std::size_t step = 0; step < states.size(); ++step) {
        std::array<double, 3> actual_resultant{};
        std::size_t sticking_count = 0, sliding_count = 0;
        for (std::size_t node = 0; node < states[step].contact.size(); ++node) {
            const auto found = std::find_if(reference.begin(), reference.end(),
                [&](const ContactReference& value) { return value.step == step + 1 && value.id == sources[node]; });
            if (found == reference.end()) throw std::invalid_argument("B4.0 contact-node mapping is incomplete");
            const fuelsim::CartesianContactNodeSummary& actual = states[step].contact[node];
            maximum_coordinate_error =
                std::max({maximum_coordinate_error, std::abs(generated.nodes()[sources[node]].x - found->point.x),
                    std::abs(generated.nodes()[sources[node]].y - found->point.y),
                    std::abs(generated.nodes()[sources[node]].z - found->point.z)});
            normal_force.add(-actual.normal_contact_force[0], found->normal_force[0]);
            tangential_y.add(-actual.tangential_contact_force[1], found->tangential_force[1]);
            tangential_z.add(-actual.tangential_contact_force[2], found->tangential_force[2]);
            slip_y.add(actual.tangential_slip[1], -found->slip_2);
            slip_z.add(actual.tangential_slip[2], found->slip_1);
            opening.add(actual.gap, found->opening);
            pressure.add(actual.pressure, found->pressure);
            maximum_zero_force = std::max(
                {maximum_zero_force, std::abs(actual.normal_contact_force[1]), std::abs(actual.normal_contact_force[2]),
                    std::abs(actual.tangential_contact_force[0]), std::abs(found->normal_force[1]),
                    std::abs(found->normal_force[2]), std::abs(found->tangential_force[0])});
            actual_resultant[0] -= actual.normal_contact_force[0];
            actual_resultant[1] -= actual.tangential_contact_force[1];
            actual_resultant[2] -= actual.tangential_contact_force[2];
            all_active = all_active && actual.projected && actual.pressure > 0.0;
            first_sticks = first_sticks && (step != 0 || !actual.sliding);
            if (step > 0) {
                const double tangential_magnitude = std::hypot(actual.tangential_contact_force[0],
                    actual.tangential_contact_force[1], actual.tangential_contact_force[2]);
                later_path_is_on_coulomb_circle = later_path_is_on_coulomb_circle &&
                                                  std::abs(tangential_magnitude - 0.3 * actual.contact_force) < 1.0e-9;
            }
            if (actual.sliding)
                ++sliding_count;
            else
                ++sticking_count;
        }
        for (std::size_t point = 0; point < states[step].histories.size(); ++point) {
            const auto& value = states[step].histories[point].cartesian_elastic_tangential_slip;
            if (step > 0) biaxial_history = biaxial_history && std::abs(value[1]) > 0.0 && std::abs(value[2]) > 0.0;
        }
        for (std::size_t component = 0; component < 3; ++component) {
            const double scale = std::max(std::abs(reactions[step].reaction[component]), 1.0);
            maximum_resultant_error = std::max(maximum_resultant_error,
                std::abs(actual_resultant[component] + reactions[step].reaction[component]) / scale);
        }
        std::cout << "b40_step_" << step + 1 << "_sticking_constraints=" << sticking_count << '\n'
                  << "b40_step_" << step + 1 << "_sliding_constraints=" << sliding_count << '\n';
        if (step > 0) later_path_reaches_sliding = later_path_reaches_sliding && sliding_count > 0;
    }
    fuelsim::test::print_relative_metrics("b40_signed_normal_force_x", normal_force);
    fuelsim::test::print_relative_metrics("b40_signed_tangential_force_y", tangential_y);
    fuelsim::test::print_relative_metrics("b40_signed_tangential_force_z", tangential_z);
    fuelsim::test::print_relative_metrics("b40_tangential_slip_y", slip_y);
    fuelsim::test::print_relative_metrics("b40_tangential_slip_z", slip_z);
    fuelsim::test::print_relative_metrics("b40_opening", opening);
    fuelsim::test::print_relative_metrics("b40_pressure", pressure);
    std::cout << "b40_zero_reference_force_count=48\n"
              << "b40_zero_reference_force_maximum_absolute_difference=" << maximum_zero_force << '\n'
              << "b40_maximum_resultant_relative_error=" << maximum_resultant_error << '\n'
              << "b40_maximum_mesh_coordinate_difference=" << maximum_coordinate_error << '\n';
    constexpr double tolerance = 1.0e-2;
    return check(maximum_coordinate_error < 1.0e-14, "B4.0 Abaqus and Fuelsim use the same tracked HEX8 mesh") &&
           check(all_active && first_sticks && later_path_reaches_sliding && later_path_is_on_coulomb_circle,
               "B4.0 keeps every constraint active, starts in sticking, and then reaches the Coulomb circle with both "
               "tangent components") &&
           check(biaxial_history, "B4.0 stores two nonzero committed elastic-slip components") &&
           check(fuelsim::test::relative_metrics_below(normal_force, tolerance) &&
                     fuelsim::test::relative_metrics_below(tangential_y, tolerance) &&
                     fuelsim::test::relative_metrics_below(tangential_z, tolerance),
               "B4.0 signed normal and two tangential nodal-force fields pass all three Abaqus metrics below 1 "
               "percent") &&
           check(fuelsim::test::relative_metrics_below(slip_y, tolerance) &&
                     fuelsim::test::relative_metrics_below(slip_z, tolerance),
               "B4.0 both total tangential-slip fields pass all three Abaqus metrics below 1 percent") &&
           check(fuelsim::test::relative_metrics_below(opening, tolerance) &&
                     fuelsim::test::relative_metrics_below(pressure, tolerance),
               "B4.0 opening and pressure pass all three Abaqus metrics below 1 percent") &&
           check(maximum_zero_force < 1.0e-9,
               "B4.0 theoretical-zero transverse force components pass their absolute check") &&
           check(maximum_resultant_error < tolerance,
               "B4.0 signed normal and two tangential resultants agree with Abaqus below 1 percent");
}

bool run(const std::string& input_path, const std::string& contact_path, const std::string& reaction_path,
    const std::string& checkpoint_path) {
    const fuelsim::FuelSimCaseDefinition input = fuelsim::read_case_input(input_path);
    const fuelsim::UnstructuredHex8Mesh generated = fuelsim::read_exodus_hex8(input.mesh_file);
    const fuelsim::SolverOptions solver = solver_options(input);
    Recorder full_recorder;
    fuelsim::TransientProblem full(input.spatial, generated);
    const fuelsim::TransientResult full_result =
        fuelsim::solve_transient(full, time_options(input, 4.0, 1.0), solver, &full_recorder);
    bool passed = check(full_result.completed && full_recorder.states.size() == 4,
        "B4.0 completes all four prescribed biaxial friction states");
    passed = compare(full, generated, full_recorder.states, contact_path, reaction_path) && passed;

    Recorder split_recorder;
    fuelsim::TransientProblem split(input.spatial, generated);
    const fuelsim::TransientResult first =
        fuelsim::solve_transient(split, time_options(input, 2.0, 1.0), solver, &split_recorder);
    passed = check(first.completed && split_recorder.states.size() == 2,
                 "B4.0 reaches its two-component sliding checkpoint") &&
             passed;
    fuelsim::write_transient_checkpoint(checkpoint_path, split, first.next_time_step);
    fuelsim::TransientProblem restarted(input.spatial, generated);
    const double restored_step = fuelsim::restore_transient_checkpoint(checkpoint_path, restarted);
    Recorder restarted_recorder;
    const fuelsim::TransientResult second =
        fuelsim::solve_transient(restarted, time_options(input, 4.0, restored_step), solver, &restarted_recorder);
    double maximum_restart_difference = 0.0;
    if (second.completed && restarted.committed_solution().size() == full.committed_solution().size())
        for (std::size_t dof = 0; dof < full.committed_solution().size(); ++dof)
            maximum_restart_difference = std::max(maximum_restart_difference,
                std::abs(restarted.committed_solution()[dof] - full.committed_solution()[dof]));
    const bool histories_match =
        histories_equal(fuelsim::cartesian::ProblemAccess::committed_contact_histories(restarted).at(0),
            fuelsim::cartesian::ProblemAccess::committed_contact_histories(full).at(0));
    passed =
        check(second.completed && restarted_recorder.states.size() == 2 && maximum_restart_difference == 0.0 &&
                  histories_match,
            "B4.0 checkpoint restart reproduces the nodal state and both elastic-slip history components exactly") &&
        passed;
    fuelsim::TransientProblem stick_problem(input.spatial, generated);
    std::vector<double> stick_trial = stick_problem.committed_solution();
    const auto& stick_spatial = fuelsim::cartesian::ProblemAccess::view(stick_problem);
    for (std::size_t local = 0; local < stick_spatial.region_mesh(1).nodes().size(); ++local) {
        const std::size_t global = stick_spatial.global_node(1, local);
        stick_trial[stick_spatial.dof(fuelsim::Field::displacement_x, global)] = -1.0e-4;
        stick_trial[stick_spatial.dof(fuelsim::Field::displacement_y, global)] = 3.0e-6;
        stick_trial[stick_spatial.dof(fuelsim::Field::displacement_z, global)] = 4.0e-6;
    }
    stick_spatial.validate_state(stick_trial);
    const double sticking_jacobian_error = jacobian_error(stick_problem, stick_trial);
    std::vector<double> sliding_trial = full.committed_solution();
    const auto& full_spatial = fuelsim::cartesian::ProblemAccess::view(full);
    for (std::size_t local = 0; local < full_spatial.region_mesh(1).nodes().size(); ++local) {
        const std::size_t global = full_spatial.global_node(1, local);
        sliding_trial[full_spatial.dof(fuelsim::Field::displacement_y, global)] += 1.0e-6;
        sliding_trial[full_spatial.dof(fuelsim::Field::displacement_z, global)] += 2.0e-6;
    }
    full_spatial.validate_state(sliding_trial);
    const double sliding_jacobian_error = jacobian_error(full, sliding_trial);
    std::cout << "b40_sticking_jacobian_directional_error=" << sticking_jacobian_error << '\n'
              << "b40_sliding_jacobian_directional_error=" << sliding_jacobian_error << '\n'
              << "b40_restart_maximum_absolute_difference=" << maximum_restart_difference << '\n';
    passed = check(sticking_jacobian_error < 1.0e-7 && sliding_jacobian_error < 1.0e-7,
                 "B4.0 strict sticking and sliding Jacobians match centered directional differences") &&
             check(std::remove(checkpoint_path.c_str()) == 0, "B4.0 removes its checkpoint artifact") && passed;
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc == 4 && std::string(argv[1]) == "--generate") {
        try {
            generate(argv[2], argv[3]);
            return 0;
        } catch (const std::exception& error) {
            std::cerr << "[FAIL] B4.0 generation raised: " << error.what() << '\n';
            return 1;
        }
    }
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_b40_hex8_sts_friction_abaqus_tests "
                     "<case.fsi> <contact.csv> <reaction.csv> <checkpoint>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim B4.0 HEX8 Abaqus friction comparison\n");
        const bool passed = run(argv[1], argv[2], argv[3], argv[4]);
        if (passed && session.rank() == 0)
            std::cout << "[PASS] B4.0 HEX8 biaxial STS friction and restart comparison\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B4.0 comparison raised: " << error.what() << '\n';
        return 1;
    }
}
