#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Matrix4 = std::array<double, 16>;

struct ProbeRow final {
    std::string step;
    std::size_t input = 0, output = 0;
    double closure_delta = 0.0, opening = 0.0, pressure = 0.0;
    std::array<double, 3> force{};
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> result;
    std::stringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

std::vector<ProbeRow> read_probe(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read the B3.8 Abaqus operator reference: " + path);
    std::string line;
    if (!std::getline(input, line) ||
        line != "step,input_local_node,input_label,output_local_node,output_label,closure_delta_m,copen_m,cpress_pa,"
                "cnormf_x_n,cnormf_y_n,cnormf_z_n")
        throw std::runtime_error("The B3.8 Abaqus operator reference header is invalid");
    std::vector<ProbeRow> result;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> values = split(line);
        if (values.size() != 11) throw std::runtime_error("The B3.8 Abaqus operator row has the wrong width");
        result.push_back({values[0], static_cast<std::size_t>(std::stoul(values[1])),
            static_cast<std::size_t>(std::stoul(values[3])), std::stod(values[5]), std::stod(values[6]),
            std::stod(values[7]), {std::stod(values[8]), std::stod(values[9]), std::stod(values[10])}});
    }
    if (result.size() != 36) throw std::runtime_error("The B3.8 Abaqus operator reference must contain 36 rows");
    return result;
}

ProbeRow row(const std::vector<ProbeRow>& rows, const std::string& step, std::size_t output) {
    const auto found = std::find_if(
        rows.begin(), rows.end(), [&](const ProbeRow& value) { return value.step == step && value.output == output; });
    if (found == rows.end()) throw std::runtime_error("The B3.8 Abaqus operator reference is incomplete");
    return *found;
}

double maximum_absolute_difference(const Matrix4& first, const Matrix4& second) {
    double result = 0.0;
    for (std::size_t entry = 0; entry < first.size(); ++entry)
        result = std::max(result, std::abs(first[entry] - second[entry]));
    return result;
}

double relative_l2_difference(const Matrix4& actual, const Matrix4& reference) {
    double difference_squared = 0.0, reference_squared = 0.0;
    for (std::size_t entry = 0; entry < actual.size(); ++entry) {
        difference_squared += std::pow(actual[entry] - reference[entry], 2);
        reference_squared += reference[entry] * reference[entry];
    }
    return std::sqrt(difference_squared / reference_squared);
}

Matrix4 identified_opening_operator(const std::vector<ProbeRow>& rows) {
    Matrix4 result{};
    for (std::size_t input = 0; input < 4; ++input)
        for (std::size_t output = 0; output < 4; ++output) {
            const ProbeRow plus = row(rows, "S" + std::to_string(input + 1) + "_PLUS", output + 1);
            const ProbeRow minus = row(rows, "S" + std::to_string(input + 1) + "_MINUS", output + 1);
            result[output * 4 + input] = -(plus.opening - minus.opening) / (2.0e-6);
        }
    return result;
}

Matrix4 identified_force_tangent(const std::vector<ProbeRow>& rows) {
    Matrix4 result{};
    for (std::size_t input = 0; input < 4; ++input)
        for (std::size_t output = 0; output < 4; ++output) {
            const ProbeRow plus = row(rows, "S" + std::to_string(input + 1) + "_PLUS", output + 1);
            const ProbeRow minus = row(rows, "S" + std::to_string(input + 1) + "_MINUS", output + 1);
            result[output * 4 + input] = (plus.force[0] - minus.force[0]) / (2.0e-6);
        }
    return result;
}

Matrix4 expected_opening_operator() {
    return {9.0 / 16.0, 3.0 / 16.0, 1.0 / 16.0, 3.0 / 16.0, 3.0 / 16.0, 9.0 / 16.0, 3.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0,
        3.0 / 16.0, 9.0 / 16.0, 3.0 / 16.0, 3.0 / 16.0, 1.0 / 16.0, 3.0 / 16.0, 9.0 / 16.0};
}

Matrix4 expected_force_tangent(const Matrix4& opening) {
    Matrix4 result{};
    constexpr double penalty = 1.0e8, area = 0.25;
    for (std::size_t row_value = 0; row_value < 4; ++row_value)
        for (std::size_t column = 0; column < 4; ++column)
            for (std::size_t constraint = 0; constraint < 4; ++constraint)
                result[row_value * 4 + column] +=
                    area * penalty * opening[constraint * 4 + row_value] * opening[constraint * 4 + column];
    return result;
}

bool check_abaqus_identification(
    const std::vector<ProbeRow>& rows, const std::string& summary_path, bool finite_sliding_probe) {
    std::ifstream input(summary_path);
    if (!input) throw std::runtime_error("Could not read the B3.8 Abaqus contact summary: " + summary_path);
    const std::string summary((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    const Matrix4 opening = identified_opening_operator(rows), expected_opening = expected_opening_operator(),
                  force_tangent = identified_force_tangent(rows), expected_tangent = expected_force_tangent(opening);
    double base_opening_error = 0.0, base_pressure_error = 0.0, base_force_error = 0.0, maximum_transverse_force = 0.0;
    for (std::size_t output = 1; output <= 4; ++output) {
        const ProbeRow base = row(rows, "BASE", output);
        base_opening_error = std::max(base_opening_error, std::abs(base.opening + 1.0e-4));
        base_pressure_error = std::max(base_pressure_error, std::abs(base.pressure - 1.0e4));
        base_force_error = std::max(base_force_error, std::abs(base.force[0] - 2500.0));
        maximum_transverse_force =
            std::max(maximum_transverse_force, std::max(std::abs(base.force[1]), std::abs(base.force[2])));
    }
    const double opening_error = maximum_absolute_difference(opening, expected_opening),
                 tangent_relative_l2 = relative_l2_difference(force_tangent, expected_tangent),
                 tangent_maximum_error = maximum_absolute_difference(force_tangent, expected_tangent);
    std::cout << "b38_abaqus_opening_operator_maximum_absolute_error=" << opening_error << '\n'
              << "b38_abaqus_force_tangent_relative_l2_error=" << tangent_relative_l2 << '\n'
              << "b38_abaqus_force_tangent_maximum_absolute_error=" << tangent_maximum_error << '\n'
              << "b38_abaqus_base_opening_maximum_absolute_error=" << base_opening_error << '\n'
              << "b38_abaqus_base_pressure_maximum_absolute_error=" << base_pressure_error << '\n'
              << "b38_abaqus_base_force_maximum_absolute_error=" << base_force_error << '\n'
              << "b38_abaqus_zero_reference_transverse_force_count=8\n"
              << "b38_abaqus_zero_reference_transverse_force_maximum_absolute_difference=" << maximum_transverse_force
              << '\n';
    const double opening_tolerance = finite_sliding_probe ? 1.0e-9 : 5.0e-13,
                 tangent_relative_tolerance = finite_sliding_probe ? 1.0e-8 : 5.0e-13,
                 tangent_absolute_tolerance = finite_sliding_probe ? 1.0e-2 : 1.0e-6,
                 base_pressure_tolerance = finite_sliding_probe ? 1.0e-8 : 1.0e-9;
    return check(summary.find("PENALTY CONSTRAINT ENFORCEMENT WITH PENALTY STIFFNESS OF   1.00000E+08") !=
                         std::string::npos &&
                     summary.find("SURFACE TO SURFACE WITH THICKNESS") != std::string::npos &&
                     summary.find("CONSTRAINT POSITION IS AT NODE") != std::string::npos &&
                     summary.find("SUPPLEMENTARY CONSTRAINTS = NO") != std::string::npos &&
                     summary.find("NUMBER OF INTERNAL ELEMENTS GENERATED FOR CONTACT         4") != std::string::npos,
               "Abaqus reports four node-positioned surface-to-surface penalty constraints") &&
           check(opening_error < opening_tolerance,
               "the identified Abaqus C3D8 opening operator is the Q4 shape matrix at parent coordinates plus or "
               "minus one half") &&
           check(tangent_relative_l2 < tangent_relative_tolerance && tangent_maximum_error < tangent_absolute_tolerance,
               "the Abaqus nodal-force tangent is penalty times A transpose W A with four equal quarter-face areas") &&
           check(base_opening_error < 1.0e-16 && base_pressure_error < base_pressure_tolerance &&
                     base_force_error < 1.0e-8 && maximum_transverse_force < 1.0e-12,
               "the Abaqus uniform-closure response has the identified opening, pressure, area, and normal direction");
}

fuelsim::UnstructuredHex8Mesh matching_mesh() {
    const std::vector<fuelsim::CartesianPoint3> nodes = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0},
        {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {1.0, 0.0, 1.0}, {2.0, 0.0, 1.0}, {2.0, 1.0, 1.0},
        {1.0, 1.0, 1.0}};
    const std::vector<fuelsim::Hex8Element> elements = {
        {{{0, 1, 2, 3, 4, 5, 6, 7}}}, {{{8, 9, 10, 11, 12, 13, 14, 15}}}};
    return fuelsim::UnstructuredHex8Mesh(nodes, elements, {1, 2}, {{1, "primary"}, {2, "secondary"}}, {},
        {{10, "primary_contact", {{0, 1}}}, {20, "secondary_contact", {{1, 3}}}});
}

fuelsim::ThermoelasticProperties material() {
    return fuelsim::test::thermoelastic(0.0, 1.0, 1.0e9, 0.0, 0.0, 300.0, 0.0, 0.0, 0.0, 1.0, 1.0);
}

fuelsim::SpatialDefinition definition() {
    fuelsim::SpatialDefinition result;
    result.regions = {
        {"primary", "primary", material(), 0.0, 300.0}, {"secondary", "secondary", material(), 0.0, 300.0}};
    fuelsim::ContactDefinition contact;
    contact.name = "matching_surface_contact";
    contact.primary = "primary_contact";
    contact.secondary = "secondary_contact";
    contact.mechanical = true;
    contact.penalty = 1.0e8;
    contact.mechanical_discretization = fuelsim::MechanicalContactDiscretization::surface_to_surface;
    contact.mechanical_sliding = fuelsim::MechanicalContactSliding::small;
    result.contacts.push_back(contact);
    return result;
}

std::size_t global_node_for_source(const fuelsim::SteadyProblem& problem, std::size_t region, std::size_t source) {
    const auto& mesh = fuelsim::cartesian::ProblemAccess::region_mesh(problem, region);
    const auto found = std::find(mesh.source_node_ids().begin(), mesh.source_node_ids().end(), source);
    if (found == mesh.source_node_ids().end()) throw std::logic_error("The B3.8 production node mapping failed");
    return fuelsim::cartesian::ProblemAccess::region_node_offset(problem, region) +
           static_cast<std::size_t>(found - mesh.source_node_ids().begin());
}

std::array<double, 4> production_openings(const fuelsim::SteadyProblem& problem, const std::vector<double>& state,
    const std::array<std::size_t, 4>& ordered_sources) {
    const std::vector<std::size_t> sources =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    const std::vector<fuelsim::CartesianContactNodeSummary> summaries =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
    std::array<double, 4> result{};
    for (std::size_t output = 0; output < 4; ++output) {
        const auto found = std::find(sources.begin(), sources.end(), ordered_sources[output]);
        if (found == sources.end()) throw std::logic_error("The B3.8 production summary mapping failed");
        result[output] = summaries[static_cast<std::size_t>(found - sources.begin())].gap;
    }
    return result;
}

std::array<double, 4> production_forces(const fuelsim::SteadyProblem& problem, const std::vector<double>& state,
    const std::array<std::size_t, 4>& ordered_sources) {
    const std::vector<std::size_t> sources =
        fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
    const std::vector<fuelsim::CartesianContactNodeSummary> summaries =
        fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
    std::array<double, 4> result{};
    for (std::size_t output = 0; output < 4; ++output) {
        const auto found = std::find(sources.begin(), sources.end(), ordered_sources[output]);
        if (found == sources.end()) throw std::logic_error("The B3.8 production force mapping failed");
        result[output] = summaries[static_cast<std::size_t>(found - sources.begin())].contact_force;
    }
    return result;
}

bool check_production_operator(const std::vector<ProbeRow>& rows, bool finite_sliding_probe) {
    const fuelsim::UnstructuredHex8Mesh mesh = matching_mesh();
    fuelsim::SteadyProblem problem(definition(), mesh);
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(problem);
    const std::array<std::size_t, 4> ordered_sources = {8, 11, 15, 12};
    std::vector<double> base = problem.initial_state();
    for (const std::size_t source : ordered_sources)
        base[dofs.dof(fuelsim::Field::displacement_x, global_node_for_source(problem, 1, source))] = -1.0e-4;
    problem.validate_state(base);

    const std::array<double, 4> base_openings = production_openings(problem, base, ordered_sources),
                                base_forces = production_forces(problem, base, ordered_sources);
    Matrix4 opening{}, force_tangent{};
    for (std::size_t input = 0; input < 4; ++input) {
        std::vector<double> plus = base, minus = base;
        const std::size_t input_dof =
            dofs.dof(fuelsim::Field::displacement_x, global_node_for_source(problem, 1, ordered_sources[input]));
        plus[input_dof] -= 1.0e-6;
        minus[input_dof] += 1.0e-6;
        problem.validate_state(plus);
        const std::array<double, 4> plus_opening = production_openings(problem, plus, ordered_sources),
                                    plus_force = production_forces(problem, plus, ordered_sources);
        problem.validate_state(minus);
        const std::array<double, 4> minus_opening = production_openings(problem, minus, ordered_sources),
                                    minus_force = production_forces(problem, minus, ordered_sources);
        for (std::size_t output = 0; output < 4; ++output) {
            opening[output * 4 + input] = -(plus_opening[output] - minus_opening[output]) / (2.0e-6);
            force_tangent[output * 4 + input] = (plus_force[output] - minus_force[output]) / (2.0e-6);
        }
    }

    double base_opening_error = 0.0, base_force_error = 0.0;
    for (std::size_t output = 0; output < 4; ++output) {
        const ProbeRow reference = row(rows, "BASE", output + 1);
        base_opening_error = std::max(base_opening_error, std::abs(base_openings[output] - reference.opening));
        base_force_error = std::max(base_force_error, std::abs(base_forces[output] - reference.force[0]));
    }
    const Matrix4 abaqus_opening = identified_opening_operator(rows), abaqus_tangent = identified_force_tangent(rows);
    const double opening_relative_l2 = relative_l2_difference(opening, abaqus_opening),
                 opening_maximum_error = maximum_absolute_difference(opening, abaqus_opening),
                 tangent_relative_l2 = relative_l2_difference(force_tangent, abaqus_tangent),
                 tangent_maximum_error = maximum_absolute_difference(force_tangent, abaqus_tangent);

    std::size_t mechanical_contributions = 0;
    std::array<double, 3> force_balance{};
    double maximum_jacobian_directional_error = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact) continue;
        ++mechanical_contributions;
        std::vector<std::size_t> local_dofs;
        spatial.contribution_dofs(contribution, local_dofs);
        std::vector<double> local_state(local_dofs.size());
        for (std::size_t local = 0; local < local_dofs.size(); ++local) local_state[local] = base[local_dofs[local]];
        std::vector<double> residual, jacobian;
        spatial.compute_contribution(contribution, local_state, nullptr, nullptr, 0.0, residual, &jacobian);
        std::vector<double> direction(local_state.size()), plus = local_state, minus = local_state;
        constexpr double perturbation = 1.0e-9;
        for (std::size_t local = 0; local < local_state.size(); ++local) {
            direction[local] = std::sin(static_cast<double>(local + 1));
            plus[local] += perturbation * direction[local];
            minus[local] -= perturbation * direction[local];
        }
        std::vector<double> plus_residual, minus_residual;
        spatial.compute_contribution(contribution, plus, nullptr, nullptr, 0.0, plus_residual, nullptr);
        spatial.compute_contribution(contribution, minus, nullptr, nullptr, 0.0, minus_residual, nullptr);
        double difference_squared = 0.0, reference_squared = 0.0;
        for (std::size_t local = 0; local < local_state.size(); ++local) {
            double analytic = 0.0;
            for (std::size_t column = 0; column < local_state.size(); ++column)
                analytic += jacobian[local * local_state.size() + column] * direction[column];
            const double reference = (plus_residual[local] - minus_residual[local]) / (2.0 * perturbation);
            difference_squared += std::pow(analytic - reference, 2);
            reference_squared += reference * reference;
            const std::size_t global = local_dofs[local];
            for (std::size_t component = 0; component < 3; ++component) {
                const auto& field = spatial.field_layout()[component + 1];
                if (global >= field.begin && global < field.end) force_balance[component] += residual[local];
            }
        }
        maximum_jacobian_directional_error =
            std::max(maximum_jacobian_directional_error, std::sqrt(difference_squared / reference_squared));
    }

    std::cout << "b38_fuelsim_abaqus_opening_operator_relative_l2_error=" << opening_relative_l2 << '\n'
              << "b38_fuelsim_abaqus_opening_operator_maximum_absolute_error=" << opening_maximum_error << '\n'
              << "b38_fuelsim_abaqus_force_tangent_relative_l2_error=" << tangent_relative_l2 << '\n'
              << "b38_fuelsim_abaqus_force_tangent_maximum_absolute_error=" << tangent_maximum_error << '\n'
              << "b38_fuelsim_abaqus_base_opening_maximum_absolute_error=" << base_opening_error << '\n'
              << "b38_fuelsim_abaqus_base_force_maximum_absolute_error=" << base_force_error << '\n'
              << "b38_fuelsim_contact_jacobian_directional_error=" << maximum_jacobian_directional_error << '\n'
              << "b38_fuelsim_contact_force_balance=" << force_balance[0] << ',' << force_balance[1] << ','
              << force_balance[2] << '\n';

    bool finite_strain_rejected = false;
    try {
        fuelsim::SpatialDefinition invalid = definition();
        invalid.regions.front().strain_formulation = fuelsim::StrainFormulation::finite;
        fuelsim::SteadyProblem invalid_problem(std::move(invalid), mesh);
        (void)invalid_problem;
    } catch (const std::invalid_argument&) { finite_strain_rejected = true; }

    const double operator_tolerance = finite_sliding_probe ? 1.0e-9 : 1.0e-12,
                 tangent_absolute_tolerance = finite_sliding_probe ? 1.0e-2 : 1.0e-6;
    return check(mechanical_contributions == 4,
               "Fuelsim constructs four averaged constraints for one C3D8 secondary face") &&
           check(opening_relative_l2 < operator_tolerance && opening_maximum_error < operator_tolerance,
               "the Fuelsim opening operator reproduces the identified Abaqus operator") &&
           check(tangent_relative_l2 < operator_tolerance && tangent_maximum_error < tangent_absolute_tolerance,
               "the Fuelsim equivalent nodal-force tangent reproduces the identified Abaqus tangent") &&
           check(base_opening_error < 1.0e-15 && base_force_error < 1.0e-8,
               "the Fuelsim uniform-closure opening and nodal force reproduce Abaqus") &&
           check(std::abs(force_balance[0]) < 1.0e-10 && std::abs(force_balance[1]) < 1.0e-10 &&
                     std::abs(force_balance[2]) < 1.0e-10,
               "the Fuelsim averaged constraints preserve exact action-reaction balance") &&
           check(maximum_jacobian_directional_error < 1.0e-7,
               "the Fuelsim averaged-constraint Jacobian matches a centered directional difference") &&
           check(finite_strain_rejected,
               "HEX8 small-sliding surface-to-surface contact explicitly rejects finite-strain regions");
}
} // namespace

int main(int argc, char** argv) {
    if ((argc != 3 && argc != 4) || (argc == 4 && std::string(argv[3]) != "finite-sliding-probe")) {
        std::cerr << "Usage: fuelsim_b38_hex8_sts_identification_tests <operator.csv> <contact_summary.txt> "
                     "[finite-sliding-probe]\n";
        return 2;
    }
    fuelsim::PetscSession session(argc, argv, "fuelsim B3.8 HEX8 small-sliding surface contact identification\n");
    try {
        const std::vector<ProbeRow> rows = read_probe(argv[1]);
        const bool finite_sliding_probe = argc == 4 && std::string(argv[3]) == "finite-sliding-probe";
        const bool passed = check_abaqus_identification(rows, argv[2], finite_sliding_probe) &&
                            check_production_operator(rows, finite_sliding_probe);
        if (passed && session.rank() == 0)
            std::cout << "[PASS] B3.8 HEX8 small-sliding surface-to-surface operator identification\n";
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
