#include "fuelsim/core/transient_problem.hpp"
#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
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
struct SourceDofs final {
    std::size_t displacement = 0, temperature = 0;
    bool present = false, temperature_active = false;
};

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> values;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ','))
        values.push_back(value);
    return values;
}

std::vector<SourceDofs> source_dofs(const fuelsim::UnstructuredHex20Mesh& mesh,
    const fuelsim::cartesian::SpatialAssembly& spatial) {
    std::vector<SourceDofs> result(mesh.nodes().size());
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.hex20_region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids().at(local);
            const std::size_t displacement = spatial.global_node(region, local);
            SourceDofs& dofs = result.at(source);
            if (dofs.present && dofs.displacement != displacement)
                throw std::runtime_error("B6.1 replay source node maps to multiple displacement nodes");
            dofs.displacement = displacement;
            dofs.present = true;
            if (region_mesh.temperature_nodes().at(local)) {
                const std::size_t temperature = spatial.global_temperature_node(region, local);
                if (dofs.temperature_active && dofs.temperature != temperature)
                    throw std::runtime_error("B6.1 replay source node maps to multiple temperature nodes");
                dofs.temperature = temperature;
                dofs.temperature_active = true;
            }
        }
    }
    return result;
}

std::vector<double> read_state(const std::string& path,
    const fuelsim::UnstructuredHex20Mesh& mesh,
    const fuelsim::cartesian::SpatialAssembly& spatial,
    const std::vector<SourceDofs>& mapping,
    const std::vector<double>& committed) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read B6.1 replay nodal file: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "id,x,y,z,temperature,displacement_x,displacement_y,displacement_z")
        throw std::runtime_error("Unexpected B6.1 replay nodal header: " + path);
    std::vector<double> result = committed;
    std::vector<bool> visited(mesh.nodes().size(), false);
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split(line);
        if (values.size() != 8)
            throw std::runtime_error("Unexpected B6.1 replay nodal row: " + path);
        const long label = std::stol(values[0]);
        if (label <= 0 || static_cast<std::size_t>(label) > mesh.nodes().size())
            throw std::runtime_error("B6.1 replay node label is outside the source mesh: " + path);
        const std::size_t source = static_cast<std::size_t>(label - 1);
        if (visited[source])
            throw std::runtime_error("Duplicate B6.1 replay node label: " + path);
        visited[source] = true;
        const auto& point = mesh.nodes().at(source);
        const std::array<double, 3> reference = {std::stod(values[1]), std::stod(values[2]), std::stod(values[3])};
        const std::array<double, 3> expected = {point.x, point.y, point.z};
        for (std::size_t component = 0; component < reference.size(); ++component)
            if (std::abs(reference[component] - expected[component]) > 1.0e-8)
                throw std::runtime_error("B6.1 replay node coordinates do not match the source mesh: " + path);
        const SourceDofs& dofs = mapping.at(source);
        if (!dofs.present)
            throw std::runtime_error("B6.1 replay node is not active in the problem: " + path);
        if (dofs.temperature_active)
            result.at(spatial.dof(fuelsim::Field::temperature, dofs.temperature)) = std::stod(values[4]);
        result.at(spatial.dof(fuelsim::Field::displacement_x, dofs.displacement)) = std::stod(values[5]);
        result.at(spatial.dof(fuelsim::Field::displacement_y, dofs.displacement)) = std::stod(values[6]);
        result.at(spatial.dof(fuelsim::Field::displacement_z, dofs.displacement)) = std::stod(values[7]);
    }
    for (std::size_t source = 0; source < visited.size(); ++source)
        if (!visited[source])
            throw std::runtime_error("B6.1 replay nodal file omits a source node: " + path);
    return result;
}

double equivalent_stress(const fuelsim::SymmetricTensor3Values& stress) {
    const double mean = (stress.xx + stress.yy + stress.zz) / 3.0;
    return std::sqrt(1.5
                     * ((stress.xx - mean) * (stress.xx - mean) + (stress.yy - mean) * (stress.yy - mean)
                         + (stress.zz - mean) * (stress.zz - mean)
                         + 2.0 * (stress.xy * stress.xy + stress.yz * stress.yz + stress.xz * stress.xz)));
}

void write_header(std::ofstream& output) {
    output << "step,time_s,element,integration_point,vonmises_stress,peeq,ceeq,"
              "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,"
              "ee11,ee22,ee33,ee12,ee13,ee23,"
              "pe11,pe22,pe33,pe12,pe13,pe23,"
              "ce11,ce22,ce33,ce12,ce13,ce23\n";
}

void write_material_step(std::ofstream& output, std::size_t step, const fuelsim::TransientProblem& problem) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& mesh = spatial.hex20_region_mesh(region);
        for (std::size_t element = 0; element < mesh.elements().size(); ++element) {
            const auto& history = fuelsim::cartesian::ProblemAccess::material_history(problem, region, element);
            for (std::size_t point = 0; point < history.size(); ++point) {
                const auto& state = history[point];
                const auto& stress = state.stress;
                output << step << ',' << (2.0 * static_cast<double>(step)) << ','
                       << (mesh.source_element_ids().at(element) + 1U) << ',' << (point + 1U) << ','
                       << equivalent_stress(stress) << ',' << state.equivalent_plastic_strain << ','
                       << state.equivalent_creep_strain << ',' << stress.xx << ',' << stress.yy << ',' << stress.zz
                       << ',' << stress.xy << ',' << stress.xz << ',' << stress.yz;
                const std::array<const std::array<double, 6>*, 3> histories = {&state.elastic_strain,
                    &state.plastic_strain,
                    &state.creep_strain};
                for (const std::array<double, 6>* tensor : histories)
                    output << ',' << (*tensor)[0] << ',' << (*tensor)[1] << ',' << (*tensor)[2] << ',' << (*tensor)[3]
                           << ',' << (*tensor)[5] << ',' << (*tensor)[4];
                output << '\n';
            }
        }
    }
}

void write_residual_step(std::ofstream& output,
    std::size_t step,
    const fuelsim::TransientProblem& problem,
    const fuelsim::cartesian::SpatialAssembly& spatial,
    const std::vector<SourceDofs>& mapping,
    const std::vector<double>& state) {
    std::vector<double> residual(problem.dof_count(), 0.0);
    fuelsim::ContributionWorkspace workspace;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        problem.evaluate_contribution(contribution, state, workspace, false);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            residual.at(workspace.dofs[local]) += workspace.residual[local];
    }
    std::vector<bool> constrained(problem.dof_count(), false);
    for (const fuelsim::DirichletCondition& condition : problem.dirichlet_conditions())
        constrained.at(condition.dof) = true;
    double free_squared = 0.0, constrained_squared = 0.0, maximum_free = 0.0;
    double free_thermal_squared = 0.0, constrained_thermal_squared = 0.0, maximum_free_thermal = 0.0;
    for (std::size_t source = 0; source < mapping.size(); ++source) {
        const SourceDofs& dofs = mapping[source];
        const std::array<std::size_t, 3> global = {spatial.dof(fuelsim::Field::displacement_x, dofs.displacement),
            spatial.dof(fuelsim::Field::displacement_y, dofs.displacement),
            spatial.dof(fuelsim::Field::displacement_z, dofs.displacement)};
        output << step << ',' << (2.0 * static_cast<double>(step)) << ',' << (source + 1U);
        for (std::size_t dof : global)
            output << ',' << residual[dof];
        for (std::size_t dof : global)
            output << ',' << (constrained[dof] ? 1 : 0);
        if (dofs.temperature_active) {
            const std::size_t temperature = spatial.dof(fuelsim::Field::temperature, dofs.temperature);
            output << ',' << residual[temperature] << ',' << (constrained[temperature] ? 1 : 0);
        } else
            output << ",0,-1";
        output << '\n';
        for (std::size_t dof : global) {
            const double squared = residual[dof] * residual[dof];
            if (constrained[dof])
                constrained_squared += squared;
            else {
                free_squared += squared;
                maximum_free = std::max(maximum_free, std::abs(residual[dof]));
            }
        }
        if (dofs.temperature_active) {
            const std::size_t temperature = spatial.dof(fuelsim::Field::temperature, dofs.temperature);
            const double squared = residual[temperature] * residual[temperature];
            if (constrained[temperature])
                constrained_thermal_squared += squared;
            else {
                free_thermal_squared += squared;
                maximum_free_thermal = std::max(maximum_free_thermal, std::abs(residual[temperature]));
            }
        }
    }
    const double free_l2 = std::sqrt(free_squared), constrained_l2 = std::sqrt(constrained_squared);
    const double free_thermal_l2 = std::sqrt(free_thermal_squared),
                 constrained_thermal_l2 = std::sqrt(constrained_thermal_squared);
    std::cout << "b61_replay_step_" << step << "_free_mechanical_residual_l2=" << free_l2 << '\n'
              << "b61_replay_step_" << step << "_constrained_mechanical_residual_l2=" << constrained_l2 << '\n'
              << "b61_replay_step_" << step
              << "_free_to_constrained_residual_ratio=" << (constrained_l2 == 0.0 ? 0.0 : free_l2 / constrained_l2)
              << '\n'
              << "b61_replay_step_" << step << "_maximum_free_mechanical_residual=" << maximum_free << '\n';
    std::cout << "b61_replay_step_" << step << "_free_thermal_residual_l2=" << free_thermal_l2 << '\n'
              << "b61_replay_step_" << step << "_constrained_thermal_residual_l2=" << constrained_thermal_l2 << '\n'
              << "b61_replay_step_" << step << "_free_to_constrained_thermal_residual_ratio="
              << (constrained_thermal_l2 == 0.0 ? 0.0 : free_thermal_l2 / constrained_thermal_l2) << '\n'
              << "b61_replay_step_" << step << "_maximum_free_thermal_residual=" << maximum_free_thermal << '\n';
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 8 && argc != 9) {
        std::cerr << "Usage: fuelsim_b61_c3d20t_material_replay <case.fsi> <step1_nodal.csv> <step2_nodal.csv> "
                     "<step3_nodal.csv> <step4_nodal.csv> <step5_nodal.csv> <material.csv> [residual.csv]\n";
        return 2;
    }
    try {
        const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(argv[1]);
        if (definition.problem != fuelsim::CaseProblem::transient || !fuelsim::exodus_uses_hex20(definition.mesh_file))
            throw std::invalid_argument("B6.1 material replay requires a transient C3D20T input");
        const fuelsim::UnstructuredHex20Mesh mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
        fuelsim::TransientProblem problem(definition.spatial, mesh);
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        const std::vector<SourceDofs> mapping = source_dofs(mesh, spatial);
        std::ofstream output(argv[7]);
        if (!output)
            throw std::runtime_error("Could not write B6.1 replay material output");
        output << std::scientific << std::setprecision(17);
        write_header(output);
        std::ofstream residual_output;
        if (argc == 9) {
            residual_output.open(argv[8]);
            if (!residual_output)
                throw std::runtime_error("Could not write B6.1 replay residual output");
            residual_output << std::scientific << std::setprecision(17)
                            << "step,time_s,node,residual_x_n,residual_y_n,residual_z_n,"
                               "constrained_x,constrained_y,constrained_z,residual_temperature_w,"
                               "constrained_temperature\n";
        }
        for (std::size_t step = 1; step <= 5; ++step) {
            const std::vector<double> state =
                read_state(argv[step + 1], mesh, spatial, mapping, problem.committed_solution());
            problem.begin_time_step({2.0 * static_cast<double>(step), 1.0, true});
            if (residual_output)
                write_residual_step(residual_output, step, problem, spatial, mapping, state);
            problem.commit_time_step(state);
            write_material_step(output, step, problem);
        }
        std::cout << "b61_replay_steps=5\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] B6.1 material replay: " << error.what() << '\n';
        return 1;
    }
}
