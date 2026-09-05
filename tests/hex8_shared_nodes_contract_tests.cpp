#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/exodus_result_reader.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    try {
        const auto definition = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_hex8(definition.mesh_file);
        const bool b35 = std::string(argv[3]) == "b35";
        if (!b35 && std::string(argv[3]) != "b36") throw std::invalid_argument("Unknown shared-node contract");
        std::unique_ptr<fuelsim::SteadyProblem> steady;
        std::unique_ptr<fuelsim::TransientProblem> transient;
        if (b35)
            steady = std::make_unique<fuelsim::SteadyProblem>(definition.spatial, mesh);
        else
            transient = std::make_unique<fuelsim::TransientProblem>(definition.spatial, mesh);
        const auto& spatial = b35 ? fuelsim::cartesian::ProblemAccess::view(*steady)
                                  : fuelsim::cartesian::ProblemAccess::view(*transient);
        const auto output = fuelsim::test::read_final_exodus_results(argv[2]);
        if (output.nodes.size() != mesh.nodes().size()) throw std::runtime_error("Shared-node result count differs");
        const std::size_t dof_count = b35 ? steady->dof_count() : transient->dof_count();
        if (spatial.node_count() != mesh.nodes().size() || dof_count != 4 * mesh.nodes().size())
            throw std::runtime_error("Shared nodes do not have exactly one global four-field node per source node");
        std::vector<std::size_t> occurrences(mesh.nodes().size(), 0), source_global(mesh.nodes().size(), 0);
        std::vector<double> state(dof_count);
        const std::array<std::string, 4> names = {"temperature", "displacement_x", "displacement_y", "displacement_z"};
        for (std::size_t region = 0; region < spatial.region_count(); ++region) {
            const auto& local_mesh = spatial.region_mesh(region);
            for (std::size_t node = 0; node < local_mesh.nodes().size(); ++node) {
                const auto source = local_mesh.source_node_ids()[node], global = spatial.global_node(region, node);
                if (occurrences[source]++ == 0)
                    source_global[source] = global;
                else if (source_global[source] != global)
                    throw std::runtime_error("Shared source maps to distinct global nodes");
                for (std::size_t f = 0; f < 4; ++f)
                    state.at(spatial.field_layout()[f].begin + global) = output.nodal(names[f]).at(source);
            }
        }
        for (std::size_t source = 0; source < mesh.nodes().size(); ++source) {
            if (occurrences[source] == 0) throw std::runtime_error("Source node was not mapped");
            if (b35 && occurrences[source] != (std::abs(mesh.nodes()[source].x - 1.0) < 1.0e-12 ? 2U : 1U))
                throw std::runtime_error("B3.5 shared interface has incorrect region ownership count");
        }
        if (b35) {
            std::vector<std::size_t> interface_dofs;
            for (std::size_t source = 0; source < mesh.nodes().size(); ++source)
                if (std::abs(mesh.nodes()[source].x - 1.0) < 1.0e-12)
                    for (const auto& field : spatial.field_layout())
                        interface_dofs.push_back(field.begin + source_global[source]);
            std::vector<double> balance(interface_dofs.size(), 0.0), scale(interface_dofs.size(), 0.0);
            fuelsim::ContributionWorkspace workspace;
            steady->validate_state(state);
            for (std::size_t contribution = 0; contribution < steady->contribution_count(); ++contribution) {
                steady->evaluate_contribution(contribution, state, workspace, false);
                for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
                    for (std::size_t item = 0; item < interface_dofs.size(); ++item)
                        if (workspace.dofs[local] == interface_dofs[item]) {
                            balance[item] += workspace.residual[local];
                            scale[item] += std::abs(workspace.residual[local]);
                        }
            }
            for (std::size_t item = 0; item < interface_dofs.size(); ++item) {
                const double relative = std::abs(balance[item]) / std::max(scale[item], 1.0);
                std::cout << "b35_interface_residual_" << item << '=' << relative << '\n';
                if (!(relative < 1.0e-8)) throw std::runtime_error("Shared interface residual is not balanced");
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
