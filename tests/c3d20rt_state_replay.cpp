#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

// Manual diagnostic only: replay a fixed-increment native nodal history.
// This does not solve equilibrium and must not replace production acceptance.
int main(int argc, char** argv) {
    try {
        if (argc != 4)
            throw std::invalid_argument("Expected input.fsi native_nodes.csv output_prefix");
        const auto input = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_hex20(input.mesh_file);
        fuelsim::TransientProblem problem(input.spatial, mesh);
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        std::ifstream stream(argv[2]);
        std::string line, token;
        if (!std::getline(stream, line))
            throw std::runtime_error("Missing node reference header");
        std::ofstream thermal(std::string(argv[3]) + "_thermal.csv"), stresses(std::string(argv[3]) + "_stress.csv"),
            residuals(std::string(argv[3]) + "_residual.csv");
        if (!thermal || !stresses || !residuals)
            throw std::runtime_error("Cannot create diagnostic output files");
        thermal << std::scientific << std::setprecision(16)
                << "time,interface_heat_rate,generated,stored,dirichlet,balance\n";
        stresses << std::scientific << std::setprecision(16) << "time,element,q,xx,yy,zz,xy,yz,xz\n";
        residuals << std::scientific << std::setprecision(16)
                  << "time,node,bulk,thermal_contact,total,native_reaction\n";
        while (std::getline(stream, line)) {
            std::vector<std::vector<double>> rows;
            for (std::size_t i = 0; i < mesh.nodes().size(); ++i) {
                if (i != 0 && !std::getline(stream, line))
                    throw std::runtime_error("Incomplete node frame");
                std::istringstream values(line);
                std::vector<double> row;
                while (std::getline(values, token, ','))
                    row.push_back(std::stod(token));
                if (row.size() != 10 || row[1] != static_cast<double>(i + 1))
                    throw std::runtime_error("Invalid node reference row");
                rows.push_back(std::move(row));
            }
            const double increment = input.transient_execution.initial_time_step;
            const double time = std::round(rows.front()[0] / increment) * increment;
            if (std::abs(time - rows.front()[0]) > 5e-8)
                throw std::runtime_error("Native frames do not lie on the input's fixed time grid");
            problem.begin_time_step({time, 1.0, true});
            auto state = problem.committed_solution();
            for (std::size_t r = 0; r < spatial.region_count(); ++r) {
                const auto& region = spatial.hex20_region_mesh(r);
                for (std::size_t i = 0; i < region.nodes().size(); ++i) {
                    const auto& row = rows.at(region.source_node_ids()[i]);
                    for (std::size_t c = 0; c < 3; ++c)
                        state[spatial.field_layout()[c + 1].begin + spatial.global_node(r, i)] = row[3 + c];
                    if (region.temperature_nodes()[i])
                        state[spatial.global_temperature_node(r, i)] = row[2];
                }
            }
            problem.validate_state(state);
            const auto interface = spatial.summarize_interface(0, state);
            const std::size_t temperatures = spatial.field_layout()[0].end;
            std::vector<double> bulk(temperatures, 0.0), contact(temperatures, 0.0), total(temperatures, 0.0);
            fuelsim::ContributionWorkspace workspace;
            for (std::size_t e = 0; e < problem.contribution_count(); ++e) {
                problem.evaluate_contribution(e, state, workspace, false);
                const auto type = spatial.contribution_type(e);
                for (std::size_t i = 0; i < workspace.dofs.size(); ++i) {
                    const auto dof = workspace.dofs[i];
                    if (dof >= temperatures)
                        continue;
                    total[dof] += workspace.residual[i];
                    if (type == fuelsim::SpatialContributionType::volume)
                        bulk[dof] += workspace.residual[i];
                    else if (type == fuelsim::SpatialContributionType::thermal_contact)
                        contact[dof] += workspace.residual[i];
                }
            }
            for (std::size_t r = 0; r < spatial.region_count(); ++r) {
                const auto& region = spatial.hex20_region_mesh(r);
                for (std::size_t i = 0; i < region.nodes().size(); ++i) {
                    if (!region.temperature_nodes()[i])
                        continue;
                    const auto dof = spatial.global_temperature_node(r, i);
                    const auto source_node = region.source_node_ids()[i];
                    residuals << rows.front()[0] << ',' << source_node + 1 << ',' << bulk[dof] << ',' << contact[dof]
                              << ',' << total[dof] << ',' << rows[source_node][9] << '\n';
                }
            }
            problem.commit_time_step(state);
            const auto& balance = problem.last_conservation_summary();
            thermal << rows.front()[0] << ',' << interface.total_heat_rate << ',' << balance.generated_heat_rate << ','
                    << balance.stored_heat_rate << ',' << balance.dirichlet_heat_input_rate << ','
                    << balance.global_thermal_balance << '\n';
            std::size_t element = 0;
            for (std::size_t r = 0; r < spatial.region_count(); ++r) {
                const auto& region = spatial.hex20_region_mesh(r);
                for (std::size_t e = 0; e < region.elements().size(); ++e, ++element) {
                    const auto& history = fuelsim::cartesian::ProblemAccess::material_history(problem, r, e);
                    for (std::size_t q = 0; q < 8; ++q) {
                        stresses << rows.front()[0] << ',' << element + 1 << ',' << q;
                        const auto& s = history[q].stress;
                        for (double component : {s.xx, s.yy, s.zz, s.xy, s.yz, s.xz})
                            stresses << ',' << component;
                        stresses << '\n';
                    }
                }
            }
            std::cerr << "replayed time=" << rows.front()[0] << '\n';
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
