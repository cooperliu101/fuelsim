#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

// Diagnostic only: prescribe native nodal histories without solving equilibrium.
int main(int argc, char** argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("Expected input.fsi native_nodes.csv");
        const auto input = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_hex20(input.mesh_file);
        fuelsim::TransientProblem problem(input.spatial, mesh);
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        const auto source = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
        std::ifstream stream(argv[2]);
        std::string line, token;
        if (!std::getline(stream, line))
            throw std::runtime_error("Missing node reference header");
        std::cout << std::scientific << std::setprecision(16)
                  << "time,node,gap,pressure,nx,ny,nz,tx,ty,tz,slipx,slipy,slipz,normalx,normaly,normalz,firstx,firsty,"
                     "firstz\n";
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
            problem.begin_time_step({rows.front()[0], 1.0, true});
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
            const auto summary = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
            problem.commit_time_step(state);
            const auto& histories = fuelsim::cartesian::ProblemAccess::committed_contact_histories(problem).at(0);
            for (std::size_t i = 0; i < source.size(); ++i) {
                const auto& value = summary[i];
                std::cout << rows.front()[0] << ',' << source[i] + 1 << ',' << value.gap << ',' << value.pressure;
                for (double v : value.normal_contact_force)
                    std::cout << ',' << -v;
                for (double v : value.tangential_contact_force)
                    std::cout << ',' << -v;
                for (double v : value.tangential_slip)
                    std::cout << ',' << v;
                for (double v : histories.at(i).cartesian_contact_normal)
                    std::cout << ',' << v;
                for (double v : histories.at(i).cartesian_contact_tangent_first)
                    std::cout << ',' << v;
                std::cout << '\n';
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
