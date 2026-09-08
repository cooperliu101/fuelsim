#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
std::vector<double> read_row(std::ifstream& input) {
    std::string line, token;
    if (!std::getline(input, line))
        throw std::runtime_error("Missing reference row");
    std::istringstream stream(line);
    std::vector<double> row;
    while (std::getline(stream, token, ','))
        row.push_back(std::stod(token));
    return row;
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5)
            throw std::invalid_argument("Expected input.fsi checkpoint nodes.csv contact.csv");
        const auto input = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_hex20(input.mesh_file);
        fuelsim::TransientProblem problem(input.spatial, mesh);
        fuelsim::restore_transient_checkpoint(argv[2], problem);
        problem.begin_time_step({problem.committed_time() + 0.025, 1.0, true});
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
        auto state = problem.committed_solution();
        const auto own = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        std::ifstream nodes(argv[3]), contacts(argv[4]);
        std::string header;
        std::getline(nodes, header);
        std::getline(contacts, header);
        std::map<std::size_t, std::vector<double>> reference;
        for (std::size_t i = 0; i < mesh.nodes().size(); ++i) {
            auto row = read_row(nodes);
            if (row.size() != 10 || std::abs(row[0] - problem.committed_time()) > 1e-6)
                throw std::runtime_error("Node reference does not match checkpoint time");
            if (!reference.emplace(static_cast<std::size_t>(row[1]) - 1, row).second)
                throw std::runtime_error("Duplicate node label");
        }
        for (std::size_t r = 0; r < spatial.region_count(); ++r) {
            const auto& region = spatial.hex20_region_mesh(r);
            for (std::size_t i = 0; i < region.nodes().size(); ++i) {
                const auto& row = reference.at(region.source_node_ids()[i]);
                for (std::size_t c = 0; c < 3; ++c)
                    state[spatial.field_layout()[c + 1].begin + spatial.global_node(r, i)] = row[3 + c];
                if (region.temperature_nodes()[i])
                    state[spatial.global_temperature_node(r, i)] = row[2];
            }
        }
        problem.validate_state(state);
        const auto replay = fuelsim::cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, state);
        const auto source = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(problem, 0);
        std::map<std::size_t, std::vector<double>> native;
        for (std::size_t i = 0; i < source.size(); ++i) {
            auto row = read_row(contacts);
            if (row.size() != 14 || std::abs(row[0] - problem.committed_time()) > 1e-6)
                throw std::runtime_error("Contact reference does not match checkpoint time");
            native.emplace(static_cast<std::size_t>(row[1]) - 1, row);
        }
        std::cout << std::scientific << std::setprecision(12)
                  << "diagnostic_only_geometry_replay\nnode,fuelsim_gap,replayed_gap,abaqus_gap,fuelsim_pressure,"
                     "replayed_pressure,abaqus_pressure,direct_nodal_gap,direct_projected\n";
        for (std::size_t i = 0; i < source.size(); ++i) {
            const auto& row = native.at(source[i]);
            static constexpr std::array<std::array<std::size_t, 8>, 6> faces = {{{{0, 1, 5, 4, 8, 13, 16, 12}},
                {{1, 2, 6, 5, 9, 14, 17, 13}},
                {{2, 3, 7, 6, 10, 15, 18, 14}},
                {{3, 0, 4, 7, 11, 12, 19, 15}},
                {{0, 3, 2, 1, 11, 10, 9, 8}},
                {{4, 5, 6, 7, 16, 17, 18, 19}}}};
            const auto current = [&](std::size_t node) {
                const auto& x = mesh.nodes()[node];
                const auto& values = reference.at(node);
                return fuelsim::CartesianPoint3{x.x + values[3], x.y + values[4], x.z + values[5]};
            };
            double direct_gap = std::numeric_limits<double>::infinity();
            for (const auto& set : mesh.side_sets()) {
                if (set.name != input.spatial.contacts[0].primary)
                    continue;
                for (const auto& side : set.sides) {
                    fuelsim::Quad8ToQuad8MechanicalGeometry geometry{};
                    geometry.secondary_coordinates.fill(current(source[i]));
                    geometry.secondary_displacement_shape[0] = 1.0;
                    geometry.normal_orientation = 1.0;
                    geometry.finite_sliding = true;
                    for (std::size_t j = 0; j < 8; ++j)
                        geometry.primary_coordinates[j] =
                            current(mesh.elements()[side.element].nodes[faces[side.local_side][j]]);
                    const auto projection = fuelsim::compute_quad8_to_quad8_contact_projection(geometry, {});
                    if (projection.projected && std::abs(projection.gap) < std::abs(direct_gap))
                        direct_gap = projection.gap;
                }
            }
            const bool direct_projected = std::isfinite(direct_gap);
            if (!direct_projected)
                direct_gap = std::numeric_limits<double>::quiet_NaN();
            std::cout << source[i] + 1 << ',' << own[i].gap << ',' << replay[i].gap << ',' << row[2] << ','
                      << own[i].pressure << ',' << replay[i].pressure << ',' << row[3] << ',' << direct_gap << ','
                      << direct_projected << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
