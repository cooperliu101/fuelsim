#include "core/problem_backend_access.hpp"
#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "pellet_mlp.hpp"
#include "solver/solve_workflows.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>

int main(int argc, char** argv) {
    using namespace fuelsim;
    try {
        PetscSession session(argc, argv, "Full FEM pellet boundary-response training sampler");
        if (argc != 4)
            throw std::invalid_argument("Usage: fuelsim_pellet_sample base.fsi samples.txt count");
        if (session.size() != 1)
            throw std::invalid_argument("The training sampler requires one MPI rank");
        auto input = read_case_input(argv[1]);
        if (input.problem != CaseProblem::steady || input.geometry != CaseGeometry::cartesian_3d
            || input.spatial.physics != Physics::thermal || input.spatial.regions.size() != 1
            || !input.spatial.contacts.empty() || input.spatial.regions[0].pellet_response != "full"
            || input.spatial.regions[0].thermal_element != ThermalElement::dc3d8)
            throw std::invalid_argument("Training requires a full single DC3D8 pellet without gap");
        const auto mesh = read_exodus_hex8(input.mesh_file);
        const SteadyProblem prototype(input.spatial, mesh);
        const auto& spatial = BackendAccess::thermal_spatial(prototype);
        std::vector<std::size_t> nodes;
        // The complete Dirichlet surface is specified by the base production card.
        for (const auto& bc : prototype.dirichlet_conditions())
            nodes.push_back(bc.dof);
        std::sort(nodes.begin(), nodes.end());
        if (nodes.empty())
            throw std::invalid_argument("Missing training surface");
        std::vector<std::size_t> source_ids;
        for (auto node : nodes) {
            const auto& sources = spatial.source_nodes(0);
            bool found = false;
            for (std::size_t i = 0; i < sources.size(); ++i)
                if (spatial.global_temperature_node(0, i) == node) {
                    source_ids.push_back(sources[i]);
                    found = true;
                    break;
                }
            if (!found)
                throw std::logic_error("Surface node mapping failed");
        }
        double volume = 0;
        for (std::size_t e = 0; e < spatial.volume_contribution_count(); ++e)
            for (const auto& point : spatial.geometry(e).points)
                volume += point.measure;
        const std::size_t count = std::stoul(argv[3]);
        if (count < 8)
            throw std::invalid_argument("At least eight samples required");
        std::ofstream output(argv[2]);
        output << std::setprecision(17) << "FUELSIM_PELLET_DATA 1\n"
               << nodes.size() << ' ' << count << ' ' << volume << ' '
               << elements::pellet_mesh_signature(spatial.coordinates(),
                      spatial.connectivity(),
                      input.spatial.regions[0].material.functions->thermal.parameters.value("conductivity"))
               << '\n';
        for (auto id : source_ids) {
            const auto& p = mesh.nodes()[id];
            output << id << ' ' << p.x << ' ' << p.y << ' ' << p.z << '\n';
        }
        std::mt19937_64 random(7319);
        std::uniform_real_distribution<double> unit(-1.0, 1.0);
        double maximum_energy = 0;
        for (std::size_t sample = 0; sample < count; ++sample) {
            auto definition = input.spatial;
            definition.boundary_conditions.clear();
            const double mean = 750 + 250 * unit(random), axial = 80 * unit(random), cosine = 80 * unit(random),
                         sine = 80 * unit(random), q = 2e8 * (1 + 0.9 * unit(random));
            std::vector<double> smooth(nodes.size());
            for (auto& v : smooth)
                v = unit(random);
            std::vector<double> temperature;
            for (auto id : source_ids) {
                const auto& p = mesh.nodes()[id];
                const double x = p.x / .004, y = p.y / .004, z = p.z / .008 - 0.5;
                double t = mean;
                if (sample % 4 == 1)
                    t += 2 * axial * z;
                if (sample % 4 == 2)
                    t += cosine * x + sine * y;
                if (sample % 4 == 3) {
                    t += 2 * axial * z + cosine * x + sine * y + 30 * x * y + 40 * z * x;
                    // Correlated smooth Gaussian fields; never independent nodal noise.
                    for (std::size_t j = 0; j < source_ids.size(); ++j) {
                        const auto& c = mesh.nodes()[source_ids[j]];
                        const double dx = (p.x - c.x) / .004, dy = (p.y - c.y) / .004, dz = (p.z - c.z) / .008;
                        t += 30 * smooth[j] * std::exp(-2 * (dx * dx + dy * dy + dz * dz));
                    }
                }
                temperature.push_back(t);
                definition.boundary_conditions.push_back({"sample_" + std::to_string(id),
                    BoundaryConditionType::dirichlet,
                    "surface_" + std::to_string(id),
                    Field::temperature,
                    t});
            }
            definition.regions[0].volumetric_heat_source = q;
            SteadyProblem problem(std::move(definition), mesh);
            const auto result = solve_steady(problem, input.steady_execution, input.solver);
            if (!result.completed || !result.solve.converged)
                throw std::runtime_error("FEM sampling failed");
            std::vector<double> residual(problem.dof_count(), 0);
            ContributionWorkspace workspace;
            for (std::size_t e = 0; e < problem.contribution_count(); ++e) {
                problem.evaluate_contribution(e, result.solve.state, workspace, false);
                for (std::size_t i = 0; i < workspace.dofs.size(); ++i)
                    residual[workspace.dofs[i]] += workspace.residual[i];
            }
            double balance = q * volume;
            for (auto node : nodes)
                balance += residual[node];
            maximum_energy = std::max(maximum_energy, std::abs(balance));
            if (std::abs(balance) > 1e-9)
                throw std::runtime_error("FEM sample violates steady energy conservation");
            for (double t : temperature)
                output << t << ' ';
            output << q;
            for (auto node : nodes)
                output << ' ' << residual[node];
            output << '\n';
        }
        if (!output)
            throw std::runtime_error("Cannot write training dataset");
        std::cout << "samples=" << count << " boundary_nodes=" << nodes.size()
                  << " maximum_energy_error_W=" << maximum_energy << '\n';
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
