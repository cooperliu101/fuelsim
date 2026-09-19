#include "core/problem_backend_access.hpp"
#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "pellet_mlp.hpp"
#include "pellet_thermal.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>

// An internal linearization/physics test, not an end-to-end replacement for fuelsim -i.
int main(int argc, char** argv) {
    using namespace fuelsim;
    try {
        if (argc != 4)
            throw std::invalid_argument("Usage: fuelsim_pellet_response_tests full.fsi model.txt evidence.txt");
        auto input = read_case_input(argv[1]);
        if (input.problem != CaseProblem::steady || input.spatial.physics != Physics::thermal
            || input.spatial.regions.size() != 1 || !input.spatial.contacts.empty()
            || input.spatial.regions[0].pellet_response != "full"
            || input.spatial.regions[0].material.functions->thermal.name != "constant_thermophysical")
            throw std::invalid_argument("Response reference requires one constant-k full FEM pellet");
        // Unit-source, zero-temperature assembly isolates the load exactly. There
        // are no gap terms and Dirichlet rows have not been substituted.
        input.spatial.regions[0].volumetric_heat_source = 1;
        const auto mesh = read_exodus_hex8(input.mesh_file);
        SteadyProblem full(input.spatial, mesh);
        full.set_load_factor(1);
        const auto& spatial = BackendAccess::thermal_spatial(full);
        const auto n = full.dof_count();
        std::vector<double> k(n * n, 0), source(n, 0), zero(n, 0);
        ContributionWorkspace workspace;
        for (std::size_t e = 0; e < full.contribution_count(); ++e) {
            full.evaluate_contribution(e, zero, workspace, true);
            const auto m = workspace.dofs.size();
            for (std::size_t i = 0; i < m; ++i) {
                const auto row = workspace.dofs[i];
                source[row] -= workspace.residual[i];
                for (std::size_t j = 0; j < m; ++j)
                    k[row * n + workspace.dofs[j]] += workspace.jacobian[i * m + j];
            }
        }
        std::vector<std::size_t> boundary, ids;
        for (const auto& condition : full.dirichlet_conditions())
            boundary.push_back(condition.dof);
        std::sort(boundary.begin(), boundary.end());
        for (auto dof : boundary) {
            bool found = false;
            for (std::size_t local = 0; local < spatial.source_nodes(0).size(); ++local)
                if (spatial.global_temperature_node(0, local) == dof) {
                    ids.push_back(spatial.source_nodes(0)[local]);
                    found = true;
                    break;
                }
            if (!found)
                throw std::logic_error("Reference surface mapping failed");
        }
        if (n != 27 || boundary.size() != 26)
            throw std::invalid_argument("Unexpected qualification mesh dimensions");
        const auto b = boundary.size();
        const elements::ExactCondensedPellet exact(k, source, boundary);
        const elements::SurrogatePellet model(argv[2]);
        model.validate_mesh(ids,
            spatial.coordinates(),
            elements::pellet_mesh_signature(spatial.coordinates(),
                spatial.connectivity(),
                input.spatial.regions[0].material.functions->thermal.parameters.value("conductivity")));
        const double volume = std::accumulate(source.begin(), source.end(), 0.0);
        const auto& tangent = exact.stiffness();
        double tangent_norm = 0;
        for (double value : tangent)
            tangent_norm = std::hypot(tangent_norm, value);
        std::ofstream evidence(argv[3]);
        if (!evidence)
            throw std::runtime_error("Cannot write response evidence");
        evidence << std::setprecision(17) << "FUELSIM_PELLET_RESPONSE 1\n" << b << ' ' << volume << '\n';
        for (auto id : ids)
            evidence << id << ' ';
        evidence << '\n';
        for (double value : tangent)
            evidence << value << ' ';
        evidence << '\n';
        std::map<std::string, double> worst;
        const auto record = [&](const std::string& name, double value) {
            if (!std::isfinite(value))
                throw std::runtime_error("Nonfinite response metric: " + name);
            worst[name] = std::max(worst[name], std::abs(value));
        };
        std::mt19937_64 generator(20260920);
        std::uniform_real_distribution<double> unit(-1, 1);
        constexpr std::size_t samples = 256;
        for (std::size_t sample = 0; sample < samples; ++sample) {
            const double mean = 700 + 200 * unit(generator), a = 60 * unit(generator), c = 60 * unit(generator),
                         d = 60 * unit(generator), e = 30 * unit(generator), q = 2e8 + 1.5e8 * unit(generator);
            std::vector<double> t;
            for (auto id : ids) {
                const auto& p = mesh.nodes()[id];
                const double x = p.x / .004, y = p.y / .004, z = p.z / .008 - .5;
                double value = mean;
                if (sample % 4 == 1)
                    value += 2 * a * z;
                if (sample % 4 == 2)
                    value += c * x + d * y;
                if (sample % 4 == 3)
                    value += 2 * a * z + c * x + d * y + e * (x * x - y * y + 2 * x * z);
                t.push_back(value);
            }
            std::vector<double> reference, residual, jacobian;
            exact.evaluate(t, q, reference);
            model.evaluate_with_jacobian(t, q, residual, jacobian);
            evidence << sample << ' ' << q << '\n';
            for (double value : t)
                evidence << value << ' ';
            evidence << '\n';
            for (double value : reference)
                evidence << value << ' ';
            evidence << '\n';
            for (double value : residual)
                evidence << value << ' ';
            evidence << '\n';
            for (double value : jacobian)
                evidence << value << ' ';
            evidence << '\n';
            double norm = 0, error = 0, derivative_error = 0;
            for (std::size_t i = 0; i < b; ++i) {
                record("residual_max_error_W", residual[i] - reference[i]);
                norm = std::hypot(norm, reference[i]);
                error = std::hypot(error, residual[i] - reference[i]);
                double row_sum = 0, column_sum = 0;
                for (std::size_t j = 0; j < b; ++j) {
                    const auto ij = i * b + j;
                    record("tangent_max_error_W_K", jacobian[ij] - tangent[ij]);
                    derivative_error = std::hypot(derivative_error, jacobian[ij] - tangent[ij]);
                    row_sum += jacobian[ij];
                    column_sum += jacobian[j * b + i];
                    record("tangent_symmetry_error_W_K", jacobian[ij] - jacobian[j * b + i]);
                }
                record("tangent_uniform_shift_error_W_K", row_sum);
                record("tangent_energy_column_error_W_K", column_sum);
            }
            record("residual_relative_l2", error / norm);
            record("tangent_relative_frobenius", derivative_error / tangent_norm);
            record("energy_error_W", std::accumulate(residual.begin(), residual.end(), 0.0) + q * volume);
            record("exact_energy_error_W", std::accumulate(reference.begin(), reference.end(), 0.0) + q * volume);
            for (double shift : {-100.0, 100.0}) {
                auto shifted = t;
                for (auto& value : shifted)
                    value += shift;
                std::vector<double> shifted_r, shifted_exact;
                model.evaluate(shifted, q, shifted_r);
                exact.evaluate(shifted, q, shifted_exact);
                for (std::size_t i = 0; i < b; ++i) {
                    record("uniform_shift_residual_error_W", shifted_r[i] - residual[i]);
                    record("exact_uniform_shift_error_W", shifted_exact[i] - reference[i]);
                }
                record("shifted_energy_error_W", std::accumulate(shifted_r.begin(), shifted_r.end(), 0.0) + q * volume);
            }
        }
        for (double temperature : {400.0, 700.0, 1000.0}) {
            std::vector<double> r;
            model.evaluate(std::vector<double>(b, temperature), 0, r);
            for (double value : r)
                record("zero_source_uniform_residual_W", value);
            record("zero_source_uniform_energy_W", std::accumulate(r.begin(), r.end(), 0.0));
        }
        if (!evidence)
            throw std::runtime_error("Failed to write response evidence");
        const std::map<std::string, double> limits{{"residual_max_error_W", 1e-3},
            {"residual_relative_l2", 1e-4},
            {"tangent_max_error_W_K", 1e-5},
            {"tangent_relative_frobenius", 1e-4},
            {"tangent_symmetry_error_W_K", 1e-5},
            {"tangent_uniform_shift_error_W_K", 1e-5},
            {"tangent_energy_column_error_W_K", 1e-5},
            {"energy_error_W", 1e-3},
            {"exact_energy_error_W", 1e-10},
            {"uniform_shift_residual_error_W", 1e-3},
            {"exact_uniform_shift_error_W", 1e-10},
            {"shifted_energy_error_W", 1e-3},
            {"zero_source_uniform_residual_W", 1e-3},
            {"zero_source_uniform_energy_W", 1e-3}};
        bool passed = true;
        std::cout << std::setprecision(17) << "samples=" << samples << '\n';
        for (const auto& metric : worst) {
            std::cout << metric.first << '=' << metric.second << '\n';
            if (metric.second > limits.at(metric.first))
                passed = false;
        }
        std::cout << "response_check_passed=" << (passed ? "true" : "false") << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
