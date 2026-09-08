#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool histories_equal(const std::vector<fuelsim::ContactPointHistory>& first,
    const std::vector<fuelsim::ContactPointHistory>& second) {
    if (first.size() != second.size())
        return false;
    for (std::size_t point = 0; point < first.size(); ++point)
        if (first[point].elastic_tangential_slip != second[point].elastic_tangential_slip
            || first[point].sliding != second[point].sliding
            || first[point].normal_multiplier != second[point].normal_multiplier
            || first[point].cartesian_elastic_tangential_slip != second[point].cartesian_elastic_tangential_slip
            || first[point].cartesian_total_tangential_slip != second[point].cartesian_total_tangential_slip
            || first[point].cartesian_tangent_basis_initialized != second[point].cartesian_tangent_basis_initialized
            || first[point].cartesian_contact_normal != second[point].cartesian_contact_normal
            || first[point].cartesian_contact_tangent_first != second[point].cartesian_contact_tangent_first)
            return false;
    return true;
}

double jacobian_error(const fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    double maximum = 0.0;
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != fuelsim::SpatialContributionType::mechanical_contact)
            continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index)
            local[index] = state[dofs[index]];
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

bool run(const std::string& input_path,
    const std::string& full_checkpoint,
    const std::string& restarted_checkpoint,
    const std::string& split_checkpoint) {
    const auto input = fuelsim::read_case_input(input_path);
    const auto generated = fuelsim::read_exodus_hex8(input.mesh_file);
    fuelsim::TransientProblem full(input.spatial, generated), restarted(input.spatial, generated),
        split(input.spatial, generated);
    const double full_step = fuelsim::restore_transient_checkpoint(full_checkpoint, full);
    const double restart_step = fuelsim::restore_transient_checkpoint(restarted_checkpoint, restarted);
    const double split_step = fuelsim::restore_transient_checkpoint(split_checkpoint, split);
    bool passed =
        check(full.committed_time() == 4.0 && restarted.committed_time() == 4.0 && split.committed_time() == 2.0
                  && full_step == 1.0 && restart_step == 1.0 && split_step == 1.0,
            "B4.0 checkpoints preserve the exact physical time and next time step");
    passed = check(full.committed_solution() == restarted.committed_solution()
                       && histories_equal(fuelsim::cartesian::ProblemAccess::committed_contact_histories(full).at(0),
                           fuelsim::cartesian::ProblemAccess::committed_contact_histories(restarted).at(0)),
                 "B4.0 production restart exactly reproduces every nodal and contact history component")
             && passed;
    const auto& split_history = fuelsim::cartesian::ProblemAccess::committed_contact_histories(split).at(0);
    passed = check(split_history.size() == 4, "B4.0 checkpoint contains four friction histories") && passed;
    for (const auto& point : split_history)
        passed = check(point.sliding && std::abs(point.cartesian_elastic_tangential_slip[1]) > 0.0
                           && std::abs(point.cartesian_elastic_tangential_slip[2]) > 0.0,
                     "B4.0 restart begins from sliding with two nonzero elastic-slip components")
                 && passed;
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
    std::cout << std::scientific << std::setprecision(12)
              << "b40_sticking_jacobian_directional_error=" << sticking_jacobian_error << '\n'
              << "b40_sliding_jacobian_directional_error=" << sliding_jacobian_error << '\n';
    return check(sticking_jacobian_error < 1.0e-7 && sliding_jacobian_error < 1.0e-7,
               "B4.0 sticking and sliding Jacobians match centered directional differences")
           && passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5)
        return 2;
    try {
        return run(argv[1], argv[2], argv[3], argv[4]) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
