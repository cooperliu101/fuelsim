#include "io/case_input.hpp"
#include "io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
bool same_point(const fuelsim::CartesianMaterialPointState& a, const fuelsim::CartesianMaterialPointState& b) {
    return a.elastic_strain == b.elastic_strain && a.plastic_strain == b.plastic_strain
           && a.creep_strain == b.creep_strain && a.equivalent_plastic_strain == b.equivalent_plastic_strain
           && a.equivalent_creep_strain == b.equivalent_creep_strain && a.stress.xx == b.stress.xx
           && a.stress.yy == b.stress.yy && a.stress.zz == b.stress.zz && a.stress.xy == b.stress.xy
           && a.stress.yz == b.stress.yz && a.stress.xz == b.stress.xz;
}

void require_same(const fuelsim::TransientCommittedState& a, const fuelsim::TransientCommittedState& b) {
    if (a.time != b.time || a.load_factor != b.load_factor || a.solution != b.solution
        || a.cartesian_material_histories.size() != b.cartesian_material_histories.size())
        throw std::runtime_error("HEX20 evaluation or rollback changed committed state");
    for (std::size_t region = 0; region < a.cartesian_material_histories.size(); ++region) {
        if (a.cartesian_material_histories[region].size() != b.cartesian_material_histories[region].size())
            throw std::runtime_error("HEX20 evaluation changed material history dimensions");
        for (std::size_t element = 0; element < a.cartesian_material_histories[region].size(); ++element)
            for (std::size_t q = 0; q < 27; ++q)
                if (!same_point(a.cartesian_material_histories[region][element][q],
                        b.cartesian_material_histories[region][element][q]))
                    throw std::runtime_error("HEX20 evaluation or rollback changed a committed material point");
    }
}

void check_history(const std::string& path) {
    const auto definition = fuelsim::read_case_input(path);
    const auto mesh = fuelsim::read_exodus_hex20(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, mesh);
    const auto initial = fuelsim::cartesian::ProblemAccess::committed_state(problem);
    problem.begin_time_step({0.1, 0.1});
    fuelsim::ContributionWorkspace workspace;
    problem.evaluate_contribution(0, problem.committed_solution(), workspace, true);
    problem.evaluate_contribution(0, problem.committed_solution(), workspace, false);
    require_same(initial, fuelsim::cartesian::ProblemAccess::committed_state(problem));
    problem.rollback_time_step();
    require_same(initial, fuelsim::cartesian::ProblemAccess::committed_state(problem));
    std::cout << "[PASS] HEX20 residual/Jacobian immutability and rollback: " << path << '\n';
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 4)
        return 2;
    try {
        for (int argument = 1; argument < argc; ++argument)
            check_history(argv[argument]);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
