#include "io/case_input.hpp"
#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/exodus_result_reader.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
bool same_point(const fuelsim::CartesianMaterialPointState& left, const fuelsim::CartesianMaterialPointState& right) {
    return left.elastic_strain == right.elastic_strain && left.plastic_strain == right.plastic_strain
           && left.creep_strain == right.creep_strain
           && left.equivalent_plastic_strain == right.equivalent_plastic_strain
           && left.equivalent_creep_strain == right.equivalent_creep_strain && left.stress.xx == right.stress.xx
           && left.stress.yy == right.stress.yy && left.stress.zz == right.stress.zz
           && left.stress.xy == right.stress.xy && left.stress.yz == right.stress.yz
           && left.stress.xz == right.stress.xz;
}

bool same_committed_state(const fuelsim::TransientCommittedState& left, const fuelsim::TransientCommittedState& right) {
    if (left.time != right.time || left.load_factor != right.load_factor || left.solution != right.solution
        || left.cartesian_material_histories.size() != right.cartesian_material_histories.size())
        return false;
    for (std::size_t region = 0; region < left.cartesian_material_histories.size(); ++region) {
        if (left.cartesian_material_histories[region].size() != right.cartesian_material_histories[region].size())
            return false;
        for (std::size_t element = 0; element < left.cartesian_material_histories[region].size(); ++element)
            for (std::size_t q = 0; q < 8; ++q)
                if (!same_point(left.cartesian_material_histories[region][element][q],
                        right.cartesian_material_histories[region][element][q]))
                    return false;
    }
    return true;
}

void verify(const std::string& card,
    const std::string& checkpoint,
    const std::string& roundtrip,
    const std::string& output) {
    const auto definition = fuelsim::read_case_input(card);
    const auto mesh = fuelsim::read_exodus_hex8(definition.mesh_file);
    fuelsim::TransientProblem problem(definition.spatial, mesh);
    const auto initial = fuelsim::cartesian::ProblemAccess::committed_state(problem);
    problem.begin_time_step({0.1, 0.1});
    fuelsim::ContributionWorkspace workspace;
    problem.evaluate_contribution(0, problem.committed_solution(), workspace, true);
    problem.evaluate_contribution(0, problem.committed_solution(), workspace, false);
    if (!same_committed_state(initial, fuelsim::cartesian::ProblemAccess::committed_state(problem)))
        throw std::runtime_error("HEX8 residual/Jacobian changed committed history");
    problem.rollback_time_step();
    if (!same_committed_state(initial, fuelsim::cartesian::ProblemAccess::committed_state(problem)))
        throw std::runtime_error("HEX8 rollback changed committed history");
    const double step = fuelsim::restore_transient_checkpoint(checkpoint, problem);
    if (!std::isfinite(step) || !(step > 0.0) || !(problem.committed_time() > 0.9))
        throw std::runtime_error("HEX8 production checkpoint did not reach its active final state");
    const auto expected = fuelsim::cartesian::ProblemAccess::committed_state(problem);
    const auto result = fuelsim::test::read_final_exodus_results(output);
    if (result.time != expected.time)
        throw std::runtime_error("Output and checkpoint times differ");
    for (const auto& field : fuelsim::transient_conservation_fields)
        if (result.global("conservation_" + std::string(field.name)) != expected.conservation.*(field.member))
            throw std::runtime_error("Output conservation diagnostic differs from committed checkpoint");
    std::vector<double> constraints(problem.dof_count(), 0.0);
    for (const auto& condition : problem.dirichlet_conditions())
        constraints.at(condition.dof) = 1.0;
    const auto& spatial = fuelsim::cartesian::ProblemAccess::view(problem);
    const std::array<std::string, 4> field_names = {"temperature",
        "displacement_x",
        "displacement_y",
        "displacement_z"};
    const std::array<std::string, 4> reaction_names = {"reaction_heat_flux",
        "reaction_force_x",
        "reaction_force_y",
        "reaction_force_z"};
    const std::array<std::string, 6> components = {"xx", "yy", "zz", "xy", "yz", "xz"};
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const auto& region_mesh = spatial.region_mesh(region);
        for (std::size_t local = 0; local < region_mesh.nodes().size(); ++local) {
            const std::size_t source = region_mesh.source_node_ids()[local],
                              global = spatial.global_node(region, local);
            for (std::size_t field = 0; field < 4; ++field) {
                const std::size_t dof = spatial.field_layout()[field].begin + global;
                if (result.nodal(field_names[field]).at(source) != expected.solution.at(dof)
                    || result.nodal(reaction_names[field]).at(source) != expected.raw_residual.at(dof)
                    || result.nodal("dirichlet_" + field_names[field]).at(source) != constraints.at(dof))
                    throw std::runtime_error("Output nodal field or pre-commit reaction differs from checkpoint");
            }
        }
        for (std::size_t element = 0; element < region_mesh.elements().size(); ++element) {
            const auto source = region_mesh.source_element_ids()[element];
            for (std::size_t q = 0; q < 8; ++q) {
                const auto& point = expected.cartesian_material_histories[region][element][q];
                for (std::size_t c = 0; c < 6; ++c) {
                    const auto suffix = components[c] + "_q" + std::to_string(q);
                    if (result.element("elastic_" + suffix).at(source) != point.elastic_strain[c]
                        || result.element("plastic_" + suffix).at(source) != point.plastic_strain[c]
                        || result.element("creep_" + suffix).at(source) != point.creep_strain[c])
                        throw std::runtime_error("Output six-component material tensors differ from checkpoint");
                }
            }
        }
    }
    bool active = false;
    for (const auto& region : expected.cartesian_material_histories)
        for (const auto& element : region)
            for (const auto& point : element)
                active = active || point.equivalent_plastic_strain > 0.0 || point.equivalent_creep_strain > 0.0;
    if (!active)
        throw std::runtime_error("HEX8 checkpoint history is not active");
    fuelsim::write_transient_checkpoint(roundtrip, problem, step);
    fuelsim::TransientProblem restored(definition.spatial, mesh);
    const double restored_step = fuelsim::restore_transient_checkpoint(roundtrip, restored);
    if (restored_step != step
        || !same_committed_state(expected, fuelsim::cartesian::ProblemAccess::committed_state(restored)))
        throw std::runtime_error("HEX8 checkpoint changed a complete six-component material history");
    if (std::remove(roundtrip.c_str()) != 0)
        throw std::runtime_error("Could not remove test checkpoint roundtrip");
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5)
        return 2;
    try {
        verify(argv[1], argv[2], argv[3], argv[4]);
        std::cout << "[PASS] HEX8 history immutability, rollback and active checkpoint roundtrip\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
