#include "contact_types.hpp"
#include "io/case_input.hpp"
#include "io/checkpoint.hpp"
#include "io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/exodus_result_reader.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
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

} // namespace

int main(int argc, char** argv) {
    if (argc != 5)
        return 2;
    try {
        const auto definition = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_hex8(definition.mesh_file);
        fuelsim::TransientProblem original(definition.spatial, mesh);
        const double step = fuelsim::restore_transient_checkpoint(argv[2], original);
        if (step != 1.0 || original.committed_time() != 1.0)
            throw std::runtime_error("B3.4 production checkpoint did not commit one physical time step");
        const auto output = fuelsim::test::read_final_exodus_results(argv[3]);
        const auto& spatial = fuelsim::cartesian::ProblemAccess::view(original);
        const std::array<std::string, 4> fields = {"temperature", "displacement_x", "displacement_y", "displacement_z"};
        for (std::size_t region = 0; region < spatial.region_count(); ++region) {
            const auto& local_mesh = spatial.region_mesh(region);
            for (std::size_t local = 0; local < local_mesh.nodes().size(); ++local)
                for (std::size_t field = 0; field < 4; ++field) {
                    const auto dof = spatial.field_layout()[field].begin + spatial.global_node(region, local);
                    if (output.nodal(fields[field]).at(local_mesh.source_node_ids()[local])
                        != original.committed_solution()[dof])
                        throw std::runtime_error("B3.4 output and committed nodal state differ");
                }
        }
        const auto& history = fuelsim::cartesian::ProblemAccess::committed_contact_histories(original).at(0);
        const auto source = fuelsim::cartesian::ProblemAccess::contact_secondary_source_nodes(original, 0);
        if (history.size() != 4 || source.size() != 4)
            throw std::runtime_error("B3.4 requires four contact history slots");
        bool sliding = false;
        for (std::size_t node = 0; node < history.size(); ++node) {
            const auto& point = history[node];
            sliding = sliding || point.sliding;
            if (output.nodal("contact_sliding_interface").at(source[node]) != (point.sliding ? 1.0 : 0.0))
                throw std::runtime_error("B3.4 output sliding state differs from checkpoint");
            for (std::size_t component = 0; component < 3; ++component) {
                const std::string suffix = std::string(1, "xyz"[component]) + "_interface";
                std::cout << std::setprecision(17) << "node=" << node << " component=" << component
                          << " output_total=" << output.nodal("contact_total_slip_" + suffix).at(source[node])
                          << " history_total=" << point.cartesian_total_tangential_slip[component]
                          << " output_elastic=" << output.nodal("contact_elastic_slip_" + suffix).at(source[node])
                          << " history_elastic=" << point.cartesian_elastic_tangential_slip[component] << '\n';
                const double actual_elastic = output.nodal("contact_elastic_slip_" + suffix).at(source[node]);
                const double committed_elastic = point.cartesian_elastic_tangential_slip[component];
                // Summary output reevaluates the Coulomb return map at the
                // committed state; this can round its elastic slip by one ULP.
                // The checkpoint roundtrip below remains exactly equal.
                const double roundoff = 4.0 * std::numeric_limits<double>::epsilon()
                                        * std::max(std::abs(actual_elastic), std::abs(committed_elastic));
                if (output.nodal("contact_total_slip_" + suffix).at(source[node])
                        != point.cartesian_total_tangential_slip[component]
                    || !std::isfinite(actual_elastic) || std::abs(actual_elastic - committed_elastic) > roundoff)
                    throw std::runtime_error("B3.4 output friction vectors differ from checkpoint");
            }
        }
        if (!sliding)
            throw std::runtime_error("B3.4 checkpoint must contain active sliding history");
        fuelsim::write_transient_checkpoint(argv[4], original, 0.25);
        fuelsim::TransientProblem restored(definition.spatial, mesh);
        const double restored_step = fuelsim::restore_transient_checkpoint(argv[4], restored);
        if (restored_step != 0.25 || restored.committed_solution() != original.committed_solution()
            || !histories_equal(history,
                fuelsim::cartesian::ProblemAccess::committed_contact_histories(restored).at(0)))
            throw std::runtime_error("B3.4 checkpoint changed the exact three-dimensional friction transaction");
        if (std::remove(argv[4]) != 0)
            throw std::runtime_error("Could not remove B3.4 roundtrip checkpoint");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
