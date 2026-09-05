#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/exodus_result_reader.hpp"
#include "support/rz_problem_access.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

// Initial candidate ownership is an internal contract; final ownership comes from the production run.
int main(int argc, char** argv) {
    if (argc != 3) return 2;
    try {
        const auto definition = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_quad4(definition.mesh_file);
        fuelsim::SteadyProblem problem(definition.spatial, mesh);
        const auto initial = fuelsim::rz::ProblemAccess::summarize_contact_nodes(problem, 0, problem.initial_state());
        const auto sources = fuelsim::rz::ProblemAccess::contact_secondary_source_nodes(problem, 0);
        const auto output = fuelsim::test::read_final_exodus_results(argv[2]);
        const auto& segments = output.nodal("contact_primary_segment_pellet_stack");
        std::size_t maximum_change = 0;
        for (std::size_t node = 0; node < initial.size(); ++node) {
            const double final_segment = segments.at(sources.at(node));
            if (!initial[node].projected || !std::isfinite(final_segment) || final_segment < 0.0 ||
                std::floor(final_segment) != final_segment)
                throw std::runtime_error("M5.2 initial or final candidate ownership is invalid");
            const auto final_index = static_cast<std::size_t>(final_segment);
            if (final_index >= initial[node].primary_segment)
                maximum_change = std::max(maximum_change, final_index - initial[node].primary_segment);
        }
        std::cout << "m52_maximum_primary_segment_change=" << maximum_change << '\n';
        if (initial.size() != 5 || maximum_change < 2)
            throw std::runtime_error("M5.2 five-node contact did not cross two primary segments");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
