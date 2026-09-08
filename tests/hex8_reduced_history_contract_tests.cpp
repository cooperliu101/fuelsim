#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/exodus_result_reader.hpp"
#include <array>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 4)
        return 2;
    try {
        const auto input = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_hex8(input.mesh_file);
        fuelsim::TransientProblem problem(input.spatial, mesh);
        for (std::size_t element = 0; element < 2; ++element)
            if (fuelsim::cartesian::ProblemAccess::material_history(problem, 0, element).size() != 1)
                throw std::runtime_error("C3D8RT must initialize exactly one material point per element");
        if (fuelsim::restore_transient_checkpoint(argv[2], problem) != 1.0 || problem.committed_time() != 1.0)
            throw std::runtime_error("C3D8RT checkpoint must preserve its final time and controller step");
        const auto output = fuelsim::test::read_final_exodus_results(argv[3]);
        for (std::size_t element = 0; element < 2; ++element) {
            const auto& history = fuelsim::cartesian::ProblemAccess::material_history(problem, 0, element);
            if (history.size() != 1)
                throw std::runtime_error("C3D8RT must commit exactly one material point per element");
            const auto& stress = history[0].stress;
            const std::array<double, 6> values = {stress.xx, stress.yy, stress.zz, stress.xy, stress.yz, stress.xz};
            const std::array<std::string, 6> names = {"xx", "yy", "zz", "xy", "yz", "xz"};
            for (std::size_t q = 0; q < 8; ++q)
                for (std::size_t field = 0; field < 6; ++field)
                    if (output.element("stress_" + names[field] + "_q" + std::to_string(q)).at(element)
                        != values[field])
                        throw std::runtime_error("C3D8RT output must repeat its single committed point exactly");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
