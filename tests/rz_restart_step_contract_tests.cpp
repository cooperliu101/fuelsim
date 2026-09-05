#include "fuelsim/io/case_input.hpp"
#include "fuelsim/io/checkpoint.hpp"
#include "fuelsim/io/results_io.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    try {
        const auto input = fuelsim::read_case_input(argv[1]);
        const auto mesh = fuelsim::read_exodus_quad4(input.mesh_file);
        fuelsim::TransientProblem problem(input.spatial, mesh);
        const double restored_step = fuelsim::restore_transient_checkpoint(argv[2], problem);
        std::cout << std::hexfloat << "restored_step=" << restored_step
                  << " committed_time=" << problem.committed_time() << '\n';
        // Ten additions of 0.05 stop one representable interval below 0.5.
        // The saved controller step retains the original exact comparison.
        if (restored_step != 0.05 ||
            std::abs(problem.committed_time() - 0.5) > 0.5 * std::numeric_limits<double>::epsilon())
            throw std::runtime_error("M5.2 checkpoint must restore the prescribed controller step and split time");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
