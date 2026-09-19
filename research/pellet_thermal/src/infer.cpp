#include "pellet_mlp.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("Usage: fuelsim_pellet_infer model.txt inputs.txt");
        const fuelsim::elements::SurrogatePellet model(argv[1]);
        std::ifstream input(argv[2]);
        std::size_t count = 0, n = 0;
        if (!(input >> count >> n) || count == 0 || n == 0 || n > 4096)
            throw std::invalid_argument("Invalid inference input dimensions");
        std::cout << std::setprecision(17);
        for (std::size_t sample = 0; sample < count; ++sample) {
            std::vector<double> temperature(n), residual, jacobian, plain;
            double source = 0;
            for (double& t : temperature)
                if (!(input >> t))
                    throw std::invalid_argument("Missing inference temperature");
            if (!(input >> source))
                throw std::invalid_argument("Missing inference heat source");
            model.evaluate_with_jacobian(temperature, source, residual, jacobian);
            model.evaluate(temperature, source, plain);
            if (plain != residual)
                throw std::logic_error("Inference residual paths differ");
            for (double v : residual)
                std::cout << v << ' ';
            for (double v : jacobian)
                std::cout << v << ' ';
            std::cout << '\n';
        }
        std::string extra;
        if (input >> extra)
            throw std::invalid_argument("Unexpected inference input");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
