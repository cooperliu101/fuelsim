#include "dcax8.hpp"
#include "thermal_test_support.hpp"

int main() {
    try {
        const std::array<fuelsim::CartesianPoint3, 8> coordinates{{{1, 0, 0},
            {2, 0, 0},
            {2, 1, 0},
            {1, 1, 0},
            {1.5, 0.0, 0.0},
            {2.0, 0.5, 0.0},
            {1.5, 1.0, 0.0},
            {1.0, 0.5, 0.0}}};
        check_thermal_element(fuelsim::elements::make_dcax8_geometry(coordinates),
            fuelsim::elements::evaluate_dcax8,
            {coordinates.begin(), coordinates.end()},
            true);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
