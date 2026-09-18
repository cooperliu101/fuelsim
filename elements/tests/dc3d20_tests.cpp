#include "dc3d20.hpp"
#include "thermal_test_support.hpp"

int main() {
    try {
        const std::array<fuelsim::CartesianPoint3, 20> coordinates{{{1, 0, 0},
            {2, 0, 0},
            {2, 1, 0},
            {1, 1, 0},
            {1, 0, 1},
            {2, 0, 1},
            {2, 1, 1},
            {1, 1, 1},
            {1.5, 0.0, 0.0},
            {2.0, 0.5, 0.0},
            {1.5, 1.0, 0.0},
            {1.0, 0.5, 0.0},
            {1.0, 0.0, 0.5},
            {2.0, 0.0, 0.5},
            {2.0, 1.0, 0.5},
            {1.0, 1.0, 0.5},
            {1.5, 0.0, 1.0},
            {2.0, 0.5, 1.0},
            {1.5, 1.0, 1.0},
            {1.0, 0.5, 1.0}}};
        check_thermal_element(fuelsim::elements::make_dc3d20_geometry(coordinates),
            fuelsim::elements::evaluate_dc3d20,
            {coordinates.begin(), coordinates.end()},
            false);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
