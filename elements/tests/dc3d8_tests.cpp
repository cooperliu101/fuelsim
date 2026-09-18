#include "dc3d8.hpp"
#include "thermal_test_support.hpp"

int main() {
    try {
        const std::array<fuelsim::CartesianPoint3, 8> coordinates{
            {{1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {1, 0, 1}, {2, 0, 1}, {2, 1, 1}, {1, 1, 1}}};
        check_thermal_element(fuelsim::elements::make_dc3d8_geometry(coordinates),
            fuelsim::elements::evaluate_dc3d8,
            {coordinates.begin(), coordinates.end()},
            false);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
