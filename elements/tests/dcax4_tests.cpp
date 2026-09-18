#include "dcax4.hpp"
#include "thermal_test_support.hpp"

int main() {
    try {
        const std::array<fuelsim::CartesianPoint3, 4> coordinates{{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}}};
        check_thermal_element(fuelsim::elements::make_dcax4_geometry(coordinates),
            fuelsim::elements::evaluate_dcax4,
            {coordinates.begin(), coordinates.end()},
            true);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
