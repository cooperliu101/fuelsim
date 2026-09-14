#pragma once
#include "core/cross_section.hpp"
#include <map>
#include <utility>

namespace fuelsim::test {
inline SectionRegion region(const std::string& name, double young, double poisson) {
    const auto registry = make_builtin_material_function_registry();
    auto functions = std::make_shared<MaterialFunctionSet>();
    functions->name = name;
    functions->thermal = registry.bind_thermal("inverse_temperature_thermophysical",
        {{"conductivity_inverse_temperature", 0.0},
            {"conductivity_constant", 1.0},
            {"density", 1.0},
            {"specific_heat", 1.0}});
    functions->elasticity = registry.bind_elasticity("linear_temperature_isotropic",
        {{"young_modulus", young},
            {"poisson_ratio", poisson},
            {"reference_temperature", 300.0},
            {"young_modulus_temperature_coefficient", 0.0},
            {"poisson_ratio_temperature_coefficient", 0.0}});
    return {name, IsotropicThermoelasticMaterial({functions, young}), 300.0};
}

inline CrossSection rectangle(const std::vector<double>& x,
    double width,
    std::vector<SectionRegion> regions,
    bool skew = false,
    int ny = 1) {
    std::vector<CartesianPoint3> nodes;
    std::vector<SectionCell> cells;
    std::map<std::pair<int, int>, std::size_t> ids;
    const std::array<std::array<int, 2>, 8> offsets{{{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 0}, {2, 1}, {1, 2}, {0, 1}}};
    for (std::size_t cell = 0; cell + 1 < x.size(); ++cell) {
        for (int row = 0; row < ny; ++row) {
            SectionCell element{};
            element.region = regions.size() == 1 ? 0 : cell;
            for (std::size_t i = 0; i < 8; ++i) {
                const auto key = std::make_pair(2 * static_cast<int>(cell) + offsets[i][0], 2 * row + offsets[i][1]);
                auto found = ids.find(key);
                if (found == ids.end()) {
                    const double px = x[cell] + 0.5 * static_cast<double>(offsets[i][0]) * (x[cell + 1] - x[cell]);
                    const double py =
                        -width / 2.0
                        + 0.5 * static_cast<double>(2 * row + offsets[i][1]) * width / static_cast<double>(ny);
                    found = ids.emplace(key, nodes.size()).first;
                    nodes.push_back({px + (skew ? 0.2 * py : 0.0), py, 0.0});
                }
                element.element.nodes[i] = found->second;
            }
            cells.push_back(element);
        }
    }
    return CrossSection(std::move(nodes), std::move(cells), std::move(regions));
}
} // namespace fuelsim::test
