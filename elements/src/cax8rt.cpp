#include "cax8rt.hpp"
#include "cax8t.hpp"
#include <stdexcept>

namespace fuelsim::elements {
Cax8Result evaluate_cax8rt(const Cax8Input& input, ElementRequest request) {
    return evaluate_cax8t(input, request, Cax8Quadrature::reduced);
}

Quad8RzGeometry make_cax8rt_geometry(const Quad8RzCoordinates& coordinates) {
    return make_cax8t_geometry(coordinates, Cax8Quadrature::reduced);
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
std::vector<double> cax8rt_creep_rates(const IsotropicThermoelasticMaterial& material,
    const Quad8RzGeometry& geometry,
    const Quad8RzValues& state,
    const Quad8MaterialHistory& history,
    double time) {
    return cax8t_creep_rates(material, geometry, state, history, time, Cax8Quadrature::reduced);
}
} // namespace fuelsim::elements
