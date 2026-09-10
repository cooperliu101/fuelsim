#include "cax8rt.hpp"
#include "detail/quad8_rz_assembly.hpp"
#include <stdexcept>

namespace fuelsim::elements {
Cax8Result evaluate_cax8rt(const Cax8Input& input, ElementRequest request) {
    if (input.geometry.point_count != 4)
        throw std::invalid_argument("CAX8RT requires its model-specific material quadrature");
    const auto& data = input;
    auto result = compute_quad8_rz(data,
        input.geometry,
        input.state,
        input.committed_state,
        input.committed_history,
        input.time_step,
        request.jacobian,
        input.include_thermal_time_term);
    if (request.stress)
        for (std::size_t q = 0; q < result.history.size(); ++q)
            result.stress[q] = result.history[q].stress;
    if (!request.residual && !request.jacobian)
        result.residual.fill(0.0);
    if (!request.history)
        result.history = {};
    return result;
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
Quad8RzGeometry make_cax8rt_geometry(const Quad8RzCoordinates& coordinates) {
    return make_quad8_rz_geometry(coordinates, 2);
}
} // namespace fuelsim::elements
