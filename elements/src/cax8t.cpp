#include "cax8t.hpp"
#include "detail/quad8_rz_assembly.hpp"
#include <stdexcept>

namespace fuelsim::elements {
Cax8Result evaluate_cax8t(const Cax8Input& input, ElementRequest request) {
    if (input.geometry.point_count != 9)
        throw std::invalid_argument("CAX8T requires its model-specific material quadrature");
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
Quad8RzGeometry make_cax8t_geometry(const Quad8RzCoordinates& coordinates) {
    return make_quad8_rz_geometry(coordinates, 3);
}
} // namespace fuelsim::elements
