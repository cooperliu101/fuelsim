#include "c3d20t.hpp"
#include "detail/hex20_assembly.hpp"
#include <stdexcept>

namespace fuelsim::elements {
C3d20Result evaluate_c3d20t(const C3d20Input& input, ElementRequest request) {
    if (input.geometry.mechanical_points.size() != 27)
        throw std::invalid_argument("C3D20T requires its model-specific material quadrature");
    const auto& data = input;
    C3d20Result result;
    if (request.residual || request.jacobian) {
        auto* tangent = request.jacobian ? &result.jacobian : nullptr;
        if (input.committed_history)
            result.residual = compute_hex20_transient(data,
                input.geometry,
                input.state,
                input.committed_state,
                *input.committed_history,
                input.time_step,
                tangent,
                input.include_thermal_time_term);
        else
            result.residual = compute_hex20_thermoelastic(data,
                input.geometry,
                input.state,
                input.time_step > 0 ? &input.committed_state : nullptr,
                input.time_step,
                tangent,
                input.include_thermal_time_term);
    }
    if (request.history && input.committed_history)
        result.history = compute_hex20_transient_update(data,
            input.geometry,
            input.state,
            input.committed_state,
            *input.committed_history,
            input.time_step);
    if (request.stress)
        result.stress = compute_hex20_stress(data, input.geometry, input.state);
    return result;
}
} // namespace fuelsim::elements

namespace fuelsim::elements {
Hex20Geometry make_c3d20t_geometry(const Hex20Coordinates& coordinates) {
    return make_hex20_geometry(coordinates, 3);
}
} // namespace fuelsim::elements
