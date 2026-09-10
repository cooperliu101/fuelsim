#include "line3_rz_boundary.hpp"
#include "boundary_types.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::elements {
Line3RzBoundaryResult compute_line3_rz_boundary(const Line3RzBoundaryData& data,
    const std::array<RzPoint, 3>& coordinates,
    const std::vector<double>& state,
    bool jacobian) {
    if (state.size() != 8)
        throw std::invalid_argument("CAX8T quadratic boundary requires eight local degrees of freedom");
    std::array<adlite::Scalar, 8> v;
    for (std::size_t i = 0; i < 8; ++i)
        v[i] = jacobian ? adlite::Scalar::independent(state[i], i, 8) : adlite::Scalar(state[i]);
    std::array<adlite::Scalar, 8> rows{};
    const bool current = data.use_displaced_geometry;
    const double load = data.load;
    const double g = std::sqrt(3.0 / 5.0);
    const std::array<double, 3> points = {-g, 0, g}, weights = {5.0 / 9, 8.0 / 9, 5.0 / 9};
    for (std::size_t q = 0; q < 3; ++q) {
        const double x = points[q];
        const std::array<double, 3> shape = {x * (x - 1) / 2, x * (x + 1) / 2, 1 - x * x},
                                    derivative = {x - .5, x + .5, -2 * x};
        const std::array<double, 2> thermal = {(1 - x) / 2, (1 + x) / 2};
        adlite::Scalar radius = 0, tr = 0, tz = 0;
        for (std::size_t n = 0; n < 3; ++n) {
            const auto& point = coordinates[n];
            const adlite::Scalar r = current ? point.r + v[2 + n] : adlite::Scalar(point.r),
                                 z = current ? point.z + v[5 + n] : adlite::Scalar(point.z);
            radius += shape[n] * r;
            tr += derivative[n] * r;
            tz += derivative[n] * z;
        }
        const auto length = adlite::hypot(tr, tz), factor = 2 * std::acos(-1.0) * weights[q] * radius,
                   measure = factor * length;
        if (!(length.value() > 0) || !(radius.value() >= 0))
            throw std::domain_error("CAX8T boundary geometry is invalid");
        if (data.kind == Line3RzBoundaryKind::pressure)
            for (std::size_t n = 0; n < 3; ++n) {
                rows[2 + n] += shape[n] * load * factor * tz;
                rows[5 + n] -= shape[n] * load * factor * tr;
            }
        else if (data.kind == Line3RzBoundaryKind::traction) {
            const auto offset = data.component == TractionComponent::radial ? 2U : 5U;
            for (std::size_t n = 0; n < 3; ++n)
                rows[offset + n] -= shape[n] * load * measure;
        } else {
            adlite::Scalar flux = -load;
            if (data.kind == Line3RzBoundaryKind::convection) {
                flux = load * (thermal[0] * v[0] + thermal[1] * v[1] - data.ambient);
            }
            for (std::size_t n = 0; n < 2; ++n)
                rows[n] += thermal[n] * flux * measure;
        }
    }
    Line3RzBoundaryResult result;
    for (std::size_t i = 0; i < 8; ++i) {
        result.residual[i] = rows[i].value();
        if (jacobian)
            rows[i].copy_derivatives(result.jacobian.data() + 8 * i, 8);
    }
    return result;
}
} // namespace fuelsim::elements
