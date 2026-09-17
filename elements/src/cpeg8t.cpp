#include "cpeg8t.hpp"
#include "c3d_common.hpp"
#include "quad8_shape.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::elements {
namespace {
// Four in-plane gradients, two point displacements, three shared controls,
// and temperature. No element-wide identity seeding is used.
using PointValues = std::array<double, 10>;
using ActivePointValues = std::array<adlite::Scalar, 10>;

PointValues point_values(const Cpeg8Point& point, const Cpeg8Values& state) {
    PointValues result{};
    for (std::size_t n = 0; n < 8; ++n) {
        result[0] += point.gradient_x[n] * state[4 + n];
        result[1] += point.gradient_y[n] * state[4 + n];
        result[2] += point.gradient_x[n] * state[12 + n];
        result[3] += point.gradient_y[n] * state[12 + n];
        result[4] += point.shape[n] * state[4 + n];
        result[5] += point.shape[n] * state[12 + n];
    }
    for (std::size_t i = 0; i < 3; ++i)
        result[6 + i] = state[20 + i];
    for (std::size_t n = 0; n < 4; ++n)
        result[9] += point.temperature_shape[n] * state[n];
    return result;
}

void add_row(Cpeg8Result& result,
    std::size_t row,
    const adlite::Scalar& value,
    const Cpeg8Point& point,
    bool jacobian) {
    result.residual[row] += value.value();
    if (!jacobian)
        return;
    PointValues derivative{};
    value.copy_derivatives(derivative.data(), derivative.size());
    for (std::size_t n = 0; n < 4; ++n)
        result.jacobian[23 * row + n] += derivative[9] * point.temperature_shape[n];
    for (std::size_t n = 0; n < 8; ++n) {
        result.jacobian[23 * row + 4 + n] +=
            derivative[0] * point.gradient_x[n] + derivative[1] * point.gradient_y[n] + derivative[4] * point.shape[n];
        result.jacobian[23 * row + 12 + n] +=
            derivative[2] * point.gradient_x[n] + derivative[3] * point.gradient_y[n] + derivative[5] * point.shape[n];
    }
    for (std::size_t i = 0; i < 3; ++i)
        result.jacobian[23 * row + 20 + i] += derivative[6 + i];
}

std::array<double, 6> tensor_values(const SymmetricTensor3& tensor) {
    return {tensor.xx.value(),
        tensor.yy.value(),
        tensor.zz.value(),
        tensor.xy.value(),
        tensor.yz.value(),
        tensor.xz.value()};
}
} // namespace

Cpeg8Geometry
make_cpeg8t_geometry(const Cpeg8Coordinates& coordinates, double thickness, std::array<double, 2> reference_point) {
    if (!std::isfinite(thickness) || !(thickness > 0.0))
        throw std::invalid_argument("CPEG8T initial thickness must be finite and positive");
    for (const auto& node : coordinates)
        for (double coordinate : node)
            if (!std::isfinite(coordinate))
                throw std::invalid_argument("CPEG8T coordinates must be finite");
    for (double coordinate : reference_point)
        if (!std::isfinite(coordinate))
            throw std::invalid_argument("CPEG8T reference point must be finite");
    Cpeg8Geometry result{coordinates, {}, reference_point, thickness};
    const double g = std::sqrt(3.0 / 5.0);
    const std::array<double, 3> locations{-g, 0.0, g}, weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    constexpr std::array<std::array<double, 2>, 4> signs{{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
    for (std::size_t j = 0; j < 3; ++j)
        for (std::size_t i = 0; i < 3; ++i) {
            const double xi = locations[i], eta = locations[j];
            auto& point = result.points[3 * j + i];
            quad8_face_detail::DoubleQuad8ShapeValues shape;
            quad8_face_detail::double_quad8_shape(xi, eta, shape);
            point.shape = shape.shape;
            double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
            for (std::size_t n = 0; n < 8; ++n) {
                a += shape.derivative_xi[n] * coordinates[n][0];
                b += shape.derivative_eta[n] * coordinates[n][0];
                c += shape.derivative_xi[n] * coordinates[n][1];
                d += shape.derivative_eta[n] * coordinates[n][1];
                point.x += shape.shape[n] * coordinates[n][0];
                point.y += shape.shape[n] * coordinates[n][1];
            }
            const double determinant = a * d - b * c;
            if (!std::isfinite(determinant) || !(determinant > 0.0))
                throw std::domain_error("CPEG8T reference mapping must have positive determinant");
            point.measure = determinant * weights[i] * weights[j] * thickness;
            for (std::size_t n = 0; n < 8; ++n) {
                point.gradient_x[n] = (d * shape.derivative_xi[n] - c * shape.derivative_eta[n]) / determinant;
                point.gradient_y[n] = (a * shape.derivative_eta[n] - b * shape.derivative_xi[n]) / determinant;
            }
            for (std::size_t n = 0; n < 4; ++n) {
                const double sx = signs[n][0], sy = signs[n][1];
                const double dx = 0.25 * sx * (1.0 + sy * eta), dy = 0.25 * sy * (1.0 + sx * xi);
                point.temperature_shape[n] = 0.25 * (1.0 + sx * xi) * (1.0 + sy * eta);
                point.temperature_gradient_x[n] = (d * dx - c * dy) / determinant;
                point.temperature_gradient_y[n] = (a * dy - b * dx) / determinant;
            }
        }
    return result;
}

Cpeg8Result evaluate_cpeg8t(const Cpeg8Input& input, ElementRequest request) {
    if (!std::isfinite(input.time_step) || input.time_step < 0.0
        || (input.include_thermal_time_term && !(input.time_step > 0.0)))
        throw std::invalid_argument("CPEG8T requires a positive time step for thermal storage");
    for (double value : input.state)
        if (!std::isfinite(value))
            throw std::domain_error("CPEG8T state must be finite");
    const bool finite = input.strain_formulation == StrainFormulation::finite;
    Cpeg8Result result;
    const CartesianMaterialPointState empty_history{};
    for (std::size_t q = 0; q < 9; ++q) {
        const auto& p = input.geometry.points[q];
        const auto values = point_values(p, input.state), old = point_values(p, input.committed_state);
        ActivePointValues v;
        for (std::size_t k = 0; k < v.size(); ++k)
            v[k] = request.jacobian ? adlite::Scalar::independent(values[k], k, v.size()) : adlite::Scalar(values[k]);
        const double h = input.geometry.thickness;
        const double rx = input.geometry.reference_point[0], ry = input.geometry.reference_point[1];
        const adlite::Scalar x = p.x + (finite ? v[4] : adlite::Scalar(0.0));
        const adlite::Scalar y = p.y + (finite ? v[5] : adlite::Scalar(0.0));
        const adlite::Scalar thickness = h + v[6] + v[7] * (y - ry) - v[8] * (x - rx);
        const adlite::Scalar a = 1.0 + v[0], b = v[1], c = v[2], d = 1.0 + v[3];
        const adlite::Scalar det = a * d - b * c;
        const adlite::Scalar measure = finite ? p.measure * det * thickness / h : adlite::Scalar(p.measure);
        SymmetricTensor3 strain{v[0], v[3], (thickness - h) / h, 0.5 * (v[1] + v[2]), 0.0, 0.0};
        CartesianRotation rotation;
        adlite::Scalar ma = 1.0, mb = 0.0, mc = 0.0, md = 1.0, mdet = 1.0;
        if (finite) {
            const double old_det = (1.0 + old[0]) * (1.0 + old[3]) - old[1] * old[2];
            const double old_h = h + old[6] + old[7] * (p.y + old[5] - ry) - old[8] * (p.x + old[4] - rx);
            ma = 1.0 + 0.5 * (v[0] + old[0]);
            mb = 0.5 * (v[1] + old[1]);
            mc = 0.5 * (v[2] + old[2]);
            md = 1.0 + 0.5 * (v[3] + old[3]);
            mdet = ma * md - mb * mc;
            if (!(det.value() > 0.0) || !(old_det > 0.0) || !(mdet.value() > 0.0) || !(thickness.value() > 0.0)
                || !(old_h > 0.0))
                throw std::domain_error("CPEG8T current, committed and midpoint configurations must be positive");
            cartesian_detail::ActiveMatrix3 increment{};
            increment[0][0] = ((v[0] - old[0]) * md - (v[1] - old[1]) * mc) / mdet;
            increment[0][1] = ((v[1] - old[1]) * ma - (v[0] - old[0]) * mb) / mdet;
            increment[1][0] = ((v[2] - old[2]) * md - (v[3] - old[3]) * mc) / mdet;
            increment[1][1] = ((v[3] - old[3]) * ma - (v[2] - old[2]) * mb) / mdet;
            const auto kinematics = cartesian_detail::evaluate_hughes_winget_increment(increment);
            strain = kinematics.strain_increment;
            strain.zz = 2.0 * (thickness - old_h) / (thickness + old_h);
            rotation = kinematics.rotation;
        }
        const MaterialFunctionContext context{input.time, p.x, p.y, 0.0};
        auto fed = tensor_values(strain);
        const auto& history = input.committed_history ? (*input.committed_history)[q] : empty_history;
        if (finite && input.committed_history) {
            auto previous_context = context;
            previous_context.time -= input.time_step;
            const auto eigen = tensor_values(input.material.eigenstrain(old[9], previous_context));
            for (std::size_t k = 0; k < 6; ++k)
                fed[k] += history.elastic_strain[k] + history.plastic_strain[k] + history.creep_strain[k] + eigen[k];
        }
        SymmetricTensor3 stress;
        if (request.jacobian) {
            const auto tangent = evaluate_stress_tangent(input.material,
                fed,
                values[9],
                input.time_step,
                input.committed_history ? &history : nullptr,
                context);
            stress = compose_cartesian_stress(tangent, strain, v[9], tangent.thermal);
        } else {
            const auto material_state = input.material.response_values({fed[0], fed[1], fed[2], fed[3], fed[4], fed[5]},
                values[9],
                input.time_step,
                history,
                context);
            const auto& s = material_state.stress;
            stress = {s.xx, s.yy, s.zz, s.xy, s.yz, s.xz};
        }
        if (finite)
            stress = rotate_cartesian_tensor(stress, rotation);
        if (request.history) {
            auto updated = input.material.response_values({fed[0], fed[1], fed[2], fed[3], fed[4], fed[5]},
                values[9],
                input.time_step,
                history,
                context);
            // Form history differences before the objective rotation, so rigid
            // rotation of accepted plastic/creep tensors adds no dissipation.
            const std::array<double, 6> previous_stress{history.stress.xx,
                history.stress.yy,
                history.stress.zz,
                history.stress.xy,
                history.stress.yz,
                history.stress.xz};
            const std::array<double, 6> updated_stress{updated.stress.xx,
                updated.stress.yy,
                updated.stress.zz,
                updated.stress.xy,
                updated.stress.yz,
                updated.stress.xz};
            for (std::size_t k = 0; k < 6; ++k) {
                const double weight = p.measure * (k < 3 ? 1.0 : 2.0);
                result.elastic_energy_change +=
                    0.5 * weight
                    * (updated_stress[k] * updated.elastic_strain[k] - previous_stress[k] * history.elastic_strain[k]);
                const double average_stress = 0.5 * (updated_stress[k] + previous_stress[k]);
                result.plastic_dissipation_increment +=
                    weight * average_stress * (updated.plastic_strain[k] - history.plastic_strain[k]);
                result.creep_dissipation_increment +=
                    weight * average_stress * (updated.creep_strain[k] - history.creep_strain[k]);
            }
            if (finite)
                for (auto* tensor : {&updated.elastic_strain, &updated.plastic_strain, &updated.creep_strain}) {
                    const auto rotated = rotate_cartesian_tensor_values(
                        {(*tensor)[0], (*tensor)[1], (*tensor)[2], (*tensor)[3], (*tensor)[4], (*tensor)[5]},
                        rotation);
                    *tensor = {rotated.xx, rotated.yy, rotated.zz, rotated.xy, rotated.yz, rotated.xz};
                }
            updated.stress = {stress.xx.value(),
                stress.yy.value(),
                stress.zz.value(),
                stress.xy.value(),
                stress.yz.value(),
                stress.xz.value()};
            result.history[q] = updated;
        }
        if (request.stress)
            result.stress[q] = {stress.xx.value(),
                stress.yy.value(),
                stress.zz.value(),
                stress.xy.value(),
                stress.yz.value(),
                stress.xz.value()};
        for (std::size_t n = 0; n < 8; ++n) {
            const adlite::Scalar gx =
                finite ? (d * p.gradient_x[n] - c * p.gradient_y[n]) / det : adlite::Scalar(p.gradient_x[n]);
            const adlite::Scalar gy =
                finite ? (a * p.gradient_y[n] - b * p.gradient_x[n]) / det : adlite::Scalar(p.gradient_y[n]);
            // Native CPEG8T uses the rotation parameters directly in the
            // in-plane virtual axial gradients. The reference-control rows
            // below instead divide by current thickness. Independently checked
            // at initial thicknesses 0.1 and 0.2 (verification/abaqus/cpeg8t).
            const adlite::Scalar bx = finite ? -v[8] * p.shape[n] : adlite::Scalar(0.0);
            const adlite::Scalar by = finite ? v[7] * p.shape[n] : adlite::Scalar(0.0);
            add_row(result, 4 + n, measure * (gx * stress.xx + gy * stress.xy + bx * stress.zz), p, request.jacobian);
            add_row(result, 12 + n, measure * (gx * stress.xy + gy * stress.yy + by * stress.zz), p, request.jacobian);
        }
        const adlite::Scalar normal_force = measure * stress.zz / (finite ? thickness : adlite::Scalar(h));
        add_row(result, 20, normal_force, p, request.jacobian);
        add_row(result, 21, (y - ry) * normal_force, p, request.jacobian);
        add_row(result, 22, -(x - rx) * normal_force, p, request.jacobian);
        std::array<adlite::Scalar, 4> gx, gy;
        adlite::Scalar tx = 0.0, ty = 0.0;
        for (std::size_t n = 0; n < 4; ++n) {
            gx[n] = (md * p.temperature_gradient_x[n] - mc * p.temperature_gradient_y[n]) / mdet;
            gy[n] = (ma * p.temperature_gradient_y[n] - mb * p.temperature_gradient_x[n]) / mdet;
            tx += gx[n] * input.state[n];
            ty += gy[n] * input.state[n];
        }
        const auto conductivity = input.material.conductivity(v[9], context);
        const adlite::Scalar storage =
            input.include_thermal_time_term
                ? input.material.reference_heat_capacity(v[9], input.initial_temperature, context) * (v[9] - old[9])
                      / input.time_step
                : adlite::Scalar(0.0);
        for (std::size_t n = 0; n < 4; ++n) {
            add_row(result,
                n,
                measure * conductivity * (gx[n] * tx + gy[n] * ty) + p.measure * p.temperature_shape[n] * storage
                    - measure * p.temperature_shape[n] * input.volumetric_heat_source,
                p,
                request.jacobian);
            if (request.jacobian)
                for (std::size_t j = 0; j < 4; ++j)
                    result.jacobian[23 * n + j] += measure.value() * conductivity.value()
                                                   * (gx[n].value() * gx[j].value() + gy[n].value() * gy[j].value());
        }
        result.generated_heat_rate += measure.value() * input.volumetric_heat_source;
        result.stored_heat_rate += p.measure * storage.value();
        if (input.body_acceleration != std::array<double, 2>{}) {
            const double mass = p.measure * input.material.initial_density(input.initial_temperature, context);
            for (std::size_t n = 0; n < 8; ++n)
                for (std::size_t k = 0; k < 2; ++k)
                    result.residual[4 + 8 * k + n] -= mass * p.shape[n] * input.body_acceleration[k];
        }
    }
    if (!request.residual && !request.jacobian)
        result.residual.fill(0.0);
    return result;
}

Cpeg8Result evaluate_cpeg8t_boundary(const Cpeg8BoundaryInput& input, bool jacobian) {
    if (input.side >= 4)
        throw std::out_of_range("CPEG8T side index");
    Cpeg8Result result;
    const double g = std::sqrt(3.0 / 5.0);
    const std::array<double, 3> locations{-g, 0.0, g}, weights{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    const std::array<std::size_t, 3> nodes{input.side, (input.side + 1) % 4, 4 + input.side};
    for (std::size_t q = 0; q < 3; ++q) {
        const double xi = locations[q];
        const std::array<double, 3> shape{0.5 * xi * (xi - 1), 0.5 * xi * (xi + 1), 1 - xi * xi};
        const std::array<double, 3> derivative{xi - 0.5, xi + 0.5, -2 * xi};
        Cpeg8Point point;
        double dx = 0.0, dy = 0.0;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto n = nodes[i];
            point.shape[n] = shape[i];
            point.gradient_x[n] = derivative[i];
            point.x += shape[i] * input.geometry.coordinates[n][0];
            point.y += shape[i] * input.geometry.coordinates[n][1];
            dx += derivative[i] * input.geometry.coordinates[n][0];
            dy += derivative[i] * input.geometry.coordinates[n][1];
        }
        point.temperature_shape[nodes[0]] = 0.5 * (1 - xi);
        point.temperature_shape[nodes[1]] = 0.5 * (1 + xi);
        const auto values = point_values(point, input.state);
        ActivePointValues active;
        for (std::size_t i = 0; i < active.size(); ++i)
            active[i] = jacobian ? adlite::Scalar::independent(values[i], i, active.size()) : adlite::Scalar(values[i]);
        adlite::Scalar tx = dx, ty = dy, thickness = input.geometry.thickness;
        if (input.current) {
            tx += active[0];
            ty += active[2];
            thickness += active[6] + active[7] * (point.y + active[5] - input.geometry.reference_point[1])
                         - active[8] * (point.x + active[4] - input.geometry.reference_point[0]);
        }
        const auto length = adlite::hypot(tx, ty);
        if (!(length.value() > 0.0) || !(thickness.value() > 0.0))
            throw std::domain_error("CPEG8T boundary length and thickness must be positive");
        const auto area = weights[q] * length * thickness;
        if (input.kind == Cpeg8BoundaryKind::heat_flux || input.kind == Cpeg8BoundaryKind::convection) {
            const auto flux = input.kind == Cpeg8BoundaryKind::heat_flux ? adlite::Scalar(-input.value)
                                                                         : input.value * (active[9] - input.ambient);
            for (std::size_t i = 0; i < 2; ++i)
                add_row(result, nodes[i], area * point.temperature_shape[nodes[i]] * flux, point, jacobian);
        } else {
            adlite::Scalar fx = 0.0, fy = 0.0;
            if (input.kind == Cpeg8BoundaryKind::pressure) {
                fx = weights[q] * input.value * thickness * ty;
                fy = -weights[q] * input.value * thickness * tx;
            } else if (input.kind == Cpeg8BoundaryKind::traction_x)
                fx = -area * input.value;
            else
                fy = -area * input.value;
            for (std::size_t i = 0; i < 3; ++i) {
                add_row(result, 4 + nodes[i], shape[i] * fx, point, jacobian);
                add_row(result, 12 + nodes[i], shape[i] * fy, point, jacobian);
            }
        }
    }
    return result;
}
} // namespace fuelsim::elements
