#include "fuelsim/hex8_thermoelastic.hpp"

#include <adlite/adlite.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

namespace fuelsim {
namespace {

constexpr double gauss = 0.577350269189625764509148780501957456;
constexpr std::array<std::array<double, 3>, 8> hex8_signs = {{{{-1.0, -1.0, -1.0}},
                                                              {{1.0, -1.0, -1.0}},
                                                              {{1.0, 1.0, -1.0}},
                                                              {{-1.0, 1.0, -1.0}},
                                                              {{-1.0, -1.0, 1.0}},
                                                              {{1.0, -1.0, 1.0}},
                                                              {{1.0, 1.0, 1.0}},
                                                              {{-1.0, 1.0, 1.0}}}};

double determinant(const std::array<std::array<double, 3>, 3>& matrix) {
    return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
           matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
           matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
}

std::array<std::array<double, 3>, 3> inverse(const std::array<std::array<double, 3>, 3>& matrix,
                                             double determinant_value) {
    return {{{{(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant_value,
               (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / determinant_value,
               (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant_value}},
             {{(matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / determinant_value,
               (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant_value,
               (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / determinant_value}},
             {{(matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant_value,
               (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / determinant_value,
               (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant_value}}}};
}

Hex8LocalAdValues passive_hex8_state(const Hex8LocalValues& state) {
    Hex8LocalAdValues result{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        result[dof] = state[dof];
    return result;
}

Hex8LocalAdValues active_hex8_state(const Hex8LocalValues& state) {
    Hex8LocalAdValues result{};
    adlite::seed_identity(state.data(), state.size(), result.data());
    return result;
}

Hex8LocalResidual hex8_residual_values(const Hex8LocalAdValues& residual) {
    Hex8LocalResidual result{};
    for (std::size_t row = 0; row < result.size(); ++row)
        result[row] = residual[row].value();
    return result;
}

Hex8LocalSystem hex8_linearized_values(const Hex8LocalAdValues& state, const Hex8LocalAdValues& residual) {
    Hex8LocalSystem result{};
    adlite::extract_jacobian(residual.data(), residual.size(), state.size(), result.residual.data(),
                             result.jacobian.data());
    return result;
}

Quad4FaceLocalAdValues passive_face_state(const Quad4FaceLocalValues& state) {
    Quad4FaceLocalAdValues result{};
    for (std::size_t dof = 0; dof < state.size(); ++dof)
        result[dof] = state[dof];
    return result;
}

Quad4FaceLocalAdValues active_face_state(const Quad4FaceLocalValues& state) {
    Quad4FaceLocalAdValues result{};
    adlite::seed_identity(state.data(), state.size(), result.data());
    return result;
}

Quad4FaceLocalResidual face_residual_values(const Quad4FaceLocalAdValues& residual) {
    Quad4FaceLocalResidual result{};
    for (std::size_t row = 0; row < result.size(); ++row)
        result[row] = residual[row].value();
    return result;
}

Quad4FaceLocalSystem face_linearized_values(const Quad4FaceLocalAdValues& state,
                                            const Quad4FaceLocalAdValues& residual) {
    Quad4FaceLocalSystem result{};
    adlite::extract_jacobian(residual.data(), residual.size(), state.size(), result.residual.data(),
                             result.jacobian.data());
    return result;
}

adlite::Scalar interpolate_hex8(const std::array<double, 8>& coefficients, const Hex8LocalAdValues& state,
                                std::size_t offset) {
    adlite::Scalar result = 0.0;
    for (std::size_t node = 0; node < 8; ++node)
        result += coefficients[node] * state[offset + node];
    return result;
}

SymmetricTensor3 strain_at(const Hex8QuadraturePoint& point, const Hex8LocalAdValues& state) {
    adlite::Scalar ux_x = 0.0;
    adlite::Scalar ux_y = 0.0;
    adlite::Scalar ux_z = 0.0;
    adlite::Scalar uy_x = 0.0;
    adlite::Scalar uy_y = 0.0;
    adlite::Scalar uy_z = 0.0;
    adlite::Scalar uz_x = 0.0;
    adlite::Scalar uz_y = 0.0;
    adlite::Scalar uz_z = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        ux_x += point.gradient[node][0] * state[8 + node];
        ux_y += point.gradient[node][1] * state[8 + node];
        ux_z += point.gradient[node][2] * state[8 + node];
        uy_x += point.gradient[node][0] * state[16 + node];
        uy_y += point.gradient[node][1] * state[16 + node];
        uy_z += point.gradient[node][2] * state[16 + node];
        uz_x += point.gradient[node][0] * state[24 + node];
        uz_y += point.gradient[node][1] * state[24 + node];
        uz_z += point.gradient[node][2] * state[24 + node];
    }
    return {ux_x, uy_y, uz_z, 0.5 * (ux_y + uy_x), 0.5 * (uy_z + uz_y), 0.5 * (ux_z + uz_x)};
}

void add_hex8_point_residual(const Hex8QuadraturePoint& point, const Hex8LocalAdValues& state,
                             const IsotropicThermoelasticMaterial& material, double time, double volumetric_heat_source,
                             double constant_heat_capacity, const Hex8LocalValues* committed_state, double time_step,
                             Hex8LocalAdValues& residual) {
    const adlite::Scalar temperature = interpolate_hex8(point.shape, state, 0);
    adlite::Scalar gradient_temperature_x = 0.0;
    adlite::Scalar gradient_temperature_y = 0.0;
    adlite::Scalar gradient_temperature_z = 0.0;
    for (std::size_t node = 0; node < 8; ++node) {
        gradient_temperature_x += point.gradient[node][0] * state[node];
        gradient_temperature_y += point.gradient[node][1] * state[node];
        gradient_temperature_z += point.gradient[node][2] * state[node];
    }
    const adlite::Scalar conductivity =
        material.conductivity_cartesian(temperature, time, point.position.x, point.position.y, point.position.z);
    const SymmetricTensor3 stress = material.stress_cartesian(strain_at(point, state), temperature, time,
                                                              point.position.x, point.position.y, point.position.z);
    adlite::Scalar temperature_rate = 0.0;
    adlite::Scalar heat_capacity = 0.0;
    if (committed_state != nullptr) {
        double old_temperature = 0.0;
        for (std::size_t node = 0; node < 8; ++node)
            old_temperature += point.shape[node] * (*committed_state)[node];
        temperature_rate = (temperature - old_temperature) / time_step;
        heat_capacity = material.properties().functions
                            ? material.heat_capacity_cartesian(temperature, time, point.position.x, point.position.y,
                                                               point.position.z)
                            : adlite::Scalar(constant_heat_capacity);
    }
    for (std::size_t node = 0; node < 8; ++node) {
        const double gradient_x = point.gradient[node][0];
        const double gradient_y = point.gradient[node][1];
        const double gradient_z = point.gradient[node][2];
        residual[node] +=
            point.weighted_measure *
            (conductivity * (gradient_x * gradient_temperature_x + gradient_y * gradient_temperature_y +
                             gradient_z * gradient_temperature_z) +
             point.shape[node] * heat_capacity * temperature_rate - point.shape[node] * volumetric_heat_source);
        residual[8 + node] +=
            point.weighted_measure * (stress.xx * gradient_x + stress.xy * gradient_y + stress.xz * gradient_z);
        residual[16 + node] +=
            point.weighted_measure * (stress.xy * gradient_x + stress.yy * gradient_y + stress.yz * gradient_z);
        residual[24 + node] +=
            point.weighted_measure * (stress.xz * gradient_x + stress.yz * gradient_y + stress.zz * gradient_z);
    }
}

std::array<SymmetricTensor3Values, 8> stress_values(const Hex8Geometry& geometry, const Hex8LocalValues& state,
                                                    const IsotropicThermoelasticMaterial& material, double time) {
    const Hex8LocalAdValues ad_state = passive_hex8_state(state);
    std::array<SymmetricTensor3Values, 8> result{};
    for (std::size_t q = 0; q < geometry.points.size(); ++q) {
        const Hex8QuadraturePoint& point = geometry.points[q];
        const adlite::Scalar temperature = interpolate_hex8(point.shape, ad_state, 0);
        const SymmetricTensor3 stress = material.stress_cartesian(strain_at(point, ad_state), temperature, time,
                                                                  point.position.x, point.position.y, point.position.z);
        result[q] = {stress.xx.value(), stress.yy.value(), stress.zz.value(),
                     stress.xy.value(), stress.yz.value(), stress.xz.value()};
    }
    return result;
}

} // namespace

Hex8Geometry make_hex8_geometry(const Hex8Coordinates& coordinates) {
    Hex8Geometry geometry{};
    std::size_t q = 0;
    for (double zeta : {-gauss, gauss}) {
        for (double eta : {-gauss, gauss}) {
            for (double xi : {-gauss, gauss}) {
                std::array<double, 8> shape{};
                std::array<std::array<double, 3>, 8> derivative{};
                for (std::size_t node = 0; node < 8; ++node) {
                    const double sx = hex8_signs[node][0];
                    const double sy = hex8_signs[node][1];
                    const double sz = hex8_signs[node][2];
                    shape[node] = 0.125 * (1.0 + sx * xi) * (1.0 + sy * eta) * (1.0 + sz * zeta);
                    derivative[node] = {{0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta),
                                         0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta),
                                         0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta)}};
                }
                std::array<std::array<double, 3>, 3> jacobian{};
                CartesianPoint3 position{0.0, 0.0, 0.0};
                for (std::size_t node = 0; node < 8; ++node) {
                    const std::array<double, 3> coordinate = {
                        coordinates[node].x,
                        coordinates[node].y,
                        coordinates[node].z,
                    };
                    position.x += shape[node] * coordinates[node].x;
                    position.y += shape[node] * coordinates[node].y;
                    position.z += shape[node] * coordinates[node].z;
                    for (std::size_t physical = 0; physical < 3; ++physical) {
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            jacobian[physical][natural] += coordinate[physical] * derivative[node][natural];
                    }
                }
                const double determinant_value = determinant(jacobian);
                if (!std::isfinite(determinant_value) || !(determinant_value > 0.0))
                    throw std::invalid_argument("Hex8Geometry requires a finite positive Jacobian determinant");
                const std::array<std::array<double, 3>, 3> inverse_jacobian = inverse(jacobian, determinant_value);
                Hex8QuadraturePoint& point = geometry.points[q++];
                point.shape = shape;
                point.position = position;
                point.weighted_measure = determinant_value;
                for (std::size_t node = 0; node < 8; ++node) {
                    for (std::size_t physical = 0; physical < 3; ++physical) {
                        for (std::size_t natural = 0; natural < 3; ++natural)
                            point.gradient[node][physical] +=
                                derivative[node][natural] * inverse_jacobian[natural][physical];
                    }
                }
            }
        }
    }
    return geometry;
}

Quad4FaceGeometry make_quad4_face_geometry(const Quad4FaceCoordinates& coordinates) {
    Quad4FaceGeometry geometry{};
    const std::array<std::array<double, 2>, 4> locations = {
        {{{-gauss, -gauss}}, {{gauss, -gauss}}, {{gauss, gauss}}, {{-gauss, gauss}}}};
    for (std::size_t q = 0; q < locations.size(); ++q) {
        const double xi = locations[q][0];
        const double eta = locations[q][1];
        const std::array<double, 4> shape = {{0.25 * (1.0 - xi) * (1.0 - eta), 0.25 * (1.0 + xi) * (1.0 - eta),
                                              0.25 * (1.0 + xi) * (1.0 + eta), 0.25 * (1.0 - xi) * (1.0 + eta)}};
        const std::array<double, 4> derivative_xi = {
            {-0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta)}};
        const std::array<double, 4> derivative_eta = {
            {-0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)}};
        CartesianPoint3 tangent_xi{0.0, 0.0, 0.0};
        CartesianPoint3 tangent_eta{0.0, 0.0, 0.0};
        for (std::size_t node = 0; node < 4; ++node) {
            tangent_xi.x += derivative_xi[node] * coordinates[node].x;
            tangent_xi.y += derivative_xi[node] * coordinates[node].y;
            tangent_xi.z += derivative_xi[node] * coordinates[node].z;
            tangent_eta.x += derivative_eta[node] * coordinates[node].x;
            tangent_eta.y += derivative_eta[node] * coordinates[node].y;
            tangent_eta.z += derivative_eta[node] * coordinates[node].z;
        }
        const CartesianPoint3 area_vector{tangent_xi.y * tangent_eta.z - tangent_xi.z * tangent_eta.y,
                                          tangent_xi.z * tangent_eta.x - tangent_xi.x * tangent_eta.z,
                                          tangent_xi.x * tangent_eta.y - tangent_xi.y * tangent_eta.x};
        const double measure =
            std::sqrt(area_vector.x * area_vector.x + area_vector.y * area_vector.y + area_vector.z * area_vector.z);
        if (!std::isfinite(measure) || !(measure > 0.0))
            throw std::invalid_argument("Quad4FaceGeometry requires a finite positive area measure");
        geometry.points[q] = {shape, area_vector, measure};
    }
    return geometry;
}

Hex8ThermoelasticKernel::Hex8ThermoelasticKernel(IsotropicThermoelasticMaterial material, double volumetric_heat_source)
    : _material(std::move(material)), _volumetric_heat_source(volumetric_heat_source), _time(0.0) {}
double Hex8ThermoelasticKernel::volumetric_heat_source() const noexcept {
    return _volumetric_heat_source;
}
void Hex8ThermoelasticKernel::set_volumetric_heat_source(double value) noexcept {
    _volumetric_heat_source = value;
}
void Hex8ThermoelasticKernel::set_time(double value) noexcept {
    _time = value;
}
void Hex8ThermoelasticKernel::residual_ad(const Hex8Geometry& geometry, const Hex8LocalAdValues& state,
                                          Hex8LocalAdValues& residual) const {
    residual.fill(adlite::Scalar(0.0));
    for (const Hex8QuadraturePoint& point : geometry.points)
        add_hex8_point_residual(point, state, _material, _time, _volumetric_heat_source, 0.0, nullptr, 0.0, residual);
}
Hex8LocalResidual Hex8ThermoelasticKernel::residual(const Hex8Geometry& geometry, const Hex8LocalValues& state) const {
    const Hex8LocalAdValues ad_state = passive_hex8_state(state);
    Hex8LocalAdValues residual{};
    residual_ad(geometry, ad_state, residual);
    return hex8_residual_values(residual);
}
Hex8LocalSystem Hex8ThermoelasticKernel::linearize(const Hex8Geometry& geometry, const Hex8LocalValues& state) const {
    const Hex8LocalAdValues ad_state = active_hex8_state(state);
    Hex8LocalAdValues residual{};
    residual_ad(geometry, ad_state, residual);
    return hex8_linearized_values(ad_state, residual);
}
std::array<SymmetricTensor3Values, 8> Hex8ThermoelasticKernel::stress_values(const Hex8Geometry& geometry,
                                                                             const Hex8LocalValues& state) const {
    return fuelsim::stress_values(geometry, state, _material, _time);
}

Hex8TransientKernel::Hex8TransientKernel(IsotropicThermoelasticMaterial material, double constant_heat_capacity,
                                         double volumetric_heat_source)
    : _material(std::move(material)), _constant_heat_capacity(constant_heat_capacity),
      _volumetric_heat_source(volumetric_heat_source), _time(0.0) {
    if (!std::isfinite(_constant_heat_capacity) || !(_constant_heat_capacity > 0.0))
        throw std::invalid_argument("Hex8TransientKernel heat capacity must be finite and positive");
}
double Hex8TransientKernel::volumetric_heat_source() const noexcept {
    return _volumetric_heat_source;
}
double Hex8TransientKernel::heat_capacity(double temperature, double x, double y, double z) const {
    return _material.properties().functions ? _material.heat_capacity_cartesian(temperature, _time, x, y, z).value()
                                            : _constant_heat_capacity;
}
void Hex8TransientKernel::set_volumetric_heat_source(double value) noexcept {
    _volumetric_heat_source = value;
}
void Hex8TransientKernel::set_time(double value) noexcept {
    _time = value;
}
void Hex8TransientKernel::residual_ad(const Hex8Geometry& geometry, const Hex8LocalAdValues& current_state,
                                      const Hex8LocalValues& committed_state, double time_step,
                                      Hex8LocalAdValues& residual) const {
    if (!std::isfinite(time_step) || !(time_step > 0.0))
        throw std::invalid_argument("Hex8TransientKernel time step must be finite and positive");
    residual.fill(adlite::Scalar(0.0));
    for (const Hex8QuadraturePoint& point : geometry.points)
        add_hex8_point_residual(point, current_state, _material, _time, _volumetric_heat_source,
                                _constant_heat_capacity, &committed_state, time_step, residual);
}
Hex8LocalResidual Hex8TransientKernel::residual(const Hex8Geometry& geometry, const Hex8LocalValues& current_state,
                                                const Hex8LocalValues& committed_state, double time_step) const {
    const Hex8LocalAdValues ad_state = passive_hex8_state(current_state);
    Hex8LocalAdValues residual{};
    residual_ad(geometry, ad_state, committed_state, time_step, residual);
    return hex8_residual_values(residual);
}
Hex8LocalSystem Hex8TransientKernel::linearize(const Hex8Geometry& geometry, const Hex8LocalValues& current_state,
                                               const Hex8LocalValues& committed_state, double time_step) const {
    const Hex8LocalAdValues ad_state = active_hex8_state(current_state);
    Hex8LocalAdValues residual{};
    residual_ad(geometry, ad_state, committed_state, time_step, residual);
    return hex8_linearized_values(ad_state, residual);
}
std::array<SymmetricTensor3Values, 8> Hex8TransientKernel::stress_values(const Hex8Geometry& geometry,
                                                                         const Hex8LocalValues& state) const {
    return fuelsim::stress_values(geometry, state, _material, _time);
}

Quad4FacePressureKernel::Quad4FacePressureKernel(double pressure) : _pressure(pressure) {}
void Quad4FacePressureKernel::set_pressure(double value) noexcept {
    _pressure = value;
}
Quad4FaceLocalResidual Quad4FacePressureKernel::residual(const Quad4FaceGeometry& geometry,
                                                         const Quad4FaceLocalValues& state) const {
    const Quad4FaceLocalAdValues ad_state = passive_face_state(state);
    Quad4FaceLocalAdValues residual{};
    for (const Quad4FaceQuadraturePoint& point : geometry.points) {
        for (std::size_t node = 0; node < 4; ++node) {
            residual[4 + node] += _pressure * point.shape[node] * point.outward_area_vector.x;
            residual[8 + node] += _pressure * point.shape[node] * point.outward_area_vector.y;
            residual[12 + node] += _pressure * point.shape[node] * point.outward_area_vector.z;
        }
    }
    static_cast<void>(ad_state);
    return face_residual_values(residual);
}
Quad4FaceLocalSystem Quad4FacePressureKernel::linearize(const Quad4FaceGeometry& geometry,
                                                        const Quad4FaceLocalValues& state) const {
    const Quad4FaceLocalAdValues ad_state = active_face_state(state);
    Quad4FaceLocalAdValues residual{};
    for (const Quad4FaceQuadraturePoint& point : geometry.points) {
        for (std::size_t node = 0; node < 4; ++node) {
            residual[4 + node] += _pressure * point.shape[node] * point.outward_area_vector.x;
            residual[8 + node] += _pressure * point.shape[node] * point.outward_area_vector.y;
            residual[12 + node] += _pressure * point.shape[node] * point.outward_area_vector.z;
        }
    }
    return face_linearized_values(ad_state, residual);
}

Quad4FaceTractionKernel::Quad4FaceTractionKernel(CartesianTractionComponent component, double traction)
    : _component(component), _traction(traction) {}
void Quad4FaceTractionKernel::set_traction(double value) noexcept {
    _traction = value;
}
Quad4FaceLocalResidual Quad4FaceTractionKernel::residual(const Quad4FaceGeometry& geometry,
                                                         const Quad4FaceLocalValues& state) const {
    const Quad4FaceLocalAdValues ad_state = passive_face_state(state);
    Quad4FaceLocalAdValues residual{};
    const std::size_t offset =
        _component == CartesianTractionComponent::x ? 4 : (_component == CartesianTractionComponent::y ? 8 : 12);
    for (const Quad4FaceQuadraturePoint& point : geometry.points) {
        for (std::size_t node = 0; node < 4; ++node)
            residual[offset + node] -= _traction * point.weighted_measure * point.shape[node];
    }
    static_cast<void>(ad_state);
    return face_residual_values(residual);
}
Quad4FaceLocalSystem Quad4FaceTractionKernel::linearize(const Quad4FaceGeometry& geometry,
                                                        const Quad4FaceLocalValues& state) const {
    const Quad4FaceLocalAdValues ad_state = active_face_state(state);
    Quad4FaceLocalAdValues residual{};
    const std::size_t offset =
        _component == CartesianTractionComponent::x ? 4 : (_component == CartesianTractionComponent::y ? 8 : 12);
    for (const Quad4FaceQuadraturePoint& point : geometry.points) {
        for (std::size_t node = 0; node < 4; ++node)
            residual[offset + node] -= _traction * point.weighted_measure * point.shape[node];
    }
    return face_linearized_values(ad_state, residual);
}

Quad4FaceConvectionKernel::Quad4FaceConvectionKernel(double heat_transfer_coefficient, double ambient_temperature)
    : _heat_transfer_coefficient(heat_transfer_coefficient), _ambient_temperature(ambient_temperature) {}
void Quad4FaceConvectionKernel::set_properties(double coefficient, double ambient) noexcept {
    _heat_transfer_coefficient = coefficient;
    _ambient_temperature = ambient;
}
Quad4FaceLocalResidual Quad4FaceConvectionKernel::residual(const Quad4FaceGeometry& geometry,
                                                           const Quad4FaceLocalValues& state) const {
    const Quad4FaceLocalAdValues ad_state = passive_face_state(state);
    Quad4FaceLocalAdValues residual{};
    for (const Quad4FaceQuadraturePoint& point : geometry.points) {
        adlite::Scalar temperature = 0.0;
        for (std::size_t node = 0; node < 4; ++node)
            temperature += point.shape[node] * ad_state[node];
        const adlite::Scalar heat_flux = _heat_transfer_coefficient * (temperature - _ambient_temperature);
        for (std::size_t node = 0; node < 4; ++node)
            residual[node] += point.weighted_measure * point.shape[node] * heat_flux;
    }
    return face_residual_values(residual);
}
Quad4FaceLocalSystem Quad4FaceConvectionKernel::linearize(const Quad4FaceGeometry& geometry,
                                                          const Quad4FaceLocalValues& state) const {
    const Quad4FaceLocalAdValues ad_state = active_face_state(state);
    Quad4FaceLocalAdValues residual{};
    for (const Quad4FaceQuadraturePoint& point : geometry.points) {
        adlite::Scalar temperature = 0.0;
        for (std::size_t node = 0; node < 4; ++node)
            temperature += point.shape[node] * ad_state[node];
        const adlite::Scalar heat_flux = _heat_transfer_coefficient * (temperature - _ambient_temperature);
        for (std::size_t node = 0; node < 4; ++node)
            residual[node] += point.weighted_measure * point.shape[node] * heat_flux;
    }
    return face_linearized_values(ad_state, residual);
}

} // namespace fuelsim
