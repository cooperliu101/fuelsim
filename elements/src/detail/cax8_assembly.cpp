#include "cax8_assembly.hpp"
#include "cax8_kinematics.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::cax8_detail {
namespace {

void add_row(elements::Cax8Result& result,
    std::size_t row,
    const adlite::Scalar& value,
    const PointKinematics& k,
    bool jacobian) {
    result.residual[row] += value.value();
    if (!jacobian)
        return;
    std::array<double, 6> d{};
    value.copy_derivatives(d.data(), 6);
    for (std::size_t j = 0; j < 20; ++j)
        for (std::size_t a = 0; a < 6; ++a)
            result.jacobian[20 * row + j] += d[a] * k.chain[a][j];
}
} // namespace

elements::Cax8Result evaluate(const elements::Cax8Input& data, elements::ElementRequest request) {
    const auto& geometry = data.geometry;
    const auto& state = data.state;
    const auto& committed = data.committed_state;
    const auto* history = data.committed_history;
    const double dt = data.time_step;
    const bool jacobian = request.jacobian;
    const bool thermal_time = data.include_thermal_time_term;
    if (history && (!(dt > 0) || !std::isfinite(dt)))
        throw std::invalid_argument("CAX8T material update requires positive finite time step");
    const bool finite = data.strain_formulation == StrainFormulation::finite;
    elements::Cax8Result result;
    for (std::size_t q = 0; q < geometry.point_count; ++q) {
        const auto& p = geometry.points[q];
        const auto k = evaluate_kinematics(p, state, committed, finite, jacobian);
        const auto& t = k.active[5];
        const MaterialFunctionContext context = {data.time, p.radius, 0, p.axial_coordinate};
        std::array<double, 5> fed = {k.strain[0].value(),
            k.strain[1].value(),
            k.strain[2].value(),
            k.strain[3].value(),
            t.value()};
        if (finite && history) {
            auto old_context = context;
            old_context.time -= dt;
            const auto eigen = data.material.eigenstrain_rz(k.old[5], old_context);
            const std::array<double, 4> imposed = {eigen.rr.value(),
                eigen.zz.value(),
                eigen.hoop.value(),
                eigen.rz.value()};
            for (std::size_t c = 0; c < 4; ++c)
                fed[c] += (*history)[q].elastic_strain[c] + (*history)[q].plastic_strain[c]
                          + (*history)[q].creep_strain[c] + imposed[c];
        }
        std::array<adlite::Scalar, 5> material;
        for (std::size_t i = 0; i < 5; ++i)
            material[i] = jacobian ? adlite::Scalar::independent(fed[i], i, 5) : adlite::Scalar(fed[i]);
        const auto raw =
            history ? data.material
                          .response(material[0],
                              material[1],
                              material[2],
                              material[3],
                              material[4],
                              dt,
                              (*history)[q],
                              context)
                          .stress
                    : data.material.stress(material[0], material[1], material[2], material[3], material[4], context);
        const std::array<adlite::Scalar, 4> components = {raw.rr, raw.zz, raw.hoop, raw.rz};
        const std::array<adlite::Scalar, 5> inputs = {k.strain[0], k.strain[1], k.strain[2], k.strain[3], t};
        std::array<adlite::Scalar, 4> composed;
        for (std::size_t c = 0; c < 4; ++c) {
            std::array<double, 5> partials{};
            if (jacobian)
                components[c].copy_derivatives(partials.data(), 5);
            composed[c] = jacobian ? adlite::compose(components[c].value(), inputs.data(), partials.data(), 5)
                                   : adlite::Scalar(components[c].value());
        }
        AxisymmetricStress stress = {composed[0], composed[1], composed[2], composed[3]};
        if (finite)
            stress = rotate_axisymmetric_tensor(stress, k.rotation);
        if (history) {
            const AxisymmetricRotation rotation = {k.rotation.rr.value(),
                k.rotation.rz.value(),
                k.rotation.zr.value(),
                k.rotation.zz.value(),
                1.0};
            const auto response =
                finite ? data.material.incremental_response(k.strain[0].value(),
                             k.strain[1].value(),
                             k.strain[2].value(),
                             k.strain[3].value(),
                             rotation,
                             t.value(),
                             k.old[5],
                             dt,
                             (*history)[q],
                             context)
                       : data.material.response(fed[0], fed[1], fed[2], fed[3], fed[4], dt, (*history)[q], context);
            result.history[q] = IsotropicThermoelasticMaterial::state_values(response.trial_state);
        }
        result.history[q].stress = {stress.rr.value(), stress.zz.value(), stress.hoop.value(), stress.rz.value()};
        for (std::size_t n = 0; n < 8; ++n) {
            const auto gradient = mechanical_gradient(p, k, n, finite);
            const auto& br = gradient.radial;
            const auto& bz = gradient.axial;
            const auto& bh = gradient.hoop;
            add_row(result, 4 + n, k.measure * (br * stress.rr + bz * stress.rz + bh * stress.hoop), k, jacobian);
            add_row(result, 12 + n, k.measure * (bz * stress.zz + br * stress.rz), k, jacobian);
        }
        // Four temperature shape functions share the full quadratic geometric map.
        const auto thermal_map = thermal_geometry(p, k, finite);
        const auto& gr = thermal_map.radial;
        const auto& gz = thermal_map.axial;
        adlite::Scalar tr = 0, tz = 0;
        for (std::size_t n = 0; n < 4; ++n) {
            tr += gr[n] * state[n];
            tz += gz[n] * state[n];
        }
        const auto conductivity = data.material.conductivity(t, context);
        const adlite::Scalar capacity =
            history && thermal_time ? data.material.heat_capacity(t, context) * (t - k.old[5]) / dt : adlite::Scalar(0);
        const auto source_map = source_geometry(p, state, finite);
        const double source_measure = source_map.measure;
        const auto& source_derivative = source_map.derivative;
        for (std::size_t n = 0; n < 4; ++n) {
            const auto row = k.measure * (conductivity * (gr[n] * tr + gz[n] * tz) + p.temperature_shape[n] * capacity);
            add_row(result, n, row, k, jacobian);
            const double source = p.temperature_shape[n] * data.volumetric_heat_source;
            result.residual[n] -= source * source_measure;
            if (jacobian)
                for (std::size_t j = 4; j < 20; ++j)
                    result.jacobian[20 * n + j] -= source * source_derivative[j];
            if (jacobian)
                for (std::size_t j = 0; j < 4; ++j)
                    result.jacobian[20 * n + j] += k.measure.value() * conductivity.value()
                                                   * (gr[n].value() * gr[j].value() + gz[n].value() * gz[j].value());
        }
        result.generated_heat_rate += source_measure * data.volumetric_heat_source;
        result.stored_heat_rate += k.measure.value() * capacity.value();
    }
    if (request.stress)
        for (std::size_t q = 0; q < result.history.size(); ++q)
            result.stress[q] = result.history[q].stress;
    if (!request.residual && !request.jacobian)
        result.residual.fill(0.0);
    if (!request.history)
        result.history = {};
    return result;
}
} // namespace fuelsim::cax8_detail
