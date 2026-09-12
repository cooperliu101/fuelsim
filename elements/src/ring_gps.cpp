#include "ring_gps.hpp"
#include "contact_common.hpp"
#include <cmath>
#include <stdexcept>

namespace fuelsim::elements {
namespace {
using adlite::Scalar;
constexpr std::size_t compact_width = 6;

void require_finite(double value, const char* message) {
    if (!std::isfinite(value))
        throw std::domain_error(message);
}

double validate(const RingGpsInput& input) {
    if (input.strain_formulation != StrainFormulation::small && input.strain_formulation != StrainFormulation::finite)
        throw std::invalid_argument("Ring GPS contact has an unknown strain formulation");
    if (input.normal.augmented_lagrangian)
        throw std::invalid_argument("Ring GPS contact does not support augmented Lagrangian enforcement");
    if (!input.thermal && !input.mechanical)
        throw std::invalid_argument("Ring GPS contact requires thermal or mechanical interaction");
    const auto& geometry = input.geometry;
    for (const double value :
        {geometry.primary_reference_radius, geometry.secondary_reference_radius, geometry.z_lower, geometry.z_upper})
        require_finite(value, "Ring GPS contact geometry must be finite");
    const double height = geometry.z_upper - geometry.z_lower;
    if (!(geometry.primary_reference_radius > 0.0) || !(geometry.secondary_reference_radius > 0.0) || !(height > 0.0)
        || !std::isfinite(height))
        throw std::domain_error("Ring GPS contact requires positive radii and finite positive axial coverage");
    const double area = 2.0 * std::acos(-1.0) * geometry.secondary_reference_radius * height;
    if (!std::isfinite(area) || !(area > 0.0))
        throw std::domain_error("Ring GPS contact has an invalid reference surface measure");
    for (const auto* state : {&input.state, &input.committed_state}) {
        for (const double value : *state)
            require_finite(value, "Ring GPS contact states must be finite");
        for (const double radius :
            {geometry.primary_reference_radius + (*state)[2], geometry.secondary_reference_radius + (*state)[3]})
            if (!(radius > 0.0) || !std::isfinite(radius))
                throw std::domain_error("Ring GPS contact requires positive current and committed radii");
        if (input.strain_formulation == StrainFormulation::finite)
            for (const double current_height : {height + (*state)[5] - (*state)[4], height + (*state)[7] - (*state)[6]})
                if (!(current_height > 0.0) || !std::isfinite(current_height))
                    throw std::domain_error("Ring GPS contact requires positive current and committed axial lengths");
    }
    for (const double value :
        {input.normal.penalty, input.normal.friction_coefficient, input.normal.maximum_elastic_slip})
        require_finite(value, "Ring GPS mechanical contact parameters must be finite");
    if ((input.mechanical && !(input.normal.penalty > 0.0)) || input.normal.penalty < 0.0
        || input.normal.friction_coefficient < 0.0 || input.normal.maximum_elastic_slip < 0.0)
        throw std::domain_error("Ring GPS contact requires a positive penalty and nonnegative friction parameters");
    for (const double value : {input.heat.gap_conductivity,
             input.heat.minimum_gap,
             input.heat.conductance,
             input.heat.clearance_derivative,
             input.heat.pressure_derivative,
             input.heat.temperature_derivative,
             input.heat.reference_temperature,
             input.heat.contact_penalty})
        require_finite(value, "Ring GPS heat contact parameters must be finite");
    if (input.heat.law != GapHeatConductanceLaw::gas_gap && input.heat.law != GapHeatConductanceLaw::affine)
        throw std::invalid_argument("Ring GPS contact has an unknown gap conductance law");
    if (input.thermal
        && ((input.heat.law == GapHeatConductanceLaw::gas_gap
                && (input.heat.gap_conductivity < 0.0 || !(input.heat.minimum_gap > 0.0)))
            || (input.heat.law == GapHeatConductanceLaw::affine && input.heat.contact_penalty < 0.0)))
        throw std::domain_error("Ring GPS contact has invalid gap conductance parameters");
    for (const auto& history : input.committed_history) {
        for (const double value :
            {history.elastic_tangential_slip, history.total_tangential_slip, history.normal_multiplier})
            require_finite(value, "Ring GPS contact history must be finite");
        for (const auto* vector : {&history.cartesian_elastic_tangential_slip,
                 &history.cartesian_total_tangential_slip,
                 &history.cartesian_contact_normal,
                 &history.cartesian_contact_tangent_first})
            for (const double value : *vector)
                require_finite(value, "Ring GPS contact history must be finite");
        if (history.normal_multiplier != 0.0)
            throw std::invalid_argument("Ring GPS contact does not accept augmented normal history");
    }
    return area;
}

Scalar local_value(double value, std::size_t component, bool jacobian) {
    require_finite(value, "Ring GPS contact has a nonfinite local kinematic value");
    return jacobian ? Scalar::independent(value, component, compact_width) : Scalar(value);
}

void accumulate(RingGpsResult& result,
    std::size_t row,
    const Scalar& value,
    const std::array<double, 2>& shape,
    ElementRequest request) {
    require_finite(value.value(), "Ring GPS contact has a nonfinite residual");
    if (request.residual)
        result.residual[row] += value.value();
    if (!request.jacobian || value.derivative_size() == 0)
        return;
    std::array<double, compact_width> derivative{};
    for (std::size_t i = 0; i < compact_width; ++i) {
        derivative[i] = value.derivative(i);
        require_finite(derivative[i], "Ring GPS contact has a nonfinite tangent");
    }
    const std::size_t offset = row * ring_gps_local_dof_count;
    result.jacobian[offset] += derivative[2];
    result.jacobian[offset + 1] += derivative[1];
    result.jacobian[offset + 2] += derivative[0];
    result.jacobian[offset + 3] -= derivative[0];
    result.jacobian[offset + 3] += derivative[4];
    for (std::size_t i = 0; i < 2; ++i) {
        result.jacobian[offset + 4 + i] -= derivative[3] * shape[i];
        result.jacobian[offset + 6 + i] += derivative[3] * shape[i];
    }
    result.jacobian[offset + 6] -= derivative[5];
    result.jacobian[offset + 7] += derivative[5];
}
} // namespace

RingGpsResult evaluate_ring_gps(const RingGpsInput& input, ElementRequest request) {
    const double reference_area = validate(input);
    RingGpsResult result;
    const auto& state = input.state;
    const auto& old_state = input.committed_state;
    Scalar point_area = 0.5 * reference_area;
    if (input.strain_formulation == StrainFormulation::finite) {
        const Scalar current_radius =
            local_value(input.geometry.secondary_reference_radius + state[3], 4, request.jacobian);
        const Scalar current_height =
            local_value(input.geometry.z_upper - input.geometry.z_lower + state[7] - state[6], 5, request.jacobian);
        point_area = std::acos(-1.0) * current_radius * current_height;
        if (!std::isfinite(point_area.value()) || !(point_area.value() > 0.0))
            throw std::domain_error("Ring GPS contact has an invalid current surface measure");
    }
    const Scalar gap = local_value(input.geometry.primary_reference_radius + state[2]
                                       - input.geometry.secondary_reference_radius - state[3],
        0,
        request.jacobian);
    const Scalar secondary_temperature = local_value(state[1], 1, request.jacobian);
    const Scalar primary_temperature = local_value(state[0], 2, request.jacobian);
    const Scalar pressure = input.mechanical ? adlite::max(-input.normal.penalty * gap, Scalar(0.0)) : Scalar(0.0);
    const Scalar heat_flux =
        input.thermal ? contact_common::gap_conductance(input.heat, gap, secondary_temperature, primary_temperature)
                            * (secondary_temperature - primary_temperature)
                      : Scalar(0.0);
    const contact_common::ActivePoint3 normal{1.0, 0.0, 0.0};
    constexpr double gauss = 0.57735026918962576451;
    for (std::size_t q = 0; q < ring_gps_quadrature_point_count; ++q) {
        const double coordinate = q == 0 ? -gauss : gauss;
        const std::array<double, 2> shape{0.5 * (1.0 - coordinate), 0.5 * (1.0 + coordinate)};
        double slip_increment = 0.0;
        for (std::size_t i = 0; i < 2; ++i)
            slip_increment += shape[i] * ((state[6 + i] - old_state[6 + i]) - (state[4 + i] - old_state[4 + i]));
        const Scalar active_increment = local_value(slip_increment, 3, request.jacobian);
        const auto& history = input.committed_history[q];
        const auto friction = contact_common::friction_return(input.normal,
            history,
            {0.0, 0.0, history.elastic_tangential_slip},
            {0.0, 0.0, history.total_tangential_slip},
            {0.0, 0.0, active_increment},
            normal,
            pressure,
            point_area);
        const Scalar tangential_force = friction.tangential_traction_vector[2] * point_area;
        accumulate(result, 0, -heat_flux * point_area, shape, request);
        accumulate(result, 1, heat_flux * point_area, shape, request);
        accumulate(result, 2, -pressure * point_area, shape, request);
        accumulate(result, 3, pressure * point_area, shape, request);
        for (std::size_t i = 0; i < 2; ++i) {
            accumulate(result, 4 + i, -tangential_force * shape[i], shape, request);
            accumulate(result, 6 + i, tangential_force * shape[i], shape, request);
        }
        auto& point = result.points[q];
        point.gap = gap.value();
        point.pressure = pressure.value();
        point.heat_flux = heat_flux.value();
        point.weighted_measure = point_area.value();
        point.signed_tangential_traction = friction.tangential_traction_vector[2].value();
        point.elastic_tangential_slip = friction.elastic_tangential_slip[2].value();
        point.total_tangential_slip = friction.tangential_slip[2].value();
        point.friction_dissipation = friction.friction_dissipation.value();
        point.sliding = friction.sliding;
        for (const double value : {point.signed_tangential_traction,
                 point.elastic_tangential_slip,
                 point.total_tangential_slip,
                 point.friction_dissipation})
            require_finite(value, "Ring GPS contact has nonfinite frictional diagnostics");
        result.heat_rate += point.heat_flux * point_area.value();
        result.normal_force += point.pressure * point_area.value();
        result.tangential_force += tangential_force.value();
        result.friction_dissipation += point.friction_dissipation;
        if (request.history) {
            auto& trial = result.history[q];
            trial.elastic_tangential_slip = point.elastic_tangential_slip;
            trial.total_tangential_slip = point.total_tangential_slip;
            trial.sliding = point.sliding;
            trial.cartesian_elastic_tangential_slip[2] = point.elastic_tangential_slip;
            trial.cartesian_total_tangential_slip[2] = point.total_tangential_slip;
            trial.cartesian_tangent_basis_initialized = true;
            trial.cartesian_contact_normal = {1.0, 0.0, 0.0};
            trial.cartesian_contact_tangent_first = {0.0, 0.0, 1.0};
        }
    }
    for (const double value :
        {result.heat_rate, result.normal_force, result.tangential_force, result.friction_dissipation})
        require_finite(value, "Ring GPS contact has nonfinite integrated diagnostics");
    for (const double value : result.residual)
        require_finite(value, "Ring GPS contact has a nonfinite integrated residual");
    for (const double value : result.jacobian)
        require_finite(value, "Ring GPS contact has a nonfinite integrated tangent");
    return result;
}
} // namespace fuelsim::elements
