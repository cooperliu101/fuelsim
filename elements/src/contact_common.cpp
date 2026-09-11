#include "contact_common.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim::contact_common {
namespace {
adlite::Scalar dot(const ActivePoint3& a, const ActivePoint3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

adlite::Scalar norm(const ActivePoint3& a) {
    return adlite::hypot(adlite::hypot(a[0], a[1]), a[2]);
}
} // namespace

FrictionResult friction_return(const NormalContactProperties& properties,
    const ContactPointHistory& history,
    const ActivePoint3& transported_history,
    const ActivePoint3& transported_total_history,
    const ActivePoint3& relative_increment,
    const ActivePoint3& normal,
    const adlite::Scalar& pressure,
    const adlite::Scalar& tributary_area) {
    FrictionResult result;
    if (properties.friction_coefficient == 0.0)
        return result;
    result.tangential_slip = transported_total_history;
    if (!(pressure.value() > 0.0))
        return result;
    const adlite::Scalar normal_increment = dot(relative_increment, normal);
    for (std::size_t component = 0; component < 3; ++component) {
        result.tangential_slip[component] =
            transported_total_history[component] + relative_increment[component] - normal_increment * normal[component];
        result.elastic_tangential_slip[component] =
            transported_history[component] + relative_increment[component] - normal_increment * normal[component];
    }
    const adlite::Scalar total_history_normal = dot(result.tangential_slip, normal);
    for (std::size_t component = 0; component < 3; ++component)
        result.tangential_slip[component] -= total_history_normal * normal[component];
    const adlite::Scalar history_normal = dot(result.elastic_tangential_slip, normal);
    for (std::size_t component = 0; component < 3; ++component)
        result.elastic_tangential_slip[component] -= history_normal * normal[component];
    const ActivePoint3 trial_elastic_tangential_slip = result.elastic_tangential_slip;
    const adlite::Scalar sliding_limit = properties.friction_coefficient * pressure,
                         stick_stiffness = properties.maximum_elastic_slip > 0.0
                                               ? sliding_limit / properties.maximum_elastic_slip
                                               : adlite::Scalar(properties.penalty);
    if (!std::isfinite(stick_stiffness.value()) || !(stick_stiffness.value() > 0.0))
        throw std::domain_error("Contact has a nonfinite or nonpositive tangential stick stiffness");
    ActivePoint3 trial_traction{};
    for (std::size_t component = 0; component < 3; ++component)
        trial_traction[component] = stick_stiffness * result.elastic_tangential_slip[component];
    const adlite::Scalar trial_magnitude = norm(trial_traction);
    if (trial_magnitude.value() < sliding_limit.value()
        || (trial_magnitude.value() == sliding_limit.value() && !history.sliding)) {
        result.tangential_traction_vector = trial_traction;
        result.tangential_traction = trial_magnitude;
    } else {
        if (!(trial_magnitude.value() > 0.0))
            throw std::domain_error("Contact has an undefined tangential direction");
        for (std::size_t component = 0; component < 3; ++component) {
            result.tangential_traction_vector[component] = sliding_limit * trial_traction[component] / trial_magnitude;
            result.elastic_tangential_slip[component] = result.tangential_traction_vector[component] / stick_stiffness;
        }
        result.tangential_traction = sliding_limit;
        result.sliding = true;
        for (std::size_t component = 0; component < 3; ++component)
            result.friction_dissipation +=
                result.tangential_traction_vector[component]
                * (trial_elastic_tangential_slip[component] - result.elastic_tangential_slip[component]);
        result.friction_dissipation *= tributary_area;
    }
    result.tangential_force = result.tangential_traction * tributary_area;

    return result;
}

bool projection_increment_converged(const adlite::Scalar& delta_xi,
    const adlite::Scalar& delta_eta,
    const adlite::Scalar& xi,
    const adlite::Scalar& eta,
    std::size_t maximum_width) {
    const std::size_t width = delta_xi.derivative_size();
    for (const auto* value : {&delta_xi, &delta_eta, &xi, &eta}) {
        if (value->derivative_size() != width || width > maximum_width)
            throw std::logic_error("Contact projection derivative widths are inconsistent");
        if (!std::isfinite(value->value()))
            throw std::domain_error("Contact projection has a nonfinite coordinate or increment");
        for (std::size_t i = 0; i < width; ++i)
            if (!std::isfinite(value->derivative(i)))
                throw std::domain_error("Contact projection has a nonfinite coordinate or increment derivative");
    }
    constexpr double value_tolerance = 64.0 * std::numeric_limits<double>::epsilon();
    const double scale = std::max({1.0, std::abs(xi.value()), std::abs(eta.value())});
    if (std::max(std::abs(delta_xi.value()), std::abs(delta_eta.value())) > value_tolerance * scale)
        return false;
    for (std::size_t i = 0; i < width; ++i) {
        const double derivative_scale = std::max({1.0, std::abs(xi.derivative(i)), std::abs(eta.derivative(i))});
        if (std::max(std::abs(delta_xi.derivative(i)), std::abs(delta_eta.derivative(i))) > 1e-12 * derivative_scale)
            return false;
    }
    return true;
}

adlite::Scalar gap_conductance(const GapHeatProperties& properties,
    const adlite::Scalar& gap,
    const adlite::Scalar& secondary_temperature,
    const adlite::Scalar& primary_temperature) {
    if (!std::isfinite(gap.value()))
        throw std::domain_error("Contact thermal gap must be finite");
    adlite::Scalar result;
    if (properties.law == GapHeatConductanceLaw::gas_gap) {
        const adlite::Scalar thermal_gap = adlite::max(gap, adlite::Scalar(properties.minimum_gap));
        result = properties.gap_conductivity / thermal_gap;
    } else {
        const adlite::Scalar pressure = adlite::max(-properties.contact_penalty * gap, adlite::Scalar(0.0));
        const adlite::Scalar average_temperature = 0.5 * (secondary_temperature + primary_temperature);
        result = properties.conductance + properties.clearance_derivative * gap
                 + properties.pressure_derivative * pressure
                 + properties.temperature_derivative * (average_temperature - properties.reference_temperature);
    }
    if (!std::isfinite(result.value()) || result.value() < 0.0)
        throw std::domain_error("Contact gap conductance must be finite and nonnegative");
    return result;
}
} // namespace fuelsim::contact_common
