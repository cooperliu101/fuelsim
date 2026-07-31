#include "fuelsim/inelastic_material.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace fuelsim {
namespace {

constexpr std::size_t component_count = 4;
constexpr int maximum_creep_iterations = 100;

struct NortonRoot final {
    double equivalent_stress;
    double trial_stress_derivative;
};

void validate_material_point_state(const MaterialPointState& state) {
    for (std::size_t component = 0; component < component_count; ++component) {
        if (!std::isfinite(state.plastic_strain[component]) ||
            !std::isfinite(state.creep_strain[component]))
            throw std::domain_error(
                "Material point inelastic strains must be finite");
    }

    if (!std::isfinite(state.equivalent_plastic_strain) ||
        !(state.equivalent_plastic_strain >= 0.0))
        throw std::domain_error(
            "Material point equivalent plastic strain must be finite and "
            "nonnegative");
    if (!std::isfinite(state.equivalent_creep_strain) ||
        !(state.equivalent_creep_strain >= 0.0))
        throw std::domain_error(
            "Material point equivalent creep strain must be finite and "
            "nonnegative");
}

MaterialPointTrialState
passive_trial_state(const MaterialPointState& committed) {
    MaterialPointTrialState trial;
    for (std::size_t component = 0; component < component_count; ++component) {
        trial.plastic_strain[component] = committed.plastic_strain[component];
        trial.creep_strain[component] = committed.creep_strain[component];
    }
    trial.equivalent_plastic_strain = committed.equivalent_plastic_strain;
    trial.equivalent_creep_strain = committed.equivalent_creep_strain;
    return trial;
}

NortonRoot solve_norton_equivalent_stress(double trial_stress,
                                          double shear_modulus,
                                          double time_step,
                                          const NortonCreepProperties& creep) {
    if (time_step == 0.0 || creep.coefficient == 0.0)
        return {trial_stress, 1.0};
    if (!(trial_stress > 0.0)) {
        if (creep.stress_exponent > 1.0)
            return {trial_stress, 1.0};

        const double log_linear_coefficient =
            std::log(3.0) + std::log(shear_modulus) + std::log(time_step) +
            std::log(creep.coefficient) - std::log(creep.reference_stress);
        const double minimum_log =
            std::log(std::numeric_limits<double>::denorm_min());
        const double maximum_log = std::log(std::numeric_limits<double>::max());
        if (log_linear_coefficient <= minimum_log)
            return {trial_stress, 1.0};
        if (log_linear_coefficient >= maximum_log)
            return {trial_stress, 0.0};
        const double linear_coefficient = std::exp(log_linear_coefficient);
        return {
            trial_stress,
            1.0 / (1.0 + linear_coefficient),
        };
    }

    const double log_coefficient =
        std::log(3.0) + std::log(shear_modulus) + std::log(time_step) +
        std::log(creep.coefficient) - std::log(trial_stress) +
        creep.stress_exponent *
            (std::log(trial_stress) - std::log(creep.reference_stress));
    if (std::isnan(log_coefficient))
        throw std::overflow_error(
            "Norton creep dimensionless local coefficient is not a number");

    const double minimum_log =
        std::log(std::numeric_limits<double>::denorm_min());
    const double maximum_log = std::log(std::numeric_limits<double>::max());
    if (log_coefficient <= minimum_log)
        return {trial_stress, 1.0};
    if (!std::isfinite(log_coefficient) || log_coefficient >= maximum_log)
        throw std::overflow_error(
            "Norton creep dimensionless local coefficient overflowed");
    const double coefficient = std::exp(log_coefficient);

    if (creep.stress_exponent == 1.0) {
        const double denominator = 1.0 + coefficient;
        if (!std::isfinite(denominator))
            throw std::overflow_error(
                "Norton creep linear local equation overflowed");
        return {
            trial_stress / denominator,
            1.0 / denominator,
        };
    }

    double lower = 0.0;
    double upper = 1.0;
    double current = coefficient > 1.0
                         ? std::exp(-log_coefficient / creep.stress_exponent)
                         : 1.0;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();

    bool converged = false;
    for (int iteration = 0; iteration < maximum_creep_iterations; ++iteration) {
        const double power = std::pow(current, creep.stress_exponent);
        const double scaled_power = coefficient * power;
        const double residual = current + scaled_power - 1.0;

        if (std::isfinite(residual) && std::fabs(residual) <= tolerance) {
            converged = true;
            break;
        }

        if (residual > 0.0 || !std::isfinite(residual))
            upper = current;
        else
            lower = current;

        if (upper - lower <= tolerance) {
            current = 0.5 * (lower + upper);
            converged = true;
            break;
        }

        const double derivative =
            1.0 + creep.stress_exponent * scaled_power / current;
        double candidate = std::isfinite(derivative)
                               ? current - residual / derivative
                               : 0.5 * (lower + upper);
        if (!std::isfinite(candidate) || !(candidate > lower) ||
            !(candidate < upper) || candidate == current)
            candidate = 0.5 * (lower + upper);
        current = candidate;
    }

    if (!converged)
        throw std::runtime_error(
            "Norton creep local Newton solve did not converge");

    const double derivative_term =
        creep.stress_exponent * ((1.0 - current) / current);
    const double trial_stress_derivative =
        std::isinf(derivative_term) ? 0.0 : 1.0 / (1.0 + derivative_term);
    const double equivalent_stress = trial_stress * current;
    if (std::isnan(derivative_term) ||
        !std::isfinite(trial_stress_derivative) ||
        !(trial_stress_derivative >= 0.0) || !std::isfinite(equivalent_stress))
        throw std::overflow_error(
            "Norton creep consistent derivative is not finite");

    return {
        equivalent_stress,
        trial_stress_derivative,
    };
}

AxisymmetricStress returned_stress(
    const adlite::Scalar& mean_stress,
    const std::array<adlite::Scalar, component_count>& deviatoric_trial,
    const adlite::Scalar& scale) {
    return {
        mean_stress + scale * deviatoric_trial[0],
        mean_stress + scale * deviatoric_trial[1],
        mean_stress + scale * deviatoric_trial[2],
        scale * deviatoric_trial[3],
    };
}

} // namespace

IsotropicInelasticMaterial::IsotropicInelasticMaterial(
    ThermoelasticProperties thermoelastic_properties,
    TransientInelasticProperties properties)
    : _thermoelastic_material(thermoelastic_properties),
      _properties(properties),
      _lame_lambda(thermoelastic_properties.young_modulus *
                   thermoelastic_properties.poisson_ratio /
                   ((1.0 + thermoelastic_properties.poisson_ratio) *
                    (1.0 - 2.0 * thermoelastic_properties.poisson_ratio))),
      _shear_modulus(thermoelastic_properties.young_modulus /
                     (2.0 * (1.0 + thermoelastic_properties.poisson_ratio))) {
    if (!std::isfinite(_properties.density) || !(_properties.density > 0.0))
        throw std::invalid_argument(
            "Inelastic material density must be finite and positive");
    if (!std::isfinite(_properties.specific_heat) ||
        !(_properties.specific_heat > 0.0))
        throw std::invalid_argument(
            "Inelastic material specific_heat must be finite and positive");

    switch (_properties.behavior) {
    case InelasticBehavior::elastic:
        break;
    case InelasticBehavior::norton_creep:
        if (!std::isfinite(_properties.creep.coefficient) ||
            !(_properties.creep.coefficient >= 0.0))
            throw std::invalid_argument(
                "Norton creep coefficient must be finite and nonnegative");
        if (!std::isfinite(_properties.creep.reference_stress) ||
            !(_properties.creep.reference_stress > 0.0))
            throw std::invalid_argument(
                "Norton creep reference_stress must be finite and positive");
        if (!std::isfinite(_properties.creep.stress_exponent) ||
            !(_properties.creep.stress_exponent >= 1.0))
            throw std::invalid_argument(
                "Norton creep stress_exponent must be finite and at least "
                "one");
        break;
    case InelasticBehavior::j2_plasticity:
        if (!std::isfinite(_properties.plasticity.yield_stress) ||
            !(_properties.plasticity.yield_stress > 0.0))
            throw std::invalid_argument(
                "J2 plasticity yield_stress must be finite and "
                "positive");
        if (!std::isfinite(
                _properties.plasticity.isotropic_hardening_modulus) ||
            !(_properties.plasticity.isotropic_hardening_modulus >= 0.0))
            throw std::invalid_argument(
                "J2 plasticity isotropic_hardening_modulus must be finite "
                "and nonnegative");
        break;
    default:
        throw std::invalid_argument(
            "Inelastic material behavior is not supported");
    }
}

const TransientInelasticProperties&
IsotropicInelasticMaterial::properties() const noexcept {
    return _properties;
}

adlite::Scalar IsotropicInelasticMaterial::conductivity(
    const adlite::Scalar& temperature) const {
    return _thermoelastic_material.conductivity(temperature);
}

InelasticStressResponse IsotropicInelasticMaterial::response(
    const adlite::Scalar& strain_rr, const adlite::Scalar& strain_zz,
    const adlite::Scalar& strain_hoop, const adlite::Scalar& strain_rz,
    const adlite::Scalar& temperature, double time_step,
    const MaterialPointState& committed) const {
    if (!std::isfinite(strain_rr.value()) ||
        !std::isfinite(strain_zz.value()) ||
        !std::isfinite(strain_hoop.value()) ||
        !std::isfinite(strain_rz.value()) ||
        !std::isfinite(temperature.value()))
        throw std::domain_error(
            "Inelastic material strain and temperature inputs must be "
            "finite");
    if (!std::isfinite(time_step) || !(time_step >= 0.0))
        throw std::domain_error(
            "Inelastic material time_step must be finite and nonnegative");
    validate_material_point_state(committed);

    MaterialPointTrialState trial_state = passive_trial_state(committed);
    const ThermoelasticProperties& thermoelastic =
        _thermoelastic_material.properties();
    const adlite::Scalar thermal_strain =
        thermoelastic.thermal_expansion *
        (temperature - thermoelastic.reference_temperature);
    const std::array<adlite::Scalar, component_count> elastic_strain = {
        strain_rr - thermal_strain - committed.plastic_strain[0] -
            committed.creep_strain[0],
        strain_zz - thermal_strain - committed.plastic_strain[1] -
            committed.creep_strain[1],
        strain_hoop - thermal_strain - committed.plastic_strain[2] -
            committed.creep_strain[2],
        strain_rz - committed.plastic_strain[3] - committed.creep_strain[3],
    };
    const adlite::Scalar elastic_trace =
        elastic_strain[0] + elastic_strain[1] + elastic_strain[2];
    const AxisymmetricStress stress_trial = {
        _lame_lambda * elastic_trace + 2.0 * _shear_modulus * elastic_strain[0],
        _lame_lambda * elastic_trace + 2.0 * _shear_modulus * elastic_strain[1],
        _lame_lambda * elastic_trace + 2.0 * _shear_modulus * elastic_strain[2],
        2.0 * _shear_modulus * elastic_strain[3],
    };
    const adlite::Scalar mean_stress =
        (stress_trial.rr + stress_trial.zz + stress_trial.hoop) / 3.0;
    const std::array<adlite::Scalar, component_count> deviatoric_trial = {
        stress_trial.rr - mean_stress,
        stress_trial.zz - mean_stress,
        stress_trial.hoop - mean_stress,
        stress_trial.rz,
    };
    const adlite::Scalar deviatoric_norm_squared =
        deviatoric_trial[0] * deviatoric_trial[0] +
        deviatoric_trial[1] * deviatoric_trial[1] +
        deviatoric_trial[2] * deviatoric_trial[2] +
        2.0 * deviatoric_trial[3] * deviatoric_trial[3];
    if (!std::isfinite(deviatoric_norm_squared.value()) ||
        !(deviatoric_norm_squared.value() >= 0.0))
        throw std::overflow_error(
            "Inelastic material trial deviatoric stress is not finite");

    const adlite::Scalar equivalent_trial_stress =
        deviatoric_norm_squared.value() == 0.0
            ? adlite::Scalar(0.0)
            : adlite::sqrt(1.5 * deviatoric_norm_squared);
    if (!std::isfinite(equivalent_trial_stress.value()))
        throw std::overflow_error(
            "Inelastic material trial equivalent stress is not finite");

    if (_properties.behavior == InelasticBehavior::elastic)
        return {stress_trial, trial_state};

    if (_properties.behavior == InelasticBehavior::j2_plasticity) {
        if (equivalent_trial_stress.value() == 0.0)
            return {stress_trial, trial_state};

        const double hardening =
            _properties.plasticity.isotropic_hardening_modulus;
        const double current_yield_stress =
            _properties.plasticity.yield_stress +
            hardening * committed.equivalent_plastic_strain;
        if (!std::isfinite(current_yield_stress))
            throw std::overflow_error(
                "J2 plasticity current yield stress is not finite");

        const double yield_function =
            equivalent_trial_stress.value() - current_yield_stress;
        if (!(yield_function > 0.0))
            return {stress_trial, trial_state};

        const double denominator = 3.0 * _shear_modulus + hardening;
        const adlite::Scalar plastic_increment =
            (equivalent_trial_stress - current_yield_stress) / denominator;
        const adlite::Scalar returned_equivalent_stress =
            equivalent_trial_stress - 3.0 * _shear_modulus * plastic_increment;
        const adlite::Scalar stress_scale =
            returned_equivalent_stress / equivalent_trial_stress;
        const AxisymmetricStress stress =
            returned_stress(mean_stress, deviatoric_trial, stress_scale);

        for (std::size_t component = 0; component < component_count;
             ++component) {
            const adlite::Scalar flow_direction =
                1.5 * deviatoric_trial[component] / equivalent_trial_stress;
            trial_state.plastic_strain[component] =
                committed.plastic_strain[component] +
                plastic_increment * flow_direction;
        }
        trial_state.equivalent_plastic_strain =
            committed.equivalent_plastic_strain + plastic_increment;
        return {stress, trial_state};
    }

    if (time_step == 0.0 || _properties.creep.coefficient == 0.0)
        return {stress_trial, trial_state};

    const NortonRoot root = solve_norton_equivalent_stress(
        equivalent_trial_stress.value(), _shear_modulus, time_step,
        _properties.creep);
    if (equivalent_trial_stress.value() == 0.0) {
        const adlite::Scalar stress_scale = root.trial_stress_derivative;
        const AxisymmetricStress stress =
            returned_stress(mean_stress, deviatoric_trial, stress_scale);
        for (std::size_t component = 0; component < component_count;
             ++component) {
            trial_state.creep_strain[component] =
                committed.creep_strain[component] +
                (1.0 - stress_scale) * deviatoric_trial[component] /
                    (2.0 * _shear_modulus);
        }
        return {stress, trial_state};
    }

    const adlite::Scalar returned_equivalent_stress =
        adlite::compose(root.equivalent_stress, equivalent_trial_stress,
                        root.trial_stress_derivative);
    const adlite::Scalar creep_increment =
        (equivalent_trial_stress - returned_equivalent_stress) /
        (3.0 * _shear_modulus);
    const adlite::Scalar stress_scale =
        returned_equivalent_stress / equivalent_trial_stress;
    const AxisymmetricStress stress =
        returned_stress(mean_stress, deviatoric_trial, stress_scale);

    for (std::size_t component = 0; component < component_count; ++component) {
        const adlite::Scalar flow_direction =
            1.5 * deviatoric_trial[component] / equivalent_trial_stress;
        trial_state.creep_strain[component] =
            committed.creep_strain[component] +
            creep_increment * flow_direction;
    }
    trial_state.equivalent_creep_strain =
        committed.equivalent_creep_strain + creep_increment;
    return {stress, trial_state};
}

MaterialPointState IsotropicInelasticMaterial::state_values(
    const MaterialPointTrialState& trial_state) {
    MaterialPointState state;
    for (std::size_t component = 0; component < component_count; ++component) {
        state.plastic_strain[component] =
            trial_state.plastic_strain[component].value();
        state.creep_strain[component] =
            trial_state.creep_strain[component].value();
    }
    state.equivalent_plastic_strain =
        trial_state.equivalent_plastic_strain.value();
    state.equivalent_creep_strain = trial_state.equivalent_creep_strain.value();
    validate_material_point_state(state);
    return state;
}

} // namespace fuelsim
