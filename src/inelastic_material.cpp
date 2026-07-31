#include "fuelsim/inelastic_material.hpp"

#include <algorithm>
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
    double relaxed_fraction;
};

struct CoupledUpdate final {
    double equivalent_stress;
    double stress_derivative;
    double plastic_increment;
    double plastic_increment_derivative;
    double creep_increment;
    double creep_increment_derivative;
};

struct CreepIncrement final {
    double value;
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

void validate_norton_properties(const NortonCreepProperties& creep) {
    if (!std::isfinite(creep.coefficient) || !(creep.coefficient >= 0.0))
        throw std::invalid_argument(
            "Norton creep coefficient must be finite and nonnegative");
    if (!std::isfinite(creep.reference_stress) ||
        !(creep.reference_stress > 0.0))
        throw std::invalid_argument(
            "Norton creep reference_stress must be finite and positive");
    if (!std::isfinite(creep.stress_exponent) ||
        !(creep.stress_exponent >= 1.0))
        throw std::invalid_argument(
            "Norton creep stress_exponent must be finite and at least one");
}

void validate_plasticity_properties(const J2PlasticityProperties& plasticity) {
    if (!std::isfinite(plasticity.yield_stress) ||
        !(plasticity.yield_stress > 0.0))
        throw std::invalid_argument(
            "J2 plasticity yield_stress must be finite and positive");
    if (!std::isfinite(plasticity.isotropic_hardening_modulus) ||
        !(plasticity.isotropic_hardening_modulus >= 0.0))
        throw std::invalid_argument(
            "J2 plasticity isotropic_hardening_modulus must be finite and "
            "nonnegative");
}

double log_add_exp(double first, double second) {
    const double maximum = std::max(first, second);
    if (std::isinf(maximum))
        return maximum;
    return maximum + std::log1p(std::exp(std::min(first, second) - maximum));
}

double second_exponential_weight(double first, double second) {
    if (second >= first)
        return 1.0 / (1.0 + std::exp(first - second));
    const double ratio = std::exp(second - first);
    return ratio / (1.0 + ratio);
}

NortonRoot solve_power_law_equivalent_stress(
    double driving_stress, double log_stress_coefficient, double time_step,
    const NortonCreepProperties& creep) {
    if (time_step == 0.0 || creep.coefficient == 0.0 ||
        log_stress_coefficient == -std::numeric_limits<double>::infinity())
        return {driving_stress, 1.0, 0.0};
    if (!(driving_stress > 0.0)) {
        if (creep.stress_exponent > 1.0)
            return {driving_stress, 1.0, 0.0};

        const double log_linear_coefficient =
            log_stress_coefficient + std::log(time_step) +
            std::log(creep.coefficient) - std::log(creep.reference_stress);
        if (std::isnan(log_linear_coefficient))
            throw std::overflow_error(
                "Norton creep linear coefficient is not a number");
        if (log_linear_coefficient == -std::numeric_limits<double>::infinity())
            return {driving_stress, 1.0, 0.0};
        if (log_linear_coefficient == std::numeric_limits<double>::infinity())
            return {driving_stress, 0.0, 0.0};
        const double log_denominator = log_add_exp(0.0, log_linear_coefficient);
        return {
            driving_stress,
            std::exp(-log_denominator),
            0.0,
        };
    }

    const double log_driving_stress = std::log(driving_stress);
    const double log_coefficient =
        log_stress_coefficient + std::log(time_step) +
        std::log(creep.coefficient) - log_driving_stress +
        creep.stress_exponent *
            (log_driving_stress - std::log(creep.reference_stress));
    if (std::isnan(log_coefficient))
        throw std::overflow_error(
            "Norton creep dimensionless local coefficient is not a number");
    if (log_coefficient == -std::numeric_limits<double>::infinity())
        return {driving_stress, 1.0, 0.0};
    if (log_coefficient == std::numeric_limits<double>::infinity())
        throw std::overflow_error(
            "Norton creep logarithmic local coefficient overflowed");

    double log_stress_ratio = 0.0;
    if (creep.stress_exponent == 1.0) {
        log_stress_ratio = log_add_exp(0.0, log_coefficient);
    } else {
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
        const double initial_residual = log_add_exp(0.0, log_coefficient);
        if (initial_residual == 0.0)
            return {driving_stress, 1.0, 0.0};

        constexpr double log_two = 0.693147180559945309417232121458176568;
        double lower = 0.0;
        double upper = std::max(log_two, (log_coefficient + log_two) /
                                             creep.stress_exponent);
        double current = log_coefficient > 0.0
                             ? log_coefficient / creep.stress_exponent
                             : std::exp(log_coefficient);
        if (!std::isfinite(upper) || !(upper > lower))
            throw std::overflow_error(
                "Norton creep logarithmic root bracket is invalid");
        if (!std::isfinite(current) || !(current > lower) || !(current < upper))
            current = 0.5 * (lower + upper);

        bool converged = false;
        for (int iteration = 0; iteration < maximum_creep_iterations;
             ++iteration) {
            const double first = -current;
            const double second =
                log_coefficient - creep.stress_exponent * current;
            const double residual = log_add_exp(first, second);
            if (!std::isfinite(residual))
                throw std::overflow_error(
                    "Norton creep logarithmic residual is not finite");
            if (std::fabs(residual) <= tolerance) {
                converged = true;
                break;
            }

            if (residual > 0.0)
                lower = current;
            else
                upper = current;
            if (upper - lower <=
                tolerance * std::max(1.0, std::fabs(current))) {
                current = 0.5 * (lower + upper);
                converged = true;
                break;
            }

            const double weight = second_exponential_weight(first, second);
            const double derivative =
                -(1.0 + (creep.stress_exponent - 1.0) * weight);
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
                "Norton creep logarithmic local Newton solve did not "
                "converge");
        log_stress_ratio = current;
    }

    const double stress_ratio = std::exp(-log_stress_ratio);
    const double relaxed_fraction = -std::expm1(-log_stress_ratio);
    const double log_equivalent_stress = log_driving_stress - log_stress_ratio;
    const double minimum_log =
        std::log(std::numeric_limits<double>::denorm_min());
    const double equivalent_stress = log_equivalent_stress <= minimum_log
                                         ? 0.0
                                         : std::exp(log_equivalent_stress);
    const double derivative_denominator =
        stress_ratio + creep.stress_exponent * relaxed_fraction;
    const double driving_stress_derivative =
        derivative_denominator == 0.0 ? 0.0
                                      : stress_ratio / derivative_denominator;
    if (!std::isfinite(equivalent_stress) ||
        !std::isfinite(driving_stress_derivative) ||
        !(driving_stress_derivative >= 0.0) ||
        !std::isfinite(relaxed_fraction) || !(relaxed_fraction >= 0.0) ||
        !(relaxed_fraction <= 1.0))
        throw std::overflow_error(
            "Norton creep logarithmic root or derivative is not finite");

    return {
        equivalent_stress,
        driving_stress_derivative,
        relaxed_fraction,
    };
}

NortonRoot solve_norton_equivalent_stress(double trial_stress,
                                          double shear_modulus,
                                          double time_step,
                                          const NortonCreepProperties& creep) {
    return solve_power_law_equivalent_stress(
        trial_stress, std::log(3.0) + std::log(shear_modulus), time_step,
        creep);
}

CreepIncrement evaluate_creep_increment(double equivalent_stress,
                                        double time_step,
                                        const NortonCreepProperties& creep) {
    if (time_step == 0.0 || creep.coefficient == 0.0)
        return {0.0};
    if (!(equivalent_stress > 0.0))
        return {0.0};

    const double log_increment =
        std::log(time_step) + std::log(creep.coefficient) +
        creep.stress_exponent *
            (std::log(equivalent_stress) - std::log(creep.reference_stress));
    if (std::isnan(log_increment))
        throw std::overflow_error(
            "Norton creep increment logarithm is not a number");

    const double minimum_log =
        std::log(std::numeric_limits<double>::denorm_min());
    const double maximum_log = std::log(std::numeric_limits<double>::max());
    if (log_increment <= minimum_log)
        return {0.0};
    if (!std::isfinite(log_increment) || log_increment >= maximum_log)
        throw std::overflow_error("Norton creep increment overflowed");

    const double increment = std::exp(log_increment);
    if (!std::isfinite(increment) || !(increment >= 0.0))
        throw std::overflow_error("Norton creep increment is not finite");
    return {increment};
}

CoupledUpdate solve_coupled_update(double trial_stress, double shear_modulus,
                                   double time_step,
                                   const NortonCreepProperties& creep,
                                   const J2PlasticityProperties& plasticity,
                                   double committed_equivalent_plastic_strain) {
    const double current_yield_stress =
        plasticity.yield_stress + plasticity.isotropic_hardening_modulus *
                                      committed_equivalent_plastic_strain;
    if (!std::isfinite(current_yield_stress))
        throw std::overflow_error(
            "Coupled plastic-creep current yield stress is not finite");

    const NortonRoot creep_only = solve_norton_equivalent_stress(
        trial_stress, shear_modulus, time_step, creep);
    const double inverse_three_shear_modulus = (1.0 / 3.0) / shear_modulus;
    if (!(creep_only.equivalent_stress > current_yield_stress)) {
        const double creep_increment = trial_stress *
                                       creep_only.relaxed_fraction *
                                       inverse_three_shear_modulus;
        const double creep_increment_derivative =
            (1.0 - creep_only.trial_stress_derivative) *
            inverse_three_shear_modulus;
        return {
            creep_only.equivalent_stress,
            creep_only.trial_stress_derivative,
            0.0,
            0.0,
            creep_increment,
            creep_increment_derivative,
        };
    }

    const double hardening = plasticity.isotropic_hardening_modulus;
    if (hardening == 0.0) {
        const CreepIncrement creep_increment =
            evaluate_creep_increment(current_yield_stress, time_step, creep);
        const double plastic_increment = (trial_stress - current_yield_stress) *
                                             inverse_three_shear_modulus -
                                         creep_increment.value;
        if (!std::isfinite(plastic_increment) || !(plastic_increment > 0.0))
            throw std::runtime_error(
                "Coupled perfect-plastic update produced a nonpositive "
                "plastic increment");
        return {
            current_yield_stress,  0.0,
            plastic_increment,     inverse_three_shear_modulus,
            creep_increment.value, 0.0,
        };
    }

    double hardening_weight = 0.0;
    double inverse_denominator = 0.0;
    double log_root_stress_coefficient = 0.0;
    if (hardening <= shear_modulus) {
        const double ratio = hardening / shear_modulus;
        const double scaled_denominator = 3.0 + ratio;
        hardening_weight = ratio / scaled_denominator;
        inverse_denominator = (1.0 / scaled_denominator) / shear_modulus;
        log_root_stress_coefficient =
            std::log(hardening) + std::log(3.0 / scaled_denominator);
    } else {
        const double ratio = shear_modulus / hardening;
        const double scaled_denominator = 1.0 + 3.0 * ratio;
        hardening_weight = 1.0 / scaled_denominator;
        inverse_denominator = (1.0 / scaled_denominator) / hardening;
        log_root_stress_coefficient =
            std::log(shear_modulus) + std::log(3.0 / scaled_denominator);
    }
    const double root_driving_stress =
        current_yield_stress +
        hardening_weight * (trial_stress - current_yield_stress);
    if (!std::isfinite(inverse_denominator) || !(inverse_denominator > 0.0) ||
        !std::isfinite(root_driving_stress) || !(root_driving_stress > 0.0) ||
        !std::isfinite(log_root_stress_coefficient))
        throw std::overflow_error(
            "Coupled plastic-creep transformed root is invalid");

    const NortonRoot root = solve_power_law_equivalent_stress(
        root_driving_stress, log_root_stress_coefficient, time_step, creep);
    const CreepIncrement creep_increment =
        evaluate_creep_increment(root.equivalent_stress, time_step, creep);
    const double plastic_increment =
        (trial_stress - current_yield_stress -
         3.0 * (shear_modulus * creep_increment.value)) *
        inverse_denominator;
    const double plastic_increment_derivative =
        root.trial_stress_derivative * inverse_denominator;
    const double stress_derivative =
        hardening_weight * root.trial_stress_derivative;
    const double creep_increment_derivative =
        (1.0 - root.trial_stress_derivative) * inverse_three_shear_modulus;
    if (!std::isfinite(root.equivalent_stress) ||
        !std::isfinite(stress_derivative) ||
        !std::isfinite(plastic_increment) || !(plastic_increment > 0.0) ||
        !std::isfinite(plastic_increment_derivative) ||
        !std::isfinite(creep_increment.value) ||
        !std::isfinite(creep_increment_derivative))
        throw std::overflow_error(
            "Coupled plastic-creep update or derivative is not finite");

    return {
        root.equivalent_stress, stress_derivative,
        plastic_increment,      plastic_increment_derivative,
        creep_increment.value,  creep_increment_derivative,
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
        validate_norton_properties(_properties.creep);
        break;
    case InelasticBehavior::j2_plasticity:
        validate_plasticity_properties(_properties.plasticity);
        break;
    case InelasticBehavior::norton_creep_j2_plasticity:
        validate_norton_properties(_properties.creep);
        validate_plasticity_properties(_properties.plasticity);
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
    const bool zero_deviatoric_stress = deviatoric_trial[0].value() == 0.0 &&
                                        deviatoric_trial[1].value() == 0.0 &&
                                        deviatoric_trial[2].value() == 0.0 &&
                                        deviatoric_trial[3].value() == 0.0;
    constexpr double square_root_two = 1.414213562373095048801688724209698079;
    constexpr double square_root_three_halves =
        1.224744871391589049098642037352945695;
    const adlite::Scalar equivalent_trial_stress =
        zero_deviatoric_stress
            ? adlite::Scalar(0.0)
            : square_root_three_halves *
                  adlite::hypot(
                      adlite::hypot(deviatoric_trial[0], deviatoric_trial[1]),
                      adlite::hypot(deviatoric_trial[2],
                                    square_root_two * deviatoric_trial[3]));
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

    if (_properties.behavior == InelasticBehavior::norton_creep_j2_plasticity) {
        const CoupledUpdate update = solve_coupled_update(
            equivalent_trial_stress.value(), _shear_modulus, time_step,
            _properties.creep, _properties.plasticity,
            committed.equivalent_plastic_strain);
        if (equivalent_trial_stress.value() == 0.0) {
            const adlite::Scalar stress_scale = update.stress_derivative;
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
            adlite::compose(update.equivalent_stress, equivalent_trial_stress,
                            update.stress_derivative);
        const adlite::Scalar plastic_increment =
            adlite::compose(update.plastic_increment, equivalent_trial_stress,
                            update.plastic_increment_derivative);
        const adlite::Scalar creep_increment =
            adlite::compose(update.creep_increment, equivalent_trial_stress,
                            update.creep_increment_derivative);
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
            trial_state.creep_strain[component] =
                committed.creep_strain[component] +
                creep_increment * flow_direction;
        }
        trial_state.equivalent_plastic_strain =
            committed.equivalent_plastic_strain + plastic_increment;
        trial_state.equivalent_creep_strain =
            committed.equivalent_creep_strain + creep_increment;
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
    const double inverse_three_shear_modulus = (1.0 / 3.0) / _shear_modulus;
    const double creep_increment_value = equivalent_trial_stress.value() *
                                         root.relaxed_fraction *
                                         inverse_three_shear_modulus;
    const double creep_increment_derivative =
        (1.0 - root.trial_stress_derivative) * inverse_three_shear_modulus;
    const adlite::Scalar creep_increment =
        adlite::compose(creep_increment_value, equivalent_trial_stress,
                        creep_increment_derivative);
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
