#include "material.hpp"
#include "material_functions.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace fuelsim {
AxisymmetricStress rotate_axisymmetric_tensor(const AxisymmetricStress& tensor, const AxisymmetricRotation& rotation) {
    return {
        rotation.rr * rotation.rr * tensor.rr + rotation.rz * rotation.rz * tensor.zz
            + 2.0 * rotation.rr * rotation.rz * tensor.rz,
        rotation.zr * rotation.zr * tensor.rr + rotation.zz * rotation.zz * tensor.zz
            + 2.0 * rotation.zr * rotation.zz * tensor.rz,
        rotation.hoop * rotation.hoop * tensor.hoop,
        rotation.rr * rotation.zr * tensor.rr + rotation.rz * rotation.zz * tensor.zz
            + (rotation.rr * rotation.zz + rotation.rz * rotation.zr) * tensor.rz,
    };
}

std::array<double, 4> rotate_axisymmetric_tensor_values(const std::array<double, 4>& tensor,
    const AxisymmetricRotation& rotation) {
    const double rr = rotation.rr.value(), rz = rotation.rz.value();
    const double zr = rotation.zr.value(), zz = rotation.zz.value(), hoop = rotation.hoop.value();
    return {
        rr * rr * tensor[0] + rz * rz * tensor[1] + 2.0 * rr * rz * tensor[3],
        zr * zr * tensor[0] + zz * zz * tensor[1] + 2.0 * zr * zz * tensor[3],
        hoop * hoop * tensor[2],
        rr * zr * tensor[0] + rz * zz * tensor[1] + (rr * zz + rz * zr) * tensor[3],
    };
}

IsotropicThermoelasticMaterial::IsotropicThermoelasticMaterial(ThermoelasticProperties properties)
    : _properties(std::move(properties)) {
    if (!_properties.functions || _properties.functions->thermal.function == nullptr)
        throw std::invalid_argument("Material function set requires a thermal function");
    if (_properties.functions->elasticity.function == nullptr)
        throw std::invalid_argument("Material function set requires an elasticity function");
    if (!std::isfinite(_properties.reference_young_modulus) || !(_properties.reference_young_modulus > 0.0))
        throw std::invalid_argument("Material reference Young modulus must be finite and positive");
}

adlite::Scalar IsotropicThermoelasticMaterial::conductivity(const adlite::Scalar& temperature,
    MaterialFunctionContext context) const {
    if (!std::isfinite(temperature.value()) || !(temperature.value() > 0.0))
        throw std::domain_error("Thermoelastic material temperature must be finite and positive");
    ThermalPropertyOutput output{};
    const ThermalFunctionInstance& instance = _properties.functions->thermal;
    instance.function({temperature, context}, output);
    if (!std::isfinite(output.conductivity.value()) || !(output.conductivity.value() > 0.0))
        throw std::domain_error("Thermal material function conductivity must be finite and positive");
    return output.conductivity;
}

double IsotropicThermoelasticMaterial::initial_density(double initial_temperature,
    MaterialFunctionContext context) const {
    if (!std::isfinite(initial_temperature))
        throw std::domain_error("Initial material temperature must be finite");
    context.time = 0.0;
    ThermalPropertyOutput output{};
    const ThermalFunctionInstance& instance = _properties.functions->thermal;
    instance.function({adlite::Scalar(initial_temperature), context}, output);
    if (!std::isfinite(output.density.value()) || !(output.density.value() > 0.0))
        throw std::domain_error("Initial material density must be finite and positive");
    return output.density.value();
}

adlite::Scalar IsotropicThermoelasticMaterial::current_density(double initial_temperature,
    const adlite::Scalar& volume_ratio,
    MaterialFunctionContext context) const {
    if (!std::isfinite(volume_ratio.value()) || !(volume_ratio.value() > 0.0))
        throw std::domain_error("Current material density requires a finite positive volume ratio");
    return initial_density(initial_temperature, context) / volume_ratio;
}

adlite::Scalar IsotropicThermoelasticMaterial::specific_heat(const adlite::Scalar& temperature,
    MaterialFunctionContext context) const {
    ThermalPropertyOutput output{};
    const ThermalFunctionInstance& instance = _properties.functions->thermal;
    instance.function({temperature, context}, output);
    if (!std::isfinite(output.specific_heat.value()) || !(output.specific_heat.value() > 0.0))
        throw std::domain_error("Material specific_heat must be finite and positive");
    return output.specific_heat;
}

adlite::Scalar IsotropicThermoelasticMaterial::reference_heat_capacity(const adlite::Scalar& temperature,
    double initial_temperature,
    MaterialFunctionContext context) const {
    return initial_density(initial_temperature, context) * specific_heat(temperature, context);
}

ActiveThermoelasticProperties IsotropicThermoelasticMaterial::active_properties(const adlite::Scalar& temperature,
    MaterialFunctionContext context) const {
    if (!std::isfinite(temperature.value()))
        throw std::domain_error("Thermoelastic material temperature must be finite");
    ElasticPropertyOutput output{};
    const ElasticFunctionInstance& instance = _properties.functions->elasticity;
    instance.function({temperature, context}, output);
    if (!std::isfinite(output.young_modulus.value()) || !(output.young_modulus.value() > 0.0))
        throw std::domain_error("Elasticity material function young_modulus must be finite and positive");
    if (!std::isfinite(output.poisson_ratio.value())
        || !(output.poisson_ratio.value() > -1.0 && output.poisson_ratio.value() < 0.5))
        throw std::domain_error("Elasticity material function poisson_ratio must lie between -1 and 0.5");
    const adlite::Scalar lame_lambda = output.young_modulus * output.poisson_ratio
                                       / ((1.0 + output.poisson_ratio) * (1.0 - 2.0 * output.poisson_ratio));
    const adlite::Scalar shear_modulus = output.young_modulus / (2.0 * (1.0 + output.poisson_ratio));
    return {lame_lambda, shear_modulus};
}

AxisymmetricStrain IsotropicThermoelasticMaterial::eigenstrain_rz(const adlite::Scalar& temperature,
    MaterialFunctionContext context) const {
    const SymmetricTensor3 value = eigenstrain(temperature, context);
    if (value.xy.value() != 0.0 || value.yz.value() != 0.0)
        throw std::domain_error("RZ material cannot represent xy or yz eigenstrain components");
    return {value.xx, value.zz, value.yy, value.xz};
}

SymmetricTensor3 IsotropicThermoelasticMaterial::eigenstrain(const adlite::Scalar& temperature,
    MaterialFunctionContext context) const {
    SymmetricTensor3 result{};
    for (const EigenstrainFunctionInstance& instance : _properties.functions->eigenstrains) {
        SymmetricTensor3 value{};
        instance.function({temperature, context}, value);
        if (!std::isfinite(value.xx.value()) || !std::isfinite(value.yy.value()) || !std::isfinite(value.zz.value())
            || !std::isfinite(value.xy.value()) || !std::isfinite(value.yz.value()) || !std::isfinite(value.xz.value()))
            throw std::domain_error("Eigenstrain material function output must be finite");
        result.xx += value.xx;
        result.yy += value.yy;
        result.zz += value.zz;
        result.xy += value.xy;
        result.yz += value.yz;
        result.xz += value.xz;
    }
    return result;
}

namespace {
AxisymmetricStress hooke_stress(const adlite::Scalar& lame_lambda,
    const adlite::Scalar& shear_modulus,
    const adlite::Scalar& strain_rr,
    const adlite::Scalar& strain_zz,
    const adlite::Scalar& strain_hoop,
    const adlite::Scalar& strain_rz) {
    const adlite::Scalar trace = strain_rr + strain_zz + strain_hoop;
    return {
        lame_lambda * trace + 2.0 * shear_modulus * strain_rr,
        lame_lambda * trace + 2.0 * shear_modulus * strain_zz,
        lame_lambda * trace + 2.0 * shear_modulus * strain_hoop,
        2.0 * shear_modulus * strain_rz,
    };
}
} // namespace

AxisymmetricStress IsotropicThermoelasticMaterial::stress(const adlite::Scalar& strain_rr,
    const adlite::Scalar& strain_zz,
    const adlite::Scalar& strain_hoop,
    const adlite::Scalar& strain_rz,
    const adlite::Scalar& temperature,
    MaterialFunctionContext context) const {
    const ActiveThermoelasticProperties active = active_properties(temperature, context);
    const AxisymmetricStrain imposed = eigenstrain_rz(temperature, context);
    return hooke_stress(active.lame_lambda,
        active.shear_modulus,
        strain_rr - imposed.rr,
        strain_zz - imposed.zz,
        strain_hoop - imposed.hoop,
        strain_rz - imposed.rz);
}

SymmetricTensor3 IsotropicThermoelasticMaterial::stress(const SymmetricTensor3& strain,
    const adlite::Scalar& temperature,
    MaterialFunctionContext context) const {
    const ActiveThermoelasticProperties active = active_properties(temperature, context);
    const SymmetricTensor3 imposed = eigenstrain(temperature, context);
    const adlite::Scalar strain_xx = strain.xx - imposed.xx, strain_yy = strain.yy - imposed.yy,
                         strain_zz = strain.zz - imposed.zz, trace = strain_xx + strain_yy + strain_zz;
    return {active.lame_lambda * trace + 2.0 * active.shear_modulus * strain_xx,
        active.lame_lambda * trace + 2.0 * active.shear_modulus * strain_yy,
        active.lame_lambda * trace + 2.0 * active.shear_modulus * strain_zz,
        2.0 * active.shear_modulus * (strain.xy - imposed.xy),
        2.0 * active.shear_modulus * (strain.yz - imposed.yz),
        2.0 * active.shear_modulus * (strain.xz - imposed.xz)};
}

SymmetricTensor3Values IsotropicThermoelasticMaterial::stress_values(const SymmetricTensor3Values& strain,
    double temperature,
    MaterialFunctionContext context) const {
    const ActiveThermoelasticProperties active = active_properties(adlite::Scalar(temperature), context);
    const SymmetricTensor3 imposed = eigenstrain(adlite::Scalar(temperature), context);
    const double lame_lambda = active.lame_lambda.value(), shear_modulus = active.shear_modulus.value();
    const double strain_xx = strain.xx - imposed.xx.value(), strain_yy = strain.yy - imposed.yy.value(),
                 strain_zz = strain.zz - imposed.zz.value(), trace = strain_xx + strain_yy + strain_zz;
    return {lame_lambda * trace + 2.0 * shear_modulus * strain_xx,
        lame_lambda * trace + 2.0 * shear_modulus * strain_yy,
        lame_lambda * trace + 2.0 * shear_modulus * strain_zz,
        2.0 * shear_modulus * (strain.xy - imposed.xy.value()),
        2.0 * shear_modulus * (strain.yz - imposed.yz.value()),
        2.0 * shear_modulus * (strain.xz - imposed.xz.value())};
}

namespace {
constexpr std::size_t component_count = 4;
constexpr int maximum_creep_iterations = 100;
const double log_denorm_min = std::log(std::numeric_limits<double>::denorm_min());

struct NortonCreepProperties final {
    double coefficient, reference_stress, stress_exponent;
};

struct J2PlasticityProperties final {
    double yield_stress, isotropic_hardening_modulus;
};

struct ActiveNortonCreepProperties final {
    adlite::Scalar coefficient, reference_stress, stress_exponent;
};

struct ActiveJ2PlasticityProperties final {
    adlite::Scalar yield_stress, isotropic_hardening_modulus;
};

struct NortonRoot final {
    double equivalent_stress, trial_stress_derivative, relaxed_fraction;
};

struct CoupledUpdate final {
    double equivalent_stress, stress_derivative, plastic_increment, plastic_increment_derivative, creep_increment,
        creep_increment_derivative;
};

void validate_material_point_state(const MaterialPointState& state) {
    for (std::size_t component = 0; component < component_count; ++component)
        if (!std::isfinite(state.elastic_strain[component]) || !std::isfinite(state.plastic_strain[component])
            || !std::isfinite(state.creep_strain[component]))
            throw std::domain_error("Material point strain histories must be finite");
    if (!std::isfinite(state.equivalent_plastic_strain) || !(state.equivalent_plastic_strain >= 0.0))
        throw std::domain_error("Material point equivalent plastic strain must be finite and nonnegative");
    if (!std::isfinite(state.equivalent_creep_strain) || !(state.equivalent_creep_strain >= 0.0))
        throw std::domain_error("Material point equivalent creep strain must be finite and nonnegative");
}

void validate_material_point_state(const CartesianMaterialPointState& state) {
    for (std::size_t component = 0; component < state.elastic_strain.size(); ++component)
        if (!std::isfinite(state.elastic_strain[component]) || !std::isfinite(state.plastic_strain[component])
            || !std::isfinite(state.creep_strain[component]))
            throw std::domain_error("Cartesian material point strain histories must be finite");
    if (!std::isfinite(state.equivalent_plastic_strain) || !(state.equivalent_plastic_strain >= 0.0)
        || !std::isfinite(state.equivalent_creep_strain) || !(state.equivalent_creep_strain >= 0.0))
        throw std::domain_error("Cartesian equivalent inelastic strains must be finite and nonnegative");
}

AxisymmetricStress tensor(const std::array<adlite::Scalar, component_count>& components) {
    return {components[0], components[1], components[2], components[3]};
}

std::array<adlite::Scalar, component_count> components(const AxisymmetricStress& tensor) {
    return {tensor.rr, tensor.zz, tensor.hoop, tensor.rz};
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

NortonRoot solve_power_law_equivalent_stress(double driving_stress,
    double log_stress_coefficient,
    double time_step,
    const NortonCreepProperties& creep) {
    if (time_step == 0.0 || creep.coefficient == 0.0
        || log_stress_coefficient == -std::numeric_limits<double>::infinity())
        return {driving_stress, 1.0, 0.0};
    if (!(driving_stress > 0.0)) {
        if (creep.stress_exponent > 1.0)
            return {driving_stress, 1.0, 0.0};
        const double log_linear_coefficient = log_stress_coefficient + std::log(time_step) + std::log(creep.coefficient)
                                              - std::log(creep.reference_stress);
        if (std::isnan(log_linear_coefficient))
            throw std::overflow_error("Norton creep linear coefficient is not a number");
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
    const double log_coefficient = log_stress_coefficient + std::log(time_step) + std::log(creep.coefficient)
                                   - log_driving_stress
                                   + creep.stress_exponent * (log_driving_stress - std::log(creep.reference_stress));
    if (std::isnan(log_coefficient))
        throw std::overflow_error("Norton creep dimensionless local coefficient is not a number");
    if (log_coefficient == -std::numeric_limits<double>::infinity())
        return {driving_stress, 1.0, 0.0};
    if (log_coefficient == std::numeric_limits<double>::infinity())
        throw std::overflow_error("Norton creep logarithmic local coefficient overflowed");
    double log_stress_ratio = 0.0;
    if (creep.stress_exponent == 1.0) {
        log_stress_ratio = log_add_exp(0.0, log_coefficient);
    } else {
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon(),
                     initial_residual = log_add_exp(0.0, log_coefficient);
        if (initial_residual == 0.0)
            return {driving_stress, 1.0, 0.0};
        constexpr double log_two = 0.693147180559945309417232121458176568;
        double lower = 0.0, upper = std::max(log_two, (log_coefficient + log_two) / creep.stress_exponent),
               current = log_coefficient > 0.0 ? log_coefficient / creep.stress_exponent : std::exp(log_coefficient);
        if (!std::isfinite(upper) || !(upper > lower))
            throw std::overflow_error("Norton creep logarithmic root bracket is invalid");
        if (!std::isfinite(current) || !(current > lower) || !(current < upper))
            current = 0.5 * (lower + upper);
        bool converged = false;
        for (int iteration = 0; iteration < maximum_creep_iterations; ++iteration) {
            const double first = -current, second = log_coefficient - creep.stress_exponent * current,
                         residual = log_add_exp(first, second);
            if (!std::isfinite(residual))
                throw std::overflow_error("Norton creep logarithmic residual is not finite");
            if (std::fabs(residual) <= tolerance) {
                converged = true;
                break;
            }
            if (residual > 0.0)
                lower = current;
            else
                upper = current;
            if (upper - lower <= tolerance * std::max(1.0, std::fabs(current))) {
                current = 0.5 * (lower + upper);
                converged = true;
                break;
            }
            const double weight = second_exponential_weight(first, second),
                         derivative = -(1.0 + (creep.stress_exponent - 1.0) * weight);
            double candidate = std::isfinite(derivative) ? current - residual / derivative : 0.5 * (lower + upper);
            if (!std::isfinite(candidate) || !(candidate > lower) || !(candidate < upper) || candidate == current)
                candidate = 0.5 * (lower + upper);
            current = candidate;
        }
        if (!converged)
            throw std::domain_error("Norton creep logarithmic local Newton solve did not converge");
        log_stress_ratio = current;
    }
    const double stress_ratio = std::exp(-log_stress_ratio), relaxed_fraction = -std::expm1(-log_stress_ratio),
                 log_equivalent_stress = log_driving_stress - log_stress_ratio,
                 equivalent_stress = log_equivalent_stress <= log_denorm_min ? 0.0 : std::exp(log_equivalent_stress),
                 derivative_denominator = stress_ratio + creep.stress_exponent * relaxed_fraction;
    const double driving_stress_derivative =
        derivative_denominator == 0.0 ? 0.0 : stress_ratio / derivative_denominator;
    if (!std::isfinite(equivalent_stress) || !std::isfinite(driving_stress_derivative)
        || !(driving_stress_derivative >= 0.0) || !std::isfinite(relaxed_fraction) || !(relaxed_fraction >= 0.0)
        || !(relaxed_fraction <= 1.0))
        throw std::overflow_error("Norton creep logarithmic root or derivative is not finite");
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
    return solve_power_law_equivalent_stress(trial_stress, std::log(3.0) + std::log(shear_modulus), time_step, creep);
}

double evaluate_creep_increment(double equivalent_stress, double time_step, const NortonCreepProperties& creep) {
    if (time_step == 0.0 || creep.coefficient == 0.0)
        return 0.0;
    if (!(equivalent_stress > 0.0))
        return 0.0;
    const double log_increment =
        std::log(time_step) + std::log(creep.coefficient)
        + creep.stress_exponent * (std::log(equivalent_stress) - std::log(creep.reference_stress));
    if (std::isnan(log_increment))
        throw std::overflow_error("Norton creep increment logarithm is not a number");
    const double maximum_log = std::log(std::numeric_limits<double>::max());
    if (log_increment <= log_denorm_min)
        return 0.0;
    if (!std::isfinite(log_increment) || log_increment >= maximum_log)
        throw std::overflow_error("Norton creep increment overflowed");
    const double increment = std::exp(log_increment);
    if (!std::isfinite(increment) || !(increment >= 0.0))
        throw std::overflow_error("Norton creep increment is not finite");
    return increment;
}

CoupledUpdate solve_coupled_update(double trial_stress,
    double shear_modulus,
    double time_step,
    const NortonCreepProperties& creep,
    const J2PlasticityProperties& plasticity,
    double committed_equivalent_plastic_strain) {
    const double current_yield_stress =
        plasticity.yield_stress + plasticity.isotropic_hardening_modulus * committed_equivalent_plastic_strain;
    if (!std::isfinite(current_yield_stress))
        throw std::overflow_error("Coupled plastic-creep current yield stress is not finite");
    const NortonRoot creep_only = solve_norton_equivalent_stress(trial_stress, shear_modulus, time_step, creep);
    const double inverse_three_shear_modulus = (1.0 / 3.0) / shear_modulus;
    if (!(creep_only.equivalent_stress > current_yield_stress)) {
        const double creep_increment = trial_stress * creep_only.relaxed_fraction * inverse_three_shear_modulus;
        const double creep_increment_derivative =
            (1.0 - creep_only.trial_stress_derivative) * inverse_three_shear_modulus;
        return {creep_only.equivalent_stress,
            creep_only.trial_stress_derivative,
            0.0,
            0.0,
            creep_increment,
            creep_increment_derivative};
    }
    const double hardening = plasticity.isotropic_hardening_modulus;
    if (hardening == 0.0) {
        const double creep_increment = evaluate_creep_increment(current_yield_stress, time_step, creep);
        const double plastic_increment =
            (trial_stress - current_yield_stress) * inverse_three_shear_modulus - creep_increment;
        if (!std::isfinite(plastic_increment) || !(plastic_increment > 0.0))
            throw std::domain_error("Coupled perfect-plastic update produced a nonpositive plastic increment");
        return {current_yield_stress, 0.0, plastic_increment, inverse_three_shear_modulus, creep_increment, 0.0};
    }
    double hardening_weight = 0.0, inverse_denominator = 0.0, log_root_stress_coefficient = 0.0;
    if (hardening <= shear_modulus) {
        const double ratio = hardening / shear_modulus, scaled_denominator = 3.0 + ratio;
        hardening_weight = ratio / scaled_denominator;
        inverse_denominator = (1.0 / scaled_denominator) / shear_modulus;
        log_root_stress_coefficient = std::log(hardening) + std::log(3.0 / scaled_denominator);
    } else {
        const double ratio = shear_modulus / hardening, scaled_denominator = 1.0 + 3.0 * ratio;
        hardening_weight = 1.0 / scaled_denominator;
        inverse_denominator = (1.0 / scaled_denominator) / hardening;
        log_root_stress_coefficient = std::log(shear_modulus) + std::log(3.0 / scaled_denominator);
    }
    const double root_driving_stress = current_yield_stress + hardening_weight * (trial_stress - current_yield_stress);
    if (!std::isfinite(inverse_denominator) || !(inverse_denominator > 0.0) || !std::isfinite(root_driving_stress)
        || !(root_driving_stress > 0.0) || !std::isfinite(log_root_stress_coefficient))
        throw std::overflow_error("Coupled plastic-creep transformed root is invalid");
    const NortonRoot root =
        solve_power_law_equivalent_stress(root_driving_stress, log_root_stress_coefficient, time_step, creep);
    const double creep_increment = evaluate_creep_increment(root.equivalent_stress, time_step, creep);
    const double plastic_increment =
        (trial_stress - current_yield_stress - 3.0 * (shear_modulus * creep_increment)) * inverse_denominator;
    const double plastic_increment_derivative = root.trial_stress_derivative * inverse_denominator,
                 stress_derivative = hardening_weight * root.trial_stress_derivative,
                 creep_increment_derivative = (1.0 - root.trial_stress_derivative) * inverse_three_shear_modulus;
    if (!std::isfinite(root.equivalent_stress) || !std::isfinite(stress_derivative) || !std::isfinite(plastic_increment)
        || !(plastic_increment > 0.0) || !std::isfinite(plastic_increment_derivative) || !std::isfinite(creep_increment)
        || !std::isfinite(creep_increment_derivative))
        throw std::overflow_error("Coupled plastic-creep update is not finite");
    return {root.equivalent_stress,
        stress_derivative,
        plastic_increment,
        plastic_increment_derivative,
        creep_increment,
        creep_increment_derivative};
}

double zero_stress_scale_shear_derivative(double scale, double time_step, const NortonCreepProperties& creep) {
    if (time_step == 0.0 || creep.coefficient == 0.0 || creep.stress_exponent > 1.0)
        return 0.0;
    return -3.0 * time_step * creep.coefficient / creep.reference_stress * scale * scale;
}

struct NortonDirectPartials final {
    double creep_stress_derivative, direct_coefficient, direct_reference, direct_exponent;
};

NortonDirectPartials
norton_direct_partials(double equivalent_stress, double creep_increment, const NortonCreepProperties& creep) {
    return {
        equivalent_stress > 0.0 ? creep.stress_exponent * creep_increment / equivalent_stress : 0.0,
        creep.coefficient > 0.0 ? creep_increment / creep.coefficient : 0.0,
        -creep.stress_exponent * creep_increment / creep.reference_stress,
        equivalent_stress > 0.0 ? creep_increment * std::log(equivalent_stress / creep.reference_stress) : 0.0,
    };
}

void fill_creep_partials(double* creep,
    std::size_t count,
    double creep_stress_derivative,
    const double* stress,
    const double* direct) {
    for (std::size_t parameter = 0; parameter < count; ++parameter)
        creep[parameter] = creep_stress_derivative * stress[parameter] + direct[parameter];
}

struct NortonPartials final {
    std::array<double, 5> stress{};
    std::array<double, 5> creep{};
};

NortonPartials norton_partials(double equivalent_trial_stress,
    double shear_modulus,
    double time_step,
    const NortonCreepProperties& creep,
    const NortonRoot& root,
    double creep_increment) {
    NortonPartials result;
    const double stress = root.equivalent_stress, stress_derivative = root.trial_stress_derivative;
    const NortonDirectPartials direct = norton_direct_partials(stress, creep_increment, creep);
    result.stress = {stress_derivative,
        -3.0 * creep_increment * stress_derivative,
        -3.0 * shear_modulus * direct.direct_coefficient * stress_derivative,
        -3.0 * shear_modulus * direct.direct_reference * stress_derivative,
        -3.0 * shear_modulus * direct.direct_exponent * stress_derivative};
    const std::array<double, 5> creep_direct = {0.0,
        0.0,
        direct.direct_coefficient,
        direct.direct_reference,
        direct.direct_exponent};
    fill_creep_partials(result.creep.data(),
        result.creep.size(),
        direct.creep_stress_derivative,
        result.stress.data(),
        creep_direct.data());
    if (time_step == 0.0 || equivalent_trial_stress == 0.0)
        result.creep.fill(0.0);
    return result;
}

struct CoupledPartials final {
    std::array<double, 7> stress{};
    std::array<double, 7> plastic{};
    std::array<double, 7> creep{};
};

CoupledPartials coupled_partials(double equivalent_trial_stress,
    double shear_modulus,
    const NortonCreepProperties& creep,
    const J2PlasticityProperties& plasticity,
    double committed_equivalent_plastic_strain,
    const CoupledUpdate& update) {
    CoupledPartials result;
    const double stress = update.equivalent_stress, plastic_increment = update.plastic_increment,
                 creep_increment = update.creep_increment;
    const NortonDirectPartials direct = norton_direct_partials(stress, creep_increment, creep);
    const std::array<double, 7> direct_creep =
        {0.0, 0.0, direct.direct_coefficient, direct.direct_reference, direct.direct_exponent, 0.0, 0.0};
    const double hardening = plasticity.isotropic_hardening_modulus;
    if (!(plastic_increment > 0.0)) {
        const NortonRoot root{stress, update.stress_derivative, 0.0};
        const NortonPartials norton =
            norton_partials(equivalent_trial_stress, shear_modulus, 1.0, creep, root, creep_increment);
        for (std::size_t parameter = 0; parameter < 5; ++parameter) {
            result.stress[parameter] = norton.stress[parameter];
            result.creep[parameter] = norton.creep[parameter];
        }
        return result;
    }
    if (hardening == 0.0) {
        result.stress[5] = 1.0;
        fill_creep_partials(result.creep.data(),
            result.creep.size(),
            direct.creep_stress_derivative,
            result.stress.data(),
            direct_creep.data());
        const double inverse_three_shear = 1.0 / (3.0 * shear_modulus);
        result.plastic[0] = inverse_three_shear;
        result.plastic[1] = -(equivalent_trial_stress - plasticity.yield_stress) / (3.0 * shear_modulus * shear_modulus)
                            - result.creep[1];
        for (std::size_t parameter = 2; parameter < result.plastic.size(); ++parameter)
            result.plastic[parameter] = -result.creep[parameter];
        result.plastic[5] -= inverse_three_shear;
        return result;
    }
    const double stress_derivative = update.stress_derivative;
    result.stress[0] = stress_derivative;
    result.stress[1] = -3.0 * (plastic_increment + creep_increment) * stress_derivative;
    result.stress[2] = -3.0 * shear_modulus * direct.direct_coefficient * stress_derivative;
    result.stress[3] = -3.0 * shear_modulus * direct.direct_reference * stress_derivative;
    result.stress[4] = -3.0 * shear_modulus * direct.direct_exponent * stress_derivative;
    result.stress[5] = 3.0 * shear_modulus / hardening * stress_derivative;
    result.stress[6] =
        3.0 * shear_modulus * (stress - plasticity.yield_stress) / (hardening * hardening) * stress_derivative;
    fill_creep_partials(result.creep.data(),
        result.creep.size(),
        direct.creep_stress_derivative,
        result.stress.data(),
        direct_creep.data());
    for (std::size_t parameter = 0; parameter < result.plastic.size(); ++parameter) {
        const double yield_derivative = parameter == 5 ? 1.0 : 0.0, hardening_derivative = parameter == 6 ? 1.0 : 0.0;
        result.plastic[parameter] =
            (result.stress[parameter] - yield_derivative
                - (committed_equivalent_plastic_strain + plastic_increment) * hardening_derivative)
            / hardening;
    }
    return result;
}

adlite::Scalar function_creep_rate(const CreepFunctionInstance& function,
    const adlite::Scalar& equivalent_stress,
    const adlite::Scalar& temperature,
    const adlite::Scalar& equivalent_creep_strain,
    MaterialFunctionContext context) {
    const adlite::Scalar rate = function.function({equivalent_stress, temperature, equivalent_creep_strain, context});
    if (std::isnan(rate.value()) || rate.value() < 0.0)
        throw std::domain_error("Creep material function rate must be nonnegative and must not be NaN");
    return rate;
}

adlite::Scalar function_flow_stress(const PlasticFunctionInstance& function,
    const adlite::Scalar& equivalent_plastic_strain,
    const adlite::Scalar& temperature,
    MaterialFunctionContext context) {
    const adlite::Scalar stress = function.function({equivalent_plastic_strain, temperature, context});
    if (!std::isfinite(stress.value()) || !(stress.value() > 0.0))
        throw std::domain_error("Plasticity material function flow stress must be finite and positive");
    return stress;
}

struct FunctionCreepRoot final {
    adlite::Scalar equivalent_stress, creep_increment;
    double zero_stress_scale;
};

double passive_creep_rate_and_derivative(const CreepFunctionInstance& function,
    double equivalent_stress,
    double temperature,
    double equivalent_creep_strain,
    MaterialFunctionContext context,
    double& derivative) {
    const adlite::Scalar active_stress = adlite::Scalar::independent(equivalent_stress, 0, 1);
    const adlite::Scalar rate =
        function_creep_rate(function, active_stress, temperature, equivalent_creep_strain, context);
    derivative = rate.is_active() ? rate.derivative(0) : 0.0;
    return rate.value();
}

struct BracketedRoot final {
    double value, derivative;
};

using ResidualFunction = adlite::Scalar (*)(double, void*);

BracketedRoot solve_bracketed_root(double upper,
    double initial,
    double residual_tolerance,
    double interval_tolerance,
    ResidualFunction function,
    void* context,
    const char* failure) {
    double lower = 0.0, current = initial;
    for (int iteration = 0; iteration < 256; ++iteration) {
        const adlite::Scalar residual = function(current, context);
        const double derivative = residual.derivative(0);
        if (std::isfinite(residual.value()) && std::isfinite(derivative) && derivative > 0.0
            && std::abs(residual.value()) <= residual_tolerance)
            return {current, derivative};
        if (!std::isfinite(residual.value()) || residual.value() > 0.0)
            upper = current;
        else
            lower = current;
        if (upper - lower <= interval_tolerance) {
            current = 0.5 * (lower + upper);
            const double final_derivative = function(current, context).derivative(0);
            if (std::isfinite(final_derivative) && final_derivative > 0.0)
                return {current, final_derivative};
            break;
        }
        double candidate = 0.5 * (lower + upper);
        if (std::isfinite(residual.value()) && std::isfinite(derivative) && derivative > 0.0) {
            const double newton = current - residual.value() / derivative;
            if (std::isfinite(newton) && newton > lower && newton < upper)
                candidate = newton;
        }
        current = candidate;
    }
    throw std::domain_error(failure);
}

struct FunctionCreepEquation final {
    const CreepFunctionInstance& function;
    double trial, shear, temperature, equivalent_creep, time_step;
    MaterialFunctionContext context;
};

adlite::Scalar function_creep_residual(double stress, void* raw) {
    auto& equation = *static_cast<FunctionCreepEquation*>(raw);
    const adlite::Scalar active = adlite::Scalar::independent(stress, 0, 1);
    const adlite::Scalar rate = function_creep_rate(equation.function,
        active,
        equation.temperature,
        equation.equivalent_creep,
        equation.context);
    if (rate.is_active() && rate.derivative(0) < 0.0)
        throw std::domain_error("Creep material function must be nondecreasing in equivalent stress");
    return active + 3.0 * equation.shear * equation.time_step * rate - equation.trial;
}

FunctionCreepRoot solve_function_creep(const CreepFunctionInstance& function,
    const adlite::Scalar& equivalent_trial_stress,
    const adlite::Scalar& shear_modulus,
    const adlite::Scalar& temperature,
    double equivalent_creep_strain,
    double time_step,
    MaterialFunctionContext context) {
    double zero_derivative = 0.0;
    const double zero_rate = passive_creep_rate_and_derivative(function,
        0.0,
        temperature.value(),
        equivalent_creep_strain,
        context,
        zero_derivative);
    if (zero_rate != 0.0)
        throw std::domain_error("Creep material function rate must be zero at zero equivalent stress");
    if (!std::isfinite(zero_derivative) || zero_derivative < 0.0)
        throw std::domain_error("Creep material function must be nondecreasing at zero equivalent stress");
    const double zero_scale = 1.0 / (1.0 + 3.0 * shear_modulus.value() * time_step * zero_derivative);
    if (!(equivalent_trial_stress.value() > 0.0))
        return {0.0, 0.0, zero_scale};
    const double trial = equivalent_trial_stress.value(), shear = shear_modulus.value();
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * trial;
    FunctionCreepEquation
        equation{function, trial, shear, temperature.value(), equivalent_creep_strain, time_step, context};
    const BracketedRoot root = solve_bracketed_root(trial,
        trial,
        tolerance,
        tolerance,
        function_creep_residual,
        &equation,
        "Registered creep material function local solve did not converge");
    const adlite::Scalar passive_stress(root.value);
    const adlite::Scalar active_rate =
        function_creep_rate(function, passive_stress, temperature, equivalent_creep_strain, context);
    const adlite::Scalar residual =
        passive_stress + 3.0 * shear_modulus * time_step * active_rate - equivalent_trial_stress;
    const adlite::Scalar returned_stress = passive_stress - residual / root.derivative;
    const adlite::Scalar increment =
        time_step * function_creep_rate(function, returned_stress, temperature, equivalent_creep_strain, context);
    return {returned_stress, increment, zero_scale};
}

struct FunctionCoupledUpdate final {
    adlite::Scalar equivalent_stress, plastic_increment, creep_increment;
};

struct FunctionCoupledEquation final {
    const MaterialFunctionSet& functions;
    double trial, shear, temperature, equivalent_plastic, equivalent_creep, time_step;
    MaterialFunctionContext context;
};

adlite::Scalar function_coupled_residual(double increment, void* raw) {
    auto& equation = *static_cast<FunctionCoupledEquation*>(raw);
    const adlite::Scalar active = adlite::Scalar::independent(increment, 0, 1);
    const adlite::Scalar flow = function_flow_stress(equation.functions.plasticity,
        equation.equivalent_plastic + active,
        equation.temperature,
        equation.context);
    if (flow.is_active() && flow.derivative(0) < 0.0)
        throw std::domain_error(
            "Plasticity material function flow stress must be nondecreasing in equivalent plastic strain");
    adlite::Scalar creep(0.0);
    if (equation.functions.has_creep())
        creep = equation.time_step
                * function_creep_rate(equation.functions.creep,
                    flow,
                    equation.temperature,
                    equation.equivalent_creep,
                    equation.context);
    return flow + 3.0 * equation.shear * (active + creep) - equation.trial;
}

FunctionCoupledUpdate solve_function_coupled(const MaterialFunctionSet& functions,
    const adlite::Scalar& equivalent_trial_stress,
    const adlite::Scalar& shear_modulus,
    const adlite::Scalar& temperature,
    double equivalent_plastic_strain,
    double equivalent_creep_strain,
    double time_step,
    MaterialFunctionContext context) {
    const adlite::Scalar flow_at_committed =
        function_flow_stress(functions.plasticity, equivalent_plastic_strain, temperature, context);
    FunctionCreepRoot creep_only{equivalent_trial_stress, 0.0, 1.0};
    if (functions.has_creep())
        creep_only = solve_function_creep(functions.creep,
            equivalent_trial_stress,
            shear_modulus,
            temperature,
            equivalent_creep_strain,
            time_step,
            context);
    if (!(creep_only.equivalent_stress.value() > flow_at_committed.value()))
        return {creep_only.equivalent_stress, 0.0, creep_only.creep_increment};
    const double trial = equivalent_trial_stress.value(), shear = shear_modulus.value(), upper = trial / (3.0 * shear);
    const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * trial;
    FunctionCoupledEquation equation{functions,
        trial,
        shear,
        temperature.value(),
        equivalent_plastic_strain,
        equivalent_creep_strain,
        time_step,
        context};
    const BracketedRoot root = solve_bracketed_root(upper,
        0.5 * upper,
        tolerance,
        tolerance / (3.0 * shear),
        function_coupled_residual,
        &equation,
        "Registered plasticity material function local solve did not converge");
    const adlite::Scalar passive_increment(root.value);
    const adlite::Scalar active_flow =
        function_flow_stress(functions.plasticity, equivalent_plastic_strain + passive_increment, temperature, context);
    adlite::Scalar active_creep(0.0);
    if (functions.has_creep())
        active_creep =
            time_step
            * function_creep_rate(functions.creep, active_flow, temperature, equivalent_creep_strain, context);
    const adlite::Scalar residual =
        active_flow + 3.0 * shear_modulus * (passive_increment + active_creep) - equivalent_trial_stress;
    const adlite::Scalar plastic_increment = passive_increment - residual / root.derivative;
    const adlite::Scalar returned_stress =
        function_flow_stress(functions.plasticity, equivalent_plastic_strain + plastic_increment, temperature, context);
    adlite::Scalar creep_increment(0.0);
    if (functions.has_creep())
        creep_increment =
            time_step
            * function_creep_rate(functions.creep, returned_stress, temperature, equivalent_creep_strain, context);
    return {returned_stress, plastic_increment, creep_increment};
}
} // namespace

namespace {
ActiveNortonCreepProperties active_creep_properties(const CreepFunctionInstance& function,
    const adlite::Scalar& temperature) {
    const CreepBuiltinParameters& parameters = function.builtin;
    adlite::Scalar coefficient = parameters.coefficient;
    adlite::Scalar reference_stress = parameters.reference_stress;
    adlite::Scalar exponent = parameters.stress_exponent;
    if (parameters.kind == CreepBuiltinParameters::Kind::linear_temperature_norton) {
        const adlite::Scalar change = temperature - parameters.reference_temperature;
        coefficient += parameters.coefficient_temperature_coefficient * change;
        reference_stress += parameters.reference_stress_temperature_coefficient * change;
        exponent += parameters.stress_exponent_temperature_coefficient * change;
    }
    ActiveNortonCreepProperties active{coefficient, reference_stress, exponent};
    if (!std::isfinite(active.coefficient.value()) || !(active.coefficient.value() >= 0.0)
        || !std::isfinite(active.reference_stress.value()) || !(active.reference_stress.value() > 0.0)
        || !std::isfinite(active.stress_exponent.value()) || !(active.stress_exponent.value() >= 1.0))
        throw std::domain_error("Active Norton creep properties violate their physical domain");
    return active;
}

ActiveJ2PlasticityProperties active_plasticity_properties(const PlasticFunctionInstance& function,
    const adlite::Scalar& temperature) {
    const PlasticBuiltinParameters& parameters = function.builtin;
    adlite::Scalar yield_stress = parameters.yield_stress;
    adlite::Scalar hardening = parameters.hardening_modulus;
    if (parameters.kind == PlasticBuiltinParameters::Kind::linear_temperature_isotropic_hardening) {
        const adlite::Scalar change = temperature - parameters.reference_temperature;
        yield_stress += parameters.yield_stress_temperature_coefficient * change;
        hardening += parameters.hardening_temperature_coefficient * change;
    }
    ActiveJ2PlasticityProperties active{yield_stress, hardening};
    if (!std::isfinite(active.yield_stress.value()) || !(active.yield_stress.value() > 0.0)
        || !std::isfinite(active.isotropic_hardening_modulus.value())
        || !(active.isotropic_hardening_modulus.value() >= 0.0))
        throw std::domain_error("Active J2 plasticity properties violate their physical domain");
    return active;
}

CartesianMaterialPointState evaluate_inelastic_tensor_values(const MaterialFunctionSet& functions,
    double lame_lambda,
    double shear_modulus,
    const std::array<double, 6>& total,
    const std::array<double, 6>& imposed,
    const CartesianMaterialPointState& committed,
    double temperature,
    double time_step) {
    CartesianMaterialPointState result = committed;
    std::array<double, 6> elastic{}, trial_stress{}, deviatoric{};
    for (std::size_t component = 0; component < 6; ++component) {
        if (!std::isfinite(total[component]))
            throw std::domain_error("Inelastic material strain is not finite");
        elastic[component] = total[component] - imposed[component] - committed.plastic_strain[component]
                             - committed.creep_strain[component];
    }
    const double trace = elastic[0] + elastic[1] + elastic[2];
    for (std::size_t component = 0; component < 3; ++component)
        trial_stress[component] = lame_lambda * trace + 2.0 * shear_modulus * elastic[component];
    for (std::size_t component = 3; component < 6; ++component)
        trial_stress[component] = 2.0 * shear_modulus * elastic[component];
    const double mean = (trial_stress[0] + trial_stress[1] + trial_stress[2]) / 3.0;
    for (std::size_t component = 0; component < 6; ++component)
        deviatoric[component] = trial_stress[component] - (component < 3 ? mean : 0.0);
    double normal_norm = 0.0, shear_norm = 0.0;
    for (std::size_t component = 0; component < 3; ++component)
        normal_norm = std::hypot(normal_norm, deviatoric[component]);
    for (std::size_t component = 3; component < 6; ++component)
        shear_norm = std::hypot(shear_norm, deviatoric[component]);
    constexpr double square_root_two = 1.414213562373095048801688724209698079;
    constexpr double square_root_three_halves = 1.224744871391589049098642037352945695;
    const double equivalent_trial_stress =
        square_root_three_halves * std::hypot(normal_norm, square_root_two * shear_norm);
    if (!std::isfinite(equivalent_trial_stress))
        throw std::overflow_error("Inelastic material trial equivalent stress is not finite");

    double equivalent_stress = equivalent_trial_stress, plastic_increment = 0.0, creep_increment = 0.0,
           zero_stress_scale = 1.0;
    if (functions.has_creep()) {
        const ActiveNortonCreepProperties active_creep =
            active_creep_properties(functions.creep, adlite::Scalar(temperature));
        const NortonCreepProperties creep = {active_creep.coefficient.value(),
            active_creep.reference_stress.value(),
            active_creep.stress_exponent.value()};
        if (functions.has_plasticity()) {
            const ActiveJ2PlasticityProperties active_plasticity =
                active_plasticity_properties(functions.plasticity, adlite::Scalar(temperature));
            const J2PlasticityProperties plasticity = {active_plasticity.yield_stress.value(),
                active_plasticity.isotropic_hardening_modulus.value()};
            const CoupledUpdate update = solve_coupled_update(equivalent_trial_stress,
                shear_modulus,
                time_step,
                creep,
                plasticity,
                committed.equivalent_plastic_strain);
            equivalent_stress = update.equivalent_stress;
            plastic_increment = update.plastic_increment;
            creep_increment = update.creep_increment;
        } else if (time_step != 0.0 && creep.coefficient != 0.0) {
            const NortonRoot root =
                solve_norton_equivalent_stress(equivalent_trial_stress, shear_modulus, time_step, creep);
            equivalent_stress = root.equivalent_stress;
            creep_increment = evaluate_creep_increment(root.equivalent_stress, time_step, creep);
            zero_stress_scale = root.trial_stress_derivative;
        }
    } else if (functions.has_plasticity() && equivalent_trial_stress != 0.0) {
        const ActiveJ2PlasticityProperties active_plasticity =
            active_plasticity_properties(functions.plasticity, adlite::Scalar(temperature));
        const double hardening = active_plasticity.isotropic_hardening_modulus.value();
        const double current_yield =
            active_plasticity.yield_stress.value() + hardening * committed.equivalent_plastic_strain;
        if (!std::isfinite(current_yield))
            throw std::overflow_error("J2 plasticity current yield stress is not finite");
        if (equivalent_trial_stress > current_yield) {
            plastic_increment = (equivalent_trial_stress - current_yield) / (3.0 * shear_modulus + hardening);
            equivalent_stress = equivalent_trial_stress - 3.0 * shear_modulus * plastic_increment;
        }
    }
    const bool zero = equivalent_trial_stress == 0.0;
    const double scale = zero ? zero_stress_scale : equivalent_stress / equivalent_trial_stress;
    for (std::size_t component = 0; component < 6; ++component) {
        if (zero)
            result.creep_strain[component] += (1.0 - scale) * deviatoric[component] / (2.0 * shear_modulus);
        else {
            const double direction = 1.5 * deviatoric[component] / equivalent_trial_stress;
            result.plastic_strain[component] += plastic_increment * direction;
            result.creep_strain[component] += creep_increment * direction;
        }
        result.elastic_strain[component] = elastic[component]
                                           - (result.plastic_strain[component] - committed.plastic_strain[component])
                                           - (result.creep_strain[component] - committed.creep_strain[component]);
    }
    result.stress = {mean + scale * deviatoric[0],
        mean + scale * deviatoric[1],
        mean + scale * deviatoric[2],
        scale * deviatoric[3],
        scale * deviatoric[4],
        scale * deviatoric[5]};
    result.equivalent_plastic_strain += plastic_increment;
    result.equivalent_creep_strain += creep_increment;
    validate_material_point_state(result);
    return result;
}

bool uses_builtin_inelastic_update(const MaterialFunctionSet& functions) {
    const bool builtin_creep =
        !functions.has_creep() || functions.creep.builtin.kind != CreepBuiltinParameters::Kind::custom;
    const bool builtin_plasticity =
        !functions.has_plasticity() || functions.plasticity.builtin.kind != PlasticBuiltinParameters::Kind::custom;
    const bool matching_reference =
        functions.creep.builtin.kind != CreepBuiltinParameters::Kind::linear_temperature_norton
        || functions.plasticity.builtin.kind != PlasticBuiltinParameters::Kind::linear_temperature_isotropic_hardening
        || functions.creep.builtin.reference_temperature == functions.plasticity.builtin.reference_temperature;
    return builtin_creep && builtin_plasticity && matching_reference;
}

struct InelasticScalarUpdate final {
    adlite::Scalar equivalent_stress, plastic_increment, creep_increment, zero_stress_scale{1.0};
};

InelasticScalarUpdate evaluate_inelastic_scalars(const MaterialFunctionSet& functions,
    const adlite::Scalar& equivalent_trial_stress,
    const adlite::Scalar& shear_modulus,
    const adlite::Scalar& temperature,
    double equivalent_plastic_strain,
    double equivalent_creep_strain,
    double time_step,
    MaterialFunctionContext context) {
    if (!uses_builtin_inelastic_update(functions)) {
        if (!functions.has_creep() && !functions.has_plasticity())
            return {equivalent_trial_stress, 0.0, 0.0, 1.0};
        if (equivalent_trial_stress.value() == 0.0 && !functions.has_creep())
            return {equivalent_trial_stress, 0.0, 0.0, 1.0};
        if (equivalent_trial_stress.value() != 0.0 && functions.has_plasticity()) {
            const FunctionCoupledUpdate update = solve_function_coupled(functions,
                equivalent_trial_stress,
                shear_modulus,
                temperature,
                equivalent_plastic_strain,
                equivalent_creep_strain,
                time_step,
                context);
            return {update.equivalent_stress, update.plastic_increment, update.creep_increment, 1.0};
        }
        const FunctionCreepRoot update = solve_function_creep(functions.creep,
            equivalent_trial_stress,
            shear_modulus,
            temperature,
            equivalent_creep_strain,
            time_step,
            context);
        return {update.equivalent_stress, 0.0, update.creep_increment, update.zero_stress_scale};
    }
    if (!functions.has_creep() && !functions.has_plasticity())
        return {equivalent_trial_stress, 0.0, 0.0, 1.0};
    if (!functions.has_creep()) {
        if (equivalent_trial_stress.value() == 0.0)
            return {equivalent_trial_stress, 0.0, 0.0, 1.0};
        const ActiveJ2PlasticityProperties plasticity = active_plasticity_properties(functions.plasticity, temperature);
        const adlite::Scalar hardening = plasticity.isotropic_hardening_modulus;
        const adlite::Scalar current_yield = plasticity.yield_stress + hardening * equivalent_plastic_strain;
        if (!std::isfinite(current_yield.value()))
            throw std::overflow_error("J2 plasticity current yield stress is not finite");
        if (!(equivalent_trial_stress.value() > current_yield.value()))
            return {equivalent_trial_stress, 0.0, 0.0, 1.0};
        const adlite::Scalar increment = (equivalent_trial_stress - current_yield) / (3.0 * shear_modulus + hardening);
        return {equivalent_trial_stress - 3.0 * shear_modulus * increment, increment, 0.0, 1.0};
    }
    const ActiveNortonCreepProperties creep = active_creep_properties(functions.creep, temperature);
    const NortonCreepProperties creep_values = {creep.coefficient.value(),
        creep.reference_stress.value(),
        creep.stress_exponent.value()};
    if (equivalent_trial_stress.value() == 0.0) {
        const NortonRoot root = solve_norton_equivalent_stress(0.0, shear_modulus.value(), time_step, creep_values);
        return {0.0,
            0.0,
            0.0,
            adlite::compose(root.trial_stress_derivative,
                shear_modulus,
                zero_stress_scale_shear_derivative(root.trial_stress_derivative, time_step, creep_values))};
    }
    if (functions.has_plasticity()) {
        const ActiveJ2PlasticityProperties plasticity = active_plasticity_properties(functions.plasticity, temperature);
        const J2PlasticityProperties plasticity_values = {plasticity.yield_stress.value(),
            plasticity.isotropic_hardening_modulus.value()};
        const CoupledUpdate update = solve_coupled_update(equivalent_trial_stress.value(),
            shear_modulus.value(),
            time_step,
            creep_values,
            plasticity_values,
            equivalent_plastic_strain);
        const CoupledPartials partials = coupled_partials(equivalent_trial_stress.value(),
            shear_modulus.value(),
            creep_values,
            plasticity_values,
            equivalent_plastic_strain,
            update);
        const std::array<adlite::Scalar, 7> inputs = {equivalent_trial_stress,
            shear_modulus,
            creep.coefficient,
            creep.reference_stress,
            creep.stress_exponent,
            plasticity.yield_stress,
            plasticity.isotropic_hardening_modulus};
        return {adlite::compose(update.equivalent_stress, inputs.data(), partials.stress.data(), inputs.size()),
            adlite::compose(update.plastic_increment, inputs.data(), partials.plastic.data(), inputs.size()),
            adlite::compose(update.creep_increment, inputs.data(), partials.creep.data(), inputs.size()),
            1.0};
    }
    if (time_step == 0.0 || creep.coefficient.value() == 0.0)
        return {equivalent_trial_stress, 0.0, 0.0, 1.0};
    const NortonRoot root =
        solve_norton_equivalent_stress(equivalent_trial_stress.value(), shear_modulus.value(), time_step, creep_values);
    const double increment = evaluate_creep_increment(root.equivalent_stress, time_step, creep_values);
    const NortonPartials partials = norton_partials(equivalent_trial_stress.value(),
        shear_modulus.value(),
        time_step,
        creep_values,
        root,
        increment);
    const std::array<adlite::Scalar, 5> inputs = {equivalent_trial_stress,
        shear_modulus,
        creep.coefficient,
        creep.reference_stress,
        creep.stress_exponent};
    return {adlite::compose(root.equivalent_stress, inputs.data(), partials.stress.data(), inputs.size()),
        0.0,
        adlite::compose(increment, inputs.data(), partials.creep.data(), inputs.size()),
        1.0};
}

struct InelasticTensorUpdate final {
    std::array<adlite::Scalar, 6> stress{}, elastic{}, plastic{}, creep{};
    adlite::Scalar equivalent_plastic{0.0}, equivalent_creep{0.0};
};

InelasticTensorUpdate evaluate_inelastic_tensor(const MaterialFunctionSet& functions,
    const ActiveThermoelasticProperties& active,
    const std::array<adlite::Scalar, 6>& total,
    const std::array<adlite::Scalar, 6>& imposed,
    const double* committed_plastic,
    const double* committed_creep,
    std::size_t count,
    double equivalent_plastic,
    double equivalent_creep,
    const adlite::Scalar& temperature,
    double time_step,
    MaterialFunctionContext context) {
    if (count < 4 || count > 6 || !std::isfinite(temperature.value()) || !std::isfinite(time_step) || time_step < 0.0)
        throw std::domain_error("Inelastic material tensor layout, temperature, or time step is invalid");
    InelasticTensorUpdate result;
    for (std::size_t component = 0; component < count; ++component) {
        if (!std::isfinite(total[component].value()))
            throw std::domain_error("Inelastic material strain is not finite");
        result.plastic[component] = committed_plastic[component];
        result.creep[component] = committed_creep[component];
        result.elastic[component] =
            total[component] - imposed[component] - committed_plastic[component] - committed_creep[component];
    }
    const adlite::Scalar trace = result.elastic[0] + result.elastic[1] + result.elastic[2];
    for (std::size_t component = 0; component < 3; ++component)
        result.stress[component] = active.lame_lambda * trace + 2.0 * active.shear_modulus * result.elastic[component];
    for (std::size_t component = 3; component < count; ++component)
        result.stress[component] = 2.0 * active.shear_modulus * result.elastic[component];
    const adlite::Scalar mean = (result.stress[0] + result.stress[1] + result.stress[2]) / 3.0;
    std::array<adlite::Scalar, 6> deviatoric = result.stress;
    bool zero = true;
    for (std::size_t component = 0; component < count; ++component) {
        if (component < 3)
            deviatoric[component] -= mean;
        zero = zero && deviatoric[component].value() == 0.0;
    }
    adlite::Scalar normal_norm(0.0), shear_norm(0.0);
    if (!zero) {
        for (std::size_t component = 0; component < 3; ++component)
            normal_norm = adlite::hypot(normal_norm, deviatoric[component]);
        for (std::size_t component = 3; component < count; ++component)
            shear_norm = adlite::hypot(shear_norm, deviatoric[component]);
    }
    constexpr double square_root_two = 1.414213562373095048801688724209698079;
    constexpr double square_root_three_halves = 1.224744871391589049098642037352945695;
    const adlite::Scalar trial_stress =
        zero ? adlite::Scalar(0.0)
             : square_root_three_halves
                   * (count == 4 ? adlite::hypot(adlite::hypot(deviatoric[0], deviatoric[1]),
                                       adlite::hypot(deviatoric[2], square_root_two * deviatoric[3]))
                                 : adlite::hypot(normal_norm, square_root_two * shear_norm));
    if (!std::isfinite(trial_stress.value()))
        throw std::overflow_error("Inelastic material trial equivalent stress is not finite");
    const InelasticScalarUpdate update = evaluate_inelastic_scalars(functions,
        trial_stress,
        active.shear_modulus,
        temperature,
        equivalent_plastic,
        equivalent_creep,
        time_step,
        context);
    const adlite::Scalar scale = zero ? update.zero_stress_scale : update.equivalent_stress / trial_stress;
    for (std::size_t component = 0; component < count; ++component) {
        if (zero)
            result.creep[component] += (1.0 - scale) * deviatoric[component] / (2.0 * active.shear_modulus);
        else {
            const adlite::Scalar direction = 1.5 * deviatoric[component] / trial_stress;
            result.plastic[component] += update.plastic_increment * direction;
            result.creep[component] += update.creep_increment * direction;
        }
        result.elastic[component] = result.elastic[component]
                                    - (result.plastic[component] - committed_plastic[component])
                                    - (result.creep[component] - committed_creep[component]);
        result.stress[component] = component < 3 ? mean + scale * deviatoric[component] : scale * deviatoric[component];
    }
    result.equivalent_plastic = equivalent_plastic + update.plastic_increment;
    result.equivalent_creep = equivalent_creep + update.creep_increment;
    return result;
}
} // namespace

double IsotropicThermoelasticMaterial::equivalent_creep_rate(const MaterialPointState& state,
    double temperature,
    MaterialFunctionContext context) const {
    CartesianMaterialPointState point;
    point.stress = {state.stress.rr, state.stress.zz, state.stress.hoop, state.stress.rz, 0.0, 0.0};
    point.equivalent_creep_strain = state.equivalent_creep_strain;
    return equivalent_creep_rate(point, temperature, context);
}

double IsotropicThermoelasticMaterial::equivalent_creep_rate(const CartesianMaterialPointState& state,
    double temperature,
    MaterialFunctionContext context) const {
    if (!functions().has_creep())
        return 0.0;
    const auto& s = state.stress;
    const double mean = s.xx / 3.0 + s.yy / 3.0 + s.zz / 3.0;
    const adlite::Scalar normal = adlite::hypot(adlite::hypot(adlite::Scalar(s.xx - mean), adlite::Scalar(s.yy - mean)),
        adlite::Scalar(s.zz - mean));
    const adlite::Scalar shear =
        adlite::hypot(adlite::hypot(adlite::Scalar(s.xy), adlite::Scalar(s.yz)), adlite::Scalar(s.xz));
    const double equivalent = (std::sqrt(1.5) * adlite::hypot(normal, std::sqrt(2.0) * shear)).value();
    if (!std::isfinite(equivalent) || !std::isfinite(temperature) || !std::isfinite(state.equivalent_creep_strain)
        || state.equivalent_creep_strain < 0.0)
        throw std::domain_error("Creep rate requires finite stress, temperature and nonnegative creep history");
    double rate;
    if (functions().creep.builtin.kind != CreepBuiltinParameters::Kind::custom) {
        const auto active = active_creep_properties(functions().creep, adlite::Scalar(temperature));
        rate = evaluate_creep_increment(equivalent,
            1.0,
            {active.coefficient.value(), active.reference_stress.value(), active.stress_exponent.value()});
    } else {
        rate = function_creep_rate(functions().creep,
            adlite::Scalar(equivalent),
            adlite::Scalar(temperature),
            adlite::Scalar(state.equivalent_creep_strain),
            context)
                   .value();
    }
    if (!std::isfinite(rate) || rate < 0.0)
        throw std::domain_error("Creep rate must be finite and nonnegative");
    return rate;
}

InelasticStressResponse IsotropicThermoelasticMaterial::response(const adlite::Scalar& strain_rr,
    const adlite::Scalar& strain_zz,
    const adlite::Scalar& strain_hoop,
    const adlite::Scalar& strain_rz,
    const adlite::Scalar& temperature,
    double time_step,
    const MaterialPointState& committed,
    MaterialFunctionContext context) const {
    validate_material_point_state(committed);
    const AxisymmetricStrain imposed = eigenstrain_rz(temperature, context);
    const auto update = evaluate_inelastic_tensor(functions(),
        active_properties(temperature, context),
        {strain_rr, strain_zz, strain_hoop, strain_rz, 0.0, 0.0},
        {imposed.rr, imposed.zz, imposed.hoop, imposed.rz, 0.0, 0.0},
        committed.plastic_strain.data(),
        committed.creep_strain.data(),
        4,
        committed.equivalent_plastic_strain,
        committed.equivalent_creep_strain,
        temperature,
        time_step,
        context);
    MaterialPointTrialState trial;
    std::copy_n(update.elastic.begin(), 4, trial.elastic_strain.begin());
    std::copy_n(update.plastic.begin(), 4, trial.plastic_strain.begin());
    std::copy_n(update.creep.begin(), 4, trial.creep_strain.begin());
    trial.equivalent_plastic_strain = update.equivalent_plastic;
    trial.equivalent_creep_strain = update.equivalent_creep;
    return {{update.stress[0], update.stress[1], update.stress[2], update.stress[3]}, trial};
}

CartesianInelasticStressResponse IsotropicThermoelasticMaterial::response(const SymmetricTensor3& strain,
    const adlite::Scalar& temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context) const {
    validate_material_point_state(committed);
    const SymmetricTensor3 imposed = eigenstrain(temperature, context);
    const auto update = evaluate_inelastic_tensor(functions(),
        active_properties(temperature, context),
        {strain.xx, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz},
        {imposed.xx, imposed.yy, imposed.zz, imposed.xy, imposed.yz, imposed.xz},
        committed.plastic_strain.data(),
        committed.creep_strain.data(),
        6,
        committed.equivalent_plastic_strain,
        committed.equivalent_creep_strain,
        temperature,
        time_step,
        context);
    CartesianMaterialPointState trial;
    for (std::size_t component = 0; component < 6; ++component) {
        trial.elastic_strain[component] = update.elastic[component].value();
        trial.plastic_strain[component] = update.plastic[component].value();
        trial.creep_strain[component] = update.creep[component].value();
    }
    trial.equivalent_plastic_strain = update.equivalent_plastic.value();
    trial.equivalent_creep_strain = update.equivalent_creep.value();
    trial.stress = {update.stress[0].value(),
        update.stress[1].value(),
        update.stress[2].value(),
        update.stress[3].value(),
        update.stress[4].value(),
        update.stress[5].value()};
    validate_material_point_state(trial);
    return {
        {update.stress[0], update.stress[1], update.stress[2], update.stress[3], update.stress[4], update.stress[5]},
        trial};
}

CartesianMaterialPointState IsotropicThermoelasticMaterial::response_values(const SymmetricTensor3Values& strain,
    double temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context) const {
    validate_material_point_state(committed);
    if (!std::isfinite(temperature) || !std::isfinite(time_step) || time_step < 0.0)
        throw std::domain_error("Inelastic material temperature or time step is invalid");
    const MaterialFunctionSet& material_functions = functions();
    if (!uses_builtin_inelastic_update(material_functions)) {
        const CartesianInelasticStressResponse response =
            this->response({strain.xx, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz},
                adlite::Scalar(temperature),
                time_step,
                committed,
                context);
        return response.trial_state;
    }
    const ActiveThermoelasticProperties active = active_properties(adlite::Scalar(temperature), context);
    const SymmetricTensor3 imposed = eigenstrain(adlite::Scalar(temperature), context);
    return evaluate_inelastic_tensor_values(material_functions,
        active.lame_lambda.value(),
        active.shear_modulus.value(),
        {strain.xx, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz},
        {imposed.xx.value(),
            imposed.yy.value(),
            imposed.zz.value(),
            imposed.xy.value(),
            imposed.yz.value(),
            imposed.xz.value()},
        committed,
        temperature,
        time_step);
}

InelasticStressResponse IsotropicThermoelasticMaterial::incremental_response(const adlite::Scalar& strain_increment_rr,
    const adlite::Scalar& strain_increment_zz,
    const adlite::Scalar& strain_increment_hoop,
    const adlite::Scalar& strain_increment_rz,
    const AxisymmetricRotation& rotation,
    const adlite::Scalar& temperature,
    double committed_temperature,
    double time_step,
    const MaterialPointState& committed,
    MaterialFunctionContext context) const {
    if (!std::isfinite(committed_temperature) || !(committed_temperature > 0.0))
        throw std::domain_error("Incremental material committed temperature must be finite and positive");
    const adlite::Scalar old_temperature(committed_temperature);
    MaterialFunctionContext old_context = context;
    old_context.time -= time_step;
    const AxisymmetricStrain old_imposed = eigenstrain_rz(old_temperature, old_context);
    const std::array<adlite::Scalar, component_count> synthetic_total = {
        committed.elastic_strain[0] + strain_increment_rr + old_imposed.rr + committed.plastic_strain[0]
            + committed.creep_strain[0],
        committed.elastic_strain[1] + strain_increment_zz + old_imposed.zz + committed.plastic_strain[1]
            + committed.creep_strain[1],
        committed.elastic_strain[2] + strain_increment_hoop + old_imposed.hoop + committed.plastic_strain[2]
            + committed.creep_strain[2],
        committed.elastic_strain[3] + strain_increment_rz + old_imposed.rz + committed.plastic_strain[3]
            + committed.creep_strain[3],
    };
    InelasticStressResponse result = response(synthetic_total[0],
        synthetic_total[1],
        synthetic_total[2],
        synthetic_total[3],
        temperature,
        time_step,
        committed,
        context);
    result.stress = rotate_axisymmetric_tensor(result.stress, rotation);
    result.trial_state.elastic_strain =
        components(rotate_axisymmetric_tensor(tensor(result.trial_state.elastic_strain), rotation));
    result.trial_state.plastic_strain =
        components(rotate_axisymmetric_tensor(tensor(result.trial_state.plastic_strain), rotation));
    result.trial_state.creep_strain =
        components(rotate_axisymmetric_tensor(tensor(result.trial_state.creep_strain), rotation));
    return result;
}

MaterialPointState IsotropicThermoelasticMaterial::state_values(const MaterialPointTrialState& trial_state) {
    MaterialPointState state;
    for (std::size_t component = 0; component < component_count; ++component) {
        state.elastic_strain[component] = trial_state.elastic_strain[component].value();
        state.plastic_strain[component] = trial_state.plastic_strain[component].value();
        state.creep_strain[component] = trial_state.creep_strain[component].value();
    }
    state.equivalent_plastic_strain = trial_state.equivalent_plastic_strain.value();
    state.equivalent_creep_strain = trial_state.equivalent_creep_strain.value();
    validate_material_point_state(state);
    return state;
}

// Evaluates the constitutive relation with AD seeded only on the four strain components and
// the temperature (width 5), returning the stress values, the consistent material tangent
// d(stress)/d(strain), the thermal coupling d(stress)/dT, and optional trial history.
AxisymmetricMaterialResponse evaluate_axisymmetric_material_response(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 4>& fed_strain,
    double temperature,
    double time_step,
    const MaterialPointState* committed_material,
    MaterialFunctionContext context,
    bool compute_tangent,
    bool compute_history) {
    const std::array<double, 5> seeds = {fed_strain[0], fed_strain[1], fed_strain[2], fed_strain[3], temperature};
    std::array<adlite::Scalar, 5> active{};
    if (compute_tangent)
        adlite::seed_identity(seeds.data(), seeds.size(), active.data());
    else
        for (std::size_t i = 0; i < 5; ++i)
            active[i] = seeds[i];
    AxisymmetricMaterialResponse result{};
    AxisymmetricStress stress;
    if (committed_material) {
        const auto response = material.response(active[0],
            active[1],
            active[2],
            active[3],
            active[4],
            time_step,
            *committed_material,
            context);
        stress = response.stress;
        if (compute_history)
            result.history = IsotropicThermoelasticMaterial::state_values(response.trial_state);
    } else
        stress = material.stress(active[0], active[1], active[2], active[3], active[4], context);
    const std::array<const adlite::Scalar*, 4> components = {&stress.rr, &stress.zz, &stress.hoop, &stress.rz};
    result.stress = {stress.rr.value(), stress.zz.value(), stress.hoop.value(), stress.rz.value()};
    if (!compute_tangent)
        return result;
    std::array<double, 5> derivatives{};
    for (std::size_t row = 0; row < 4; ++row) {
        components[row]->copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t column = 0; column < 4; ++column)
            result.tangent[row][column] = derivatives[column];
        result.thermal[row] = derivatives[4];
    }
    return result;
}

void rotate_axisymmetric_strain_history(MaterialPointState& history, const AxisymmetricRotation& rotation) {
    for (auto* strain : {&history.elastic_strain, &history.plastic_strain, &history.creep_strain})
        *strain = rotate_axisymmetric_tensor_values(*strain, rotation);
    validate_material_point_state(history);
}

AxisymmetricStress compose_axisymmetric_stress(const AxisymmetricMaterialResponse& response,
    const std::array<adlite::Scalar, 5>& inputs,
    bool compute_tangent) {
    const std::array<double, 4> values = {response.stress.rr,
        response.stress.zz,
        response.stress.hoop,
        response.stress.rz};
    std::array<adlite::Scalar, 4> result;
    for (std::size_t i = 0; i < 4; ++i) {
        std::array<double, 5> partials{};
        for (std::size_t j = 0; j < 4; ++j)
            partials[j] = response.tangent[i][j];
        partials[4] = response.thermal[i];
        result[i] =
            compute_tangent ? adlite::compose(values[i], inputs.data(), partials.data(), 5) : adlite::Scalar(values[i]);
    }
    return {result[0], result[1], result[2], result[3]};
}

namespace {
template <typename Scalar>
std::array<Scalar, 6> rotate_cartesian_tensor_impl(const std::array<Scalar, 6>& tensor,
    const std::array<std::array<Scalar, 3>, 3>& rotation) {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, adlite::Scalar>,
        "Material tensor rotation only supports double and adlite::Scalar");
    const std::array<std::array<Scalar, 3>, 3> value = {{{{tensor[0], tensor[3], tensor[5]}},
        {{tensor[3], tensor[1], tensor[4]}},
        {{tensor[5], tensor[4], tensor[2]}}}};
    std::array<std::array<Scalar, 3>, 3> left{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                left[i][j] += rotation[i][k] * value[k][j];
    std::array<std::array<Scalar, 3>, 3> rotated{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = i; j < 3; ++j)
            for (std::size_t k = 0; k < 3; ++k)
                rotated[i][j] += left[i][k] * rotation[j][k];
    return {rotated[0][0], rotated[1][1], rotated[2][2], rotated[0][1], rotated[1][2], rotated[0][2]};
}
} // namespace

SymmetricTensor3 compose_cartesian_stress(const CartesianStressTangent& response,
    const SymmetricTensor3& strain,
    const adlite::Scalar& temperature,
    const std::array<double, 6>& thermal) {
    const std::array<adlite::Scalar, 7> inputs =
        {strain.xx, strain.yy, strain.zz, strain.xy, strain.yz, strain.xz, temperature};
    const std::array<double, 6> values = {response.stress.xx,
        response.stress.yy,
        response.stress.zz,
        response.stress.xy,
        response.stress.yz,
        response.stress.xz};
    std::array<double, 7> partials{};
    std::array<adlite::Scalar, 6> stress{};
    for (std::size_t row = 0; row < 6; ++row) {
        for (std::size_t column = 0; column < 6; ++column)
            partials[column] = response.tangent[row][column];
        partials[6] = thermal[row];
        stress[row] = adlite::compose(values[row], inputs.data(), partials.data(), inputs.size());
    }
    return {stress[0], stress[1], stress[2], stress[3], stress[4], stress[5]};
}

SymmetricTensor3Values rotate_cartesian_tensor_values(const SymmetricTensor3Values& tensor,
    const std::array<std::array<double, 3>, 3>& rotation) {
    const auto rotated =
        rotate_cartesian_tensor_impl<double>({tensor.xx, tensor.yy, tensor.zz, tensor.xy, tensor.yz, tensor.xz},
            rotation);
    return {rotated[0], rotated[1], rotated[2], rotated[3], rotated[4], rotated[5]};
}

CartesianInelasticStressResponse evaluate_incremental_cartesian_response(const IsotropicThermoelasticMaterial& material,
    const SymmetricTensor3& strain_increment,
    const adlite::Scalar& temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context) {
    if (!std::isfinite(committed_temperature) || !(committed_temperature > 0.0))
        throw std::domain_error("Incremental Cartesian material committed temperature must be finite and positive");
    MaterialFunctionContext old_context = context;
    old_context.time -= time_step;
    const SymmetricTensor3 old_imposed = material.eigenstrain(adlite::Scalar(committed_temperature), old_context);
    const SymmetricTensor3 synthetic_total{committed.elastic_strain[0] + strain_increment.xx + old_imposed.xx
                                               + committed.plastic_strain[0] + committed.creep_strain[0],
        committed.elastic_strain[1] + strain_increment.yy + old_imposed.yy + committed.plastic_strain[1]
            + committed.creep_strain[1],
        committed.elastic_strain[2] + strain_increment.zz + old_imposed.zz + committed.plastic_strain[2]
            + committed.creep_strain[2],
        committed.elastic_strain[3] + strain_increment.xy + old_imposed.xy + committed.plastic_strain[3]
            + committed.creep_strain[3],
        committed.elastic_strain[4] + strain_increment.yz + old_imposed.yz + committed.plastic_strain[4]
            + committed.creep_strain[4],
        committed.elastic_strain[5] + strain_increment.xz + old_imposed.xz + committed.plastic_strain[5]
            + committed.creep_strain[5]};
    return material.response(synthetic_total, temperature, time_step, committed, context);
}

SymmetricTensor3 rotate_cartesian_tensor(const SymmetricTensor3& tensor, const CartesianRotation& rotation) {
    const std::array<std::array<adlite::Scalar, 3>, 3> r = {{{rotation.xx, rotation.xy, rotation.xz},
        {rotation.yx, rotation.yy, rotation.yz},
        {rotation.zx, rotation.zy, rotation.zz}}};
    const auto rotated =
        rotate_cartesian_tensor_impl<adlite::Scalar>({tensor.xx, tensor.yy, tensor.zz, tensor.xy, tensor.yz, tensor.xz},
            r);
    return {rotated[0], rotated[1], rotated[2], rotated[3], rotated[4], rotated[5]};
}

CartesianInelasticStressResponse IsotropicThermoelasticMaterial::incremental_response(
    const SymmetricTensor3& strain_increment,
    const CartesianRotation& rotation,
    const adlite::Scalar& temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context) const {
    CartesianInelasticStressResponse result = evaluate_incremental_cartesian_response(*this,
        strain_increment,
        temperature,
        committed_temperature,
        time_step,
        committed,
        context);
    result.stress = rotate_cartesian_tensor(result.stress, rotation);
    std::array<double, 6>* histories[3] = {&result.trial_state.elastic_strain,
        &result.trial_state.plastic_strain,
        &result.trial_state.creep_strain};
    for (std::array<double, 6>* history : histories) {
        const SymmetricTensor3 rotated = rotate_cartesian_tensor(
            {(*history)[0], (*history)[1], (*history)[2], (*history)[3], (*history)[4], (*history)[5]},
            rotation);
        *history = {rotated.xx.value(),
            rotated.yy.value(),
            rotated.zz.value(),
            rotated.xy.value(),
            rotated.yz.value(),
            rotated.xz.value()};
        for (double component : *history)
            if (!std::isfinite(component))
                throw std::domain_error("Rotated Cartesian material history must contain only finite values");
    }
    result.trial_state.stress = {result.stress.xx.value(),
        result.stress.yy.value(),
        result.stress.zz.value(),
        result.stress.xy.value(),
        result.stress.yz.value(),
        result.stress.xz.value()};
    return result;
}

CartesianMaterialPointState IsotropicThermoelasticMaterial::incremental_response_values(
    const SymmetricTensor3Values& strain_increment,
    const CartesianRotation& rotation,
    double temperature,
    double committed_temperature,
    double time_step,
    const CartesianMaterialPointState& committed,
    MaterialFunctionContext context) const {
    if (!std::isfinite(committed_temperature) || !(committed_temperature > 0.0))
        throw std::domain_error("Incremental Cartesian material committed temperature must be finite and positive");
    MaterialFunctionContext old_context = context;
    old_context.time -= time_step;
    const SymmetricTensor3 old_imposed = eigenstrain(adlite::Scalar(committed_temperature), old_context);
    const SymmetricTensor3Values synthetic_total{committed.elastic_strain[0] + strain_increment.xx
                                                     + old_imposed.xx.value() + committed.plastic_strain[0]
                                                     + committed.creep_strain[0],
        committed.elastic_strain[1] + strain_increment.yy + old_imposed.yy.value() + committed.plastic_strain[1]
            + committed.creep_strain[1],
        committed.elastic_strain[2] + strain_increment.zz + old_imposed.zz.value() + committed.plastic_strain[2]
            + committed.creep_strain[2],
        committed.elastic_strain[3] + strain_increment.xy + old_imposed.xy.value() + committed.plastic_strain[3]
            + committed.creep_strain[3],
        committed.elastic_strain[4] + strain_increment.yz + old_imposed.yz.value() + committed.plastic_strain[4]
            + committed.creep_strain[4],
        committed.elastic_strain[5] + strain_increment.xz + old_imposed.xz.value() + committed.plastic_strain[5]
            + committed.creep_strain[5]};
    CartesianMaterialPointState result = response_values(synthetic_total, temperature, time_step, committed, context);
    const std::array<std::array<double, 3>, 3> rotation_values = {
        {{{rotation.xx.value(), rotation.xy.value(), rotation.xz.value()}},
            {{rotation.yx.value(), rotation.yy.value(), rotation.yz.value()}},
            {{rotation.zx.value(), rotation.zy.value(), rotation.zz.value()}}}};
    result.stress = rotate_cartesian_tensor_values(result.stress, rotation_values);
    std::array<double, 6>* histories[3] = {&result.elastic_strain, &result.plastic_strain, &result.creep_strain};
    for (std::array<double, 6>* history : histories) {
        const SymmetricTensor3Values rotated = rotate_cartesian_tensor_values(
            {(*history)[0], (*history)[1], (*history)[2], (*history)[3], (*history)[4], (*history)[5]},
            rotation_values);
        *history = {rotated.xx, rotated.yy, rotated.zz, rotated.xy, rotated.yz, rotated.xz};
        for (double component : *history)
            if (!std::isfinite(component))
                throw std::domain_error("Rotated Cartesian material history must contain only finite values");
    }
    return result;
}

MaterialFunctionContext material_context(double time, const CartesianPoint3& point) {
    return {time, point.x, point.y, point.z};
}

CartesianStressTangent evaluate_stress_tangent(const IsotropicThermoelasticMaterial& material,
    const std::array<double, 6>& fed_strain,
    double temperature,
    double time_step,
    const CartesianMaterialPointState* committed_material,
    MaterialFunctionContext context) {
    std::array<double, 7> seeds{};
    for (std::size_t component = 0; component < 6; ++component)
        seeds[component] = fed_strain[component];
    seeds[6] = temperature;
    std::array<adlite::Scalar, 7> active{};
    adlite::seed_identity(seeds.data(), seeds.size(), active.data());
    const SymmetricTensor3 strain{active[0], active[1], active[2], active[3], active[4], active[5]};
    const SymmetricTensor3 stress =
        committed_material == nullptr
            ? material.stress(strain, active[6], context)
            : material.response(strain, active[6], time_step, *committed_material, context).stress;
    const std::array<const adlite::Scalar*, 6> components =
        {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz};
    CartesianStressTangent result{};
    result.stress = {stress.xx.value(),
        stress.yy.value(),
        stress.zz.value(),
        stress.xy.value(),
        stress.yz.value(),
        stress.xz.value()};
    std::array<double, 7> derivatives{};
    for (std::size_t row = 0; row < 6; ++row) {
        components[row]->copy_derivatives(derivatives.data(), derivatives.size());
        for (std::size_t column = 0; column < 6; ++column)
            result.tangent[row][column] = derivatives[column];
        result.thermal[row] = derivatives[6];
    }
    return result;
}

SymmetricTensor3Values rotate_cartesian_tensor_values(const SymmetricTensor3Values& tensor,
    const CartesianRotation& rotation) {
    const std::array<std::array<double, 3>, 3> values = {
        {{{rotation.xx.value(), rotation.xy.value(), rotation.xz.value()}},
            {{rotation.yx.value(), rotation.yy.value(), rotation.yz.value()}},
            {{rotation.zx.value(), rotation.zy.value(), rotation.zz.value()}}}};
    return rotate_cartesian_tensor_values(tensor, values);
}
} // namespace fuelsim
