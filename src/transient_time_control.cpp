#include "transient_time_control.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fuelsim::time_control {
namespace {

struct ErrorAccumulator final {
    double difference_squared = 0.0;
    double solution_squared = 0.0;
    std::size_t count = 0;
};

void accumulate_error(ErrorAccumulator& accumulator, double full_step,
                      double two_half_steps) {
    const double difference = two_half_steps - full_step;
    accumulator.difference_squared += difference * difference;
    accumulator.solution_squared += two_half_steps * two_half_steps;
    ++accumulator.count;
}

double normalized_error(const ErrorAccumulator& accumulator,
                        double absolute_tolerance,
                        double relative_tolerance) {
    if (accumulator.count == 0)
        return 0.0;
    const double denominator =
        absolute_tolerance *
            std::sqrt(static_cast<double>(accumulator.count)) +
        relative_tolerance * std::sqrt(accumulator.solution_squared);
    return std::sqrt(accumulator.difference_squared) / denominator;
}

} // namespace

TransientConservationSummary combine_half_step_conservation(
    const TransientConservationSummary& first,
    const TransientConservationSummary& second) {
    TransientConservationSummary result;
    const auto average = [](double left, double right) {
        return 0.5 * (left + right);
    };
    result.generated_heat_rate =
        average(first.generated_heat_rate, second.generated_heat_rate);
    result.stored_heat_rate =
        average(first.stored_heat_rate, second.stored_heat_rate);
    result.convection_heat_rate =
        average(first.convection_heat_rate, second.convection_heat_rate);
    result.interface_heat_imbalance = average(
        first.interface_heat_imbalance, second.interface_heat_imbalance);
    result.dirichlet_heat_input_rate = average(
        first.dirichlet_heat_input_rate, second.dirichlet_heat_input_rate);
    result.global_thermal_balance =
        result.stored_heat_rate + result.convection_heat_rate +
        result.interface_heat_imbalance - result.generated_heat_rate -
        result.dirichlet_heat_input_rate;
    const double thermal_scale =
        std::abs(result.generated_heat_rate) +
        std::abs(result.stored_heat_rate) +
        std::abs(result.convection_heat_rate) +
        std::abs(result.interface_heat_imbalance) +
        std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0
            ? std::abs(result.global_thermal_balance) / thermal_scale
            : 0.0;
    result.unconstrained_thermal_residual_l2 = std::max(
        first.unconstrained_thermal_residual_l2,
        second.unconstrained_thermal_residual_l2);

    result.internal_mechanical_work_increment =
        first.internal_mechanical_work_increment +
        second.internal_mechanical_work_increment;
    result.pressure_traction_work_increment =
        first.pressure_traction_work_increment +
        second.pressure_traction_work_increment;
    result.dirichlet_reaction_work_increment =
        first.dirichlet_reaction_work_increment +
        second.dirichlet_reaction_work_increment;
    result.contact_work_increment = first.contact_work_increment +
                                    second.contact_work_increment;
    result.mechanical_work_balance =
        result.internal_mechanical_work_increment +
        result.contact_work_increment -
        result.pressure_traction_work_increment -
        result.dirichlet_reaction_work_increment;
    const double mechanical_scale =
        std::abs(result.internal_mechanical_work_increment) +
        std::abs(result.pressure_traction_work_increment) +
        std::abs(result.dirichlet_reaction_work_increment) +
        std::abs(result.contact_work_increment);
    result.relative_mechanical_work_balance =
        mechanical_scale > 0.0
            ? std::abs(result.mechanical_work_balance) / mechanical_scale
            : 0.0;
    result.unconstrained_mechanical_residual_l2 = std::max(
        first.unconstrained_mechanical_residual_l2,
        second.unconstrained_mechanical_residual_l2);
    result.elastic_energy_change = first.elastic_energy_change +
                                   second.elastic_energy_change;
    result.plastic_dissipation_increment =
        first.plastic_dissipation_increment +
        second.plastic_dissipation_increment;
    result.creep_dissipation_increment =
        first.creep_dissipation_increment +
        second.creep_dissipation_increment;
    return result;
}

TransientTimeErrorEstimate step_doubling_error(
    const TransientCommittedState& full_step,
    const TransientCommittedState& two_half_steps,
    const TransientTimeOptions& options) {
    if (full_step.solution.size() != two_half_steps.solution.size() ||
        full_step.solution.size() % 3 != 0)
        throw std::logic_error(
            "step-doubling states must share the [T, ur, uz] layout");
    if (full_step.material_histories.size() !=
            two_half_steps.material_histories.size() ||
        full_step.material_stresses.size() !=
            two_half_steps.material_stresses.size())
        throw std::logic_error(
            "step-doubling material-state region layouts differ");

    const std::size_t node_count = full_step.solution.size() / 3;
    std::array<ErrorAccumulator, 3> nodal{};
    for (std::size_t field = 0; field < 3; ++field) {
        const std::size_t begin = field * node_count;
        const std::size_t end = begin + node_count;
        for (std::size_t dof = begin; dof < end; ++dof)
            accumulate_error(nodal[field], full_step.solution[dof],
                             two_half_steps.solution[dof]);
    }

    ErrorAccumulator elastic;
    ErrorAccumulator plastic;
    ErrorAccumulator creep;
    ErrorAccumulator equivalent_plastic;
    ErrorAccumulator equivalent_creep;
    ErrorAccumulator stress;
    ErrorAccumulator contact_friction;
    ErrorAccumulator contact_normal_multiplier;
    bool contact_state_mismatch = false;
    for (std::size_t region = 0;
         region < full_step.material_histories.size(); ++region) {
        const auto& full_history = full_step.material_histories[region];
        const auto& half_history = two_half_steps.material_histories[region];
        const auto& full_stress = full_step.material_stresses[region];
        const auto& half_stress = two_half_steps.material_stresses[region];
        if (full_history.size() != half_history.size() ||
            full_stress.size() != half_stress.size() ||
            full_history.size() != full_stress.size())
            throw std::logic_error(
                "step-doubling material-state element layouts differ");
        for (std::size_t element = 0; element < full_history.size();
             ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const MaterialPointState& full_point =
                    full_history[element][q];
                const MaterialPointState& half_point =
                    half_history[element][q];
                for (std::size_t component = 0; component < 4; ++component) {
                    accumulate_error(elastic,
                                     full_point.elastic_strain[component],
                                     half_point.elastic_strain[component]);
                    accumulate_error(plastic,
                                     full_point.plastic_strain[component],
                                     half_point.plastic_strain[component]);
                    accumulate_error(creep,
                                     full_point.creep_strain[component],
                                     half_point.creep_strain[component]);
                }
                accumulate_error(equivalent_plastic,
                                 full_point.equivalent_plastic_strain,
                                 half_point.equivalent_plastic_strain);
                accumulate_error(equivalent_creep,
                                 full_point.equivalent_creep_strain,
                                 half_point.equivalent_creep_strain);

                const AxisymmetricStressValues& full_value =
                    full_stress[element][q];
                const AxisymmetricStressValues& half_value =
                    half_stress[element][q];
                accumulate_error(stress, full_value.rr, half_value.rr);
                accumulate_error(stress, full_value.zz, half_value.zz);
                accumulate_error(stress, full_value.hoop, half_value.hoop);
                accumulate_error(stress, full_value.rz, half_value.rz);
            }
        }
    }
    if (full_step.contact_histories.size() !=
        two_half_steps.contact_histories.size())
        throw std::logic_error(
            "step-doubling contact-history layouts differ");
    for (std::size_t contact = 0;
         contact < full_step.contact_histories.size(); ++contact) {
        if (full_step.contact_histories[contact].size() !=
            two_half_steps.contact_histories[contact].size())
            throw std::logic_error(
                "step-doubling contact-node history layouts differ");
        for (std::size_t node = 0;
             node < full_step.contact_histories[contact].size(); ++node) {
            const ContactPointHistory& full =
                full_step.contact_histories[contact][node];
            const ContactPointHistory& half =
                two_half_steps.contact_histories[contact][node];
            accumulate_error(contact_friction,
                             full.elastic_tangential_slip,
                             half.elastic_tangential_slip);
            accumulate_error(contact_normal_multiplier,
                             full.normal_multiplier,
                             half.normal_multiplier);
            contact_state_mismatch =
                contact_state_mismatch || full.sliding != half.sliding;
        }
    }

    TransientTimeErrorEstimate result;
    result.temperature = normalized_error(
        nodal[0], options.temperature_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.radial_displacement = normalized_error(
        nodal[1], options.displacement_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.axial_displacement = normalized_error(
        nodal[2], options.displacement_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.elastic_strain = normalized_error(
        elastic, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.plastic_strain = normalized_error(
        plastic, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.creep_strain = normalized_error(
        creep, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.equivalent_plastic_strain = normalized_error(
        equivalent_plastic, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.equivalent_creep_strain = normalized_error(
        equivalent_creep, options.strain_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.stress = normalized_error(
        stress, options.stress_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.contact_friction =
        contact_state_mismatch
            ? std::numeric_limits<double>::infinity()
            : normalized_error(
                  contact_friction,
                  options.displacement_time_absolute_tolerance,
                  options.time_error_relative_tolerance);
    result.contact_normal_multiplier = normalized_error(
        contact_normal_multiplier,
        options.stress_history_time_absolute_tolerance,
        options.time_error_relative_tolerance);
    result.maximum = std::max(
        {result.temperature, result.radial_displacement,
         result.axial_displacement, result.elastic_strain,
         result.plastic_strain, result.creep_strain,
         result.equivalent_plastic_strain,
         result.equivalent_creep_strain, result.stress,
         result.contact_friction, result.contact_normal_multiplier});
    return result;
}

double step_factor(const TransientTimeOptions& options, double error) {
    if (!(error > 0.0))
        return options.growth_factor;
    return std::clamp(options.time_error_safety_factor / std::sqrt(error),
                      0.1, options.growth_factor);
}

} // namespace fuelsim::time_control
