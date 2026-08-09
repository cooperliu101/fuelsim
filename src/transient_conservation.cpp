#include "fuelsim/transient_problem.hpp"

#include "spatial_assembly.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace fuelsim {
namespace {

double stress_strain_inner_product(
    const AxisymmetricStressValues& stress,
    const std::array<double, 4>& strain) noexcept {
    return stress.rr * strain[0] + stress.zz * strain[1] +
           stress.hoop * strain[2] + 2.0 * stress.rz * strain[3];
}

std::array<double, 4> strain_difference(const std::array<double, 4>& current,
                                        const std::array<double, 4>& old) {
    std::array<double, 4> result{};
    for (std::size_t component = 0; component < result.size(); ++component)
        result[component] = current[component] - old[component];
    return result;
}

} // namespace

TransientConservationSummary TransientProblem::summarize_active_step(
    const std::vector<double>& converged_solution,
    const std::vector<std::vector<Quad4MaterialHistory>>& staged_histories,
    const std::vector<std::vector<std::array<AxisymmetricStressValues, 4>>>&
        staged_stresses) const {
    require_active_time_step();
    TransientConservationSummary result;
    std::vector<double> raw_residual(dof_count(), 0.0);

    for (std::size_t contribution = 0; contribution < contribution_count();
         ++contribution) {
        const LocalDofs dofs = contribution_dofs(contribution);
        const LocalValues state =
            contribution_state(contribution, converged_solution);
        const LocalResidual residual =
            contribution_residual(contribution, state);
        const SpatialContributionType type =
            _spatial->contribution_type(contribution);
        for (std::size_t local = 0; local < local_dof_count; ++local) {
            raw_residual[dofs[local]] += residual[local];
            const double increment =
                converged_solution[dofs[local]] -
                _committed_solution[dofs[local]];
            if (local < 4) {
                if (type == SpatialContributionType::thermal_contact)
                    result.interface_heat_imbalance += residual[local];
                else if (type == SpatialContributionType::convection)
                    result.convection_heat_rate += residual[local];
                continue;
            }
            const double work = residual[local] * increment;
            if (type == SpatialContributionType::volume)
                result.internal_mechanical_work_increment += work;
            else if (type == SpatialContributionType::mechanical_contact)
                result.contact_work_increment += work;
            else if (type == SpatialContributionType::pressure ||
                     type == SpatialContributionType::traction)
                result.pressure_traction_work_increment -= work;
        }
    }
    for (std::size_t region_value = 0; region_value < region_count();
         ++region_value) {
        const Quad4RzTransientKernel& kernel =
            _region_kernels[region_value];
        const TransientInelasticProperties& properties = kernel.properties();
        const double heat_capacity =
            properties.density * properties.specific_heat;
        const std::size_t offset =
            _spatial->region_element_offset(region_value);
        for (std::size_t element = 0;
             element < staged_histories[region_value].size(); ++element) {
            const LocalValues current =
                contribution_state(offset + element, converged_solution);
            const LocalValues old =
                contribution_state(offset + element, _committed_solution);
            const Quad4RzGeometry& geometry =
                region_element_geometry(region_value, element);
            for (std::size_t q = 0; q < geometry.points.size(); ++q) {
                const RzQuadraturePoint& point = geometry.points[q];
                double current_temperature = 0.0;
                double old_temperature = 0.0;
                for (std::size_t node = 0; node < quad4_node_count; ++node) {
                    current_temperature += point.shape[node] * current[node];
                    old_temperature += point.shape[node] * old[node];
                }
                result.stored_heat_rate +=
                    point.weighted_measure * heat_capacity *
                    (current_temperature - old_temperature) /
                    _active_time_step;
                result.generated_heat_rate +=
                    point.weighted_measure * kernel.volumetric_heat_source();

                const MaterialPointState& old_history =
                    _material_histories[region_value][element][q];
                const MaterialPointState& new_history =
                    staged_histories[region_value][element][q];
                const AxisymmetricStressValues& old_stress =
                    _material_stresses[region_value][element][q];
                const AxisymmetricStressValues& new_stress =
                    staged_stresses[region_value][element][q];
                result.elastic_energy_change +=
                    0.5 * point.weighted_measure *
                    (stress_strain_inner_product(
                         new_stress, new_history.elastic_strain) -
                     stress_strain_inner_product(
                         old_stress, old_history.elastic_strain));
                result.plastic_dissipation_increment +=
                    point.weighted_measure * stress_strain_inner_product(
                        new_stress,
                        strain_difference(new_history.plastic_strain,
                                          old_history.plastic_strain));
                result.creep_dissipation_increment +=
                    point.weighted_measure * stress_strain_inner_product(
                        new_stress,
                        strain_difference(new_history.creep_strain,
                                          old_history.creep_strain));
            }
        }
    }

    std::vector<bool> constrained(dof_count(), false);
    for (const DirichletCondition& condition : dirichlet_conditions()) {
        constrained[condition.dof] = true;
        const double increment = converged_solution[condition.dof] -
                                 _committed_solution[condition.dof];
        if (condition.dof < dof_map().node_count())
            result.dirichlet_heat_input_rate += raw_residual[condition.dof];
        else
            result.dirichlet_reaction_work_increment +=
                raw_residual[condition.dof] * increment;
    }
    double thermal_residual_squared = 0.0;
    double mechanical_residual_squared = 0.0;
    for (std::size_t dof = 0; dof < dof_count(); ++dof) {
        if (constrained[dof])
            continue;
        const double value = raw_residual[dof];
        if (dof < dof_map().node_count())
            thermal_residual_squared += value * value;
        else
            mechanical_residual_squared += value * value;
    }
    result.unconstrained_thermal_residual_l2 =
        std::sqrt(thermal_residual_squared);
    result.unconstrained_mechanical_residual_l2 =
        std::sqrt(mechanical_residual_squared);
    result.global_thermal_balance =
        result.stored_heat_rate + result.convection_heat_rate +
        result.interface_heat_imbalance - result.generated_heat_rate -
        result.dirichlet_heat_input_rate;
    const double thermal_scale =
        std::abs(result.stored_heat_rate) +
        std::abs(result.convection_heat_rate) +
        std::abs(result.interface_heat_imbalance) +
        std::abs(result.generated_heat_rate) +
        std::abs(result.dirichlet_heat_input_rate);
    result.relative_thermal_balance =
        thermal_scale > 0.0
            ? std::abs(result.global_thermal_balance) / thermal_scale
            : 0.0;
    result.mechanical_work_balance =
        result.internal_mechanical_work_increment +
        result.contact_work_increment -
        result.pressure_traction_work_increment -
        result.dirichlet_reaction_work_increment;
    const double mechanical_scale =
        std::abs(result.internal_mechanical_work_increment) +
        std::abs(result.contact_work_increment) +
        std::abs(result.pressure_traction_work_increment) +
        std::abs(result.dirichlet_reaction_work_increment);
    result.relative_mechanical_work_balance =
        mechanical_scale > 0.0
            ? std::abs(result.mechanical_work_balance) / mechanical_scale
            : 0.0;
    return result;
}

} // namespace fuelsim
