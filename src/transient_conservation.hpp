#ifndef FUELSIM_TRANSIENT_CONSERVATION_HPP
#define FUELSIM_TRANSIENT_CONSERVATION_HPP

#include "fuelsim/transient_problem.hpp"

#include <array>

namespace fuelsim {

struct TransientConservationField final {
    const char* name;
    double TransientConservationSummary::*member;
};

inline constexpr std::array<TransientConservationField, 18>
    transient_conservation_fields = {{
        {"generated_heat_rate",
         &TransientConservationSummary::generated_heat_rate},
        {"stored_heat_rate", &TransientConservationSummary::stored_heat_rate},
        {"convection_heat_rate",
         &TransientConservationSummary::convection_heat_rate},
        {"interface_heat_imbalance",
         &TransientConservationSummary::interface_heat_imbalance},
        {"dirichlet_heat_input_rate",
         &TransientConservationSummary::dirichlet_heat_input_rate},
        {"global_thermal_balance",
         &TransientConservationSummary::global_thermal_balance},
        {"relative_thermal_balance",
         &TransientConservationSummary::relative_thermal_balance},
        {"unconstrained_thermal_residual_l2",
         &TransientConservationSummary::unconstrained_thermal_residual_l2},
        {"internal_mechanical_work_increment",
         &TransientConservationSummary::internal_mechanical_work_increment},
        {"pressure_traction_work_increment",
         &TransientConservationSummary::pressure_traction_work_increment},
        {"dirichlet_reaction_work_increment",
         &TransientConservationSummary::dirichlet_reaction_work_increment},
        {"contact_work_increment",
         &TransientConservationSummary::contact_work_increment},
        {"mechanical_work_balance",
         &TransientConservationSummary::mechanical_work_balance},
        {"relative_mechanical_work_balance",
         &TransientConservationSummary::relative_mechanical_work_balance},
        {"unconstrained_mechanical_residual_l2",
         &TransientConservationSummary::unconstrained_mechanical_residual_l2},
        {"elastic_energy_change",
         &TransientConservationSummary::elastic_energy_change},
        {"plastic_dissipation_increment",
         &TransientConservationSummary::plastic_dissipation_increment},
        {"creep_dissipation_increment",
         &TransientConservationSummary::creep_dissipation_increment},
    }};

class TransientConservationCalculator final {
  public:
    static TransientConservationSummary summarize(
        const TransientProblem& problem,
        const std::vector<double>& converged_solution,
        const std::vector<std::vector<Quad4MaterialHistory>>& staged_histories,
        const std::vector<
            std::vector<std::array<AxisymmetricStressValues, 4>>>&
            staged_stresses);
};

} // namespace fuelsim

#endif
