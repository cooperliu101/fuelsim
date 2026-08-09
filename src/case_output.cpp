#include "case_output.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace fuelsim::app {

CaseOutput::CaseOutput(bool console) : _console(console) {
    if (_console)
        std::cout << std::boolalpha << std::scientific
                  << std::setprecision(12);
}

CaseOutput::CaseOutput(const CaseOutputInput& options, bool force_console,
                       bool active)
    : CaseOutput(active && (options.console || force_console)) {
    if (!active || options.csv_file.empty())
        return;
    _csv.open(options.csv_file, std::ios::out | std::ios::trunc);
    if (!_csv)
        throw std::runtime_error("Could not open CSV output file '" +
                                 options.csv_file + "'");
    _csv << "metric,value\n" << std::scientific << std::setprecision(12);
}

void CaseOutput::value(const std::string& key, const std::string& data) {
    if (_console)
        std::cout << key << '=' << data << '\n';
    if (_csv)
        _csv << key << ',' << data << '\n';
}

void CaseOutput::value(const std::string& key, const char* data) {
    value(key, std::string(data));
}

void CaseOutput::value(const std::string& key, double data) {
    if (_console)
        std::cout << key << '=' << data << '\n';
    if (_csv)
        _csv << key << ',' << data << '\n';
}

void CaseOutput::value(const std::string& key, std::size_t data) {
    if (_console)
        std::cout << key << '=' << data << '\n';
    if (_csv)
        _csv << key << ',' << data << '\n';
}

void CaseOutput::value(const std::string& key, int data) {
    if (_console)
        std::cout << key << '=' << data << '\n';
    if (_csv)
        _csv << key << ',' << data << '\n';
}

void CaseOutput::value(const std::string& key, bool data) {
    if (_console)
        std::cout << key << '=' << std::boolalpha << data << '\n';
    if (_csv)
        _csv << key << ',' << std::boolalpha << data << '\n';
}

void write_conservation_summary(
    const std::string& prefix,
    const TransientConservationSummary& summary,
    CaseOutput& output) {
    output.value(prefix + "generated_heat_rate", summary.generated_heat_rate);
    output.value(prefix + "stored_heat_rate", summary.stored_heat_rate);
    output.value(prefix + "convection_heat_rate", summary.convection_heat_rate);
    output.value(prefix + "interface_heat_imbalance",
                 summary.interface_heat_imbalance);
    output.value(prefix + "dirichlet_heat_input_rate",
                 summary.dirichlet_heat_input_rate);
    output.value(prefix + "global_thermal_balance",
                 summary.global_thermal_balance);
    output.value(prefix + "relative_thermal_balance",
                 summary.relative_thermal_balance);
    output.value(prefix + "unconstrained_thermal_residual_l2",
                 summary.unconstrained_thermal_residual_l2);
    output.value(prefix + "internal_mechanical_work_increment",
                 summary.internal_mechanical_work_increment);
    output.value(prefix + "pressure_traction_work_increment",
                 summary.pressure_traction_work_increment);
    output.value(prefix + "dirichlet_reaction_work_increment",
                 summary.dirichlet_reaction_work_increment);
    output.value(prefix + "contact_work_increment",
                 summary.contact_work_increment);
    output.value(prefix + "mechanical_work_balance",
                 summary.mechanical_work_balance);
    output.value(prefix + "relative_mechanical_work_balance",
                 summary.relative_mechanical_work_balance);
    output.value(prefix + "unconstrained_mechanical_residual_l2",
                 summary.unconstrained_mechanical_residual_l2);
    output.value(prefix + "elastic_energy_change",
                 summary.elastic_energy_change);
    output.value(prefix + "plastic_dissipation_increment",
                 summary.plastic_dissipation_increment);
    output.value(prefix + "creep_dissipation_increment",
                 summary.creep_dissipation_increment);
}

void write_time_error_components(
    const std::string& prefix,
    const TransientTimeErrorEstimate& estimate,
    CaseOutput& output) {
    output.value(prefix + "temperature", estimate.temperature);
    output.value(prefix + "radial_displacement",
                 estimate.radial_displacement);
    output.value(prefix + "axial_displacement", estimate.axial_displacement);
    output.value(prefix + "elastic_strain", estimate.elastic_strain);
    output.value(prefix + "plastic_strain", estimate.plastic_strain);
    output.value(prefix + "creep_strain", estimate.creep_strain);
    output.value(prefix + "equivalent_plastic_strain",
                 estimate.equivalent_plastic_strain);
    output.value(prefix + "equivalent_creep_strain",
                 estimate.equivalent_creep_strain);
    output.value(prefix + "stress", estimate.stress);
    output.value(prefix + "contact_friction", estimate.contact_friction);
    output.value(prefix + "contact_normal_multiplier",
                 estimate.contact_normal_multiplier);
}

} // namespace fuelsim::app
