#include "case_output.hpp"
#include "transient_conservation.hpp"

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
    for (const TransientConservationField& field :
         transient_conservation_fields)
        output.value(prefix + field.name, summary.*field.member);
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
