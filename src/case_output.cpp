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

} // namespace fuelsim::app
