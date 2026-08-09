#ifndef FUELSIM_CASE_OUTPUT_HPP
#define FUELSIM_CASE_OUTPUT_HPP

#include "fuelsim/case_input.hpp"
#include "fuelsim/problem_solver.hpp"

#include <cstddef>
#include <fstream>
#include <string>

namespace fuelsim::app {

class CaseOutput final {
  public:
    explicit CaseOutput(bool console);
    CaseOutput(const CaseOutputInput& options, bool force_console,
               bool active);

    void value(const std::string& key, const std::string& data);
    void value(const std::string& key, const char* data);
    void value(const std::string& key, double data);
    void value(const std::string& key, std::size_t data);
    void value(const std::string& key, int data);
    void value(const std::string& key, bool data);

  private:
    bool _console;
    std::ofstream _csv;
};

void write_conservation_summary(
    const std::string& prefix,
    const TransientConservationSummary& summary,
    CaseOutput& output);
void write_time_error_components(
    const std::string& prefix,
    const TransientTimeErrorEstimate& estimate,
    CaseOutput& output);

} // namespace fuelsim::app

#endif
