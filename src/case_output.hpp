#ifndef FUELSIM_CASE_OUTPUT_HPP
#define FUELSIM_CASE_OUTPUT_HPP

#include "fuelsim/case_input.hpp"

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

} // namespace fuelsim::app

#endif
