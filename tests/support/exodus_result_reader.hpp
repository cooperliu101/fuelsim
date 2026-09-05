#ifndef FUELSIM_TEST_EXODUS_RESULT_READER_HPP
#define FUELSIM_TEST_EXODUS_RESULT_READER_HPP
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace fuelsim::test {
struct ExodusResults final {
    std::vector<std::array<double, 3>> nodes;
    std::vector<std::string> nodal_variable_names;
    std::vector<std::vector<double>> nodal_variables;
    std::vector<std::string> element_variable_names;
    std::vector<std::vector<double>> element_variables;
    std::size_t step_count = 0;
    double time = 0.0;

    const std::vector<double>& nodal(const std::string& name) const;
    const std::vector<double>& element(const std::string& name) const;
};

ExodusResults read_final_exodus_results(const std::string& path);
} // namespace fuelsim::test
#endif
