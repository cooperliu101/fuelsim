#ifndef FUELSIM_TEST_EXODUS_RESULT_READER_HPP
#define FUELSIM_TEST_EXODUS_RESULT_READER_HPP
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace fuelsim::test {
struct ExodusResults final {
    std::vector<std::array<double, 3>> nodes;
    std::vector<std::string> block_names;
    std::vector<std::size_t> block_element_counts;
    std::vector<std::string> side_set_names;
    std::vector<std::size_t> side_set_sizes;
    std::vector<std::vector<std::vector<std::size_t>>> side_set_face_nodes;
    std::vector<std::string> nodal_variable_names;
    std::vector<std::vector<double>> nodal_variables;
    std::vector<std::string> element_variable_names;
    std::vector<std::vector<double>> element_variables;
    std::vector<std::string> global_variable_names;
    std::vector<double> global_variables;
    std::size_t step_count = 0;
    double time = 0.0;

    const std::vector<double>& nodal(const std::string& name) const;
    const std::vector<double>& element(const std::string& name) const;
    double global(const std::string& name) const;
};

ExodusResults read_final_exodus_results(const std::string& path);
// Exodus step numbers are one-based; zero selects the final frame.
ExodusResults read_exodus_results(const std::string& path, std::size_t step);
// Includes the initial frame and keeps the file open for the complete read.
std::vector<ExodusResults> read_exodus_history(const std::string& path);
// Reads coordinates and nodal variables only; intended for contact and nodal-history comparisons.
std::vector<ExodusResults> read_exodus_nodal_history(const std::string& path);
} // namespace fuelsim::test
#endif
