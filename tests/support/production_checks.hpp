#pragma once
#include <string>

namespace fuelsim::test {
bool check_hex8_multi_contact(const std::string& results, const std::string& nodal_reference,
    const std::string& first_reference, const std::string& second_reference);
bool check_hex20_nonmatching(const std::string& output, const std::string& displacement_path,
    const std::string& reaction_path, const std::string& pressure_path, const std::string& mortar_path,
    const std::string& mortar_reaction_path);
bool check_rz_multi_contact(const std::string& results, const std::string& first_reference,
    const std::string& second_reference, double friction);
} // namespace fuelsim::test
