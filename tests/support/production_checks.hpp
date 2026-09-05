#pragma once
#include <string>

namespace fuelsim::test {
bool check_hex8_sliding(const std::string& output_path, const std::string& nodes_path, const std::string& contact_path,
    const std::string& reaction_path);
bool check_hex20_curved_friction(
    const std::string& output_path, const std::string& displacement_path, const std::string& contact_path);
bool check_hex20_curved(const std::string& output_path, const std::string& name, const std::string& displacement_path,
    const std::string& force_path, const std::string& resultant_path);
bool check_hex8_sticking(const std::string& output_path, const std::string& thermal_path,
    const std::string& mechanical_path, const std::string& contact_path, const std::string& baseline_path);
bool check_hex8_shared_plate(
    const std::string& output_path, const std::string& node_path, const std::string& element_path);
bool check_rz_noncoaxial(
    const std::string& output_path, const std::string& variant_name, const std::string& history_reference_path);
bool check_rz_inelastic(const std::string& output_path, const std::string& branch, const std::string& history_path);
bool check_rz_pcmi(const std::string& output_path, const std::string& coordinate_path, const std::string& value_path,
    const std::string& scalar_path, const std::string& contact_path, double tolerance);
bool check_rz_integrated(const std::string& output_path, const std::string& coordinate_path,
    const std::string& value_path, const std::string& scalar_path);
bool check_hex8_multi_contact(const std::string& results, const std::string& nodal_reference,
    const std::string& first_reference, const std::string& second_reference);
bool check_hex20_nonmatching(const std::string& output, const std::string& displacement_path,
    const std::string& reaction_path, const std::string& pressure_path, const std::string& mortar_path,
    const std::string& mortar_reaction_path);
bool check_rz_multi_contact(const std::string& results, const std::string& first_reference,
    const std::string& second_reference, double friction);
} // namespace fuelsim::test
