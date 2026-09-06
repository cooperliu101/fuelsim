#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::vector<std::string> split(const std::string& text, char delimiter) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t end = text.find(delimiter, begin);
        fields.push_back(text.substr(begin, end - begin));
        if (end == std::string::npos) return fields;
        begin = end + 1;
    }
}

std::set<std::string> read_registered_tests(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read registered CTest list: " + path);
    std::set<std::string> tests;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && !tests.insert(line).second)
            throw std::runtime_error("Duplicate registered CTest: " + line);
    }
    return tests;
}

void check_evidence(
    const std::filesystem::path& repository, const std::string& row_id, const std::string& evidence_text) {
    if (evidence_text == "-") throw std::runtime_error(row_id + " has no evidence");
    for (const std::string& relative : split(evidence_text, ';')) {
        const std::filesystem::path path(relative);
        if (relative.empty() || path.is_absolute() || relative.find("..") != std::string::npos)
            throw std::runtime_error(row_id + " has an invalid evidence path: " + relative);
        const std::filesystem::path absolute = repository / path;
        if (!std::filesystem::is_regular_file(absolute))
            throw std::runtime_error(row_id + " evidence is not a regular file: " + relative);
        if (std::filesystem::file_size(absolute) == 0)
            throw std::runtime_error(row_id + " evidence is empty: " + relative);
    }
}

void check_c3d8t_contract(const std::string& path, const std::set<std::string>& registered_tests) {
    const std::string nodal_contract = "T,U1,U2,U3,RFL,RF1,RF2,RF3";
    const std::string integration_contract =
        "HFL1,HFL2,HFL3,S11,S22,S33,S12,S13,S23,E11,E22,E33,E12,E13,E23,active_history";
    const std::string contact_contract = "opening,pressure,slip1,slip2,contact_force,resultant,moment,force_center";
    const std::string zero_contract = "separate_absolute_no_denominator_floor";
    const std::set<std::string> existing_scoped_tests = {
        "fuelsim_b49_hex8_c3d8t_operator_abaqus_tests",
        "fuelsim_b50_hex8_c3d8t_capacity_abaqus_tests",
        "fuelsim_b51_hex8_c3d8t_finite_heat_abaqus_tests",
        "fuelsim_b52_hex8_c3d8t_thermal_contact_abaqus_tests",
        "fuelsim_b53_hex8_c3d8t_thermal_load_abaqus_tests",
        "fuelsim_b54_hex8_c3d8t_finite_thermal_load_abaqus_tests",
        "fuelsim_b56_hex8_c3d8t_temperature_operator_abaqus_tests",
        "fuelsim_b57_hex8_c3d8t_temperature_capacity_abaqus_tests",
        "fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests",
        "fuelsim_b511_hex8_c3d8t_small_norton_abaqus_tests",
        "fuelsim_b512_hex8_c3d8t_small_coupled_abaqus_tests",
        "fuelsim_b513_hex8_c3d8t_small_noncoaxial_abaqus_tests",
        "fuelsim_b514_hex8_c3d8t_finite_elastic_abaqus_tests",
        "fuelsim_b515_hex8_c3d8t_finite_j2_abaqus_tests",
        "fuelsim_b516_hex8_c3d8t_finite_norton_abaqus_tests",
        "fuelsim_b517_hex8_c3d8t_finite_coupled_abaqus_tests",
        "fuelsim_b518_hex8_c3d8t_finite_noncoaxial_abaqus_tests",
        "fuelsim_b519_hex8_c3d8t_finite_selective_abaqus_tests",
        "fuelsim_b520_hex8_c3d8t_gap_conductance_abaqus_tests",
        "fuelsim_b521_hex8_c3d8t_thermal_contact_path_abaqus_tests",
        "fuelsim_b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests",
        "fuelsim_b552_hex8_c3d8t_sts_cross_face_abaqus_tests",
    };

    std::set<std::string> registered_c3d8t;
    for (const std::string& test : registered_tests)
        if (test.find("hex8_c3d8t") != std::string::npos && test.find("_abaqus_tests") != std::string::npos)
            registered_c3d8t.insert(test);

    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read C3D8T validation contract: " + path);
    std::string line;
    if (!std::getline(input, line) || line != "ctest\tscope\tall_time_steps\tnodal_fields\tintegration_fields\tcontact_"
                                              "fields\tacceptance_class\tzero_reference_policy\tknown_boundary")
        throw std::runtime_error("C3D8T validation contract header does not match schema");
    std::set<std::string> found_tests;
    std::size_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) throw std::runtime_error("Blank C3D8T contract row at line " + std::to_string(line_number));
        const std::vector<std::string> fields = split(line, '\t');
        if (fields.size() != 9)
            throw std::runtime_error("C3D8T contract row must have nine fields at line " + std::to_string(line_number));
        for (const std::string& field : fields)
            if (field.empty())
                throw std::runtime_error("Empty C3D8T contract field at line " + std::to_string(line_number));
        const std::string& test = fields[0];
        const std::string& scope = fields[1];
        if (registered_c3d8t.count(test) == 0)
            throw std::runtime_error("C3D8T contract names an unregistered test: " + test);
        if (!found_tests.insert(test).second) throw std::runtime_error("Duplicate C3D8T contract test: " + test);
        if (fields[7] != zero_contract)
            throw std::runtime_error(test + " does not use the required zero-reference policy");

        const bool complete = scope == "full_field" || scope == "contact_full_field";
        if (!complete && existing_scoped_tests.count(test) == 0)
            throw std::runtime_error(test + " is a new C3D8T test without complete full-field coverage");
        if (!complete) continue;
        if (fields[2] != "yes" || fields[3] != nodal_contract || fields[4] != integration_contract)
            throw std::runtime_error(
                test + " does not satisfy the complete C3D8T time, nodal, and integration contract");
        if ((scope == "contact_full_field" && fields[5] != contact_contract) ||
            (scope == "full_field" && fields[5] != "none"))
            throw std::runtime_error(test + " does not satisfy its C3D8T contact-field contract");
        if (fields[6] != "base_0.1_percent" && fields[6] != "integrated_finite_0.5_percent" &&
            fields[6] != "integrated_contact_0.5_percent" && fields[6] != "contact_1_percent_complete_force" &&
            fields[6] != "qualified_contact")
            throw std::runtime_error(test + " has an invalid C3D8T acceptance class");
        if (fields[6] == "qualified_contact" && (scope != "contact_full_field" || fields[8] == "none"))
            throw std::runtime_error(test + " lacks a documented contact qualification boundary");
        if (fields[6] == "contact_1_percent_complete_force" &&
            (test != "fuelsim_b522_hex8_c3d8t_faceted_thermal_contact_abaqus_tests" || scope != "contact_path" ||
                fields[8] == "none"))
            throw std::runtime_error(test + " misuses the B5.22 complete-force acceptance class");
    }
    if (found_tests != registered_c3d8t)
        throw std::runtime_error("C3D8T validation contract does not cover every registered Abaqus test");
    std::cout << "c3d8t_validation_contract_rows=" << found_tests.size() << '\n';
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_verification_matrix_tests "
                     "<matrix.tsv> <repository> <registered-tests> <c3d8t-contract.tsv>\n";
        return 2;
    }
    try {
        const std::string matrix_path = argv[1];
        const std::filesystem::path repository = argv[2];
        const std::set<std::string> registered_tests = read_registered_tests(argv[3]);
        std::ifstream matrix(matrix_path);
        if (!matrix) throw std::runtime_error("Could not read verification matrix: " + matrix_path);
        std::string line;
        if (!std::getline(matrix, line) || line != "id\tstatus\tcapability\tctest\tevidence\tacceptance")
            throw std::runtime_error("Verification matrix header does not match schema");
        const std::set<std::string> required_ids = {
            "b100.cax8t_small_probe",
            "b101.cax8t_finite_probe",
            "b102.cax8t_small_contact",
            "b103.cax8t_finite_contact",
            "b104.cax8t_small_coupled",
            "b105.cax8t_finite_coupled",
            "b106.cax8t_small_creep",
            "b107.cax8t_finite_creep",
            "b108.cax8t_small_plastic",
            "b109.cax8t_finite_plastic",
            "b110.cax8t_small_thermal",
            "b111.cax8t_finite_thermal",
            "b112.cax8t_small_friction",
            "b113.cax8t_finite_sliding",
            "b114.cax8t_recovery",
            "b120.cax8rt_small_probe",
            "b121.cax8rt_finite_probe",
            "b122.cax8rt_small_contact",
            "b123.cax8rt_finite_contact",
            "b124.cax8rt_small_coupled",
            "b125.cax8rt_finite_coupled",
            "b126.cax8rt_small_creep",
            "b127.cax8rt_finite_creep",
            "b128.cax8rt_small_plastic",
            "b129.cax8rt_finite_plastic",
            "b1210.cax8rt_small_thermal",
            "b1211.cax8rt_finite_thermal",
            "b1212.cax8rt_small_friction",
            "b1213.cax8rt_finite_sliding",
            "b1214.cax8rt_recovery",
            "b90.cax4rt_small_probe",
            "b91.cax4rt_finite_probe",
            "b92.cax4rt_rectangle_probe",
            "b93.cax4rt_small_coupled",
            "b94.cax4rt_finite_coupled",
            "b95.cax4rt_small_creep",
            "b96.cax4rt_finite_creep",
            "b97.cax4rt_small_plastic",
            "b98.cax4rt_finite_plastic",
            "b99.cax4rt_small_thermal",
            "b910.cax4rt_finite_thermal",
            "b911.cax4rt_small_contact",
            "b912.cax4rt_finite_contact",
            "b913.cax4rt_small_friction",
            "b914.cax4rt_finite_sliding",
            "b915.cax4rt_thermal_identification",
            "b70.rz_small_coupled",
            "b71.rz_finite_coupled",
            "b72.rz_small_creep",
            "b73.rz_finite_creep",
            "b74.rz_small_plastic",
            "b75.rz_finite_plastic",
            "b76.rz_small_thermal",
            "b77.rz_finite_thermal",
            "b78.rz_small_contact",
            "b79.rz_finite_contact",
            "b80.rz_small_friction",
            "b81.rz_finite_sliding",
            "input.v3",
            "build.reproducibility",
            "architecture.geometry_independence",
            "b3.hex8_thermomechanics",
            "b6.hex20_u2_t1",
            "h20.hex20_contact",
            "h20.abaqus_sts_identification",
            "h20.nonmatching_contact",
            "h20.node_to_face_three_way",
            "h20.sts_finite_sliding_planar",
            "h20.sts_finite_sliding_curved",
            "h20.sts_finite_sliding_partial_contact",
            "h20.sts_multicase",
            "h20.abaqus_curved_sts_identification",
            "h20.sts_curved",
            "h20.sts_nonmatching_friction_path",
            "h20.sts_curved_friction",
            "h20.sts_curved_friction_path",
            "b35.hex8_shared_material_interface",
            "b31.hex8_inelastic",
            "b32.hex8_finite_strain",
            "b33.hex8_contact",
            "b37.hex8_multi_contact_friction",
            "b38.hex8_sts_identification",
            "b39.hex8_sts_multicase",
            "b40.hex8_sts_friction",
            "b41.hex8_sts_finite_sliding",
            "b42.hex8_sts_friction_objectivity",
            "b43.hex8_sts_finite_strain",
            "b44.hex8_sts_nonmatching_partial",
            "b45.hex8_sts_nonplanar_smoothing",
            "b46.hex8_sts_release_recontact",
            "b47.hex8_sts_scale_parameter",
            "b48.hex8_sts_multi_contact_path",
            "b49.hex8_c3d8t_operator",
            "b50.hex8_c3d8t_capacity",
            "b51.hex8_c3d8t_finite_heat",
            "b52.hex8_c3d8t_thermal_contact",
            "b53.hex8_c3d8t_thermal_load",
            "b54.hex8_c3d8t_finite_thermal_load",
            "b55.hex8_c3d8t_transient_full_field",
            "b56.hex8_c3d8t_temperature_operator",
            "b57.hex8_c3d8t_temperature_capacity",
            "b58.hex8_c3d8t_multistep",
            "b59.hex8_c3d8t_multimaterial",
            "b510.hex8_c3d8t_small_j2",
            "b511.hex8_c3d8t_small_norton",
            "b512.hex8_c3d8t_small_coupled",
            "b513.hex8_c3d8t_small_noncoaxial",
            "b514.hex8_c3d8t_finite_elastic",
            "b515.hex8_c3d8t_finite_j2",
            "b516.hex8_c3d8t_finite_norton",
            "b517.hex8_c3d8t_finite_coupled",
            "b518.hex8_c3d8t_finite_noncoaxial",
            "b519.hex8_c3d8t_finite_selective",
            "b520.hex8_c3d8t_gap_conductance",
            "b521.hex8_c3d8t_thermal_contact_path",
            "b522.hex8_c3d8t_faceted_thermal_contact",
            "b523.hex8_c3d8t_integrated_path",
            "b524.hex8_load_parameter_independence",
            "b525.hex8_distorted_nearly_incompressible",
            "b526.hex8_nonlinear_transitions",
            "b527.hex8_engineering_fuel_clad",
            "b528.hex8_c3d8rt_operator",
            "b529.hex8_c3d8rt_capacity",
            "b530.hex8_c3d8rt_thermal_load",
            "b531.hex8_c3d8rt_transient_full_field",
            "b532.hex8_c3d8rt_finite_operator",
            "b533.hex8_c3d8rt_finite_hourglass",
            "b534.hex8_c3d8rt_finite_transient_full_field",
            "b535.hex8_c3d8rt_finite_j2",
            "b536.hex8_c3d8rt_finite_norton",
            "b537.hex8_c3d8rt_finite_coupled",
            "b538.hex8_c3d8rt_contact_cycle",
            "b539.hex8_c3d8rt_friction_reversal",
            "b540.hex8_c3d8rt_nonmatching_contact_cycle",
            "b541.hex8_c3d8rt_small_j2",
            "b542.hex8_c3d8rt_small_norton",
            "b543.hex8_c3d8rt_small_coupled",
            "b544.hex8_c3d8rt_distorted_bending",
            "b545.hex8_c3d8rt_finite_noncoaxial",
            "b546.m58_c3d8t_abaqus",
            "b547.m58_c3d8rt_abaqus",
            "b548.m58_c3d20t_timing",
            "b549.small_c3d20t_volume",
            "b550.small_c3d20t_contact",
            "b551.m58_c3d20t_sts_high_heat",
            "b552.hex8_c3d8t_sts_cross_face",
            "b553.m58_c3d8t_sts",
            "b554.m58_c3d8rt_sts",
            "b555.small_c3d20t_friction",
            "b556.m58_c3d20t_sts_friction",
            "b60.c3d8rt_finite_steady_bending",
            "b60.c3d8rt_finite_ramped_bending",
            "b60.c3d20t_finite_steady_bending",
            "b60.c3d20t_finite_ramped_bending",
            "b61.c3d20t_finite_inelastic_bending",
            "b61.c3d8rt_finite_inelastic_bending",
            "io.exodus",
            "m0.steady",
            "m1.contact",
            "m51.friction",
            "m52.large_sliding",
            "m53.time_integration",
            "m54.augmented_contact",
            "m55.shadow_state",
            "m56.iterative_solver",
            "m57.integrated",
            "m58.integrated_hex8",
            "m59.c3d8rt_parallel",
            "m21.transient",
            "m22.inelastic",
            "m23.pcmi",
            "m30.restart",
            "m31.loads",
            "m32.diagnostics",
            "m33.contact",
            "m34.parallel",
            "m40.foundation",
            "m41.finite_strain",
            "m42.follower_pressure",
            "m43.noncoaxial_finite_strain",
            "performance.m34",
            "performance.c3d8rt",
            "scope.boundary",
        };
        const std::set<std::string> manually_qualified_ids = {
            "b60.c3d20t_finite_steady_bending",
            "b60.c3d20t_finite_ramped_bending",
        };
        std::set<std::string> found_ids;
        std::size_t verified = 0;
        std::size_t qualified = 0;
        std::size_t measured = 0;
        std::size_t limitations = 0;
        std::size_t line_number = 1;
        while (std::getline(matrix, line)) {
            ++line_number;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) throw std::runtime_error("Blank matrix row at line " + std::to_string(line_number));
            const std::vector<std::string> fields = split(line, '\t');
            if (fields.size() != 6)
                throw std::runtime_error("Matrix row must have six fields at "
                                         "line " +
                                         std::to_string(line_number));
            for (const std::string& field : fields)
                if (field.empty())
                    throw std::runtime_error("Empty matrix field at line " + std::to_string(line_number));
            const std::string& id = fields[0];
            const std::string& status = fields[1];
            if (!found_ids.insert(id).second) throw std::runtime_error("Duplicate matrix id: " + id);
            if (required_ids.count(id) == 0) throw std::runtime_error("Unexpected matrix id: " + id);
            if (status == "verified" || status == "qualified") {
                if (status == "verified")
                    ++verified;
                else
                    ++qualified;
                if (fields[3] == "-") {
                    if (status != "qualified" || manually_qualified_ids.count(id) == 0 ||
                        fields[5].find("manual case is not registered in CTest") == std::string::npos)
                        throw std::runtime_error(id + " has no qualifying CTest");
                } else {
                    for (const std::string& test : split(fields[3], ';'))
                        if (registered_tests.count(test) == 0)
                            throw std::runtime_error(id + " names unknown CTest: " + test);
                }
            } else if (status == "measured") {
                ++measured;
                if (fields[3] != "-") throw std::runtime_error(id + " measured evidence must not masquerade as CTest");
            } else if (status == "limitation") {
                ++limitations;
                if (fields[3] != "-") throw std::runtime_error(id + " limitation must not masquerade as CTest");
            } else {
                throw std::runtime_error(id + " has invalid status: " + status);
            }
            check_evidence(repository, id, fields[4]);
        }
        if (found_ids != required_ids)
            throw std::runtime_error("Verification matrix is missing one or more required rows");
        if (verified != 160 || qualified != 8 || measured != 12 || limitations != 3)
            throw std::runtime_error("Verification matrix status counts differ from release schema");
        check_c3d8t_contract(argv[4], registered_tests);
        std::cout << "verification_matrix_rows=" << found_ids.size() << '\n'
                  << "verification_matrix_verified=" << verified << '\n'
                  << "verification_matrix_qualified=" << qualified << '\n'
                  << "verification_matrix_measured=" << measured << '\n'
                  << "verification_matrix_limitations=" << limitations << '\n'
                  << "[PASS] engineering verification matrix is complete\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] verification matrix audit raised: " << error.what() << '\n';
        return 1;
    }
}
