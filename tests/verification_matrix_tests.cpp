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
} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_verification_matrix_tests "
                     "<matrix.tsv> <repository> <registered-tests>\n";
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
            "scope.boundary",
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
                if (fields[3] == "-") throw std::runtime_error(id + " has no qualifying CTest");
                for (const std::string& test : split(fields[3], ';'))
                    if (registered_tests.count(test) == 0)
                        throw std::runtime_error(id + " names unknown CTest: " + test);
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
        if (verified != 45 || qualified != 9 || measured != 3 || limitations != 1)
            throw std::runtime_error("Verification matrix status counts differ from release schema");
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
