#include "fuelsim/case_input.hpp"
#include "fuelsim/input_file.hpp"

#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool expect_parse_failure(const std::string& path, const std::string& contents,
                          const std::string& expected_message) {
    {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        if (!output)
            return check(false, "could not create malformed input fixture");
        output << contents;
    }

    bool failed_as_expected = false;
    try {
        (void)fuelsim::InputParser::parse_file(path);
    } catch (const std::invalid_argument& error) {
        failed_as_expected = std::string(error.what()).find(expected_message) !=
                             std::string::npos;
    }
    const int remove_status = std::remove(path.c_str());
    return check(remove_status == 0 && failed_as_expected,
                 "malformed input reports '" + expected_message + "'");
}

std::string read_text(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read input fixture '" + path + "'");
    std::ostringstream contents;
    contents << input.rdbuf();
    return contents.str();
}

bool expect_case_failure(const std::string& path, const std::string& contents,
                         const std::string& expected_message) {
    {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        if (!output)
            return check(false, "could not create invalid case fixture");
        output << contents;
    }

    bool failed_as_expected = false;
    try {
        (void)fuelsim::CaseInputReader::read(path);
    } catch (const std::invalid_argument& error) {
        failed_as_expected = std::string(error.what()).find(expected_message) !=
                             std::string::npos;
    }
    const int remove_status = std::remove(path.c_str());
    return check(remove_status == 0 && failed_as_expected,
                 "invalid case reports '" + expected_message + "'");
}

bool run_tests(const std::string& steady_path,
               const std::string& transient_path,
               const std::string& malformed_path) {
    const fuelsim::FuelSimCaseDefinition steady =
        fuelsim::CaseInputReader::read(steady_path);
    const fuelsim::FuelSimCaseDefinition transient =
        fuelsim::CaseInputReader::read(transient_path);

    bool passed =
        check(steady.version == 1 &&
                  steady.problem == fuelsim::CaseProblem::steady,
              "steady input selects the physical steady problem") &&
        check(steady.regions.size() == 2 &&
                  steady.regions[0].spatial.block == "fuel" &&
                  steady.regions[1].spatial.block == "clad",
              "steady input preserves arbitrary Exodus regions") &&
        check(steady.contacts.size() == 1 &&
                  steady.contacts[0].primary == "clad_left" &&
                  steady.contacts[0].secondary == "fuel_right" &&
                  steady.contacts[0].thermal && steady.contacts[0].mechanical,
              "contact is defined only by primary and secondary side sets") &&
        check(steady.steady_execution.load_steps == 20 &&
                  steady.solver.maximum_iterations == 50,
              "steady execution and solver fields are parsed") &&
        check(transient.problem == fuelsim::CaseProblem::transient,
              "transient input selects the physical transient problem") &&
        check(transient.regions[0].transient_material.behavior ==
                      fuelsim::InelasticBehavior::elastic &&
                  transient.regions[1].transient_material.behavior ==
                      fuelsim::InelasticBehavior::norton_creep_j2_plasticity,
              "transient material behaviors are parsed") &&
        check(transient.transient_execution.end_time == 20.0 &&
                  transient.transient_execution.heat_source_ramp_time == 20.0 &&
                  transient.solver.maximum_iterations == 80,
              "transient execution, ramp, and solver fields are parsed");

    passed = expect_parse_failure(malformed_path,
                                  "[Case]\n  version = 1\n  version = 1\n[]\n",
                                  "duplicate key 'version'") &&
             passed;
    passed = expect_parse_failure(malformed_path, "[Case]\n  version = 1\n",
                                  "missing its closing []") &&
             passed;

    std::string unknown_key_case = read_text(steady_path);
    const std::string console = "console = true";
    const std::size_t console_position = unknown_key_case.find(console);
    if (console_position == std::string::npos)
        return check(false, "steady fixture has the expected console key");
    unknown_key_case.replace(console_position, console.size(),
                             "mystery = true");
    passed = expect_case_failure(malformed_path, unknown_key_case,
                                 "unknown key 'mystery'") &&
             passed;

    std::string invalid_material_case = read_text(transient_path);
    const std::string elastic_model = "inelastic_model = elastic";
    const std::size_t model_position =
        invalid_material_case.find(elastic_model);
    if (model_position == std::string::npos)
        return check(false, "transient fixture has the elastic fuel model");
    invalid_material_case.insert(model_position + elastic_model.size(),
                                 "\n    yield_stress = 1e8");
    passed = expect_case_failure(malformed_path, invalid_material_case,
                                 "not valid for inelastic_model='elastic'") &&
             passed;

    std::string block_contact_case = read_text(steady_path);
    const std::string primary = "primary = clad_left";
    const std::size_t primary_position = block_contact_case.find(primary);
    if (primary_position == std::string::npos)
        return check(false, "steady fixture has the expected primary key");
    block_contact_case.insert(primary_position, "primary_block = clad\n    ");
    passed = expect_case_failure(malformed_path, block_contact_case,
                                 "unknown key 'primary_block'") &&
             passed;
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: fuelsim_input_tests <steady.fsi> "
                     "<transient.fsi> <malformed.fsi>\n";
        return 2;
    }

    try {
        if (!run_tests(argv[1], argv[2], argv[3]))
            return 1;
        std::cout << "[PASS] fuelsim strict input-card tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] input tests raised: " << error.what() << '\n';
        return 1;
    }
}
