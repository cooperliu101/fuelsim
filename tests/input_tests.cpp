#include "fuelsim/case_input.hpp"
#include "fuelsim/input_file.hpp"

#include <cstdio>
#include <cstdint>
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

bool verify_m3_output_input(const std::string& path,
                            const std::string& contents) {
    {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        if (!output)
            return check(false, "could not create M3 input fixture");
        output << contents;
    }
    const fuelsim::FuelSimCaseDefinition definition =
        fuelsim::CaseInputReader::read(path);
    const int remove_status = std::remove(path.c_str());
    return check(
        remove_status == 0 &&
            definition.transient_execution.restart_file.find("restart.bin") !=
                std::string::npos &&
            definition.outputs.exodus_file.find("results.e") !=
                std::string::npos &&
            definition.outputs.exodus_interval == 3 &&
            definition.outputs.history_file.find("history.csv") !=
                std::string::npos &&
            definition.outputs.history_interval == 2 &&
            definition.outputs.progress_interval == 4 &&
            definition.outputs.checkpoint_file.find("checkpoint.bin") !=
                std::string::npos &&
            definition.outputs.checkpoint_interval == 5 &&
            definition.time_tables.size() == 1 &&
            definition.time_tables[0].value(1.0) == 0.5 &&
            definition.regions[0].spatial.heat_source_function == "power" &&
            definition.boundary_conditions.back().type ==
                fuelsim::BoundaryConditionType::convection &&
            definition.boundary_conditions.back().coefficient_function ==
                "power" &&
            definition.transient_execution.target_nonlinear_iterations == 6 &&
            definition.transient_execution.iteration_window == 2 &&
            definition.transient_execution.time_error_relative_tolerance ==
                2.0e-4 &&
            definition.transient_execution
                    .temperature_time_absolute_tolerance == 1.0e-3 &&
            definition.transient_execution
                    .displacement_time_absolute_tolerance == 1.0e-10 &&
            definition.transient_execution.time_error_safety_factor == 0.85 &&
            definition.solver.linear_solver == "gmres" &&
            definition.solver.preconditioner == "field_split" &&
            definition.solver.linear_relative_tolerance == 1.0e-7 &&
            definition.solver.maximum_linear_iterations == 700 &&
            definition.solver.backtracking_fallback &&
            definition.solver.field_residual_scaling &&
            definition.solver.residual_reduction_tolerance == 2.0e-6 &&
            definition.solver.temperature_residual_absolute_tolerance ==
                3.0e-8 &&
            definition.solver.mechanical_residual_absolute_tolerance ==
                4.0e-6 &&
            definition.regions[1]
                    .transient_material.creep
                    .coefficient_temperature_coefficient == 1.0e-8 &&
            definition.regions[1]
                    .transient_material.plasticity
                    .yield_stress_temperature_coefficient == -100.0,
        "restart, time functions, convection, solver and outputs are parsed");
}

bool fuzz_input_parser(const std::string& seed, const std::string& path) {
    std::uint64_t generator = 0x6a09e667f3bcc909ULL;
    for (std::size_t iteration = 0; iteration < 512; ++iteration) {
        std::string mutation = seed;
        generator = generator * 6364136223846793005ULL + 1ULL;
        const std::size_t edits = 1 + static_cast<std::size_t>(generator % 8);
        for (std::size_t edit = 0; edit < edits && !mutation.empty(); ++edit) {
            generator = generator * 6364136223846793005ULL + 1ULL;
            const std::size_t position =
                static_cast<std::size_t>(generator % mutation.size());
            generator = generator * 6364136223846793005ULL + 1ULL;
            mutation[position] = static_cast<char>(generator & 0x7fU);
        }
        if (iteration % 7 == 0 && !mutation.empty()) {
            generator = generator * 6364136223846793005ULL + 1ULL;
            mutation.resize(static_cast<std::size_t>(generator %
                                                     mutation.size()));
        }
        {
            std::ofstream output(path, std::ios::out | std::ios::trunc);
            if (!output)
                return check(false, "could not create parser fuzz fixture");
            output << mutation;
        }
        try {
            (void)fuelsim::InputParser::parse_file(path);
        } catch (const std::exception&) {
        } catch (...) {
            std::remove(path.c_str());
            return check(false,
                         "parser fuzz mutation raised a non-standard error");
        }
    }
    const int remove_status = std::remove(path.c_str());
    return check(remove_status == 0,
                 "deterministic parser fuzz mutations complete safely");
}

bool run_tests(const std::string& steady_path,
               const std::string& transient_path,
               const std::string& finite_strain_path,
               const std::string& scaled_displacement_path,
               const std::string& traction_path,
               const std::string& malformed_path) {
    const fuelsim::FuelSimCaseDefinition steady =
        fuelsim::CaseInputReader::read(steady_path);
    const fuelsim::FuelSimCaseDefinition transient =
        fuelsim::CaseInputReader::read(transient_path);
    const fuelsim::FuelSimCaseDefinition finite_strain =
        fuelsim::CaseInputReader::read(finite_strain_path);
    const fuelsim::FuelSimCaseDefinition scaled_displacement =
        fuelsim::CaseInputReader::read(scaled_displacement_path);
    const fuelsim::FuelSimCaseDefinition traction =
        fuelsim::CaseInputReader::read(traction_path);

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
                  steady.steady_execution.cutback_factor == 0.5 &&
                  steady.steady_execution.maximum_cutbacks == 12 &&
                  steady.steady_execution.minimum_load_increment == 1.0e-6 &&
                  steady.solver.maximum_iterations == 50 &&
                  steady.solver.linear_solver == "automatic" &&
                  steady.solver.preconditioner == "automatic" &&
                  steady.solver.linear_relative_tolerance == 1.0e-8 &&
                  steady.solver.maximum_linear_iterations == 500,
              "steady execution and solver fields are parsed") &&
        check(transient.problem == fuelsim::CaseProblem::transient,
              "transient input selects the physical transient problem") &&
        check(finite_strain.regions.size() == 2 &&
                  finite_strain.regions[0].spatial.strain_formulation ==
                      fuelsim::StrainFormulation::finite &&
                  finite_strain.regions[1].spatial.strain_formulation ==
                      fuelsim::StrainFormulation::finite,
              "finite strain is parsed independently for every region") &&
        check(transient.regions[0].transient_material.behavior ==
                      fuelsim::InelasticBehavior::elastic &&
                  transient.regions[1].transient_material.behavior ==
                      fuelsim::InelasticBehavior::norton_creep_j2_plasticity,
              "transient material behaviors are parsed") &&
        check(transient.transient_execution.end_time == 20.0 &&
                  transient.transient_execution.load_ramp_time == 20.0 &&
                  transient.solver.maximum_iterations == 80,
              "transient execution, ramp, and solver fields are parsed") &&
        check(
            scaled_displacement.regions.size() == 1 &&
                scaled_displacement.regions[0].spatial.block.empty() &&
                scaled_displacement.regions[0].spatial.block_id == 0 &&
                scaled_displacement.boundary_conditions.back().scale_with_load,
            "block ID and scaled displacement are parsed") &&
        check(traction.boundary_conditions.back().type ==
                      fuelsim::BoundaryConditionType::traction &&
                  traction.boundary_conditions.back().field ==
                      fuelsim::Field::axial_displacement &&
                  traction.boundary_conditions.back().scale_with_load,
              "scaled axial traction is parsed");

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

    std::string invalid_strain_case = read_text(finite_strain_path);
    const std::string finite_strain_key = "strain = finite";
    const std::size_t strain_position =
        invalid_strain_case.find(finite_strain_key);
    if (strain_position == std::string::npos)
        return check(false, "finite fixture has the expected strain key");
    invalid_strain_case.replace(strain_position, finite_strain_key.size(),
                                "strain = large");
    passed = expect_case_failure(malformed_path, invalid_strain_case,
                                 "unknown strain formulation 'large'") &&
             passed;

    std::string missing_strain_case = read_text(steady_path);
    const std::string small_strain_line = "    strain = small\n";
    const std::size_t small_strain_position =
        missing_strain_case.find(small_strain_line);
    if (small_strain_position == std::string::npos)
        return check(false, "steady fixture has the expected strain key");
    missing_strain_case.erase(small_strain_position,
                              small_strain_line.size());
    passed = expect_case_failure(malformed_path, missing_strain_case,
                                 "missing required key 'strain'") &&
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

    std::string m3_case = read_text(transient_path);
    m3_case.insert(0, "[TimeFunctions]\n  [power]\n    type = "
                      "piecewise_linear\n    times = 0 2 5\n    values = 0 1 "
                      "0.4\n  []\n[]\n\n");
    const std::string heat_source = "volumetric_heat_source = 2e8";
    const std::size_t heat_source_position = m3_case.find(heat_source);
    if (heat_source_position == std::string::npos)
        return check(false, "transient fixture has the expected heat source");
    m3_case.insert(heat_source_position + heat_source.size(),
                   "\n    heat_source_function = power");
    const std::string executioner_start = "\n[Executioner]";
    const std::size_t executioner_start_position =
        m3_case.find(executioner_start);
    if (executioner_start_position == std::string::npos)
        return check(false, "transient fixture has an executioner section");
    const std::size_t boundary_close =
        m3_case.rfind("[]", executioner_start_position);
    if (boundary_close == std::string::npos)
        return check(false, "transient fixture closes boundary conditions");
    m3_case.insert(
        boundary_close,
        "  [coolant]\n    type = convection\n    boundary = clad_outer\n"
        "    heat_transfer_coefficient = 1000\n"
        "    ambient_temperature = 600\n"
        "    coefficient_function = power\n  []\n");
    const std::string executioner_type = "type = transient";
    const std::size_t executioner_position = m3_case.find(executioner_type);
    if (executioner_position == std::string::npos)
        return check(false, "transient fixture has an executioner type");
    m3_case.insert(executioner_position + executioner_type.size(),
                   "\n  restart = restart.bin"
                   "\n  target_nonlinear_iterations = 6"
                   "\n  iteration_window = 2"
                   "\n  time_error_relative_tolerance = 2e-4"
                   "\n  temperature_time_absolute_tolerance = 1e-3"
                   "\n  displacement_time_absolute_tolerance = 1e-10"
                   "\n  time_error_safety_factor = 0.85");
    const std::string solver_start = "[Solver]";
    const std::size_t solver_position = m3_case.find(solver_start);
    if (solver_position == std::string::npos)
        return check(false, "transient fixture has a solver section");
    m3_case.insert(solver_position + solver_start.size(),
                   "\n  linear_solver = gmres"
                   "\n  preconditioner = field_split"
                   "\n  linear_relative_tolerance = 1e-7"
                   "\n  maximum_linear_iterations = 700"
                   "\n  backtracking_fallback = true"
                   "\n  field_residual_scaling = true"
                   "\n  residual_reduction_tolerance = 2e-6"
                   "\n  temperature_residual_absolute_tolerance = 3e-8"
                   "\n  mechanical_residual_absolute_tolerance = 4e-6");
    const std::string creep_exponent = "creep_exponent = 3";
    const std::size_t creep_position = m3_case.find(creep_exponent);
    if (creep_position == std::string::npos)
        return check(false, "transient fixture has coupled creep properties");
    m3_case.insert(creep_position + creep_exponent.size(),
                   "\n    creep_coefficient_temperature_coefficient = 1e-8"
                   "\n    creep_reference_stress_temperature_coefficient = 10"
                   "\n    creep_exponent_temperature_coefficient = 1e-4"
                   "\n    yield_stress_temperature_coefficient = -100"
                   "\n    hardening_temperature_coefficient = -1000");
    const std::size_t output_position = m3_case.find(console);
    if (output_position == std::string::npos)
        return check(false, "transient fixture has the expected console key");
    m3_case.insert(output_position + console.size(),
                   "\n  exodus = results.e\n  exodus_interval = 3"
                   "\n  history = history.csv\n  history_interval = 2"
                   "\n  progress_interval = 4"
                   "\n  checkpoint = checkpoint.bin"
                   "\n  checkpoint_interval = 5");
    passed = verify_m3_output_input(malformed_path, m3_case) && passed;

    std::string invalid_preconditioner = m3_case;
    const std::string valid_preconditioner = "preconditioner = field_split";
    const std::size_t preconditioner_position =
        invalid_preconditioner.find(valid_preconditioner);
    invalid_preconditioner.replace(preconditioner_position,
                                   valid_preconditioner.size(),
                                   "preconditioner = magic");
    passed = expect_case_failure(malformed_path, invalid_preconditioner,
                                 "preconditioner must be") &&
             passed;

    std::string unknown_function = read_text(transient_path);
    const std::size_t unknown_heat_position =
        unknown_function.find(heat_source);
    unknown_function.insert(unknown_heat_position + heat_source.size(),
                            "\n    heat_source_function = missing");
    passed = expect_case_failure(malformed_path, unknown_function,
                                 "unknown time function 'missing'") &&
             passed;

    std::string invalid_table = m3_case;
    const std::string valid_times = "times = 0 2 5";
    const std::size_t valid_times_position = invalid_table.find(valid_times);
    if (valid_times_position == std::string::npos)
        return check(false, "M3 fixture has time-table nodes");
    invalid_table.replace(valid_times_position, valid_times.size(),
                          "times = 0 0 5");
    passed = expect_case_failure(malformed_path, invalid_table,
                                 "strictly increasing") &&
             passed;

    std::string invalid_window = m3_case;
    const std::string valid_window = "iteration_window = 2";
    const std::size_t valid_window_position = invalid_window.find(valid_window);
    invalid_window.replace(valid_window_position, valid_window.size(),
                           "iteration_window = 6");
    passed = expect_case_failure(malformed_path, invalid_window,
                                 "must be smaller") &&
             passed;

    std::string orphan_interval = read_text(transient_path);
    const std::size_t orphan_output = orphan_interval.find(console);
    orphan_interval.insert(orphan_output + console.size(),
                           "\n  checkpoint_interval = 2");
    passed = expect_case_failure(malformed_path, orphan_interval,
                                 "checkpoint_interval requires checkpoint") &&
             passed;

    std::string steady_checkpoint = read_text(steady_path);
    const std::size_t steady_output = steady_checkpoint.find(console);
    steady_checkpoint.insert(steady_output + console.size(),
                             "\n  checkpoint = checkpoint.bin");
    passed = expect_case_failure(malformed_path, steady_checkpoint,
                                 "only valid for transient cases") &&
             passed;
    std::string mesh_overwrite = read_text(steady_path);
    const std::string mesh_file =
        "file = ../moose/m1_fuel_cladding_gap_rz_mesh.e";
    const std::size_t mesh_output = mesh_overwrite.find(console);
    if (mesh_overwrite.find(mesh_file) == std::string::npos ||
        mesh_output == std::string::npos)
        return check(false, "steady fixture has expected mesh and output");
    mesh_overwrite.insert(mesh_output + console.size(),
                          "\n  csv = ../moose/"
                          "m1_fuel_cladding_gap_rz_mesh.e");
    passed = expect_case_failure(malformed_path, mesh_overwrite,
                                 "must not overwrite the input mesh") &&
             passed;

    std::string output_collision = read_text(transient_path);
    const std::size_t collision_output = output_collision.find(console);
    output_collision.insert(collision_output + console.size(),
                            "\n  csv = collision.dat"
                            "\n  checkpoint = collision.dat");
    passed = expect_case_failure(malformed_path, output_collision,
                                 "paths must differ") &&
             passed;
    passed = fuzz_input_parser(read_text(transient_path), malformed_path) &&
             passed;
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "Usage: fuelsim_input_tests <steady.fsi> "
                     "<transient.fsi> <finite-strain.fsi> "
                     "<scaled-displacement.fsi> "
                     "<traction.fsi> <malformed.fsi>\n";
        return 2;
    }

    try {
        if (!run_tests(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]))
            return 1;
        std::cout << "[PASS] fuelsim strict input-card tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] input tests raised: " << error.what() << '\n';
        return 1;
    }
}
