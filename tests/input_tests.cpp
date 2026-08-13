#include "fuelsim/case_input.hpp"
#include <cstdint>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
namespace {
bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}
void registered_test_thermal(const fuelsim::ThermoelasticFunctionInput& input, fuelsim::ThermalPropertyOutput& output) {
    output.conductivity = input.parameters->value("inverse_coefficient") / input.temperature +
                          input.parameters->value("constant_coefficient");
    output.density = input.parameters->value("density");
    output.specific_heat = input.parameters->value("specific_heat");
}
bool expect_parse_failure(const std::string& path, const std::string& contents, const std::string& expected_message) {
    {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create malformed input fixture");
        output << contents;
    }
    bool failed_as_expected = false;
    try {
        (void)fuelsim::parse_input_file(path);
    } catch (const std::invalid_argument& error) {
        failed_as_expected = std::string(error.what()).find(expected_message) != std::string::npos;
    }
    const int remove_status = std::remove(path.c_str());
    return check(remove_status == 0 && failed_as_expected, "malformed input reports '" + expected_message + "'");
}
std::string read_text(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read input fixture '" + path + "'");
    std::ostringstream contents;
    contents << input.rdbuf();
    return contents.str();
}
bool expect_case_failure(const std::string& path, const std::string& contents, const std::string& expected_message) {
    {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create invalid case fixture");
        output << contents;
    }
    bool failed_as_expected = false;
    try {
        (void)fuelsim::read_case_input(path);
    } catch (const std::invalid_argument& error) {
        failed_as_expected = std::string(error.what()).find(expected_message) != std::string::npos;
    }
    const int remove_status = std::remove(path.c_str());
    return check(remove_status == 0 && failed_as_expected, "invalid case reports '" + expected_message + "'");
}
bool expect_registered_case_failure(const std::string& path, const std::string& contents,
    const fuelsim::MaterialFunctionRegistry& registry, const std::string& expected_message) {
    {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create invalid registered-material fixture");
        output << contents;
    }
    bool failed_as_expected = false;
    try {
        (void)fuelsim::read_case_input(path, registry);
    } catch (const std::invalid_argument& error) {
        failed_as_expected = std::string(error.what()).find(expected_message) != std::string::npos;
    }
    const int remove_status = std::remove(path.c_str());
    return check(remove_status == 0 && failed_as_expected,
        "invalid registered-material case reports '" + expected_message + "'");
}
bool verify_m3_output_input(const std::string& path, const std::string& contents) {
    {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create M3 input fixture");
        output << contents;
    }
    const fuelsim::FuelSimCaseDefinition definition = fuelsim::read_case_input(path);
    const int remove_status = std::remove(path.c_str());
    return check(remove_status == 0 && definition.restart_file.find("restart.bin") != std::string::npos &&
                     definition.outputs.exodus_file.find("results.e") != std::string::npos &&
                     definition.outputs.exodus_interval == 3 &&
                     definition.outputs.history_file.find("history.csv") != std::string::npos &&
                     definition.outputs.history_interval == 2 && definition.outputs.progress_interval == 4 &&
                     definition.outputs.checkpoint_file.find("checkpoint.bin") != std::string::npos &&
                     definition.outputs.checkpoint_interval == 5 && definition.spatial.time_tables.size() == 1 &&
                     definition.spatial.time_tables[0].value(1.0) == 0.5 &&
                     definition.spatial.regions[0].heat_source_function == "power" &&
                     definition.spatial.boundary_conditions.back().type == fuelsim::BoundaryConditionType::convection &&
                     definition.spatial.boundary_conditions.back().coefficient_function == "power" &&
                     definition.transient_execution.target_nonlinear_iterations == 6 &&
                     definition.transient_execution.iteration_window == 2 &&
                     definition.transient_execution.time_error_relative_tolerance == 2.0e-4 &&
                     definition.transient_execution.temperature_time_absolute_tolerance == 1.0e-3 &&
                     definition.transient_execution.displacement_time_absolute_tolerance == 1.0e-10 &&
                     definition.transient_execution.strain_history_time_absolute_tolerance == 2.0e-9 &&
                     definition.transient_execution.stress_history_time_absolute_tolerance == 5.0 &&
                     definition.transient_execution.time_error_safety_factor == 0.85 &&
                     definition.solver.linear_solver == fuelsim::SolverOptions::LinearSolver::gmres &&
                     definition.solver.preconditioner == fuelsim::SolverOptions::Preconditioner::field_split &&
                     definition.solver.linear_relative_tolerance == 1.0e-7 &&
                     definition.solver.maximum_linear_iterations == 700 && definition.solver.backtracking_fallback &&
                     definition.solver.field_residual_scaling &&
                     definition.solver.residual_reduction_tolerance == 2.0e-6 &&
                     definition.solver.temperature_residual_absolute_tolerance == 3.0e-8 &&
                     definition.solver.mechanical_residual_absolute_tolerance == 4.0e-6 &&
                     definition.spatial.regions[1].material.functions->creep.parameters.value(
                         "coefficient_temperature_coefficient") == 1.0e-8 &&
                     definition.spatial.regions[1].material.functions->plasticity.parameters.value(
                         "yield_stress_temperature_coefficient") == -100.0,
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
            const std::size_t position = static_cast<std::size_t>(generator % mutation.size());
            generator = generator * 6364136223846793005ULL + 1ULL;
            mutation[position] = static_cast<char>(generator & 0x7fU);
        }
        if (iteration % 7 == 0 && !mutation.empty()) {
            generator = generator * 6364136223846793005ULL + 1ULL;
            mutation.resize(static_cast<std::size_t>(generator % mutation.size()));
        }
        {
            std::ofstream output(path, std::ios::out | std::ios::trunc);
            if (!output) return check(false, "could not create parser fuzz fixture");
            output << mutation;
        }
        try {
            (void)fuelsim::parse_input_file(path);
        } catch (const std::exception&) {
        } catch (...) {
            std::remove(path.c_str());
            return check(false, "parser fuzz mutation raised a non-standard error");
        }
    }
    const int remove_status = std::remove(path.c_str());
    return check(remove_status == 0, "deterministic parser fuzz mutations complete safely");
}
bool run_tests(const std::string& steady_path, const std::string& transient_path, const std::string& finite_strain_path,
    const std::string& scaled_displacement_path, const std::string& traction_path, const std::string& malformed_path) {
    const fuelsim::FuelSimCaseDefinition steady = fuelsim::read_case_input(steady_path);
    const fuelsim::FuelSimCaseDefinition transient = fuelsim::read_case_input(transient_path);
    const fuelsim::FuelSimCaseDefinition finite_strain = fuelsim::read_case_input(finite_strain_path);
    const fuelsim::FuelSimCaseDefinition scaled_displacement = fuelsim::read_case_input(scaled_displacement_path);
    const fuelsim::FuelSimCaseDefinition traction = fuelsim::read_case_input(traction_path);
    bool passed =
        check(steady.version == 3 && steady.problem == fuelsim::CaseProblem::steady &&
                  steady.geometry == fuelsim::CaseGeometry::axisymmetric_rz,
            "steady input selects the physical steady problem") &&
        check(steady.spatial.regions.size() == 2 && steady.spatial.regions[0].block == "fuel" &&
                  steady.spatial.regions[1].block == "clad",
            "steady input preserves arbitrary Exodus regions") &&
        check(steady.spatial.contacts.size() == 1 && steady.spatial.contacts[0].primary == "clad_left" &&
                  steady.spatial.contacts[0].secondary == "fuel_right" && steady.spatial.contacts[0].thermal &&
                  steady.spatial.contacts[0].mechanical && steady.spatial.contacts[0].friction_coefficient == 0.0,
            "contact is defined only by primary and secondary side sets") &&
        check(steady.steady_execution.load_steps == 20 && steady.steady_execution.cutback_factor == 0.5 &&
                  steady.steady_execution.maximum_cutbacks_per_step == 12 &&
                  steady.steady_execution.minimum_load_increment == 1.0e-6 && steady.solver.maximum_iterations == 50 &&
                  steady.solver.linear_solver == fuelsim::SolverOptions::LinearSolver::automatic &&
                  steady.solver.preconditioner == fuelsim::SolverOptions::Preconditioner::automatic &&
                  steady.solver.linear_relative_tolerance == 1.0e-8 && steady.solver.maximum_linear_iterations == 500,
            "steady execution and solver fields are parsed") &&
        check(transient.problem == fuelsim::CaseProblem::transient,
            "transient input selects the physical transient problem") &&
        check(finite_strain.spatial.regions.size() == 2 &&
                  finite_strain.spatial.regions[0].strain_formulation == fuelsim::StrainFormulation::finite &&
                  finite_strain.spatial.regions[1].strain_formulation == fuelsim::StrainFormulation::finite,
            "finite strain is parsed independently for every region") &&
        check(!transient.spatial.regions[0].material.functions->has_creep() &&
                  !transient.spatial.regions[0].material.functions->has_plasticity() &&
                  transient.spatial.regions[1].material.functions->has_creep() &&
                  transient.spatial.regions[1].material.functions->has_plasticity(),
            "transient material behaviors are parsed") &&
        check(transient.transient_execution.end_time == 20.0 && transient.transient_execution.load_ramp_time == 20.0 &&
                  transient.transient_execution.time_error_relative_tolerance == 0.0 &&
                  transient.solver.maximum_iterations == 80,
            "transient execution keeps time-error control opt-in and "
            "parses ramp and solver fields") &&
        check(scaled_displacement.spatial.regions.size() == 1 && scaled_displacement.spatial.regions[0].block.empty() &&
                  scaled_displacement.spatial.regions[0].block_id == 0 &&
                  scaled_displacement.spatial.boundary_conditions.back().scale_with_load,
            "block ID and scaled displacement are parsed") &&
        check(traction.spatial.boundary_conditions.back().type == fuelsim::BoundaryConditionType::traction &&
                  traction.spatial.boundary_conditions.back().field == fuelsim::Field::axial_displacement &&
                  traction.spatial.boundary_conditions.back().scale_with_load &&
                  !traction.spatial.boundary_conditions.back().use_displaced_geometry,
            "scaled axial traction is parsed");
    passed =
        expect_parse_failure(malformed_path, "[Case]\n  version = 1\n  version = 1\n[]\n", "duplicate key 'version'") &&
        passed;
    passed = expect_parse_failure(malformed_path, "[Case]\n  version = 1\n", "missing its closing []") && passed;
    fuelsim::MaterialFunctionRegistry custom_registry = fuelsim::make_builtin_material_function_registry();
    custom_registry.add_thermal("registered_test_thermal",
        {{"inverse_coefficient", "W/m"}, {"constant_coefficient", "W/(m*K)"}, {"density", "kg/m^3"},
            {"specific_heat", "J/(kg*K)"}},
        &registered_test_thermal);
    std::string registered_case = read_text(steady_path);
    const std::string builtin_thermal = "function = inverse_temperature_thermophysical";
    const std::size_t builtin_thermal_position = registered_case.find(builtin_thermal);
    if (builtin_thermal_position == std::string::npos)
        return check(false, "steady fixture has an inverse-temperature thermal function");
    registered_case.replace(builtin_thermal_position, builtin_thermal.size(), "function = registered_test_thermal");
    const std::string inverse_name = "conductivity_inverse_temperature";
    const std::string constant_name = "conductivity_constant";
    registered_case.replace(registered_case.find(inverse_name), inverse_name.size(), "inverse_coefficient");
    registered_case.replace(registered_case.find(constant_name), constant_name.size(), "constant_coefficient");
    {
        std::ofstream output(malformed_path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create registered material input fixture");
        output << registered_case;
    }
    const fuelsim::FuelSimCaseDefinition registered = fuelsim::read_case_input(malformed_path, custom_registry);
    passed = check(registered.spatial.regions[0].material.functions->thermal.name == "registered_test_thermal" &&
                       registered.spatial.regions[0].material.functions->thermal.parameters.value(
                           "inverse_coefficient") == 3824.0,
                 "a custom registered function receives strictly named input parameters") &&
             passed;
    if (std::remove(malformed_path.c_str()) != 0)
        return check(false, "could not remove registered material input fixture");
    std::string missing_material_parameter = registered_case;
    const std::string density_line = "      density = 1\n";
    const std::size_t density_position = missing_material_parameter.find(density_line);
    if (density_position == std::string::npos)
        return check(false, "registered material fixture has a named density parameter");
    missing_material_parameter.erase(density_position, density_line.size());
    passed = expect_registered_case_failure(
                 malformed_path, missing_material_parameter, custom_registry, "missing required key 'density'") &&
             passed;
    std::string unknown_key_case = read_text(steady_path);
    const std::string console = "console = true";
    const std::size_t console_position = unknown_key_case.find(console);
    if (console_position == std::string::npos) return check(false, "steady fixture has the expected console key");
    unknown_key_case.replace(console_position, console.size(), "mystery = true");
    passed = expect_case_failure(malformed_path, unknown_key_case, "unknown key 'mystery'") && passed;
    std::string friction_case = read_text(steady_path);
    const std::string contact_penalty = "penalty = 1e14";
    const std::size_t contact_penalty_position = friction_case.find(contact_penalty);
    if (contact_penalty_position == std::string::npos)
        return check(false, "steady fixture has the contact penalty key");
    friction_case.insert(contact_penalty_position + contact_penalty.size(), "\n      mu = 0.25");
    {
        std::ofstream output(malformed_path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create friction input fixture");
        output << friction_case;
    }
    const fuelsim::FuelSimCaseDefinition friction = fuelsim::read_case_input(malformed_path);
    passed = check(friction.spatial.contacts.size() == 1 && friction.spatial.contacts[0].friction_coefficient == 0.25,
                 "optional contact mu is parsed as the Coulomb friction "
                 "coefficient") &&
             passed;
    std::string automatic_penalty_case = read_text(steady_path);
    const std::string contact_penalty_line = "      penalty = 1e14\n";
    const std::size_t automatic_penalty_position = automatic_penalty_case.find(contact_penalty_line);
    if (automatic_penalty_position == std::string::npos)
        return check(false, "steady fixture has the contact penalty line");
    automatic_penalty_case.erase(automatic_penalty_position, contact_penalty_line.size());
    {
        std::ofstream output(malformed_path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create automatic-penalty input fixture");
        output << automatic_penalty_case;
    }
    const fuelsim::FuelSimCaseDefinition automatic_penalty = fuelsim::read_case_input(malformed_path);
    passed = check(automatic_penalty.spatial.contacts[0].automatic_penalty &&
                       automatic_penalty.spatial.contacts[0].penalty_factor == 1.0,
                 "omitting penalty selects the documented automatic "
                 "contact factor") &&
             passed;
    std::string augmented_case = read_text(steady_path);
    const std::string penalty_formulation = "formulation = penalty";
    const std::size_t formulation_position = augmented_case.find(penalty_formulation);
    if (formulation_position == std::string::npos) return check(false, "steady fixture has penalty formulation");
    augmented_case.replace(formulation_position, penalty_formulation.size(), "formulation = augmented_lagrangian");
    const std::size_t augmented_penalty_position = augmented_case.find(contact_penalty);
    augmented_case.insert(augmented_penalty_position + contact_penalty.size(),
        "\n      penetration_tolerance = 2e-9"
        "\n      maximum_augmented_iterations = 15");
    {
        std::ofstream output(malformed_path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create augmented input fixture");
        output << augmented_case;
    }
    const fuelsim::FuelSimCaseDefinition augmented = fuelsim::read_case_input(malformed_path);
    passed = check(augmented.spatial.contacts[0].mechanical_formulation ==
                           fuelsim::MechanicalContactFormulation::augmented_lagrangian &&
                       augmented.spatial.contacts[0].penetration_tolerance == 2.0e-9 &&
                       augmented.spatial.contacts[0].maximum_augmented_iterations == 15,
                 "augmented contact tolerance and iteration limit are "
                 "parsed") &&
             passed;
    std::string ambiguous_penalty_case = read_text(steady_path);
    ambiguous_penalty_case.insert(
        ambiguous_penalty_case.find(contact_penalty) + contact_penalty.size(), "\n      penalty_factor = 10");
    passed = expect_case_failure(malformed_path, ambiguous_penalty_case, "mutually exclusive") && passed;
    std::string negative_friction_case = read_text(steady_path);
    negative_friction_case.insert(
        negative_friction_case.find(contact_penalty) + contact_penalty.size(), "\n      mu = -0.1");
    {
        std::ofstream output(malformed_path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create trusted-value fixture");
        output << negative_friction_case;
    }
    const fuelsim::FuelSimCaseDefinition trusted_values = fuelsim::read_case_input(malformed_path);
    passed = check(trusted_values.spatial.contacts[0].friction_coefficient == -0.1,
                 "input preserves user-provided physical values") &&
             passed;
    if (std::remove(malformed_path.c_str()) != 0) return check(false, "could not remove trusted-value fixture");
    std::string invalid_material_case = read_text(transient_path);
    const std::string elastic_model = "function = constant_isotropic";
    const std::size_t model_position = invalid_material_case.find(elastic_model);
    if (model_position == std::string::npos) return check(false, "transient fixture has the fuel elasticity function");
    invalid_material_case.insert(model_position + elastic_model.size(), "\n      yield_stress = 1e8");
    passed = expect_case_failure(malformed_path, invalid_material_case, "unknown key 'yield_stress'") && passed;
    std::string invalid_strain_case = read_text(finite_strain_path);
    const std::string finite_strain_key = "strain = finite";
    const std::size_t strain_position = invalid_strain_case.find(finite_strain_key);
    if (strain_position == std::string::npos) return check(false, "finite fixture has the expected strain key");
    invalid_strain_case.replace(strain_position, finite_strain_key.size(), "strain = large");
    passed = expect_case_failure(malformed_path, invalid_strain_case, "unknown strain formulation 'large'") && passed;
    std::string missing_strain_case = read_text(steady_path);
    const std::string small_strain_line = "    strain = small\n";
    const std::size_t small_strain_position = missing_strain_case.find(small_strain_line);
    if (small_strain_position == std::string::npos) return check(false, "steady fixture has the expected strain key");
    missing_strain_case.erase(small_strain_position, small_strain_line.size());
    passed = expect_case_failure(malformed_path, missing_strain_case, "missing required key 'strain'") && passed;
    std::string block_contact_case = read_text(steady_path);
    const std::string primary = "primary = clad_left";
    const std::size_t primary_position = block_contact_case.find(primary);
    if (primary_position == std::string::npos) return check(false, "steady fixture has the expected primary key");
    block_contact_case.insert(primary_position, "primary_block = clad\n    ");
    passed = expect_case_failure(malformed_path, block_contact_case, "unknown key 'primary_block'") && passed;
    std::string m3_case = read_text(transient_path);
    m3_case.insert(0, "[TimeFunctions]\n  [power]\n    type = "
                      "piecewise_linear\n    times = 0 2 5\n    values = 0 1 "
                      "0.4\n  []\n[]\n\n");
    const std::string heat_source = "volumetric_heat_source = 2e8";
    const std::size_t heat_source_position = m3_case.find(heat_source);
    if (heat_source_position == std::string::npos)
        return check(false, "transient fixture has the expected heat source");
    m3_case.insert(heat_source_position + heat_source.size(), "\n    heat_source_function = power");
    const std::string executioner_start = "\n[Executioner]";
    const std::size_t executioner_start_position = m3_case.find(executioner_start);
    if (executioner_start_position == std::string::npos)
        return check(false, "transient fixture has an executioner section");
    const std::size_t boundary_close = m3_case.rfind("[]", executioner_start_position);
    if (boundary_close == std::string::npos) return check(false, "transient fixture closes boundary conditions");
    m3_case.insert(boundary_close, "  [coolant]\n    type = convection\n    boundary = clad_outer\n"
                                   "    heat_transfer_coefficient = 1000\n"
                                   "    ambient_temperature = 600\n"
                                   "    coefficient_function = power\n  []\n");
    const std::string executioner_type = "type = transient";
    const std::size_t executioner_position = m3_case.find(executioner_type);
    if (executioner_position == std::string::npos) return check(false, "transient fixture has an executioner type");
    m3_case.insert(executioner_position + executioner_type.size(), "\n  restart = restart.bin"
                                                                   "\n  target_nonlinear_iterations = 6"
                                                                   "\n  iteration_window = 2"
                                                                   "\n  time_error_relative_tolerance = 2e-4"
                                                                   "\n  temperature_time_absolute_tolerance = 1e-3"
                                                                   "\n  displacement_time_absolute_tolerance = 1e-10"
                                                                   "\n  strain_history_time_absolute_tolerance = 2e-9"
                                                                   "\n  stress_history_time_absolute_tolerance = 5"
                                                                   "\n  time_error_safety_factor = 0.85");
    const std::string solver_start = "[Solver]";
    const std::size_t solver_position = m3_case.find(solver_start);
    if (solver_position == std::string::npos) return check(false, "transient fixture has a solver section");
    m3_case.insert(solver_position + solver_start.size(), "\n  linear_solver = gmres"
                                                          "\n  preconditioner = field_split"
                                                          "\n  linear_relative_tolerance = 1e-7"
                                                          "\n  maximum_linear_iterations = 700"
                                                          "\n  backtracking_fallback = true"
                                                          "\n  field_residual_scaling = true"
                                                          "\n  residual_reduction_tolerance = 2e-6"
                                                          "\n  temperature_residual_absolute_tolerance = 3e-8"
                                                          "\n  mechanical_residual_absolute_tolerance = 4e-6");
    const std::string creep_function = "function = norton";
    const std::size_t creep_function_position = m3_case.rfind(creep_function);
    if (creep_function_position == std::string::npos)
        return check(false, "transient fixture has coupled creep function");
    m3_case.replace(creep_function_position, creep_function.size(), "function = linear_temperature_norton");
    const std::string creep_exponent = "stress_exponent = 3";
    const std::size_t creep_position = m3_case.find(creep_exponent);
    if (creep_position == std::string::npos) return check(false, "transient fixture has coupled creep properties");
    m3_case.insert(creep_position + creep_exponent.size(), "\n      reference_temperature = 600"
                                                           "\n      coefficient_temperature_coefficient = 1e-8"
                                                           "\n      reference_stress_temperature_coefficient = 10"
                                                           "\n      stress_exponent_temperature_coefficient = 1e-4");
    const std::string plastic_function = "function = linear_isotropic_hardening";
    const std::size_t plastic_function_position = m3_case.rfind(plastic_function);
    if (plastic_function_position == std::string::npos)
        return check(false, "transient fixture has coupled plasticity function");
    m3_case.replace(
        plastic_function_position, plastic_function.size(), "function = linear_temperature_isotropic_hardening");
    const std::string hardening = "hardening_modulus = 2e9";
    const std::size_t hardening_position = m3_case.rfind(hardening);
    if (hardening_position == std::string::npos) return check(false, "transient fixture has coupled hardening modulus");
    m3_case.insert(hardening_position + hardening.size(), "\n      reference_temperature = 600"
                                                          "\n      yield_stress_temperature_coefficient = -100"
                                                          "\n      hardening_temperature_coefficient = -1000");
    const std::size_t output_position = m3_case.find(console);
    if (output_position == std::string::npos) return check(false, "transient fixture has the expected console key");
    m3_case.insert(output_position + console.size(), "\n  exodus = results.e\n  exodus_interval = 3"
                                                     "\n  history = history.csv\n  history_interval = 2"
                                                     "\n  progress_interval = 4"
                                                     "\n  checkpoint = checkpoint.bin"
                                                     "\n  checkpoint_interval = 5");
    passed = verify_m3_output_input(malformed_path, m3_case) && passed;
    std::string invalid_preconditioner = m3_case;
    const std::string valid_preconditioner = "preconditioner = field_split";
    const std::size_t preconditioner_position = invalid_preconditioner.find(valid_preconditioner);
    invalid_preconditioner.replace(preconditioner_position, valid_preconditioner.size(), "preconditioner = magic");
    passed = expect_case_failure(malformed_path, invalid_preconditioner, "preconditioner must be") && passed;
    std::string fixed_scaling_case = read_text(transient_path);
    const std::size_t fixed_scaling_solver = fixed_scaling_case.find(solver_start);
    if (fixed_scaling_solver == std::string::npos) return check(false, "transient fixture has a solver section");
    fixed_scaling_case.insert(fixed_scaling_solver + solver_start.size(), "\n  temperature_residual_scale = 1e4"
                                                                          "\n  mechanical_residual_scale = 1e3");
    {
        std::ofstream output(malformed_path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create fixed-scale input fixture");
        output << fixed_scaling_case;
    }
    const fuelsim::FuelSimCaseDefinition fixed_scaling = fuelsim::read_case_input(malformed_path);
    passed = check(fixed_scaling.solver.temperature_residual_scale == 1.0e4 &&
                       fixed_scaling.solver.mechanical_residual_scale == 1.0e3,
                 "fixed physical residual scales are parsed") &&
             passed;
    if (std::remove(malformed_path.c_str()) != 0) return check(false, "could not remove fixed-scale input fixture");
    std::string current_traction_case = read_text(traction_path);
    const std::string traction_type = "type = traction";
    const std::size_t traction_type_position = current_traction_case.find(traction_type);
    if (traction_type_position == std::string::npos) return check(false, "traction fixture has a traction condition");
    current_traction_case.insert(traction_type_position + traction_type.size(), "\n    configuration = current");
    {
        std::ofstream output(malformed_path, std::ios::out | std::ios::trunc);
        if (!output) return check(false, "could not create current-traction input fixture");
        output << current_traction_case;
    }
    const fuelsim::FuelSimCaseDefinition current_traction = fuelsim::read_case_input(malformed_path);
    passed = check(current_traction.spatial.boundary_conditions.back().use_displaced_geometry,
                 "current-configuration traction is parsed") &&
             passed;
    if (std::remove(malformed_path.c_str()) != 0)
        return check(false, "could not remove current-traction input fixture");
    std::string invalid_configuration = read_text(traction_path);
    const std::size_t invalid_configuration_position = invalid_configuration.find(traction_type);
    invalid_configuration.insert(
        invalid_configuration_position + traction_type.size(), "\n    configuration = rotating");
    passed = expect_case_failure(malformed_path, invalid_configuration, "must be reference or current") && passed;
    std::string dirichlet_configuration = read_text(transient_path);
    const std::string dirichlet_type = "type = dirichlet";
    const std::size_t dirichlet_position = dirichlet_configuration.find(dirichlet_type);
    if (dirichlet_position == std::string::npos) return check(false, "transient fixture has a Dirichlet condition");
    dirichlet_configuration.insert(dirichlet_position + dirichlet_type.size(), "\n    configuration = current");
    passed = expect_case_failure(malformed_path, dirichlet_configuration, "not valid for type='dirichlet'") && passed;
    std::string pressure_configuration = current_traction_case;
    pressure_configuration.replace(pressure_configuration.find(traction_type), traction_type.size(), "type = pressure");
    passed = expect_case_failure(malformed_path, pressure_configuration, "not valid for type='pressure'") && passed;
    std::string convection_configuration = m3_case;
    const std::string convection_type = "type = convection";
    const std::size_t convection_position = convection_configuration.find(convection_type);
    if (convection_position == std::string::npos) return check(false, "M3 fixture has a convection condition");
    convection_configuration.insert(convection_position + convection_type.size(), "\n    configuration = current");
    passed = expect_case_failure(malformed_path, convection_configuration, "not valid for type='convection'") && passed;
    std::string unknown_function = read_text(transient_path);
    const std::size_t unknown_heat_position = unknown_function.find(heat_source);
    unknown_function.insert(unknown_heat_position + heat_source.size(), "\n    heat_source_function = missing");
    passed = expect_case_failure(malformed_path, unknown_function, "unknown time function 'missing'") && passed;
    std::string invalid_table = m3_case;
    const std::string valid_times = "times = 0 2 5";
    const std::size_t valid_times_position = invalid_table.find(valid_times);
    if (valid_times_position == std::string::npos) return check(false, "M3 fixture has time-table nodes");
    invalid_table.replace(valid_times_position, valid_times.size(), "times = 0 0 5");
    passed = expect_case_failure(malformed_path, invalid_table, "strictly increasing") && passed;
    std::string orphan_interval = read_text(transient_path);
    const std::size_t orphan_output = orphan_interval.find(console);
    orphan_interval.insert(orphan_output + console.size(), "\n  checkpoint_interval = 2");
    passed = expect_case_failure(malformed_path, orphan_interval, "checkpoint_interval requires checkpoint") && passed;
    std::string steady_checkpoint = read_text(steady_path);
    const std::size_t steady_output = steady_checkpoint.find(console);
    steady_checkpoint.insert(steady_output + console.size(), "\n  checkpoint = checkpoint.bin");
    passed = expect_case_failure(malformed_path, steady_checkpoint, "only valid for transient cases") && passed;
    std::string mesh_overwrite = read_text(steady_path);
    const std::string mesh_file = "file = ../moose/m1_fuel_cladding_gap_rz_mesh.e";
    const std::size_t mesh_output = mesh_overwrite.find(console);
    if (mesh_overwrite.find(mesh_file) == std::string::npos || mesh_output == std::string::npos)
        return check(false, "steady fixture has expected mesh and output");
    mesh_overwrite.insert(mesh_output + console.size(), "\n  csv = ../moose/"
                                                        "m1_fuel_cladding_gap_rz_mesh.e");
    passed = expect_case_failure(malformed_path, mesh_overwrite, "must not overwrite the input mesh") && passed;
    std::string output_collision = read_text(transient_path);
    const std::size_t collision_output = output_collision.find(console);
    output_collision.insert(collision_output + console.size(), "\n  csv = collision.dat"
                                                               "\n  checkpoint = collision.dat");
    passed = expect_case_failure(malformed_path, output_collision, "paths must differ") && passed;
    passed = fuzz_input_parser(read_text(transient_path), malformed_path) && passed;
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
        if (!run_tests(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6])) return 1;
        std::cout << "[PASS] fuelsim strict input-card tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] input tests raised: " << error.what() << '\n';
        return 1;
    }
}
