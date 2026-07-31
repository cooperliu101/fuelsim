#include "fuelsim/case_input.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fuelsim {
namespace {

[[noreturn]] void value_error(const InputDocument& document,
                              const InputEntry& entry,
                              const std::string& message) {
    throw std::invalid_argument(document.source_path() + ":" +
                                std::to_string(entry.line) + ": " + message);
}

const InputEntry* find_entry(const InputSection& section,
                             const std::string& key) {
    const auto found = std::find_if(
        section.entries().begin(), section.entries().end(),
        [&key](const InputEntry& entry) { return entry.key == key; });
    return found == section.entries().end() ? nullptr : &*found;
}

const InputEntry& required_entry(const InputDocument& document,
                                 const InputSection& section,
                                 const std::string& key) {
    const InputEntry* entry = find_entry(section, key);
    if (entry == nullptr)
        throw std::invalid_argument(document.source_path() + ":" +
                                    std::to_string(section.line()) +
                                    ": section [" + section.path() +
                                    "] is missing required key '" + key + "'");
    return *entry;
}

void validate_sections(const InputDocument& document,
                       const std::vector<std::string>& allowed) {
    for (const InputSection& section : document.sections()) {
        if (std::find(allowed.begin(), allowed.end(), section.path()) ==
            allowed.end())
            throw std::invalid_argument(
                document.source_path() + ":" + std::to_string(section.line()) +
                ": unknown section [" + section.path() + "]");
    }
}

void validate_keys(const InputDocument& document, const InputSection& section,
                   const std::vector<std::string>& allowed) {
    for (const InputEntry& entry : section.entries()) {
        if (std::find(allowed.begin(), allowed.end(), entry.key) ==
            allowed.end())
            value_error(document, entry,
                        "unknown key '" + entry.key + "' in [" +
                            section.path() + "]");
    }
}

std::string read_string(const InputDocument& document,
                        const InputSection& section, const std::string& key) {
    return required_entry(document, section, key).value;
}

double parse_double(const InputDocument& document, const InputEntry& entry) {
    errno = 0;
    char* end = nullptr;
    const double result = std::strtod(entry.value.c_str(), &end);
    if (errno == ERANGE || end == entry.value.c_str() || *end != '\0' ||
        !std::isfinite(result))
        value_error(document, entry,
                    "key '" + entry.key + "' requires a finite real value");
    return result;
}

double read_double(const InputDocument& document, const InputSection& section,
                   const std::string& key) {
    return parse_double(document, required_entry(document, section, key));
}

double read_optional_double(const InputDocument& document,
                            const InputSection& section, const std::string& key,
                            double fallback) {
    const InputEntry* entry = find_entry(section, key);
    return entry == nullptr ? fallback : parse_double(document, *entry);
}

std::size_t parse_size(const InputDocument& document, const InputEntry& entry) {
    if (!entry.value.empty() && entry.value.front() == '-')
        value_error(document, entry,
                    "key '" + entry.key +
                        "' requires a nonnegative integer value");
    errno = 0;
    char* end = nullptr;
    const unsigned long long result =
        std::strtoull(entry.value.c_str(), &end, 10);
    if (errno == ERANGE || end == entry.value.c_str() || *end != '\0' ||
        result > std::numeric_limits<std::size_t>::max())
        value_error(document, entry,
                    "key '" + entry.key +
                        "' requires a nonnegative integer value");
    return static_cast<std::size_t>(result);
}

std::size_t read_size(const InputDocument& document,
                      const InputSection& section, const std::string& key) {
    return parse_size(document, required_entry(document, section, key));
}

int read_optional_int(const InputDocument& document,
                      const InputSection& section, const std::string& key,
                      int fallback) {
    const InputEntry* entry = find_entry(section, key);
    if (entry == nullptr)
        return fallback;
    if (!entry->value.empty() && entry->value.front() == '-')
        value_error(document, *entry,
                    "key '" + key + "' requires a nonnegative integer value");
    const std::size_t value = parse_size(document, *entry);
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        value_error(document, *entry,
                    "key '" + key + "' is outside the supported range");
    return static_cast<int>(value);
}

bool read_optional_bool(const InputDocument& document,
                        const InputSection& section, const std::string& key,
                        bool fallback) {
    const InputEntry* entry = find_entry(section, key);
    if (entry == nullptr)
        return fallback;
    if (entry->value == "true")
        return true;
    if (entry->value == "false")
        return false;
    value_error(document, *entry,
                "key '" + key + "' requires 'true' or 'false'");
}

std::string read_optional_string(const InputSection& section,
                                 const std::string& key,
                                 const std::string& fallback) {
    const InputEntry* entry = find_entry(section, key);
    return entry == nullptr ? fallback : entry->value;
}

void forbid_key(const InputDocument& document, const InputSection& section,
                const std::string& key, const std::string& model) {
    const InputEntry* entry = find_entry(section, key);
    if (entry != nullptr)
        value_error(document, *entry,
                    "key '" + key + "' is not valid for inelastic_model='" +
                        model + "'");
}

ThermoelasticProperties read_thermoelastic(const InputDocument& document,
                                           const InputSection& section) {
    return {
        read_double(document, section, "conductivity_inverse_temperature"),
        read_double(document, section, "conductivity_constant"),
        read_double(document, section, "young_modulus"),
        read_double(document, section, "poisson_ratio"),
        read_double(document, section, "thermal_expansion"),
        read_double(document, section, "reference_temperature"),
    };
}

TransientInelasticProperties read_transient(const InputDocument& document,
                                            const InputSection& section) {
    const std::string model = read_string(document, section, "inelastic_model");
    InelasticBehavior behavior = InelasticBehavior::elastic;
    NortonCreepProperties creep{0.0, 1.0, 1.0};
    J2PlasticityProperties plasticity{1.0, 0.0};

    if (model == "elastic") {
        forbid_key(document, section, "creep_coefficient", model);
        forbid_key(document, section, "creep_reference_stress", model);
        forbid_key(document, section, "creep_exponent", model);
        forbid_key(document, section, "yield_stress", model);
        forbid_key(document, section, "hardening_modulus", model);
    } else if (model == "norton_creep") {
        behavior = InelasticBehavior::norton_creep;
        creep = {
            read_double(document, section, "creep_coefficient"),
            read_double(document, section, "creep_reference_stress"),
            read_double(document, section, "creep_exponent"),
        };
        forbid_key(document, section, "yield_stress", model);
        forbid_key(document, section, "hardening_modulus", model);
    } else if (model == "j2_plasticity") {
        behavior = InelasticBehavior::j2_plasticity;
        plasticity = {
            read_double(document, section, "yield_stress"),
            read_double(document, section, "hardening_modulus"),
        };
        forbid_key(document, section, "creep_coefficient", model);
        forbid_key(document, section, "creep_reference_stress", model);
        forbid_key(document, section, "creep_exponent", model);
    } else if (model == "norton_creep_j2_plasticity") {
        behavior = InelasticBehavior::norton_creep_j2_plasticity;
        creep = {
            read_double(document, section, "creep_coefficient"),
            read_double(document, section, "creep_reference_stress"),
            read_double(document, section, "creep_exponent"),
        };
        plasticity = {
            read_double(document, section, "yield_stress"),
            read_double(document, section, "hardening_modulus"),
        };
    } else {
        value_error(document,
                    required_entry(document, section, "inelastic_model"),
                    "unknown inelastic_model '" + model + "'");
    }

    return {
        read_double(document, section, "density"),
        read_double(document, section, "specific_heat"),
        behavior,
        creep,
        plasticity,
    };
}

MeshRegionInput read_region(const InputDocument& document,
                            const InputSection& section) {
    validate_keys(document, section,
                  {"block", "radial_inner", "radial_outer", "bottom", "top"});
    return {
        read_string(document, section, "block"),
        {
            read_string(document, section, "radial_inner"),
            read_string(document, section, "radial_outer"),
            read_string(document, section, "bottom"),
            read_string(document, section, "top"),
        },
    };
}

std::string resolved_path(const std::string& input_path,
                          const std::string& configured_path) {
    const std::filesystem::path configured(configured_path);
    if (configured.is_absolute())
        return configured.lexically_normal().string();
    const std::filesystem::path parent =
        std::filesystem::absolute(input_path).parent_path();
    return (parent / configured).lexically_normal().string();
}

} // namespace

SteadyFuelCladdingParameters FuelSimCaseDefinition::steady_parameters(
    const StructuredRzMesh& fuel_mesh,
    const StructuredRzMesh& cladding_mesh) const {
    return {
        fuel_mesh.outer_radius(),
        cladding_mesh.inner_radius(),
        cladding_mesh.outer_radius(),
        fuel_mesh.length(),
        cladding_mesh.length(),
        fuel_mesh.radial_elements(),
        cladding_mesh.radial_elements(),
        fuel_mesh.axial_elements(),
        fuel_thermoelastic,
        cladding_thermoelastic,
        physics.final_heat_source,
        physics.outer_temperature,
        physics.initial_temperature,
        physics.gap_conductivity,
        physics.minimum_gap,
        physics.contact_penalty,
    };
}

TransientFuelCladdingParameters FuelSimCaseDefinition::transient_parameters(
    const StructuredRzMesh& fuel_mesh,
    const StructuredRzMesh& cladding_mesh) const {
    return {
        steady_parameters(fuel_mesh, cladding_mesh),
        fuel_transient,
        cladding_transient,
    };
}

FuelSimCaseDefinition CaseInputReader::read(const std::string& path) {
    const InputDocument document = InputParser::parse_file(path);
    validate_sections(document,
                      {"Case", "Mesh", "Mesh/fuel", "Mesh/cladding",
                       "Materials", "Materials/fuel", "Materials/cladding",
                       "Physics", "Executioner", "Solver", "Outputs"});

    FuelSimCaseDefinition result{};
    const InputSection& case_section = document.section("Case");
    validate_keys(document, case_section, {"version", "problem"});
    const std::size_t version = read_size(document, case_section, "version");
    if (version != 1)
        value_error(document, case_section.entry("version"),
                    "unsupported fuelsim input version '" +
                        std::to_string(version) + "'");
    result.version = 1;
    const std::string problem = read_string(document, case_section, "problem");
    if (problem == "steady_fuel_cladding")
        result.problem = CaseProblem::steady_fuel_cladding;
    else if (problem == "transient_fuel_cladding")
        result.problem = CaseProblem::transient_fuel_cladding;
    else
        value_error(document, case_section.entry("problem"),
                    "unknown problem '" + problem + "'");

    const InputSection& mesh_section = document.section("Mesh");
    validate_keys(document, mesh_section, {"type", "file"});
    const std::string mesh_type = read_string(document, mesh_section, "type");
    if (mesh_type != "exodus")
        value_error(document, mesh_section.entry("type"),
                    "only mesh type 'exodus' is supported");
    result.mesh.file =
        resolved_path(path, read_string(document, mesh_section, "file"));
    result.mesh.fuel = read_region(document, document.section("Mesh/fuel"));
    result.mesh.cladding =
        read_region(document, document.section("Mesh/cladding"));

    const InputSection& materials = document.section("Materials");
    validate_keys(document, materials, {});
    const InputSection& fuel_material = document.section("Materials/fuel");
    const InputSection& cladding_material =
        document.section("Materials/cladding");
    const std::vector<std::string> steady_material_keys = {
        "conductivity_inverse_temperature",
        "conductivity_constant",
        "young_modulus",
        "poisson_ratio",
        "thermal_expansion",
        "reference_temperature"};
    std::vector<std::string> material_keys = steady_material_keys;
    if (result.problem == CaseProblem::transient_fuel_cladding) {
        material_keys.insert(material_keys.end(),
                             {"density", "specific_heat", "inelastic_model",
                              "creep_coefficient", "creep_reference_stress",
                              "creep_exponent", "yield_stress",
                              "hardening_modulus"});
    }
    validate_keys(document, fuel_material, material_keys);
    validate_keys(document, cladding_material, material_keys);
    result.fuel_thermoelastic = read_thermoelastic(document, fuel_material);
    result.cladding_thermoelastic =
        read_thermoelastic(document, cladding_material);
    if (result.problem == CaseProblem::transient_fuel_cladding) {
        result.fuel_transient = read_transient(document, fuel_material);
        result.cladding_transient = read_transient(document, cladding_material);
    }

    const InputSection& physics = document.section("Physics");
    validate_keys(document, physics,
                  {"initial_temperature", "outer_temperature",
                   "final_heat_source", "heat_source_ramp_time",
                   "gap_conductivity", "minimum_gap", "contact_penalty"});
    result.physics = {
        read_double(document, physics, "initial_temperature"),
        read_double(document, physics, "outer_temperature"),
        read_double(document, physics, "final_heat_source"),
        read_optional_double(document, physics, "heat_source_ramp_time", 0.0),
        read_double(document, physics, "gap_conductivity"),
        read_double(document, physics, "minimum_gap"),
        read_double(document, physics, "contact_penalty"),
    };

    const InputSection& executioner = document.section("Executioner");
    const std::string executioner_type =
        read_string(document, executioner, "type");
    if (result.problem == CaseProblem::steady_fuel_cladding) {
        validate_keys(document, executioner, {"type", "load_steps"});
        if (executioner_type != "steady")
            value_error(document, executioner.entry("type"),
                        "steady_fuel_cladding requires type='steady'");
        result.steady_execution.load_steps =
            read_size(document, executioner, "load_steps");
        if (result.steady_execution.load_steps == 0)
            value_error(document, executioner.entry("load_steps"),
                        "load_steps must be positive");
    } else {
        validate_keys(document, executioner,
                      {"type", "end_time", "initial_time_step",
                       "minimum_time_step", "maximum_time_step",
                       "growth_factor", "cutback_factor", "maximum_cutbacks"});
        if (executioner_type != "transient")
            value_error(document, executioner.entry("type"),
                        "transient_fuel_cladding requires type='transient'");
        result.transient_execution = {
            read_double(document, executioner, "end_time"),
            read_double(document, executioner, "initial_time_step"),
            read_double(document, executioner, "minimum_time_step"),
            read_double(document, executioner, "maximum_time_step"),
            read_double(document, executioner, "growth_factor"),
            read_double(document, executioner, "cutback_factor"),
            read_size(document, executioner, "maximum_cutbacks"),
        };
    }

    const InputSection& solver = document.section("Solver");
    validate_keys(document, solver,
                  {"absolute_tolerance", "relative_tolerance", "step_tolerance",
                   "maximum_iterations"});
    result.solver = {
        read_optional_double(document, solver, "absolute_tolerance", 1.0e-8),
        read_optional_double(document, solver, "relative_tolerance", 1.0e-10),
        read_optional_double(document, solver, "step_tolerance", 1.0e-12),
        read_optional_int(document, solver, "maximum_iterations", 40),
    };
    if (!(result.solver.absolute_tolerance > 0.0))
        value_error(document, solver.entry("absolute_tolerance"),
                    "absolute_tolerance must be positive");
    if (!(result.solver.relative_tolerance > 0.0))
        value_error(document, solver.entry("relative_tolerance"),
                    "relative_tolerance must be positive");
    if (!(result.solver.step_tolerance > 0.0))
        value_error(document, solver.entry("step_tolerance"),
                    "step_tolerance must be positive");
    if (result.solver.maximum_iterations <= 0)
        value_error(document, solver.entry("maximum_iterations"),
                    "maximum_iterations must be positive");

    const InputSection& outputs = document.section("Outputs");
    validate_keys(document, outputs, {"console", "csv"});
    result.outputs.console =
        read_optional_bool(document, outputs, "console", true);
    const std::string csv = read_optional_string(outputs, "csv", {});
    result.outputs.csv_file =
        csv.empty() ? std::string{} : resolved_path(path, csv);
    return result;
}

} // namespace fuelsim
