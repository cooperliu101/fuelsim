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

const InputSection* find_section(const InputDocument& document,
                                 const std::string& path) {
    const auto found =
        std::find_if(document.sections().begin(), document.sections().end(),
                     [&path](const InputSection& section) {
                         return section.path() == path;
                     });
    return found == document.sections().end() ? nullptr : &*found;
}

std::vector<const InputSection*> direct_children(const InputDocument& document,
                                                 const std::string& parent) {
    std::vector<const InputSection*> result;
    const std::string prefix = parent + "/";
    for (const InputSection& section : document.sections()) {
        if (section.path().compare(0, prefix.size(), prefix) != 0)
            continue;
        const std::string suffix = section.path().substr(prefix.size());
        if (suffix.find('/') == std::string::npos)
            result.push_back(&section);
    }
    return result;
}

std::string leaf_name(const InputSection& section) {
    const std::size_t separator = section.path().find_last_of('/');
    return separator == std::string::npos
               ? section.path()
               : section.path().substr(separator + 1);
}

void validate_sections(const InputDocument& document) {
    const std::vector<std::string> fixed = {
        "Case",        "Mesh",   "Regions", "Contact", "BoundaryConditions",
        "Executioner", "Solver", "Outputs"};
    for (const InputSection& section : document.sections()) {
        if (std::find(fixed.begin(), fixed.end(), section.path()) !=
            fixed.end())
            continue;
        const std::string& path = section.path();
        const bool region = path.compare(0, 8, "Regions/") == 0 &&
                            path.find('/', 8) == std::string::npos;
        const bool boundary = path.compare(0, 19, "BoundaryConditions/") == 0 &&
                              path.find('/', 19) == std::string::npos;
        bool contact = false;
        if (path.compare(0, 8, "Contact/") == 0) {
            const std::string suffix = path.substr(8);
            const std::size_t separator = suffix.find('/');
            contact = separator == std::string::npos;
            if (separator != std::string::npos) {
                const std::string child = suffix.substr(separator + 1);
                contact =
                    suffix.find('/', separator + 1) == std::string::npos &&
                    (child == "thermal" || child == "mechanical");
            }
        }
        if (!region && !boundary && !contact)
            throw std::invalid_argument(document.source_path() + ":" +
                                        std::to_string(section.line()) +
                                        ": unknown section [" + path + "]");
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
    const std::size_t result = parse_size(document, *entry);
    if (result > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        value_error(document, *entry,
                    "key '" + key + "' is outside the supported range");
    return static_cast<int>(result);
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
                const std::string& key, const std::string& context) {
    const InputEntry* entry = find_entry(section, key);
    if (entry != nullptr)
        value_error(document, *entry,
                    "key '" + key + "' is not valid for " + context);
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
        forbid_key(document, section, "creep_coefficient",
                   "inelastic_model='elastic'");
        forbid_key(document, section, "creep_reference_stress",
                   "inelastic_model='elastic'");
        forbid_key(document, section, "creep_exponent",
                   "inelastic_model='elastic'");
        forbid_key(document, section, "yield_stress",
                   "inelastic_model='elastic'");
        forbid_key(document, section, "hardening_modulus",
                   "inelastic_model='elastic'");
    } else if (model == "norton_creep") {
        behavior = InelasticBehavior::norton_creep;
        creep = {read_double(document, section, "creep_coefficient"),
                 read_double(document, section, "creep_reference_stress"),
                 read_double(document, section, "creep_exponent")};
        forbid_key(document, section, "yield_stress", model);
        forbid_key(document, section, "hardening_modulus", model);
    } else if (model == "j2_plasticity") {
        behavior = InelasticBehavior::j2_plasticity;
        plasticity = {read_double(document, section, "yield_stress"),
                      read_double(document, section, "hardening_modulus")};
        forbid_key(document, section, "creep_coefficient", model);
        forbid_key(document, section, "creep_reference_stress", model);
        forbid_key(document, section, "creep_exponent", model);
    } else if (model == "norton_creep_j2_plasticity") {
        behavior = InelasticBehavior::norton_creep_j2_plasticity;
        creep = {read_double(document, section, "creep_coefficient"),
                 read_double(document, section, "creep_reference_stress"),
                 read_double(document, section, "creep_exponent")};
        plasticity = {read_double(document, section, "yield_stress"),
                      read_double(document, section, "hardening_modulus")};
    } else {
        value_error(document,
                    required_entry(document, section, "inelastic_model"),
                    "unknown inelastic_model '" + model + "'");
    }
    return {read_double(document, section, "density"),
            read_double(document, section, "specific_heat"), behavior, creep,
            plasticity};
}

Field parse_field(const InputDocument& document, const InputEntry& entry) {
    if (entry.value == "temperature")
        return Field::temperature;
    if (entry.value == "radial_displacement")
        return Field::radial_displacement;
    if (entry.value == "axial_displacement")
        return Field::axial_displacement;
    value_error(document, entry, "unknown field '" + entry.value + "'");
}

std::string resolved_path(const std::string& input_path,
                          const std::string& configured_path) {
    const std::filesystem::path configured(configured_path);
    if (configured.is_absolute())
        return configured.lexically_normal().string();
    return (std::filesystem::absolute(input_path).parent_path() / configured)
        .lexically_normal()
        .string();
}

CaseRegionDefinition read_region(const InputDocument& document,
                                 const InputSection& section,
                                 CaseProblem problem) {
    std::vector<std::string> keys = {
        "block",
        "conductivity_inverse_temperature",
        "conductivity_constant",
        "young_modulus",
        "poisson_ratio",
        "thermal_expansion",
        "reference_temperature",
        "initial_temperature",
        "volumetric_heat_source",
    };
    if (problem == CaseProblem::transient) {
        keys.insert(keys.end(),
                    {"density", "specific_heat", "inelastic_model",
                     "creep_coefficient", "creep_reference_stress",
                     "creep_exponent", "yield_stress", "hardening_modulus"});
    }
    validate_keys(document, section, keys);
    CaseRegionDefinition result{};
    result.spatial = {leaf_name(section),
                      read_string(document, section, "block"),
                      read_thermoelastic(document, section),
                      read_double(document, section, "volumetric_heat_source"),
                      read_double(document, section, "initial_temperature")};
    if (problem == CaseProblem::transient)
        result.transient_material = read_transient(document, section);
    return result;
}

ContactDefinition read_contact(const InputDocument& document,
                               const InputSection& section) {
    validate_keys(document, section, {"primary", "secondary"});
    const std::string base = section.path();
    const InputSection* thermal = find_section(document, base + "/thermal");
    const InputSection* mechanical =
        find_section(document, base + "/mechanical");
    if (thermal == nullptr && mechanical == nullptr)
        throw std::invalid_argument(
            document.source_path() + ":" + std::to_string(section.line()) +
            ": contact [" + base + "] requires [thermal] or [mechanical]");
    ContactDefinition result{leaf_name(section),
                             read_string(document, section, "primary"),
                             read_string(document, section, "secondary"),
                             thermal != nullptr,
                             mechanical != nullptr,
                             1.0,
                             1.0,
                             1.0};
    if (thermal != nullptr) {
        validate_keys(document, *thermal, {"gap_conductivity", "minimum_gap"});
        result.gap_conductivity =
            read_double(document, *thermal, "gap_conductivity");
        result.minimum_gap = read_double(document, *thermal, "minimum_gap");
    }
    if (mechanical != nullptr) {
        validate_keys(document, *mechanical, {"formulation", "penalty"});
        const std::string formulation =
            read_string(document, *mechanical, "formulation");
        if (formulation != "penalty")
            value_error(document, mechanical->entry("formulation"),
                        "only mechanical formulation 'penalty' is supported");
        result.penalty = read_double(document, *mechanical, "penalty");
    }
    return result;
}

BoundaryConditionDefinition
read_boundary_condition(const InputDocument& document,
                        const InputSection& section) {
    validate_keys(document, section, {"type", "boundary", "field", "value"});
    const std::string type = read_string(document, section, "type");
    if (type == "dirichlet") {
        return {
            leaf_name(section), BoundaryConditionType::dirichlet,
            read_string(document, section, "boundary"),
            parse_field(document, required_entry(document, section, "field")),
            read_double(document, section, "value")};
    }
    if (type == "pressure") {
        forbid_key(document, section, "field", "type='pressure'");
        return {leaf_name(section), BoundaryConditionType::pressure,
                read_string(document, section, "boundary"),
                Field::radial_displacement,
                read_double(document, section, "value")};
    }
    value_error(document, section.entry("type"),
                "unknown boundary-condition type '" + type + "'");
}

} // namespace

SteadyProblemDefinition FuelSimCaseDefinition::steady_definition() const {
    SteadyProblemDefinition result;
    result.contacts = contacts;
    result.boundary_conditions = boundary_conditions;
    result.regions.reserve(regions.size());
    for (const CaseRegionDefinition& region : regions)
        result.regions.push_back(region.spatial);
    return result;
}

TransientProblemDefinition FuelSimCaseDefinition::transient_definition() const {
    TransientProblemDefinition result;
    result.spatial = steady_definition();
    result.regions.reserve(regions.size());
    for (const CaseRegionDefinition& region : regions)
        result.regions.push_back(
            {region.spatial.name, region.transient_material});
    return result;
}

FuelSimCaseDefinition CaseInputReader::read(const std::string& path) {
    const InputDocument document = InputParser::parse_file(path);
    validate_sections(document);
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
    if (problem == "steady")
        result.problem = CaseProblem::steady;
    else if (problem == "transient")
        result.problem = CaseProblem::transient;
    else
        value_error(document, case_section.entry("problem"),
                    "unknown problem '" + problem + "'");

    const InputSection& mesh = document.section("Mesh");
    validate_keys(document, mesh, {"type", "file"});
    if (read_string(document, mesh, "type") != "exodus")
        value_error(document, mesh.entry("type"),
                    "only mesh type 'exodus' is supported");
    result.mesh_file = resolved_path(path, read_string(document, mesh, "file"));

    const InputSection& regions = document.section("Regions");
    validate_keys(document, regions, {});
    for (const InputSection* section : direct_children(document, "Regions"))
        result.regions.push_back(
            read_region(document, *section, result.problem));
    if (result.regions.empty())
        throw std::invalid_argument(path +
                                    ": [Regions] requires a child region");

    const InputSection& contacts = document.section("Contact");
    validate_keys(document, contacts, {});
    for (const InputSection* section : direct_children(document, "Contact"))
        result.contacts.push_back(read_contact(document, *section));

    const InputSection& boundary_conditions =
        document.section("BoundaryConditions");
    validate_keys(document, boundary_conditions, {});
    for (const InputSection* section :
         direct_children(document, "BoundaryConditions"))
        result.boundary_conditions.push_back(
            read_boundary_condition(document, *section));

    const InputSection& executioner = document.section("Executioner");
    const std::string executioner_type =
        read_string(document, executioner, "type");
    if (result.problem == CaseProblem::steady) {
        validate_keys(document, executioner, {"type", "load_steps"});
        if (executioner_type != "steady")
            value_error(document, executioner.entry("type"),
                        "problem='steady' requires type='steady'");
        result.steady_execution.load_steps =
            read_size(document, executioner, "load_steps");
        if (result.steady_execution.load_steps == 0)
            value_error(document, executioner.entry("load_steps"),
                        "load_steps must be positive");
    } else {
        validate_keys(document, executioner,
                      {"type", "end_time", "initial_time_step",
                       "minimum_time_step", "maximum_time_step",
                       "growth_factor", "cutback_factor", "maximum_cutbacks",
                       "heat_source_ramp_time"});
        if (executioner_type != "transient")
            value_error(document, executioner.entry("type"),
                        "problem='transient' requires type='transient'");
        result.transient_execution = {
            read_double(document, executioner, "end_time"),
            read_double(document, executioner, "initial_time_step"),
            read_double(document, executioner, "minimum_time_step"),
            read_double(document, executioner, "maximum_time_step"),
            read_double(document, executioner, "growth_factor"),
            read_double(document, executioner, "cutback_factor"),
            read_size(document, executioner, "maximum_cutbacks"),
            read_double(document, executioner, "heat_source_ramp_time")};
    }

    const InputSection& solver = document.section("Solver");
    validate_keys(document, solver,
                  {"absolute_tolerance", "relative_tolerance", "step_tolerance",
                   "maximum_iterations"});
    result.solver = {
        read_optional_double(document, solver, "absolute_tolerance", 1.0e-8),
        read_optional_double(document, solver, "relative_tolerance", 1.0e-10),
        read_optional_double(document, solver, "step_tolerance", 1.0e-12),
        read_optional_int(document, solver, "maximum_iterations", 40)};
    if (!(result.solver.absolute_tolerance > 0.0) ||
        !(result.solver.relative_tolerance > 0.0) ||
        !(result.solver.step_tolerance > 0.0) ||
        result.solver.maximum_iterations <= 0)
        throw std::invalid_argument(
            path + ": solver tolerances and iteration limit must be positive");

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
