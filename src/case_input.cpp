#include "fuelsim/case_input.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <sstream>
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
        "Case",        "Mesh",    "TimeFunctions",
        "Regions",     "Contact", "BoundaryConditions",
        "Executioner", "Solver",  "Outputs"};
    for (const InputSection& section : document.sections()) {
        if (std::find(fixed.begin(), fixed.end(), section.path()) !=
            fixed.end())
            continue;
        const std::string& path = section.path();
        const bool region = path.compare(0, 8, "Regions/") == 0 &&
                            path.find('/', 8) == std::string::npos;
        const bool boundary = path.compare(0, 19, "BoundaryConditions/") == 0 &&
                              path.find('/', 19) == std::string::npos;
        const bool function = path.compare(0, 14, "TimeFunctions/") == 0 &&
                              path.find('/', 14) == std::string::npos;
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
        if (!region && !boundary && !function && !contact)
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

std::vector<double> parse_double_list(const InputDocument& document,
                                      const InputEntry& entry) {
    std::istringstream input(entry.value);
    std::vector<double> result;
    double value = 0.0;
    while (input >> value) {
        if (!std::isfinite(value))
            value_error(document, entry,
                        "key '" + entry.key + "' requires finite real values");
        result.push_back(value);
    }
    if (!input.eof() || result.empty())
        value_error(document, entry,
                    "key '" + entry.key +
                        "' requires a whitespace-separated real list");
    return result;
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

std::size_t read_optional_size(const InputDocument& document,
                               const InputSection& section,
                               const std::string& key, std::size_t fallback) {
    const InputEntry* entry = find_entry(section, key);
    return entry == nullptr ? fallback : parse_size(document, *entry);
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

void forbid_convection_keys(const InputDocument& document,
                            const InputSection& section,
                            const std::string& context) {
    forbid_key(document, section, "heat_transfer_coefficient", context);
    forbid_key(document, section, "ambient_temperature", context);
    forbid_key(document, section, "coefficient_function", context);
    forbid_key(document, section, "ambient_temperature_function", context);
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
        read_optional_double(document, section,
                             "young_modulus_temperature_coefficient", 0.0),
        read_optional_double(document, section,
                             "poisson_ratio_temperature_coefficient", 0.0),
        read_optional_double(
            document, section,
            "thermal_expansion_temperature_coefficient", 0.0),
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
        forbid_key(document, section,
                   "creep_coefficient_temperature_coefficient",
                   "inelastic_model='elastic'");
        forbid_key(document, section,
                   "creep_reference_stress_temperature_coefficient",
                   "inelastic_model='elastic'");
        forbid_key(document, section,
                   "creep_exponent_temperature_coefficient",
                   "inelastic_model='elastic'");
        forbid_key(document, section, "yield_stress",
                   "inelastic_model='elastic'");
        forbid_key(document, section, "hardening_modulus",
                   "inelastic_model='elastic'");
        forbid_key(document, section,
                   "yield_stress_temperature_coefficient",
                   "inelastic_model='elastic'");
        forbid_key(document, section, "hardening_temperature_coefficient",
                   "inelastic_model='elastic'");
    } else if (model == "norton_creep") {
        behavior = InelasticBehavior::norton_creep;
        creep = {read_double(document, section, "creep_coefficient"),
                 read_double(document, section, "creep_reference_stress"),
                 read_double(document, section, "creep_exponent"),
                 read_optional_double(
                     document, section,
                     "creep_coefficient_temperature_coefficient", 0.0),
                 read_optional_double(
                     document, section,
                     "creep_reference_stress_temperature_coefficient", 0.0),
                 read_optional_double(
                     document, section,
                     "creep_exponent_temperature_coefficient", 0.0)};
        forbid_key(document, section, "yield_stress", model);
        forbid_key(document, section, "hardening_modulus", model);
        forbid_key(document, section,
                   "yield_stress_temperature_coefficient", model);
        forbid_key(document, section, "hardening_temperature_coefficient",
                   model);
    } else if (model == "j2_plasticity") {
        behavior = InelasticBehavior::j2_plasticity;
        plasticity = {read_double(document, section, "yield_stress"),
                      read_double(document, section, "hardening_modulus"),
                      read_optional_double(
                          document, section,
                          "yield_stress_temperature_coefficient", 0.0),
                      read_optional_double(
                          document, section,
                          "hardening_temperature_coefficient", 0.0)};
        forbid_key(document, section, "creep_coefficient", model);
        forbid_key(document, section, "creep_reference_stress", model);
        forbid_key(document, section, "creep_exponent", model);
        forbid_key(document, section,
                   "creep_coefficient_temperature_coefficient", model);
        forbid_key(document, section,
                   "creep_reference_stress_temperature_coefficient", model);
        forbid_key(document, section,
                   "creep_exponent_temperature_coefficient", model);
    } else if (model == "norton_creep_j2_plasticity") {
        behavior = InelasticBehavior::norton_creep_j2_plasticity;
        creep = {read_double(document, section, "creep_coefficient"),
                 read_double(document, section, "creep_reference_stress"),
                 read_double(document, section, "creep_exponent"),
                 read_optional_double(
                     document, section,
                     "creep_coefficient_temperature_coefficient", 0.0),
                 read_optional_double(
                     document, section,
                     "creep_reference_stress_temperature_coefficient", 0.0),
                 read_optional_double(
                     document, section,
                     "creep_exponent_temperature_coefficient", 0.0)};
        plasticity = {read_double(document, section, "yield_stress"),
                      read_double(document, section, "hardening_modulus"),
                      read_optional_double(
                          document, section,
                          "yield_stress_temperature_coefficient", 0.0),
                      read_optional_double(
                          document, section,
                          "hardening_temperature_coefficient", 0.0)};
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
        "block_id",
        "conductivity_inverse_temperature",
        "conductivity_constant",
        "young_modulus",
        "poisson_ratio",
        "thermal_expansion",
        "reference_temperature",
        "young_modulus_temperature_coefficient",
        "poisson_ratio_temperature_coefficient",
        "thermal_expansion_temperature_coefficient",
        "initial_temperature",
        "volumetric_heat_source",
        "heat_source_function",
    };
    if (problem == CaseProblem::transient) {
        keys.insert(keys.end(),
                    {"density", "specific_heat", "inelastic_model",
                     "creep_coefficient", "creep_reference_stress",
                     "creep_exponent", "yield_stress", "hardening_modulus",
                     "creep_coefficient_temperature_coefficient",
                     "creep_reference_stress_temperature_coefficient",
                     "creep_exponent_temperature_coefficient",
                     "yield_stress_temperature_coefficient",
                     "hardening_temperature_coefficient"});
    }
    validate_keys(document, section, keys);
    const InputEntry* block = find_entry(section, "block");
    const InputEntry* block_id = find_entry(section, "block_id");
    if ((block == nullptr) == (block_id == nullptr))
        throw std::invalid_argument(
            document.source_path() + ":" + std::to_string(section.line()) +
            ": region [" + section.path() +
            "] requires exactly one of 'block' or 'block_id'");
    std::int64_t resolved_block_id = -1;
    if (block_id != nullptr) {
        const std::size_t value = parse_size(document, *block_id);
        if (value >
            static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
            value_error(document, *block_id,
                        "key 'block_id' is outside the supported range");
        resolved_block_id = static_cast<std::int64_t>(value);
    }
    CaseRegionDefinition result{};
    result.spatial = {leaf_name(section),
                      block == nullptr ? std::string{} : block->value,
                      read_thermoelastic(document, section),
                      read_double(document, section, "volumetric_heat_source"),
                      read_double(document, section, "initial_temperature"),
                      resolved_block_id};
    result.spatial.heat_source_function =
        read_optional_string(section, "heat_source_function", {});
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
    validate_keys(document, section,
                  {"type", "boundary", "field", "value", "scale_with_load",
                   "function", "heat_transfer_coefficient",
                   "ambient_temperature", "coefficient_function",
                   "ambient_temperature_function"});
    const std::string type = read_string(document, section, "type");
    const bool scale_with_load =
        read_optional_bool(document, section, "scale_with_load", false);
    const std::string function = read_optional_string(section, "function", {});
    if (scale_with_load && !function.empty())
        value_error(document, section.entry("scale_with_load"),
                    "scale_with_load cannot be combined with function");
    if (type == "dirichlet") {
        forbid_convection_keys(document, section, "type='dirichlet'");
        BoundaryConditionDefinition result{
            leaf_name(section),
            BoundaryConditionType::dirichlet,
            read_string(document, section, "boundary"),
            parse_field(document, required_entry(document, section, "field")),
            read_double(document, section, "value"),
            scale_with_load};
        result.function = function;
        return result;
    }
    if (type == "pressure") {
        forbid_key(document, section, "field", "type='pressure'");
        forbid_convection_keys(document, section, "type='pressure'");
        BoundaryConditionDefinition result{
            leaf_name(section),
            BoundaryConditionType::pressure,
            read_string(document, section, "boundary"),
            Field::radial_displacement,
            read_double(document, section, "value"),
            scale_with_load};
        result.function = function;
        return result;
    }
    if (type == "traction") {
        forbid_convection_keys(document, section, "type='traction'");
        const Field field =
            parse_field(document, required_entry(document, section, "field"));
        if (field == Field::temperature)
            value_error(document, section.entry("field"),
                        "traction requires a displacement field");
        BoundaryConditionDefinition result{
            leaf_name(section),
            BoundaryConditionType::traction,
            read_string(document, section, "boundary"),
            field,
            read_double(document, section, "value"),
            scale_with_load};
        result.function = function;
        return result;
    }
    if (type == "convection") {
        forbid_key(document, section, "field", "type='convection'");
        forbid_key(document, section, "value", "type='convection'");
        forbid_key(document, section, "scale_with_load", "type='convection'");
        forbid_key(document, section, "function", "type='convection'");
        BoundaryConditionDefinition result{
            leaf_name(section),
            BoundaryConditionType::convection,
            read_string(document, section, "boundary"),
            Field::temperature,
            0.0,
            false};
        result.heat_transfer_coefficient =
            read_double(document, section, "heat_transfer_coefficient");
        result.ambient_temperature =
            read_double(document, section, "ambient_temperature");
        result.coefficient_function =
            read_optional_string(section, "coefficient_function", {});
        result.ambient_temperature_function =
            read_optional_string(section, "ambient_temperature_function", {});
        return result;
    }
    value_error(document, section.entry("type"),
                "unknown boundary-condition type '" + type + "'");
}

} // namespace

SteadyProblemDefinition FuelSimCaseDefinition::steady_definition() const {
    SteadyProblemDefinition result;
    result.contacts = contacts;
    result.boundary_conditions = boundary_conditions;
    result.time_tables = time_tables;
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

    const InputSection* time_functions =
        find_section(document, "TimeFunctions");
    if (time_functions != nullptr) {
        validate_keys(document, *time_functions, {});
        for (const InputSection* section :
             direct_children(document, "TimeFunctions")) {
            validate_keys(document, *section, {"type", "times", "values"});
            if (read_string(document, *section, "type") != "piecewise_linear")
                value_error(document, section->entry("type"),
                            "only time-function type 'piecewise_linear' is "
                            "supported");
            result.time_tables.emplace_back(
                leaf_name(*section),
                parse_double_list(document, section->entry("times")),
                parse_double_list(document, section->entry("values")));
        }
    }

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

    const auto require_function = [&](const std::string& name,
                                      const std::string& owner) {
        if (name.empty())
            return;
        const auto found =
            std::find_if(result.time_tables.begin(), result.time_tables.end(),
                         [&name](const PiecewiseLinearTimeTable& table) {
                             return table.name() == name;
                         });
        if (found == result.time_tables.end())
            throw std::invalid_argument(path + ": " + owner +
                                        " references unknown time function '" +
                                        name + "'");
    };
    for (const CaseRegionDefinition& region : result.regions)
        require_function(region.spatial.heat_source_function,
                         "region '" + region.spatial.name + "'");
    for (const BoundaryConditionDefinition& boundary :
         result.boundary_conditions) {
        require_function(boundary.function,
                         "boundary condition '" + boundary.name + "'");
        require_function(boundary.coefficient_function,
                         "boundary condition '" + boundary.name + "'");
        require_function(boundary.ambient_temperature_function,
                         "boundary condition '" + boundary.name + "'");
    }
    if (result.problem == CaseProblem::steady && !result.time_tables.empty())
        throw std::invalid_argument(
            path + ": time functions are only valid for transient cases");

    const InputSection& executioner = document.section("Executioner");
    const std::string executioner_type =
        read_string(document, executioner, "type");
    if (result.problem == CaseProblem::steady) {
        validate_keys(document, executioner,
                      {"type", "load_steps", "cutback_factor",
                       "maximum_cutbacks", "minimum_load_increment"});
        if (executioner_type != "steady")
            value_error(document, executioner.entry("type"),
                        "problem='steady' requires type='steady'");
        result.steady_execution.load_steps =
            read_size(document, executioner, "load_steps");
        if (result.steady_execution.load_steps == 0)
            value_error(document, executioner.entry("load_steps"),
                        "load_steps must be positive");
        result.steady_execution.cutback_factor = read_optional_double(
            document, executioner, "cutback_factor", 0.5);
        result.steady_execution.maximum_cutbacks = read_optional_size(
            document, executioner, "maximum_cutbacks", 12);
        result.steady_execution.minimum_load_increment = read_optional_double(
            document, executioner, "minimum_load_increment", 1.0e-6);
        if (!(result.steady_execution.cutback_factor > 0.0 &&
              result.steady_execution.cutback_factor < 1.0))
            value_error(document, executioner.entry("cutback_factor"),
                        "cutback_factor must lie between zero and one");
        if (!(result.steady_execution.minimum_load_increment > 0.0 &&
              result.steady_execution.minimum_load_increment <= 1.0))
            value_error(document,
                        executioner.entry("minimum_load_increment"),
                        "minimum_load_increment must lie in (0, 1]");
    } else {
        validate_keys(document, executioner,
                      {"type", "end_time", "initial_time_step",
                       "minimum_time_step", "maximum_time_step",
                       "growth_factor", "cutback_factor", "maximum_cutbacks",
                       "load_ramp_time", "restart",
                       "target_nonlinear_iterations", "iteration_window",
                       "time_error_relative_tolerance",
                       "temperature_time_absolute_tolerance",
                       "displacement_time_absolute_tolerance",
                       "time_error_safety_factor"});
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
            read_double(document, executioner, "load_ramp_time"),
            {},
            read_optional_size(document, executioner,
                               "target_nonlinear_iterations", 0),
            read_optional_size(document, executioner, "iteration_window", 0),
            read_optional_double(document, executioner,
                                 "time_error_relative_tolerance", 0.0),
            read_optional_double(document, executioner,
                                 "temperature_time_absolute_tolerance",
                                 1.0e-3),
            read_optional_double(document, executioner,
                                 "displacement_time_absolute_tolerance",
                                 1.0e-10),
            read_optional_double(document, executioner,
                                 "time_error_safety_factor", 0.9)};
        const std::string restart =
            read_optional_string(executioner, "restart", {});
        result.transient_execution.restart_file =
            restart.empty() ? std::string{} : resolved_path(path, restart);
        if (result.transient_execution.target_nonlinear_iterations == 0 &&
            result.transient_execution.iteration_window != 0)
            value_error(document, executioner.entry("iteration_window"),
                        "iteration_window requires positive "
                        "target_nonlinear_iterations");
        if (result.transient_execution.target_nonlinear_iterations > 0 &&
            result.transient_execution.iteration_window >=
                result.transient_execution.target_nonlinear_iterations)
            value_error(document, executioner.entry("iteration_window"),
                        "iteration_window must be smaller than "
                        "target_nonlinear_iterations");
        if (!(result.transient_execution.time_error_relative_tolerance >=
                  0.0) ||
            !(result.transient_execution
                      .temperature_time_absolute_tolerance > 0.0) ||
            !(result.transient_execution
                      .displacement_time_absolute_tolerance > 0.0) ||
            !(result.transient_execution.time_error_safety_factor > 0.0 &&
              result.transient_execution.time_error_safety_factor < 1.0))
            throw std::invalid_argument(
                path + ": transient time-error tolerances must be finite and "
                       "nonnegative/positive, and the safety factor must lie "
                       "in (0, 1)");
    }

    const InputSection& solver = document.section("Solver");
    validate_keys(document, solver,
                  {"absolute_tolerance", "relative_tolerance", "step_tolerance",
                   "maximum_iterations", "linear_solver", "preconditioner",
                   "linear_relative_tolerance",
                   "maximum_linear_iterations", "backtracking_fallback",
                   "field_residual_scaling",
                   "residual_reduction_tolerance",
                   "temperature_residual_absolute_tolerance",
                   "mechanical_residual_absolute_tolerance"});
    result.solver = {
        read_optional_double(document, solver, "absolute_tolerance", 1.0e-8),
        read_optional_double(document, solver, "relative_tolerance", 1.0e-10),
        read_optional_double(document, solver, "step_tolerance", 1.0e-12),
        read_optional_int(document, solver, "maximum_iterations", 40),
        read_optional_string(solver, "linear_solver", "automatic"),
        read_optional_string(solver, "preconditioner", "automatic"),
        read_optional_double(document, solver, "linear_relative_tolerance",
                             1.0e-8),
        read_optional_int(document, solver, "maximum_linear_iterations",
                          500),
        read_optional_bool(document, solver, "backtracking_fallback", true),
        read_optional_bool(document, solver, "field_residual_scaling", false),
        read_optional_double(document, solver,
                             "residual_reduction_tolerance", 1.0e-6),
        read_optional_double(document, solver,
                             "temperature_residual_absolute_tolerance",
                             1.0e-8),
        read_optional_double(document, solver,
                             "mechanical_residual_absolute_tolerance",
                             1.0e-4)};
    if (!(result.solver.absolute_tolerance > 0.0) ||
        !(result.solver.relative_tolerance > 0.0) ||
        !(result.solver.step_tolerance > 0.0) ||
        result.solver.maximum_iterations <= 0 ||
        !(result.solver.linear_relative_tolerance > 0.0) ||
        !(result.solver.residual_reduction_tolerance > 0.0) ||
        !(result.solver.temperature_residual_absolute_tolerance > 0.0) ||
        !(result.solver.mechanical_residual_absolute_tolerance > 0.0) ||
        result.solver.maximum_linear_iterations <= 0)
        throw std::invalid_argument(
            path + ": solver tolerances and iteration limit must be positive");
    if (result.solver.linear_solver != "automatic" &&
        result.solver.linear_solver != "direct" &&
        result.solver.linear_solver != "gmres")
        throw std::invalid_argument(
            path + ": linear_solver must be automatic, direct, or gmres");
    if (result.solver.preconditioner != "automatic" &&
        result.solver.preconditioner != "lu" &&
        result.solver.preconditioner != "block_jacobi" &&
        result.solver.preconditioner != "field_split" &&
        result.solver.preconditioner != "hypre")
        throw std::invalid_argument(
            path + ": preconditioner must be automatic, lu, block_jacobi, "
                   "field_split, or hypre");

    const InputSection& outputs = document.section("Outputs");
    validate_keys(
        document, outputs,
        {"console", "csv", "exodus", "exodus_interval", "history",
         "history_interval", "progress_interval", "checkpoint",
         "checkpoint_interval"});
    result.outputs.console =
        read_optional_bool(document, outputs, "console", true);
    const std::string csv = read_optional_string(outputs, "csv", {});
    result.outputs.csv_file =
        csv.empty() ? std::string{} : resolved_path(path, csv);
    const std::string exodus = read_optional_string(outputs, "exodus", {});
    result.outputs.exodus_file =
        exodus.empty() ? std::string{} : resolved_path(path, exodus);
    result.outputs.exodus_interval =
        read_optional_size(document, outputs, "exodus_interval", 1);
    const std::string history = read_optional_string(outputs, "history", {});
    result.outputs.history_file =
        history.empty() ? std::string{} : resolved_path(path, history);
    result.outputs.history_interval =
        read_optional_size(document, outputs, "history_interval", 1);
    result.outputs.progress_interval =
        read_optional_size(document, outputs, "progress_interval", 1);
    const std::string checkpoint =
        read_optional_string(outputs, "checkpoint", {});
    result.outputs.checkpoint_file =
        checkpoint.empty() ? std::string{} : resolved_path(path, checkpoint);
    result.outputs.checkpoint_interval =
        read_optional_size(document, outputs, "checkpoint_interval", 1);
    const std::string input_file =
        std::filesystem::absolute(path).lexically_normal().string();
    const std::vector<std::pair<std::string, std::string>> protected_inputs = {
        {"input card", input_file},
        {"input mesh", result.mesh_file},
        {"restart checkpoint", result.transient_execution.restart_file},
    };
    const std::vector<std::pair<std::string, std::string>> output_paths = {
        {"CSV output", result.outputs.csv_file},
        {"Exodus results", result.outputs.exodus_file},
        {"engineering history", result.outputs.history_file},
        {"checkpoint output", result.outputs.checkpoint_file},
    };
    for (const auto& output_path : output_paths) {
        if (output_path.second.empty())
            continue;
        for (const auto& protected_input : protected_inputs) {
            if (!protected_input.second.empty() &&
                output_path.second == protected_input.second)
                throw std::invalid_argument(
                    path + ": " + output_path.first +
                    " must not overwrite the " + protected_input.first);
        }
    }
    for (std::size_t first = 0; first < output_paths.size(); ++first) {
        if (output_paths[first].second.empty())
            continue;
        for (std::size_t second = first + 1; second < output_paths.size();
             ++second) {
            if (!output_paths[second].second.empty() &&
                output_paths[first].second == output_paths[second].second)
                throw std::invalid_argument(
                    path + ": " + output_paths[first].first + " and " +
                    output_paths[second].first + " paths must differ");
        }
    }
    if (result.problem == CaseProblem::steady &&
        (!result.outputs.history_file.empty() ||
         find_entry(outputs, "history_interval") != nullptr ||
         find_entry(outputs, "progress_interval") != nullptr ||
         find_entry(outputs, "exodus_interval") != nullptr ||
         !result.outputs.checkpoint_file.empty() ||
         find_entry(outputs, "checkpoint_interval") != nullptr))
        throw std::invalid_argument(
            path + ": transient output controls are only valid for transient "
                   "cases");
    if (result.outputs.exodus_file.empty() &&
        find_entry(outputs, "exodus_interval") != nullptr)
        value_error(document, outputs.entry("exodus_interval"),
                    "exodus_interval requires exodus");
    if (result.outputs.history_file.empty() &&
        find_entry(outputs, "history_interval") != nullptr)
        value_error(document, outputs.entry("history_interval"),
                    "history_interval requires history");
    if (result.outputs.checkpoint_file.empty() &&
        find_entry(outputs, "checkpoint_interval") != nullptr)
        value_error(document, outputs.entry("checkpoint_interval"),
                    "checkpoint_interval requires checkpoint");
    if (!result.outputs.checkpoint_file.empty() &&
        result.outputs.checkpoint_interval == 0)
        value_error(document, outputs.entry("checkpoint_interval"),
                    "checkpoint_interval must be positive");
    if (result.outputs.exodus_interval == 0 ||
        result.outputs.history_interval == 0 ||
        result.outputs.progress_interval == 0)
        throw std::invalid_argument(
            path + ": output intervals must be positive");
    return result;
}

} // namespace fuelsim
