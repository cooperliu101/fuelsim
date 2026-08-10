#include "fuelsim/case_input.hpp"
#include "fuelsim/input_file.hpp"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {
namespace {

std::string trim(const std::string& value) {
    const auto first = std::find_if_not(
        value.begin(), value.end(),
        [](unsigned char character) { return std::isspace(character) != 0; });
    if (first == value.end())
        return {};
    const auto last = std::find_if_not(
        value.rbegin(), value.rend(),
        [](unsigned char character) { return std::isspace(character) != 0; });
    return std::string(first, last.base());
}

[[noreturn]] void input_error(const std::string& path, std::size_t line,
                              const std::string& message) {
    throw std::invalid_argument(path + ":" + std::to_string(line) + ": " +
                                message);
}

std::string strip_comment(const std::string& line, const std::string& path,
                          std::size_t line_number) {
    char quote = '\0';
    for (std::size_t index = 0; index < line.size(); ++index) {
        const char character = line[index];
        if (quote == '\0' && (character == '\'' || character == '"')) {
            quote = character;
            continue;
        }
        if (quote != '\0' && character == quote) {
            quote = '\0';
            continue;
        }
        if (quote == '\0' && character == '#')
            return line.substr(0, index);
    }
    if (quote != '\0')
        input_error(path, line_number, "unterminated quoted value");
    return line;
}

bool valid_name(const std::string& name) {
    if (name.empty() ||
        !(std::isalpha(static_cast<unsigned char>(name.front())) != 0 ||
          name.front() == '_'))
        return false;
    return std::all_of(name.begin() + 1, name.end(), [](unsigned char value) {
        return std::isalnum(value) != 0 || value == '_';
    });
}

std::string section_path(const std::vector<std::string>& stack) {
    std::string result;
    for (const std::string& name : stack) {
        if (!result.empty())
            result += '/';
        result += name;
    }
    return result;
}

std::string normalized_value(const std::string& raw, const std::string& path,
                             std::size_t line_number) {
    std::string value = trim(raw);
    if (value.empty())
        input_error(path, line_number, "input value must not be empty");
    const bool starts_quoted = value.front() == '\'' || value.front() == '"';
    const bool ends_quoted = value.back() == '\'' || value.back() == '"';
    if (starts_quoted || ends_quoted) {
        if (value.size() < 2 || value.front() != value.back())
            input_error(path, line_number, "mismatched value quotes");
        value = value.substr(1, value.size() - 2);
    }
    return value;
}

} // namespace

const std::string& InputSection::path() const noexcept {
    return _path;
}

std::size_t InputSection::line() const noexcept {
    return _line;
}

const std::vector<InputEntry>& InputSection::entries() const noexcept {
    return _entries;
}

const InputEntry& InputSection::entry(const std::string& key) const {
    const auto found = std::find_if(
        _entries.begin(), _entries.end(),
        [&key](const InputEntry& candidate) { return candidate.key == key; });
    if (found == _entries.end())
        throw std::invalid_argument("Input section [" + _path +
                                    "] is missing required key '" + key + "'");
    return *found;
}

const std::string& InputDocument::source_path() const noexcept {
    return _source_path;
}

const std::vector<InputSection>& InputDocument::sections() const noexcept {
    return _sections;
}

const InputSection& InputDocument::section(const std::string& path) const {
    const auto found = std::find_if(_sections.begin(), _sections.end(),
                                    [&path](const InputSection& candidate) {
                                        return candidate.path() == path;
                                    });
    if (found == _sections.end())
        throw std::invalid_argument(
            _source_path + ": missing required section [" + path + "]");
    return *found;
}

InputDocument InputParser::parse_file(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not open fuelsim input file '" + path +
                                 "'");

    InputDocument document;
    document._source_path = path;
    std::vector<std::string> stack;
    InputSection* current = nullptr;
    std::string raw_line;
    std::size_t line_number = 0;
    while (std::getline(input, raw_line)) {
        ++line_number;
        const std::string line =
            trim(strip_comment(raw_line, path, line_number));
        if (line.empty())
            continue;

        if (line.front() == '[') {
            if (line.back() != ']')
                input_error(path, line_number,
                            "section header must end with ']'");
            const std::string name = trim(line.substr(1, line.size() - 2));
            if (name.empty()) {
                if (stack.empty())
                    input_error(path, line_number,
                                "unexpected section terminator []");
                stack.pop_back();
                if (stack.empty()) {
                    current = nullptr;
                } else {
                    const std::string parent_path = section_path(stack);
                    current = nullptr;
                    for (InputSection& section : document._sections) {
                        if (section.path() == parent_path) {
                            current = &section;
                            break;
                        }
                    }
                    if (current == nullptr)
                        throw std::logic_error(
                            "InputParser lost its parent section");
                }
                continue;
            }
            if (!valid_name(name))
                input_error(path, line_number,
                            "invalid section name '" + name + "'");
            stack.push_back(name);
            const std::string full_path = section_path(stack);
            const bool duplicate = std::any_of(
                document._sections.begin(), document._sections.end(),
                [&full_path](const InputSection& section) {
                    return section.path() == full_path;
                });
            if (duplicate)
                input_error(path, line_number,
                            "duplicate section [" + full_path + "]");
            document._sections.push_back({});
            current = &document._sections.back();
            current->_path = full_path;
            current->_line = line_number;
            continue;
        }

        if (current == nullptr)
            input_error(path, line_number,
                        "key-value entry must be inside a section");
        const std::size_t equals = line.find('=');
        if (equals == std::string::npos)
            input_error(path, line_number,
                        "expected a key-value entry containing '='");
        const std::string key = trim(line.substr(0, equals));
        if (!valid_name(key))
            input_error(path, line_number, "invalid key '" + key + "'");
        const bool duplicate = std::any_of(
            current->_entries.begin(), current->_entries.end(),
            [&key](const InputEntry& entry) { return entry.key == key; });
        if (duplicate)
            input_error(path, line_number,
                        "duplicate key '" + key + "' in [" + current->_path +
                            "]");
        current->_entries.push_back(
            {key, normalized_value(line.substr(equals + 1), path, line_number),
             line_number});
    }

    if (!input.eof())
        throw std::runtime_error("Could not read fuelsim input file '" + path +
                                 "'");
    if (!stack.empty())
        input_error(path, line_number,
                    "section [" + section_path(stack) +
                        "] is missing its closing []");
    if (document._sections.empty())
        throw std::invalid_argument(path + ": input file contains no sections");
    return document;
}

// Case-definition validation and translation.
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

void forbid_keys(const InputDocument& document, const InputSection& section,
                 const std::vector<std::string>& keys,
                 const std::string& context) {
    for (const std::string& key : keys)
        forbid_key(document, section, key, context);
}

void forbid_convection_keys(const InputDocument& document,
                            const InputSection& section,
                            const std::string& context) {
    forbid_keys(document, section,
                {"heat_transfer_coefficient", "ambient_temperature",
                 "coefficient_function", "ambient_temperature_function"},
                context);
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

NortonCreepProperties read_creep(const InputDocument& document,
                                 const InputSection& section) {
    return {read_double(document, section, "creep_coefficient"),
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
}

J2PlasticityProperties read_plasticity(const InputDocument& document,
                                       const InputSection& section) {
    return {read_double(document, section, "yield_stress"),
            read_double(document, section, "hardening_modulus"),
            read_optional_double(document, section,
                                 "yield_stress_temperature_coefficient", 0.0),
            read_optional_double(document, section,
                                 "hardening_temperature_coefficient", 0.0)};
}

TransientInelasticProperties read_transient(const InputDocument& document,
                                            const InputSection& section) {
    const std::string model = read_string(document, section, "inelastic_model");
    InelasticBehavior behavior;
    if (model == "elastic")
        behavior = InelasticBehavior::elastic;
    else if (model == "norton_creep")
        behavior = InelasticBehavior::norton_creep;
    else if (model == "j2_plasticity")
        behavior = InelasticBehavior::j2_plasticity;
    else if (model == "norton_creep_j2_plasticity")
        behavior = InelasticBehavior::norton_creep_j2_plasticity;
    else
        value_error(document,
                    required_entry(document, section, "inelastic_model"),
                    "unknown inelastic_model '" + model + "'");

    const bool uses_creep =
        behavior == InelasticBehavior::norton_creep ||
        behavior == InelasticBehavior::norton_creep_j2_plasticity;
    const bool uses_plasticity =
        behavior == InelasticBehavior::j2_plasticity ||
        behavior == InelasticBehavior::norton_creep_j2_plasticity;
    const std::string context =
        model == "elastic" ? "inelastic_model='elastic'" : model;
    NortonCreepProperties creep{0.0, 1.0, 1.0};
    J2PlasticityProperties plasticity{1.0, 0.0};
    if (uses_creep)
        creep = read_creep(document, section);
    else
        forbid_keys(document, section,
                    {"creep_coefficient", "creep_reference_stress",
                     "creep_exponent",
                     "creep_coefficient_temperature_coefficient",
                     "creep_reference_stress_temperature_coefficient",
                     "creep_exponent_temperature_coefficient"},
                    context);
    if (uses_plasticity)
        plasticity = read_plasticity(document, section);
    else
        forbid_keys(document, section,
                    {"yield_stress", "hardening_modulus",
                     "yield_stress_temperature_coefficient",
                     "hardening_temperature_coefficient"},
                    context);
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

std::string read_optional_path(const std::string& input_path,
                               const InputSection& section,
                               const std::string& key) {
    const std::string value = read_optional_string(section, key, {});
    return value.empty() ? std::string{} : resolved_path(input_path, value);
}

CaseRegionDefinition read_region(const InputDocument& document,
                                 const InputSection& section,
                                 CaseProblem problem) {
    std::vector<std::string> keys = {
        "block",
        "block_id",
        "strain",
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
    const InputEntry strain = required_entry(document, section, "strain");
    if (strain.value == "small")
        result.spatial.strain_formulation = StrainFormulation::small;
    else if (strain.value == "finite")
        result.spatial.strain_formulation = StrainFormulation::finite;
    else
        value_error(document, strain,
                    "unknown strain formulation '" + strain.value + "'");
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
                             1.0,
                             0.0};
    if (thermal != nullptr) {
        validate_keys(document, *thermal, {"gap_conductivity", "minimum_gap"});
        result.gap_conductivity =
            read_double(document, *thermal, "gap_conductivity");
        result.minimum_gap = read_double(document, *thermal, "minimum_gap");
    }
    if (mechanical != nullptr) {
        validate_keys(document, *mechanical,
                      {"formulation", "penalty", "penalty_factor", "mu",
                       "penetration_tolerance",
                       "maximum_augmented_iterations"});
        const std::string formulation =
            read_string(document, *mechanical, "formulation");
        if (formulation == "penalty")
            result.mechanical_formulation =
                MechanicalContactFormulation::penalty;
        else if (formulation == "augmented_lagrangian")
            result.mechanical_formulation =
                MechanicalContactFormulation::augmented_lagrangian;
        else
            value_error(document, mechanical->entry("formulation"),
                        "mechanical formulation must be 'penalty' or "
                        "'augmented_lagrangian'");
        const InputEntry* penalty = find_entry(*mechanical, "penalty");
        const InputEntry* penalty_factor =
            find_entry(*mechanical, "penalty_factor");
        if (penalty != nullptr && penalty_factor != nullptr)
            value_error(document, *penalty_factor,
                        "penalty and penalty_factor are mutually exclusive");
        result.automatic_penalty = penalty == nullptr;
        result.penalty =
            penalty == nullptr ? 0.0 : parse_double(document, *penalty);
        result.penalty_factor =
            penalty_factor == nullptr
                ? 1.0
                : parse_double(document, *penalty_factor);
        result.friction_coefficient =
            read_optional_double(document, *mechanical, "mu", 0.0);
        if (result.friction_coefficient < 0.0)
            value_error(document, mechanical->entry("mu"),
                        "friction coefficient mu must be nonnegative");
        const InputEntry* penetration_tolerance =
            find_entry(*mechanical, "penetration_tolerance");
        const InputEntry* maximum_augmented_iterations =
            find_entry(*mechanical, "maximum_augmented_iterations");
        if (result.mechanical_formulation ==
            MechanicalContactFormulation::penalty) {
            if (penetration_tolerance != nullptr)
                value_error(document, *penetration_tolerance,
                            "penetration_tolerance requires "
                            "formulation = augmented_lagrangian");
            if (maximum_augmented_iterations != nullptr)
                value_error(document, *maximum_augmented_iterations,
                            "maximum_augmented_iterations requires "
                            "formulation = augmented_lagrangian");
        } else {
            result.penetration_tolerance = read_optional_double(
                document, *mechanical, "penetration_tolerance", 1.0e-8);
            result.maximum_augmented_iterations = read_optional_size(
                document, *mechanical, "maximum_augmented_iterations", 20);
        }
    }
    return result;
}

BoundaryConditionDefinition make_boundary_condition(
    const InputDocument& document, const InputSection& section,
    BoundaryConditionType type, Field field, double value,
    bool scale_with_load, const std::string& function) {
    BoundaryConditionDefinition result{
        leaf_name(section), type,
        read_string(document, section, "boundary"), field, value,
        scale_with_load};
    result.function = function;
    return result;
}

BoundaryConditionDefinition
read_boundary_condition(const InputDocument& document,
                        const InputSection& section) {
    validate_keys(document, section,
                  {"type", "boundary", "field", "value", "scale_with_load",
                   "function", "heat_transfer_coefficient",
                   "ambient_temperature", "coefficient_function",
                   "ambient_temperature_function", "configuration"});
    const std::string type = read_string(document, section, "type");
    const bool scale_with_load =
        read_optional_bool(document, section, "scale_with_load", false);
    const std::string function = read_optional_string(section, "function", {});
    if (scale_with_load && !function.empty())
        value_error(document, section.entry("scale_with_load"),
                    "scale_with_load cannot be combined with function");
    if (type == "dirichlet") {
        forbid_key(document, section, "configuration", "type='dirichlet'");
        forbid_convection_keys(document, section, "type='dirichlet'");
        return make_boundary_condition(
            document, section, BoundaryConditionType::dirichlet,
            parse_field(document, required_entry(document, section, "field")),
            read_double(document, section, "value"), scale_with_load,
            function);
    }
    if (type == "pressure") {
        forbid_key(document, section, "configuration", "type='pressure'");
        forbid_key(document, section, "field", "type='pressure'");
        forbid_convection_keys(document, section, "type='pressure'");
        return make_boundary_condition(
            document, section, BoundaryConditionType::pressure,
            Field::radial_displacement,
            read_double(document, section, "value"), scale_with_load,
            function);
    }
    if (type == "traction") {
        forbid_convection_keys(document, section, "type='traction'");
        const Field field =
            parse_field(document, required_entry(document, section, "field"));
        if (field == Field::temperature)
            value_error(document, section.entry("field"),
                        "traction requires a displacement field");
        BoundaryConditionDefinition result = make_boundary_condition(
            document, section, BoundaryConditionType::traction, field,
            read_double(document, section, "value"), scale_with_load,
            function);
        const std::string configuration =
            read_optional_string(section, "configuration", "reference");
        if (configuration != "reference" && configuration != "current")
            value_error(document, section.entry("configuration"),
                        "traction configuration must be reference or current");
        result.use_displaced_geometry = configuration == "current";
        return result;
    }
    if (type == "convection") {
        forbid_key(document, section, "configuration", "type='convection'");
        forbid_key(document, section, "field", "type='convection'");
        forbid_key(document, section, "value", "type='convection'");
        forbid_key(document, section, "scale_with_load", "type='convection'");
        forbid_key(document, section, "function", "type='convection'");
        BoundaryConditionDefinition result = make_boundary_condition(
            document, section, BoundaryConditionType::convection,
            Field::temperature, 0.0, false, {});
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

SpatialDefinition FuelSimCaseDefinition::spatial_definition() const {
    SpatialDefinition result;
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
    result.spatial = spatial_definition();
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
                       "time_error_safety_factor",
                       "strain_history_time_absolute_tolerance",
                       "stress_history_time_absolute_tolerance"});
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
                                 "time_error_safety_factor", 0.9),
            read_optional_double(
                document, executioner,
                "strain_history_time_absolute_tolerance", 1.0e-10),
            read_optional_double(
                document, executioner,
                "stress_history_time_absolute_tolerance", 1.0)};
        result.transient_execution.restart_file =
            read_optional_path(path, executioner, "restart");
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
            !(result.transient_execution
                      .strain_history_time_absolute_tolerance > 0.0) ||
            !(result.transient_execution
                      .stress_history_time_absolute_tolerance > 0.0) ||
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
                   "mechanical_residual_absolute_tolerance",
                   "temperature_residual_scale",
                   "mechanical_residual_scale"});
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
                             1.0e-4),
        read_optional_double(document, solver, "temperature_residual_scale",
                             0.0),
        read_optional_double(document, solver, "mechanical_residual_scale",
                             0.0)};
    if (!(result.solver.absolute_tolerance > 0.0) ||
        !(result.solver.relative_tolerance > 0.0) ||
        !(result.solver.step_tolerance > 0.0) ||
        result.solver.maximum_iterations <= 0 ||
        !(result.solver.linear_relative_tolerance > 0.0) ||
        !(result.solver.residual_reduction_tolerance > 0.0) ||
        !(result.solver.temperature_residual_absolute_tolerance > 0.0) ||
        !(result.solver.mechanical_residual_absolute_tolerance > 0.0) ||
        result.solver.temperature_residual_scale < 0.0 ||
        result.solver.mechanical_residual_scale < 0.0 ||
        result.solver.maximum_linear_iterations <= 0)
        throw std::invalid_argument(
            path + ": solver tolerances and iteration limit must be positive");
    const bool fixed_temperature_scale =
        result.solver.temperature_residual_scale > 0.0;
    const bool fixed_mechanical_scale =
        result.solver.mechanical_residual_scale > 0.0;
    if (fixed_temperature_scale != fixed_mechanical_scale)
        throw std::invalid_argument(
            path + ": temperature_residual_scale and "
                   "mechanical_residual_scale must both be zero or positive");
    if (fixed_temperature_scale && result.solver.field_residual_scaling)
        throw std::invalid_argument(
            path + ": fixed residual scales cannot be combined with "
                   "field_residual_scaling");
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
    result.outputs.csv_file = read_optional_path(path, outputs, "csv");
    result.outputs.exodus_file = read_optional_path(path, outputs, "exodus");
    result.outputs.exodus_interval =
        read_optional_size(document, outputs, "exodus_interval", 1);
    result.outputs.history_file = read_optional_path(path, outputs, "history");
    result.outputs.history_interval =
        read_optional_size(document, outputs, "history_interval", 1);
    result.outputs.progress_interval =
        read_optional_size(document, outputs, "progress_interval", 1);
    result.outputs.checkpoint_file =
        read_optional_path(path, outputs, "checkpoint");
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
