#include "fuelsim/io/case_input.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
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
struct InputEntry final {
    std::string key, value;
    std::size_t line;
};

struct InputSection final {
    std::string name, path, parent;
    std::size_t line = 0;
    std::vector<InputEntry> entries;
};

struct InputDocument final {
    std::string source_path;
    std::vector<InputSection> sections;
};

std::string trim(const std::string& value) {
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char character) {
        return std::isspace(character) != 0;
    });
    if (first == value.end())
        return {};
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char character) {
        return std::isspace(character) != 0;
    });
    return std::string(first, last.base());
}

[[noreturn]] void input_error(const std::string& path, std::size_t line, const std::string& message) {
    throw std::invalid_argument(path + ":" + std::to_string(line) + ": " + message);
}

std::string strip_comment(const std::string& line, const std::string& path, std::size_t line_number) {
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
    if (name.empty() || !(std::isalpha(static_cast<unsigned char>(name.front())) != 0 || name.front() == '_'))
        return false;
    return std::all_of(name.begin() + 1, name.end(), [](unsigned char value) {
        return std::isalnum(value) != 0 || value == '_';
    });
}

std::string normalized_value(const std::string& raw, const std::string& path, std::size_t line_number) {
    std::string value = trim(raw);
    if (value.empty())
        input_error(path, line_number, "input value must not be empty");
    const bool starts_quoted = value.front() == '\'' || value.front() == '"',
               ends_quoted = value.back() == '\'' || value.back() == '"';
    if (starts_quoted || ends_quoted) {
        if (value.size() < 2 || value.front() != value.back())
            input_error(path, line_number, "mismatched value quotes");
        value = value.substr(1, value.size() - 2);
    }
    return value;
}

InputDocument parse_input_file(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not open fuelsim input file '" + path + "'");
    InputDocument document;
    document.source_path = path;
    std::vector<std::size_t> stack;
    InputSection* current = nullptr;
    std::string raw_line;
    std::size_t line_number = 0;
    while (std::getline(input, raw_line)) {
        ++line_number;
        const std::string line = trim(strip_comment(raw_line, path, line_number));
        if (line.empty())
            continue;
        if (line.front() == '[') {
            if (line.back() != ']')
                input_error(path, line_number, "section header must end with ']'");
            const std::string name = trim(line.substr(1, line.size() - 2));
            if (name.empty()) {
                if (stack.empty())
                    input_error(path, line_number, "unexpected section terminator []");
                stack.pop_back();
                current = stack.empty() ? nullptr : &document.sections[stack.back()];
                continue;
            }
            if (!valid_name(name))
                input_error(path, line_number, "invalid section name '" + name + "'");
            const std::string full_path = current == nullptr ? name : current->path + '/' + name;
            const bool duplicate = std::any_of(document.sections.begin(),
                document.sections.end(),
                [&full_path](const InputSection& section) { return section.path == full_path; });
            if (duplicate)
                input_error(path, line_number, "duplicate section [" + full_path + "]");
            document.sections.push_back({});
            stack.push_back(document.sections.size() - 1);
            current = &document.sections.back();
            current->name = name;
            current->path = full_path;
            if (stack.size() > 1)
                current->parent = document.sections[stack[stack.size() - 2]].path;
            current->line = line_number;
            continue;
        }
        if (current == nullptr)
            input_error(path, line_number, "key-value entry must be inside a section");
        const std::size_t equals = line.find('=');
        if (equals == std::string::npos)
            input_error(path, line_number, "expected a key-value entry containing '='");
        const std::string key = trim(line.substr(0, equals));
        if (!valid_name(key))
            input_error(path, line_number, "invalid key '" + key + "'");
        const bool duplicate = std::any_of(current->entries.begin(),
            current->entries.end(),
            [&key](const InputEntry& entry) { return entry.key == key; });
        if (duplicate)
            input_error(path, line_number, "duplicate key '" + key + "' in [" + current->path + "]");
        current->entries.push_back({key, normalized_value(line.substr(equals + 1), path, line_number), line_number});
    }
    if (!input.eof())
        throw std::runtime_error("Could not read fuelsim input file '" + path + "'");
    if (!stack.empty())
        input_error(path,
            line_number,
            "section [" + document.sections[stack.back()].path + "] is missing its closing []");
    if (document.sections.empty())
        throw std::invalid_argument(path + ": input file contains no sections");
    return document;
}

[[noreturn]] void value_error(const InputDocument& document, const InputEntry& entry, const std::string& message) {
    throw std::invalid_argument(document.source_path + ":" + std::to_string(entry.line) + ": " + message);
}

const InputEntry* find_entry(const InputSection& section, const std::string& key) {
    const auto found = std::find_if(section.entries.begin(), section.entries.end(), [&key](const InputEntry& entry) {
        return entry.key == key;
    });
    return found == section.entries.end() ? nullptr : &*found;
}

const InputEntry& required_entry(const InputDocument& document, const InputSection& section, const std::string& key) {
    const InputEntry* entry = find_entry(section, key);
    if (entry == nullptr)
        throw std::invalid_argument(document.source_path + ":" + std::to_string(section.line) + ": section ["
                                    + section.path + "] is missing required key '" + key + "'");
    return *entry;
}

const InputSection* find_section(const InputDocument& document, const std::string& path) {
    const auto found = std::find_if(document.sections.begin(),
        document.sections.end(),
        [&path](const InputSection& section) { return section.path == path; });
    return found == document.sections.end() ? nullptr : &*found;
}

const InputSection& required_section(const InputDocument& document, const char* path) {
    const InputSection* section = find_section(document, path);
    if (section == nullptr)
        throw std::invalid_argument(document.source_path + ": missing required section [" + path + "]");
    return *section;
}

std::vector<const InputSection*> direct_children(const InputDocument& document, const std::string& parent) {
    std::vector<const InputSection*> result;
    for (const InputSection& section : document.sections)
        if (section.parent == parent)
            result.push_back(&section);
    return result;
}

void validate_sections(const InputDocument& document) {
    const std::vector<std::string> fixed = {"Case",
        "Mesh",
        "TimeFunctions",
        "Materials",
        "Regions",
        "Contact",
        "BoundaryConditions",
        "Executioner",
        "Solver",
        "Outputs"};
    for (const InputSection& section : document.sections) {
        bool known = section.parent.empty() && std::find(fixed.begin(), fixed.end(), section.name) != fixed.end();
        known = known || section.parent == "Regions" || section.parent == "BoundaryConditions"
                || section.parent == "TimeFunctions" || section.parent == "Materials" || section.parent == "Contact";
        if (section.parent.compare(0, 10, "Materials/") == 0) {
            const std::string material_path = section.parent.substr(10);
            if (material_path.find('/') == std::string::npos)
                known = section.name == "thermal" || section.name == "elasticity" || section.name == "eigenstrains"
                        || section.name == "creep" || section.name == "plasticity";
            else
                known = section.parent.compare(section.parent.size() - 13, 13, "/eigenstrains") == 0
                        && material_path.find('/') == material_path.rfind('/');
        } else if (section.parent.compare(0, 8, "Contact/") == 0) {
            known = section.parent.find('/', 8) == std::string::npos
                    && (section.name == "thermal" || section.name == "mechanical");
        }
        if (!known)
            throw std::invalid_argument(
                document.source_path + ":" + std::to_string(section.line) + ": unknown section [" + section.path + "]");
    }
}

void validate_keys(const InputDocument& document,
    const InputSection& section,
    const std::vector<std::string>& allowed) {
    for (const InputEntry& entry : section.entries)
        if (std::find(allowed.begin(), allowed.end(), entry.key) == allowed.end())
            value_error(document, entry, "unknown key '" + entry.key + "' in [" + section.path + "]");
}

std::string read_string(const InputDocument& document, const InputSection& section, const std::string& key) {
    return required_entry(document, section, key).value;
}

double parse_double(const InputDocument& document, const InputEntry& entry) {
    errno = 0;
    char* end = nullptr;
    const double result = std::strtod(entry.value.c_str(), &end);
    if (errno == ERANGE || end == entry.value.c_str() || *end != '\0' || !std::isfinite(result))
        value_error(document, entry, "key '" + entry.key + "' requires a finite real value");
    return result;
}

double read_double(const InputDocument& document, const InputSection& section, const std::string& key) {
    return parse_double(document, required_entry(document, section, key));
}

double read_optional_double(const InputDocument& document,
    const InputSection& section,
    const std::string& key,
    double fallback) {
    const InputEntry* entry = find_entry(section, key);
    return entry == nullptr ? fallback : parse_double(document, *entry);
}

std::vector<double> parse_double_list(const InputDocument& document, const InputEntry& entry) {
    std::istringstream input(entry.value);
    std::vector<double> result;
    double value = 0.0;
    while (input >> value) {
        if (!std::isfinite(value))
            value_error(document, entry, "key '" + entry.key + "' requires finite real values");
        result.push_back(value);
    }
    if (!input.eof() || result.empty())
        value_error(document, entry, "key '" + entry.key + "' requires a whitespace-separated real list");
    return result;
}

std::size_t parse_size(const InputDocument& document, const InputEntry& entry) {
    if (!entry.value.empty() && entry.value.front() == '-')
        value_error(document, entry, "key '" + entry.key + "' requires a nonnegative integer value");
    errno = 0;
    char* end = nullptr;
    const unsigned long long result = std::strtoull(entry.value.c_str(), &end, 10);
    if (errno == ERANGE || end == entry.value.c_str() || *end != '\0'
        || result > std::numeric_limits<std::size_t>::max())
        value_error(document, entry, "key '" + entry.key + "' requires a nonnegative integer value");
    return static_cast<std::size_t>(result);
}

std::size_t read_size(const InputDocument& document, const InputSection& section, const std::string& key) {
    return parse_size(document, required_entry(document, section, key));
}

std::size_t read_optional_size(const InputDocument& document,
    const InputSection& section,
    const std::string& key,
    std::size_t fallback) {
    const InputEntry* entry = find_entry(section, key);
    return entry == nullptr ? fallback : parse_size(document, *entry);
}

int read_optional_int(const InputDocument& document,
    const InputSection& section,
    const std::string& key,
    int fallback) {
    const InputEntry* entry = find_entry(section, key);
    if (entry == nullptr)
        return fallback;
    const std::size_t result = parse_size(document, *entry);
    if (result > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        value_error(document, *entry, "key '" + key + "' is outside the supported range");
    return static_cast<int>(result);
}

bool read_optional_bool(const InputDocument& document,
    const InputSection& section,
    const std::string& key,
    bool fallback) {
    const InputEntry* entry = find_entry(section, key);
    if (entry == nullptr)
        return fallback;
    if (entry->value == "true")
        return true;
    if (entry->value == "false")
        return false;
    value_error(document, *entry, "key '" + key + "' requires 'true' or 'false'");
}

std::string read_optional_string(const InputSection& section, const std::string& key, const std::string& fallback) {
    const InputEntry* entry = find_entry(section, key);
    return entry == nullptr ? fallback : entry->value;
}

void forbid_key(const InputDocument& document,
    const InputSection& section,
    const std::string& key,
    const std::string& context) {
    const InputEntry* entry = find_entry(section, key);
    if (entry != nullptr)
        value_error(document, *entry, "key '" + key + "' is not valid for " + context);
}

void forbid_keys(const InputDocument& document,
    const InputSection& section,
    const std::vector<std::string>& keys,
    const std::string& context) {
    for (const std::string& key : keys)
        forbid_key(document, section, key, context);
}

void forbid_convection_keys(const InputDocument& document, const InputSection& section, const std::string& context) {
    forbid_keys(document,
        section,
        {"heat_transfer_coefficient", "ambient_temperature", "coefficient_function", "ambient_temperature_function"},
        context);
}

struct ParsedMaterial final {
    std::string name;
    std::shared_ptr<const MaterialFunctionSet> functions;
};

std::vector<MaterialParameterValue> read_material_parameters(const InputDocument& document,
    const InputSection& section) {
    std::vector<MaterialParameterValue> values;
    values.reserve(section.entries.size() - 1);
    for (const InputEntry& entry : section.entries)
        if (entry.key != "function")
            values.push_back({entry.key, parse_double(document, entry)});
    return values;
}

std::vector<ParsedMaterial> read_materials(const InputDocument& document, const MaterialFunctionRegistry& registry) {
    const InputSection& materials = required_section(document, "Materials");
    validate_keys(document, materials, {});
    std::vector<ParsedMaterial> result;
    for (const InputSection* material_section : direct_children(document, "Materials")) {
        validate_keys(document, *material_section, {});
        const std::string material_name = material_section->name;
        const std::string base = material_section->path;
        const InputSection* thermal = find_section(document, base + "/thermal");
        const InputSection* elasticity = find_section(document, base + "/elasticity");
        if (thermal == nullptr || elasticity == nullptr)
            throw std::invalid_argument(document.source_path + ":" + std::to_string(material_section->line)
                                        + ": material [" + base + "] requires [thermal] and [elasticity]");
        auto functions = std::make_shared<MaterialFunctionSet>();
        functions->name = material_name;
        try {
            const std::string thermal_name = read_string(document, *thermal, "function");
            functions->thermal = registry.bind_thermal(thermal_name, read_material_parameters(document, *thermal));
            const std::string elasticity_name = read_string(document, *elasticity, "function");
            functions->elasticity =
                registry.bind_elasticity(elasticity_name, read_material_parameters(document, *elasticity));
            const InputSection* eigenstrains = find_section(document, base + "/eigenstrains");
            if (eigenstrains != nullptr) {
                validate_keys(document, *eigenstrains, {});
                for (const InputSection* section : direct_children(document, eigenstrains->path)) {
                    const std::string function_name = read_string(document, *section, "function");
                    functions->eigenstrains.push_back(registry.bind_eigenstrain(section->name,
                        function_name,
                        read_material_parameters(document, *section)));
                }
            }
            const InputSection* creep = find_section(document, base + "/creep");
            if (creep != nullptr) {
                const std::string function_name = read_string(document, *creep, "function");
                functions->creep = registry.bind_creep(function_name, read_material_parameters(document, *creep));
            }
            const InputSection* plasticity = find_section(document, base + "/plasticity");
            if (plasticity != nullptr) {
                const std::string function_name = read_string(document, *plasticity, "function");
                functions->plasticity =
                    registry.bind_plasticity(function_name, read_material_parameters(document, *plasticity));
            }
        } catch (const std::invalid_argument& error) {
            throw std::invalid_argument(document.source_path + ": material '" + material_name + "': " + error.what());
        }
        result.push_back({material_name, std::move(functions)});
    }
    if (result.empty())
        throw std::invalid_argument(document.source_path + ": [Materials] requires a child material");
    return result;
}

Field parse_field(const InputDocument& document, const InputEntry& entry) {
    if (entry.value == "temperature")
        return Field::temperature;
    if (entry.value == "radial_displacement")
        return Field::radial_displacement;
    if (entry.value == "axial_displacement")
        return Field::axial_displacement;
    if (entry.value == "displacement_x")
        return Field::displacement_x;
    if (entry.value == "displacement_y")
        return Field::displacement_y;
    if (entry.value == "displacement_z")
        return Field::displacement_z;
    value_error(document, entry, "unknown field '" + entry.value + "'");
}

std::string resolved_path(const std::string& input_path, const std::string& configured_path) {
    const std::filesystem::path configured(configured_path);
    if (configured.is_absolute())
        return configured.lexically_normal().string();
    return (std::filesystem::absolute(input_path).parent_path() / configured).lexically_normal().string();
}

std::string read_optional_path(const std::string& input_path, const InputSection& section, const std::string& key) {
    const std::string value = read_optional_string(section, key, {});
    return value.empty() ? std::string{} : resolved_path(input_path, value);
}

RegionDefinition read_region(const InputDocument& document,
    const InputSection& section,
    const std::vector<ParsedMaterial>& materials,
    CaseGeometry geometry) {
    validate_keys(document,
        section,
        {"block",
            "block_id",
            "material",
            "strain",
            "element",
            "initial_temperature",
            "volumetric_heat_source",
            "heat_source_function",
            "heat_source_time_evaluation"});
    const InputEntry* block = find_entry(section, "block");
    const InputEntry* block_id = find_entry(section, "block_id");
    if ((block == nullptr) == (block_id == nullptr))
        throw std::invalid_argument(document.source_path + ":" + std::to_string(section.line) + ": region ["
                                    + section.path + "] requires exactly one of 'block' or 'block_id'");
    std::int64_t resolved_block_id = -1;
    if (block_id != nullptr) {
        const std::size_t value = parse_size(document, *block_id);
        if (value > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
            value_error(document, *block_id, "key 'block_id' is outside the supported range");
        resolved_block_id = static_cast<std::int64_t>(value);
    }
    const std::string material_name = read_string(document, section, "material");
    const auto material = std::find_if(materials.begin(), materials.end(), [&](const ParsedMaterial& candidate) {
        return candidate.name == material_name;
    });
    if (material == materials.end())
        value_error(document,
            required_entry(document, section, "material"),
            "unknown material '" + material_name + "'");
    const double initial_temperature = read_double(document, section, "initial_temperature");
    ElasticPropertyOutput initial_elasticity{};
    const ElasticFunctionInstance& elasticity = material->functions->elasticity;
    elasticity.function({initial_temperature, {}}, initial_elasticity);
    if (!std::isfinite(initial_elasticity.young_modulus.value()) || !(initial_elasticity.young_modulus.value() > 0.0))
        value_error(document,
            required_entry(document, section, "material"),
            "material elasticity must produce positive young_modulus at initial_temperature");
    ThermoelasticProperties thermoelastic{material->functions, initial_elasticity.young_modulus.value()};
    RegionDefinition result = {section.name,
        block == nullptr ? std::string{} : block->value,
        std::move(thermoelastic),
        read_double(document, section, "volumetric_heat_source"),
        initial_temperature,
        resolved_block_id};
    result.heat_source_function = read_optional_string(section, "heat_source_function", {});
    const std::string heat_source_time_evaluation =
        read_optional_string(section, "heat_source_time_evaluation", "end_time");
    if (heat_source_time_evaluation == "end_time")
        result.heat_source_time_evaluation = HeatSourceTimeEvaluation::end_time;
    else if (heat_source_time_evaluation == "interval_average")
        result.heat_source_time_evaluation = HeatSourceTimeEvaluation::interval_average;
    else
        value_error(document,
            required_entry(document, section, "heat_source_time_evaluation"),
            "heat_source_time_evaluation must be end_time or interval_average");
    if (result.heat_source_time_evaluation == HeatSourceTimeEvaluation::interval_average
        && result.heat_source_function.empty())
        value_error(document,
            required_entry(document, section, "heat_source_time_evaluation"),
            "heat_source_time_evaluation=interval_average requires heat_source_function");
    const InputEntry strain = required_entry(document, section, "strain");
    if (strain.value == "small")
        result.strain_formulation = StrainFormulation::small;
    else if (strain.value == "finite")
        result.strain_formulation = StrainFormulation::finite;
    else
        value_error(document, strain, "unknown strain formulation '" + strain.value + "'");
    const std::string element =
        read_optional_string(section, "element", geometry == CaseGeometry::axisymmetric_rz ? "quad4" : "c3d8t");
    if (geometry == CaseGeometry::axisymmetric_rz) {
        if (element == "quad4")
            result.rz_element_formulation = RzElementFormulation::quad4;
        else if (element == "cax4t")
            result.rz_element_formulation = RzElementFormulation::cax4t;
        else if (element == "cax4rt")
            result.rz_element_formulation = RzElementFormulation::cax4rt;
        else if (element == "cax8t")
            result.rz_element_formulation = RzElementFormulation::cax8t;
        else if (element == "cax8rt")
            result.rz_element_formulation = RzElementFormulation::cax8rt;
        else
            value_error(document, required_entry(document, section, "element"), "unknown RZ element '" + element + "'");
        return result;
    }
    if (element == "c3d8t")
        result.hex8_element_formulation = Hex8ElementFormulation::c3d8t;
    else if (element == "c3d8rt")
        result.hex8_element_formulation = Hex8ElementFormulation::c3d8rt;
    else if (element == "c3d20rt")
        result.hex20_element_formulation = Hex20ElementFormulation::c3d20rt;
    else
        value_error(document,
            required_entry(document, section, "element"),
            "unknown Cartesian element '" + element + "'");
    return result;
}

ContactDefinition read_contact(const InputDocument& document, const InputSection& section, CaseGeometry geometry) {
    validate_keys(document, section, {"primary", "secondary"});
    const std::string base = section.path;
    const InputSection* thermal = find_section(document, base + "/thermal");
    const InputSection* mechanical = find_section(document, base + "/mechanical");
    if (thermal == nullptr && mechanical == nullptr)
        throw std::invalid_argument(document.source_path + ":" + std::to_string(section.line) + ": contact [" + base
                                    + "] requires [thermal] or [mechanical]");
    ContactDefinition result{section.name,
        read_string(document, section, "primary"),
        read_string(document, section, "secondary"),
        thermal != nullptr,
        mechanical != nullptr,
        1.0,
        1.0,
        1.0,
        0.0};
    if (thermal != nullptr) {
        const std::string discretization = read_optional_string(*thermal, "discretization", "surface_to_surface");
        if (discretization == "node_to_surface") {
            if (geometry != CaseGeometry::axisymmetric_rz)
                value_error(document,
                    required_entry(document, *thermal, "discretization"),
                    "node_to_surface thermal contact currently requires axisymmetric_rz");
            result.thermal_discretization = ThermalContactDiscretization::node_to_surface;
        } else if (discretization != "surface_to_surface")
            value_error(document,
                required_entry(document, *thermal, "discretization"),
                "thermal discretization must be 'node_to_surface' or 'surface_to_surface'");
        const std::string law = read_optional_string(*thermal, "law", "gas_gap");
        if (law == "gas_gap") {
            validate_keys(document, *thermal, {"law", "gap_conductivity", "minimum_gap", "discretization"});
            result.gap_conductivity = read_double(document, *thermal, "gap_conductivity");
            result.minimum_gap = read_double(document, *thermal, "minimum_gap");
        } else if (law == "affine") {
            validate_keys(document,
                *thermal,
                {"law",
                    "conductance",
                    "clearance_derivative",
                    "pressure_derivative",
                    "temperature_derivative",
                    "reference_temperature",
                    "discretization"});
            result.gap_heat_conductance_law = GapHeatConductanceLaw::affine;
            result.gap_conductance = read_double(document, *thermal, "conductance");
            result.gap_conductance_clearance_derivative =
                read_optional_double(document, *thermal, "clearance_derivative", 0.0);
            result.gap_conductance_pressure_derivative =
                read_optional_double(document, *thermal, "pressure_derivative", 0.0);
            result.gap_conductance_temperature_derivative =
                read_optional_double(document, *thermal, "temperature_derivative", 0.0);
            result.gap_conductance_reference_temperature =
                read_optional_double(document, *thermal, "reference_temperature", 0.0);
            if (!(result.gap_conductance >= 0.0))
                value_error(document,
                    required_entry(document, *thermal, "conductance"),
                    "affine thermal contact conductance must be nonnegative at the reference state");
        } else {
            value_error(document,
                required_entry(document, *thermal, "law"),
                "thermal contact law must be 'gas_gap' or 'affine'");
        }
    }
    if (mechanical != nullptr) {
        validate_keys(document,
            *mechanical,
            {"formulation",
                "penalty",
                "penalty_factor",
                "mu",
                "penetration_tolerance",
                "maximum_augmented_iterations",
                "discretization",
                "sliding",
                "quad8_nodal_area_rule",
                "slip_tolerance"});
        const std::string formulation = read_string(document, *mechanical, "formulation");
        if (formulation == "penalty")
            result.mechanical_formulation = MechanicalContactFormulation::penalty;
        else if (formulation == "augmented_lagrangian")
            result.mechanical_formulation = MechanicalContactFormulation::augmented_lagrangian;
        else
            value_error(document,
                required_entry(document, *mechanical, "formulation"),
                "mechanical formulation must be 'penalty' or 'augmented_lagrangian'");
        const std::string discretization = read_optional_string(*mechanical, "discretization", "automatic");
        if (discretization == "automatic")
            result.mechanical_discretization = MechanicalContactDiscretization::automatic;
        else if (discretization == "node_to_surface")
            result.mechanical_discretization = MechanicalContactDiscretization::node_to_surface;
        else if (discretization == "surface_to_surface")
            result.mechanical_discretization = MechanicalContactDiscretization::surface_to_surface;
        else
            value_error(document,
                required_entry(document, *mechanical, "discretization"),
                "mechanical discretization must be 'automatic', 'node_to_surface', or 'surface_to_surface'");
        const std::string sliding = read_optional_string(*mechanical, "sliding", "small");
        if (sliding == "small")
            result.mechanical_sliding = MechanicalContactSliding::small;
        else if (sliding == "finite")
            result.mechanical_sliding = MechanicalContactSliding::finite;
        else
            value_error(document,
                required_entry(document, *mechanical, "sliding"),
                "mechanical sliding must be 'small' or 'finite'");
        if (find_entry(*mechanical, "sliding") != nullptr
            && result.mechanical_discretization != MechanicalContactDiscretization::surface_to_surface
            && !(geometry == CaseGeometry::axisymmetric_rz && sliding == "finite"))
            value_error(document,
                required_entry(document, *mechanical, "sliding"),
                "sliding applies only to discretization = surface_to_surface");
        const InputEntry* penalty = find_entry(*mechanical, "penalty");
        const InputEntry* penalty_factor = find_entry(*mechanical, "penalty_factor");
        if (penalty != nullptr && penalty_factor != nullptr)
            value_error(document, *penalty_factor, "penalty and penalty_factor are mutually exclusive");
        result.automatic_penalty = penalty == nullptr;
        result.penalty = penalty == nullptr ? 0.0 : parse_double(document, *penalty);
        result.penalty_factor = penalty_factor == nullptr ? 1.0 : parse_double(document, *penalty_factor);
        result.friction_coefficient = read_optional_double(document, *mechanical, "mu", 0.0);
        result.friction_slip_tolerance = read_optional_double(document, *mechanical, "slip_tolerance", 0.0);
        if (result.friction_slip_tolerance < 0.0)
            value_error(document,
                required_entry(document, *mechanical, "slip_tolerance"),
                "slip_tolerance must be nonnegative");
        if (result.friction_slip_tolerance > 0.0 && !(result.friction_coefficient > 0.0))
            value_error(document,
                required_entry(document, *mechanical, "slip_tolerance"),
                "slip_tolerance requires a positive mu");
        const std::string area_rule = read_optional_string(*mechanical, "quad8_nodal_area_rule", "positive_lumped");
        if (area_rule == "positive_lumped")
            result.quad8_nodal_area_rule = Quad8NodalAreaRule::positive_lumped;
        else if (area_rule == "consistent_shape")
            result.quad8_nodal_area_rule = Quad8NodalAreaRule::consistent_shape;
        else
            value_error(document,
                required_entry(document, *mechanical, "quad8_nodal_area_rule"),
                "quad8_nodal_area_rule must be 'positive_lumped' or 'consistent_shape'");
        if (result.mechanical_discretization == MechanicalContactDiscretization::surface_to_surface
            && find_entry(*mechanical, "quad8_nodal_area_rule") != nullptr)
            value_error(document,
                required_entry(document, *mechanical, "quad8_nodal_area_rule"),
                "quad8_nodal_area_rule applies only to discretization = node_to_surface");
        const InputEntry* penetration_tolerance = find_entry(*mechanical, "penetration_tolerance");
        const InputEntry* maximum_augmented_iterations = find_entry(*mechanical, "maximum_augmented_iterations");
        if (result.mechanical_formulation == MechanicalContactFormulation::penalty) {
            if (penetration_tolerance != nullptr)
                value_error(document,
                    *penetration_tolerance,
                    "penetration_tolerance requires formulation = augmented_lagrangian");
            if (maximum_augmented_iterations != nullptr)
                value_error(document,
                    *maximum_augmented_iterations,
                    "maximum_augmented_iterations requires formulation = augmented_lagrangian");
        } else {
            result.penetration_tolerance = read_optional_double(document, *mechanical, "penetration_tolerance", 1.0e-8);
            result.maximum_augmented_iterations =
                read_optional_size(document, *mechanical, "maximum_augmented_iterations", 20);
        }
    }
    if (thermal != nullptr && result.gap_heat_conductance_law == GapHeatConductanceLaw::affine
        && result.gap_conductance_pressure_derivative != 0.0 && mechanical == nullptr)
        value_error(document,
            required_entry(document, *thermal, "pressure_derivative"),
            "pressure-dependent thermal contact requires a mechanical contact definition");
    if (thermal != nullptr && result.gap_heat_conductance_law == GapHeatConductanceLaw::affine
        && result.gap_conductance_pressure_derivative != 0.0
        && result.mechanical_formulation != MechanicalContactFormulation::penalty)
        value_error(document,
            required_entry(document, *thermal, "pressure_derivative"),
            "pressure-dependent thermal contact requires formulation = penalty");
    return result;
}

BoundaryConditionDefinition make_boundary_condition(const InputDocument& document,
    const InputSection& section,
    BoundaryConditionType type,
    Field field,
    double value,
    bool scale_with_load,
    const std::string& function) {
    BoundaryConditionDefinition result{section.name,
        type,
        read_string(document, section, "boundary"),
        field,
        value,
        scale_with_load};
    result.function = function;
    return result;
}

BoundaryConditionDefinition read_boundary_condition(const InputDocument& document, const InputSection& section) {
    validate_keys(document,
        section,
        {"type",
            "boundary",
            "field",
            "value",
            "scale_with_load",
            "function",
            "heat_transfer_coefficient",
            "ambient_temperature",
            "coefficient_function",
            "ambient_temperature_function",
            "configuration"});
    const std::string type = read_string(document, section, "type");
    const bool scale_with_load = read_optional_bool(document, section, "scale_with_load", false);
    const std::string function = read_optional_string(section, "function", {});
    if (scale_with_load && !function.empty())
        value_error(document,
            required_entry(document, section, "scale_with_load"),
            "scale_with_load cannot be combined with function");
    if (type == "dirichlet") {
        forbid_key(document, section, "configuration", "type='dirichlet'");
        forbid_convection_keys(document, section, "type='dirichlet'");
        return make_boundary_condition(document,
            section,
            BoundaryConditionType::dirichlet,
            parse_field(document, required_entry(document, section, "field")),
            read_double(document, section, "value"),
            scale_with_load,
            function);
    }
    if (type == "pressure") {
        forbid_key(document, section, "field", "type='pressure'");
        forbid_convection_keys(document, section, "type='pressure'");
        BoundaryConditionDefinition result = make_boundary_condition(document,
            section,
            BoundaryConditionType::pressure,
            Field::radial_displacement,
            read_double(document, section, "value"),
            scale_with_load,
            function);
        const InputEntry* configuration_entry = find_entry(section, "configuration");
        const std::string configuration = read_optional_string(section, "configuration", "reference");
        if (configuration != "reference" && configuration != "current")
            value_error(document,
                required_entry(document, section, "configuration"),
                "pressure configuration must be reference or current");
        result.use_displaced_geometry = configuration == "current";
        result.configuration_explicit = configuration_entry != nullptr;
        return result;
    }
    if (type == "traction") {
        forbid_convection_keys(document, section, "type='traction'");
        const Field field = parse_field(document, required_entry(document, section, "field"));
        if (field == Field::temperature)
            value_error(document, required_entry(document, section, "field"), "traction requires a displacement field");
        BoundaryConditionDefinition result = make_boundary_condition(document,
            section,
            BoundaryConditionType::traction,
            field,
            read_double(document, section, "value"),
            scale_with_load,
            function);
        const InputEntry* configuration_entry = find_entry(section, "configuration");
        const std::string configuration = read_optional_string(section, "configuration", "reference");
        if (configuration != "reference" && configuration != "current")
            value_error(document,
                required_entry(document, section, "configuration"),
                "traction configuration must be reference or current");
        result.use_displaced_geometry = configuration == "current";
        result.configuration_explicit = configuration_entry != nullptr;
        return result;
    }
    if (type == "heat_flux") {
        forbid_key(document, section, "field", "type='heat_flux'");
        forbid_convection_keys(document, section, "type='heat_flux'");
        BoundaryConditionDefinition result = make_boundary_condition(document,
            section,
            BoundaryConditionType::heat_flux,
            Field::temperature,
            read_double(document, section, "value"),
            scale_with_load,
            function);
        const InputEntry* configuration_entry = find_entry(section, "configuration");
        const std::string configuration = read_optional_string(section, "configuration", "reference");
        if (configuration != "reference" && configuration != "current")
            value_error(document,
                required_entry(document, section, "configuration"),
                "heat_flux configuration must be reference or current");
        result.use_displaced_geometry = configuration == "current";
        result.configuration_explicit = configuration_entry != nullptr;
        return result;
    }
    if (type == "convection") {
        forbid_key(document, section, "field", "type='convection'");
        forbid_key(document, section, "value", "type='convection'");
        forbid_key(document, section, "scale_with_load", "type='convection'");
        forbid_key(document, section, "function", "type='convection'");
        BoundaryConditionDefinition result = make_boundary_condition(document,
            section,
            BoundaryConditionType::convection,
            Field::temperature,
            0.0,
            false,
            {});
        result.heat_transfer_coefficient = read_double(document, section, "heat_transfer_coefficient");
        result.ambient_temperature = read_double(document, section, "ambient_temperature");
        result.coefficient_function = read_optional_string(section, "coefficient_function", {});
        result.ambient_temperature_function = read_optional_string(section, "ambient_temperature_function", {});
        const InputEntry* configuration_entry = find_entry(section, "configuration");
        const std::string configuration = read_optional_string(section, "configuration", "reference");
        if (configuration != "reference" && configuration != "current")
            value_error(document,
                required_entry(document, section, "configuration"),
                "convection configuration must be reference or current");
        result.use_displaced_geometry = configuration == "current";
        result.configuration_explicit = configuration_entry != nullptr;
        return result;
    }
    value_error(document, required_entry(document, section, "type"), "unknown boundary-condition type '" + type + "'");
}

void read_case(const InputDocument& document, FuelSimCaseDefinition& result) {
    const InputSection& case_section = required_section(document, "Case");
    validate_keys(document, case_section, {"version", "problem", "geometry"});
    const std::size_t version = read_size(document, case_section, "version");
    if (version != 3)
        value_error(document,
            required_entry(document, case_section, "version"),
            "unsupported fuelsim input version '" + std::to_string(version) + "'");
    result.version = 3;
    const std::string problem = read_string(document, case_section, "problem");
    if (problem == "steady")
        result.problem = CaseProblem::steady;
    else if (problem == "transient")
        result.problem = CaseProblem::transient;
    else
        value_error(document, required_entry(document, case_section, "problem"), "unknown problem '" + problem + "'");
    const std::string geometry = read_string(document, case_section, "geometry");
    if (geometry == "axisymmetric_rz")
        result.geometry = CaseGeometry::axisymmetric_rz;
    else if (geometry == "cartesian_3d")
        result.geometry = CaseGeometry::cartesian_3d;
    else
        value_error(document,
            required_entry(document, case_section, "geometry"),
            "unknown geometry '" + geometry + "'");
}

void read_mesh(const InputDocument& document, const std::string& path, FuelSimCaseDefinition& result) {
    const InputSection& mesh = required_section(document, "Mesh");
    validate_keys(document, mesh, {"type", "file"});
    if (read_string(document, mesh, "type") != "exodus")
        value_error(document, required_entry(document, mesh, "type"), "only mesh type 'exodus' is supported");
    result.mesh_file = resolved_path(path, read_string(document, mesh, "file"));
}

void read_time_functions(const InputDocument& document, FuelSimCaseDefinition& result) {
    const InputSection* time_functions = find_section(document, "TimeFunctions");
    if (time_functions != nullptr) {
        validate_keys(document, *time_functions, {});
        for (const InputSection* section : direct_children(document, "TimeFunctions")) {
            validate_keys(document, *section, {"type", "times", "values"});
            if (read_string(document, *section, "type") != "piecewise_linear")
                value_error(document,
                    required_entry(document, *section, "type"),
                    "only time-function type 'piecewise_linear' is supported");
            result.spatial.time_tables.emplace_back(section->name,
                parse_double_list(document, required_entry(document, *section, "times")),
                parse_double_list(document, required_entry(document, *section, "values")));
        }
    }
}

void read_regions(const InputDocument& document,
    const std::string& path,
    const std::vector<ParsedMaterial>& materials,
    FuelSimCaseDefinition& result) {
    const InputSection& regions = required_section(document, "Regions");
    validate_keys(document, regions, {});
    for (const InputSection* section : direct_children(document, "Regions"))
        result.spatial.regions.push_back(read_region(document, *section, materials, result.geometry));
    if (result.spatial.regions.empty())
        throw std::invalid_argument(path + ": [Regions] requires a child region");
}

void read_contacts(const InputDocument& document, FuelSimCaseDefinition& result) {
    if (const InputSection* contacts = find_section(document, "Contact"); contacts != nullptr)
        validate_keys(document, *contacts, {});
    for (const InputSection* section : direct_children(document, "Contact"))
        result.spatial.contacts.push_back(read_contact(document, *section, result.geometry));
}

void read_boundary_conditions(const InputDocument& document, FuelSimCaseDefinition& result) {
    if (const InputSection* boundary_conditions = find_section(document, "BoundaryConditions");
        boundary_conditions != nullptr)
        validate_keys(document, *boundary_conditions, {});
    for (const InputSection* section : direct_children(document, "BoundaryConditions"))
        result.spatial.boundary_conditions.push_back(read_boundary_condition(document, *section));
}

void validate_time_function_references(const std::string& path, const FuelSimCaseDefinition& result) {
    const auto require_function = [&](const std::string& name, const std::string& owner) {
        if (name.empty())
            return;
        const auto found = std::find_if(result.spatial.time_tables.begin(),
            result.spatial.time_tables.end(),
            [&name](const PiecewiseLinearTimeTable& table) { return table.name() == name; });
        if (found == result.spatial.time_tables.end())
            throw std::invalid_argument(path + ": " + owner + " references unknown time function '" + name + "'");
    };
    for (const RegionDefinition& region : result.spatial.regions)
        require_function(region.heat_source_function, "region '" + region.name + "'");
    for (const BoundaryConditionDefinition& boundary : result.spatial.boundary_conditions) {
        require_function(boundary.function, "boundary condition '" + boundary.name + "'");
        require_function(boundary.coefficient_function, "boundary condition '" + boundary.name + "'");
        require_function(boundary.ambient_temperature_function, "boundary condition '" + boundary.name + "'");
    }
    if (result.problem == CaseProblem::steady && !result.spatial.time_tables.empty())
        throw std::invalid_argument(path + ": time functions are only valid for transient cases");
}

void read_executioner(const InputDocument& document, const std::string& path, FuelSimCaseDefinition& result) {
    const InputSection& executioner = required_section(document, "Executioner");
    const std::string executioner_type = read_string(document, executioner, "type");
    if (result.problem == CaseProblem::steady) {
        validate_keys(document,
            executioner,
            {"type",
                "load_steps",
                "cutback_factor",
                "maximum_cutbacks",
                "minimum_load_increment",
                "use_small_strain_predictor",
                "use_linear_load_predictor"});
        if (executioner_type != "steady")
            value_error(document,
                required_entry(document, executioner, "type"),
                "problem='steady' requires type='steady'");
        result.steady_execution.load_steps = read_size(document, executioner, "load_steps");
        result.steady_execution.cutback_factor = read_optional_double(document, executioner, "cutback_factor", 0.5);
        result.steady_execution.maximum_cutbacks_per_step =
            read_optional_size(document, executioner, "maximum_cutbacks", 12);
        result.steady_execution.minimum_load_increment =
            read_optional_double(document, executioner, "minimum_load_increment", 1.0e-6);
        result.steady_execution.use_small_strain_predictor =
            read_optional_bool(document, executioner, "use_small_strain_predictor", false);
        result.steady_execution.use_linear_load_predictor =
            read_optional_bool(document, executioner, "use_linear_load_predictor", false);
    } else {
        validate_keys(document,
            executioner,
            {"type",
                "end_time",
                "initial_time_step",
                "minimum_time_step",
                "maximum_time_step",
                "growth_factor",
                "cutback_factor",
                "maximum_cutbacks",
                "load_ramp_time",
                "restart",
                "target_nonlinear_iterations",
                "iteration_window",
                "time_error_relative_tolerance",
                "temperature_time_absolute_tolerance",
                "displacement_time_absolute_tolerance",
                "time_error_safety_factor",
                "strain_history_time_absolute_tolerance",
                "stress_history_time_absolute_tolerance",
                "include_thermal_time_term",
                "use_linear_time_predictor"});
        if (executioner_type != "transient")
            value_error(document,
                required_entry(document, executioner, "type"),
                "problem='transient' requires type='transient'");
        result.transient_execution = {read_double(document, executioner, "end_time"),
            read_double(document, executioner, "initial_time_step"),
            read_double(document, executioner, "minimum_time_step"),
            read_double(document, executioner, "maximum_time_step"),
            read_double(document, executioner, "growth_factor"),
            read_double(document, executioner, "cutback_factor"),
            read_size(document, executioner, "maximum_cutbacks"),
            read_double(document, executioner, "load_ramp_time"),
            read_optional_size(document, executioner, "target_nonlinear_iterations", 0),
            read_optional_size(document, executioner, "iteration_window", 0),
            read_optional_double(document, executioner, "time_error_relative_tolerance", 0.0),
            read_optional_double(document, executioner, "temperature_time_absolute_tolerance", 1.0e-3),
            read_optional_double(document, executioner, "displacement_time_absolute_tolerance", 1.0e-10),
            read_optional_double(document, executioner, "time_error_safety_factor", 0.9),
            read_optional_double(document, executioner, "strain_history_time_absolute_tolerance", 1.0e-10),
            read_optional_double(document, executioner, "stress_history_time_absolute_tolerance", 1.0),
            read_optional_bool(document, executioner, "include_thermal_time_term", true),
            read_optional_bool(document, executioner, "use_linear_time_predictor", false)};
        result.restart_file = read_optional_path(path, executioner, "restart");
    }
}

void read_solver(const InputDocument& document, const std::string& path, FuelSimCaseDefinition& result) {
    const InputSection* solver_section = find_section(document, "Solver");
    if (solver_section == nullptr)
        return;
    const InputSection& solver = *solver_section;
    validate_keys(document,
        solver,
        {"absolute_tolerance",
            "relative_tolerance",
            "step_tolerance",
            "maximum_iterations",
            "linear_solver",
            "preconditioner",
            "direct_factorization",
            "mumps_ordering",
            "linear_relative_tolerance",
            "maximum_linear_iterations",
            "backtracking_fallback",
            "field_residual_scaling",
            "field_residual_convergence",
            "residual_reduction_tolerance",
            "temperature_residual_absolute_tolerance",
            "mechanical_residual_absolute_tolerance",
            "temperature_residual_scale",
            "mechanical_residual_scale",
            "jacobian_lag",
            "predictor_jacobian_lag",
            "line_search"});
    result.solver.absolute_tolerance = read_optional_double(document, solver, "absolute_tolerance", 1.0e-8);
    result.solver.relative_tolerance = read_optional_double(document, solver, "relative_tolerance", 1.0e-10);
    result.solver.step_tolerance = read_optional_double(document, solver, "step_tolerance", 1.0e-12);
    result.solver.maximum_iterations = read_optional_int(document, solver, "maximum_iterations", 40);
    const std::string linear_solver = read_optional_string(solver, "linear_solver", "automatic");
    if (linear_solver == "direct")
        result.solver.linear_solver = SolverOptions::LinearSolver::direct;
    else if (linear_solver == "gmres")
        result.solver.linear_solver = SolverOptions::LinearSolver::gmres;
    else if (linear_solver != "automatic")
        throw std::invalid_argument(path + ": linear_solver must be automatic, direct, or gmres");
    const std::string preconditioner = read_optional_string(solver, "preconditioner", "automatic");
    if (preconditioner == "lu")
        result.solver.preconditioner = SolverOptions::Preconditioner::lu;
    else if (preconditioner == "block_jacobi")
        result.solver.preconditioner = SolverOptions::Preconditioner::block_jacobi;
    else if (preconditioner == "field_split")
        result.solver.preconditioner = SolverOptions::Preconditioner::field_split;
    else if (preconditioner == "hypre")
        result.solver.preconditioner = SolverOptions::Preconditioner::hypre;
    else if (preconditioner != "automatic")
        throw std::invalid_argument(path
                                    + ": preconditioner must be automatic, lu, block_jacobi, "
                                      "field_split, or hypre");
    const std::string direct_factorization = read_optional_string(solver, "direct_factorization", "automatic");
    if (direct_factorization == "mumps")
        result.solver.direct_factorization = SolverOptions::DirectFactorization::mumps;
    else if (direct_factorization != "automatic")
        throw std::invalid_argument(path + ": direct_factorization must be automatic or mumps");
    const std::string mumps_ordering = read_optional_string(solver, "mumps_ordering", "automatic");
    if (mumps_ordering == "scotch")
        result.solver.mumps_ordering = SolverOptions::MumpsOrdering::scotch;
    else if (mumps_ordering == "pord")
        result.solver.mumps_ordering = SolverOptions::MumpsOrdering::pord;
    else if (mumps_ordering != "automatic")
        throw std::invalid_argument(path + ": mumps_ordering must be automatic, scotch, or pord");
    if (result.solver.direct_factorization == SolverOptions::DirectFactorization::mumps
        && (result.solver.linear_solver == SolverOptions::LinearSolver::gmres
            || (result.solver.preconditioner != SolverOptions::Preconditioner::automatic
                && result.solver.preconditioner != SolverOptions::Preconditioner::lu)))
        throw std::invalid_argument(path + ": direct_factorization=mumps requires a direct LU solve");
    result.solver.linear_relative_tolerance =
        read_optional_double(document, solver, "linear_relative_tolerance", 1.0e-8);
    result.solver.maximum_linear_iterations = read_optional_int(document, solver, "maximum_linear_iterations", 500);
    result.solver.jacobian_lag = read_optional_int(document, solver, "jacobian_lag", 1);
    if (result.solver.jacobian_lag < 1)
        throw std::invalid_argument(path + ": jacobian_lag must be at least one");
    result.solver.predictor_jacobian_lag = read_optional_int(document, solver, "predictor_jacobian_lag", 0);
    const std::string line_search = read_optional_string(solver, "line_search", "basic");
    if (line_search == "basic")
        result.solver.line_search = SolverOptions::LineSearch::basic;
    else if (line_search == "backtracking")
        result.solver.line_search = SolverOptions::LineSearch::backtracking;
    else if (line_search == "critical_point")
        result.solver.line_search = SolverOptions::LineSearch::critical_point;
    else
        throw std::invalid_argument(path + ": line_search must be 'basic', 'backtracking', or 'critical_point'");
    result.solver.backtracking_fallback = read_optional_bool(document, solver, "backtracking_fallback", true);
    result.solver.field_residual_scaling = read_optional_bool(document, solver, "field_residual_scaling", false);
    result.solver.field_residual_convergence =
        read_optional_bool(document, solver, "field_residual_convergence", false);
    result.solver.residual_reduction_tolerance =
        read_optional_double(document, solver, "residual_reduction_tolerance", 1.0e-6);
    result.solver.temperature_residual_absolute_tolerance =
        read_optional_double(document, solver, "temperature_residual_absolute_tolerance", 1.0e-8);
    result.solver.mechanical_residual_absolute_tolerance =
        read_optional_double(document, solver, "mechanical_residual_absolute_tolerance", 1.0e-4);
    result.solver.temperature_residual_scale =
        read_optional_double(document, solver, "temperature_residual_scale", 0.0);
    result.solver.mechanical_residual_scale = read_optional_double(document, solver, "mechanical_residual_scale", 0.0);
}

void read_outputs(const InputDocument& document, const std::string& path, FuelSimCaseDefinition& result) {
    const InputSection* output_section = find_section(document, "Outputs");
    if (output_section == nullptr)
        return;
    const InputSection& outputs = *output_section;
    validate_keys(document,
        outputs,
        {"console",
            "csv",
            "exodus",
            "exodus_interval",
            "history",
            "history_interval",
            "progress_interval",
            "checkpoint",
            "checkpoint_interval"});
    result.outputs.console = read_optional_bool(document, outputs, "console", true);
    result.outputs.csv_file = read_optional_path(path, outputs, "csv");
    result.outputs.exodus_file = read_optional_path(path, outputs, "exodus");
    result.outputs.exodus_interval = read_optional_size(document, outputs, "exodus_interval", 1);
    result.outputs.history_file = read_optional_path(path, outputs, "history");
    result.outputs.history_interval = read_optional_size(document, outputs, "history_interval", 1);
    result.outputs.progress_interval = read_optional_size(document, outputs, "progress_interval", 1);
    result.outputs.checkpoint_file = read_optional_path(path, outputs, "checkpoint");
    result.outputs.checkpoint_interval = read_optional_size(document, outputs, "checkpoint_interval", 1);
    const std::string input_file = std::filesystem::absolute(path).lexically_normal().string();
    const std::array<std::string, 4> output_paths = {result.outputs.csv_file,
        result.outputs.exodus_file,
        result.outputs.history_file,
        result.outputs.checkpoint_file};
    for (const std::string& output_path : output_paths) {
        if (output_path.empty())
            continue;
        if (output_path == input_file)
            throw std::invalid_argument(path + ": an output file must not overwrite the input card");
        if (!result.mesh_file.empty() && output_path == result.mesh_file)
            throw std::invalid_argument(path + ": an output file must not overwrite the input mesh");
        if (!result.restart_file.empty() && output_path == result.restart_file)
            throw std::invalid_argument(path + ": an output file must not overwrite the restart checkpoint");
    }
    for (std::size_t first = 0; first < output_paths.size(); ++first) {
        if (output_paths[first].empty())
            continue;
        for (std::size_t second = first + 1; second < output_paths.size(); ++second)
            if (!output_paths[second].empty() && output_paths[first] == output_paths[second])
                throw std::invalid_argument(path + ": output file paths must differ");
    }
    if (result.problem == CaseProblem::steady
        && (!result.outputs.history_file.empty() || find_entry(outputs, "history_interval") != nullptr
            || find_entry(outputs, "progress_interval") != nullptr || find_entry(outputs, "exodus_interval") != nullptr
            || !result.outputs.checkpoint_file.empty() || find_entry(outputs, "checkpoint_interval") != nullptr))
        throw std::invalid_argument(path + ": transient output controls are only valid for transient cases");
    const auto require_output_file =
        [&](const std::string& file, const std::string& interval_key, const std::string& output_key) {
            if (file.empty() && find_entry(outputs, interval_key) != nullptr)
                value_error(document,
                    required_entry(document, outputs, interval_key),
                    interval_key + " requires " + output_key);
        };
    require_output_file(result.outputs.exodus_file, "exodus_interval", "exodus");
    require_output_file(result.outputs.history_file, "history_interval", "history");
    require_output_file(result.outputs.checkpoint_file, "checkpoint_interval", "checkpoint");
}
} // namespace

FuelSimCaseDefinition read_case_input(const std::string& path) {
    const MaterialFunctionRegistry registry = make_builtin_material_function_registry();
    return read_case_input(path, registry);
}

FuelSimCaseDefinition read_case_input(const std::string& path, const MaterialFunctionRegistry& registry) {
    const InputDocument document = parse_input_file(path);
    validate_sections(document);
    FuelSimCaseDefinition result{};
    read_case(document, result);
    read_mesh(document, path, result);
    read_time_functions(document, result);
    const std::vector<ParsedMaterial> materials = read_materials(document, registry);
    read_regions(document, path, materials, result);
    if (result.geometry != CaseGeometry::cartesian_3d)
        for (const RegionDefinition& region : result.spatial.regions)
            if (region.hex8_element_formulation != Hex8ElementFormulation::c3d8t)
                throw std::invalid_argument(path + ": region '" + region.name
                                            + "' selects c3d8rt outside Cartesian three-dimensional geometry");
    read_contacts(document, result);
    read_boundary_conditions(document, result);
    validate_time_function_references(path, result);
    read_executioner(document, path, result);
    read_solver(document, path, result);
    read_outputs(document, path, result);
    return result;
}
} // namespace fuelsim
