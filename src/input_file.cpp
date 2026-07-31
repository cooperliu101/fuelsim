#include "fuelsim/input_file.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <utility>

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

} // namespace fuelsim
