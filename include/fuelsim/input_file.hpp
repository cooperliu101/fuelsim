#ifndef FUELSIM_INPUT_FILE_HPP
#define FUELSIM_INPUT_FILE_HPP
#include <cstddef>
#include <string>
#include <vector>
namespace fuelsim {
struct InputEntry final {
    std::string key;
    std::string value;
    std::size_t line;
};
class InputSection final {
  public:
    const std::string& path() const noexcept;
    std::size_t line() const noexcept;
    const std::vector<InputEntry>& entries() const noexcept;
    const InputEntry& entry(const std::string& key) const;

  private:
    friend class InputParser;
    std::string _path;
    std::size_t _line = 0;
    std::vector<InputEntry> _entries;
};
class InputDocument final {
  public:
    const std::string& source_path() const noexcept;
    const std::vector<InputSection>& sections() const noexcept;
    const InputSection& section(const std::string& path) const;

  private:
    friend class InputParser;
    std::string _source_path;
    std::vector<InputSection> _sections;
};
class InputParser final {
  public:
    static InputDocument parse_file(const std::string& path);
};
} // namespace fuelsim
#endif
