#include "support/exodus_result_reader.hpp"
#include <algorithm>
#include <cstdint>
#include <exodusII.h>
#include <limits>
#include <stdexcept>

namespace fuelsim::test {
namespace {
class ExodusFile final {
  public:
    explicit ExodusFile(const std::string& path) : _id(-1) {
        int cpu_word_size = static_cast<int>(sizeof(double));
        int io_word_size = 0;
        float version = 0.0F;
        _id = ex_open(path.c_str(), EX_READ, &cpu_word_size, &io_word_size, &version);
        if (_id < 0) throw std::runtime_error("Could not read fuelsim Exodus results: " + path);
        ex_set_int64_status(_id, EX_ALL_INT64_API);
    }

    ~ExodusFile() {
        if (_id >= 0) ex_close(_id);
    }

    ExodusFile(const ExodusFile&) = delete;
    ExodusFile& operator=(const ExodusFile&) = delete;

    int id() const noexcept { return _id; }

  private:
    int _id;
};

void check_exodus(int status, const std::string& message) {
    if (status < 0) throw std::runtime_error(message);
}

std::size_t count(std::int64_t value, const std::string& name) {
    if (value < 0 || static_cast<std::uint64_t>(value) > std::numeric_limits<std::size_t>::max())
        throw std::runtime_error("Invalid " + name + " in fuelsim Exodus results");
    return static_cast<std::size_t>(value);
}

std::vector<std::string> variable_names(int exoid, ex_entity_type type) {
    int variable_count = 0;
    check_exodus(ex_get_variable_param(exoid, type, &variable_count),
        "Could not read variable count from fuelsim Exodus results");
    const std::size_t maximum_name_length =
        count(ex_inquire_int(exoid, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH), "maximum variable-name length");
    std::vector<char> storage(maximum_name_length + 1, '\0');
    std::vector<std::string> result;
    result.reserve(static_cast<std::size_t>(variable_count));
    for (int variable = 1; variable <= variable_count; ++variable) {
        std::fill(storage.begin(), storage.end(), '\0');
        check_exodus(ex_get_variable_name(exoid, type, variable, storage.data()),
            "Could not read a variable name from fuelsim Exodus results");
        result.emplace_back(storage.data());
    }
    return result;
}

std::vector<std::vector<double>> read_nodal_variables(
    int exoid, int step, std::size_t node_count, std::size_t variable_count) {
    std::vector<std::vector<double>> result(variable_count, std::vector<double>(node_count, 0.0));
    for (std::size_t variable = 0; variable < variable_count; ++variable)
        check_exodus(ex_get_var(exoid, step, EX_NODAL, static_cast<int>(variable + 1), 1,
                         static_cast<std::int64_t>(node_count), result[variable].data()),
            "Could not read a nodal variable from fuelsim Exodus results");
    return result;
}

std::vector<std::vector<double>> read_element_variables(
    int exoid, int step, std::size_t element_count, std::size_t variable_count) {
    const std::size_t block_count = count(ex_inquire_int(exoid, EX_INQ_ELEM_BLK), "element-block count");
    std::vector<std::int64_t> block_ids(block_count, 0);
    if (block_count != 0)
        check_exodus(ex_get_ids(exoid, EX_ELEM_BLOCK, block_ids.data()),
            "Could not read element-block IDs from fuelsim Exodus results");
    std::vector<std::vector<double>> result(variable_count);
    for (std::vector<double>& values : result) values.reserve(element_count);
    for (const std::int64_t block_id : block_ids) {
        ex_block block{};
        block.type = EX_ELEM_BLOCK;
        block.id = block_id;
        check_exodus(ex_get_block_param(exoid, &block), "Could not read an element block from fuelsim Exodus results");
        const std::size_t block_size = count(block.num_entry, "element-block size");
        for (std::size_t variable = 0; variable < variable_count; ++variable) {
            std::vector<double> values(block_size, 0.0);
            if (block_size != 0)
                check_exodus(ex_get_var(exoid, step, EX_ELEM_BLOCK, static_cast<int>(variable + 1), block_id,
                                 static_cast<std::int64_t>(block_size), values.data()),
                    "Could not read an element variable from fuelsim Exodus results");
            result[variable].insert(result[variable].end(), values.begin(), values.end());
        }
    }
    for (const std::vector<double>& values : result)
        if (values.size() != element_count)
            throw std::runtime_error("Element count changed while reading fuelsim Exodus results");
    return result;
}

const std::vector<double>& named_variable(
    const std::vector<std::string>& names, const std::vector<std::vector<double>>& values, const std::string& name) {
    const auto found = std::find(names.begin(), names.end(), name);
    if (found == names.end()) throw std::invalid_argument("fuelsim Exodus results are missing variable '" + name + "'");
    return values.at(static_cast<std::size_t>(found - names.begin()));
}
} // namespace

const std::vector<double>& ExodusResults::nodal(const std::string& name) const {
    return named_variable(nodal_variable_names, nodal_variables, name);
}

const std::vector<double>& ExodusResults::element(const std::string& name) const {
    return named_variable(element_variable_names, element_variables, name);
}

ExodusResults read_final_exodus_results(const std::string& path) {
    ExodusFile file(path);
    ExodusResults result;
    const std::size_t dimension = count(ex_inquire_int(file.id(), EX_INQ_DIM), "coordinate dimension");
    const std::size_t node_count = count(ex_inquire_int(file.id(), EX_INQ_NODES), "node count");
    const std::size_t element_count = count(ex_inquire_int(file.id(), EX_INQ_ELEM), "element count");
    result.step_count = count(ex_inquire_int(file.id(), EX_INQ_TIME), "time-step count");
    if ((dimension != 2 && dimension != 3) || node_count == 0 || result.step_count == 0)
        throw std::invalid_argument("fuelsim Exodus results must contain a two- or three-dimensional solved model");
    if (result.step_count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Too many time steps in fuelsim Exodus results");

    std::array<std::vector<double>, 3> coordinates;
    for (std::vector<double>& values : coordinates) values.resize(node_count, 0.0);
    check_exodus(ex_get_coord(file.id(), coordinates[0].data(), coordinates[1].data(),
                     dimension == 3 ? coordinates[2].data() : nullptr),
        "Could not read coordinates from fuelsim Exodus results");
    result.nodes.resize(node_count);
    for (std::size_t node = 0; node < node_count; ++node)
        result.nodes[node] = {coordinates[0][node], coordinates[1][node], coordinates[2][node]};

    const int final_step = static_cast<int>(result.step_count);
    check_exodus(ex_get_time(file.id(), final_step, &result.time), "Could not read result time");
    result.nodal_variable_names = variable_names(file.id(), EX_NODAL);
    result.nodal_variables =
        read_nodal_variables(file.id(), final_step, node_count, result.nodal_variable_names.size());
    result.element_variable_names = variable_names(file.id(), EX_ELEM_BLOCK);
    result.element_variables =
        read_element_variables(file.id(), final_step, element_count, result.element_variable_names.size());
    return result;
}
} // namespace fuelsim::test
