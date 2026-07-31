#include "fuelsim/exodus_mesh_io.hpp"

#include <exodusII.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fuelsim {
namespace {

void check_exodus(int status, const std::string& operation) {
    if (status < 0)
        throw std::runtime_error(operation + ": " + ex_strerror(status));
}

std::size_t checked_size(std::int64_t value, const char* description) {
    if (value < 0 || static_cast<std::uint64_t>(value) >
                         std::numeric_limits<std::size_t>::max())
        throw std::length_error(std::string(description) + " is out of range");
    return static_cast<std::size_t>(value);
}

std::int64_t checked_count(std::size_t value, const char* description) {
    if (value >
        static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error(std::string(description) + " is out of range");
    return static_cast<std::int64_t>(value);
}

class ExodusFile final {
  public:
    explicit ExodusFile(int id) : _id(id) {}

    ExodusFile(const ExodusFile&) = delete;
    ExodusFile& operator=(const ExodusFile&) = delete;

    ~ExodusFile() {
        if (_id >= 0)
            ex_close(_id);
    }

    int id() const noexcept {
        return _id;
    }

    void close() {
        const int id = _id;
        _id = -1;
        check_exodus(ex_close(id), "Could not close Exodus file");
    }

  private:
    int _id;
};

std::string normalized_topology(const char* topology) {
    std::string result(topology);
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char character) {
                       return static_cast<char>(std::toupper(character));
                   });
    return result;
}

struct BlockConnectivity final {
    std::int64_t _id;
    std::vector<std::int64_t> _nodes;
};

} // namespace

UnstructuredQuad4Mesh ExodusMeshIo::read_quad4(const std::string& path) {
    int cpu_word_size = static_cast<int>(sizeof(double));
    int io_word_size = 0;
    float version = 0.0F;
    const int exoid =
        ex_open(path.c_str(), EX_READ, &cpu_word_size, &io_word_size, &version);
    if (exoid < 0)
        throw std::runtime_error("Could not open Exodus file '" + path +
                                 "': " + ex_strerror(exoid));
    ExodusFile file(exoid);
    ex_set_int64_status(file.id(), EX_ALL_INT64_API);

    ex_init_params parameters{};
    check_exodus(ex_get_init_ext(file.id(), &parameters),
                 "Could not read Exodus model parameters");
    if (parameters.num_dim != 2)
        throw std::runtime_error("Exodus Quad4 mesh must be two-dimensional");

    const std::size_t node_count =
        checked_size(parameters.num_nodes, "Exodus node count");
    const std::size_t block_count =
        checked_size(parameters.num_elem_blk, "Exodus element block count");
    if (node_count == 0 || block_count == 0)
        throw std::runtime_error(
            "Exodus Quad4 mesh must contain nodes and element blocks");

    std::vector<double> radial_coordinates(node_count);
    std::vector<double> axial_coordinates(node_count);
    check_exodus(ex_get_coord(file.id(), radial_coordinates.data(),
                              axial_coordinates.data(), nullptr),
                 "Could not read Exodus coordinates");

    std::vector<RzPoint> nodes;
    nodes.reserve(node_count);
    for (std::size_t node = 0; node < node_count; ++node)
        nodes.push_back({radial_coordinates[node], axial_coordinates[node]});

    std::vector<std::int64_t> block_ids(block_count);
    check_exodus(ex_get_ids(file.id(), EX_ELEM_BLOCK, block_ids.data()),
                 "Could not read Exodus element block IDs");

    std::vector<Quad4Element> elements;
    std::vector<std::int64_t> element_block_ids;
    elements.reserve(checked_size(parameters.num_elem, "Exodus element count"));
    element_block_ids.reserve(elements.capacity());

    for (const std::int64_t block_id : block_ids) {
        ex_block block{};
        block.id = block_id;
        block.type = EX_ELEM_BLOCK;
        check_exodus(ex_get_block_param(file.id(), &block),
                     "Could not read Exodus element block");
        const std::string topology = normalized_topology(block.topology);
        if (block.num_nodes_per_entry != 4 ||
            (topology != "QUAD4" && topology != "QUAD"))
            throw std::runtime_error(
                "Exodus element blocks must contain only Quad4 elements");

        const std::size_t block_elements =
            checked_size(block.num_entry, "Exodus block element count");
        if (block_elements > std::numeric_limits<std::size_t>::max() / 4U)
            throw std::length_error("Exodus connectivity is too large");
        std::vector<std::int64_t> connectivity(block_elements * 4U);
        check_exodus(ex_get_conn(file.id(), EX_ELEM_BLOCK, block_id,
                                 connectivity.data(), nullptr, nullptr),
                     "Could not read Exodus Quad4 connectivity");

        for (std::size_t element = 0; element < block_elements; ++element) {
            Quad4Element quad{};
            for (std::size_t local_node = 0; local_node < 4U; ++local_node) {
                const std::int64_t exodus_node =
                    connectivity[element * 4U + local_node];
                if (exodus_node <= 0 ||
                    static_cast<std::uint64_t>(exodus_node) > node_count)
                    throw std::runtime_error(
                        "Exodus connectivity contains an invalid node ID");
                quad.nodes[local_node] =
                    static_cast<std::size_t>(exodus_node - 1);
            }
            elements.push_back(quad);
            element_block_ids.push_back(block_id);
        }
    }

    if (elements.size() !=
        checked_size(parameters.num_elem, "Exodus element count"))
        throw std::runtime_error(
            "Exodus element count does not match its element blocks");

    file.close();
    return UnstructuredQuad4Mesh(std::move(nodes), std::move(elements),
                                 std::move(element_block_ids));
}

void ExodusMeshIo::write_quad4(const std::string& path,
                               const UnstructuredQuad4Mesh& mesh) {
    int cpu_word_size = static_cast<int>(sizeof(double));
    int io_word_size = static_cast<int>(sizeof(double));
    const int mode = EX_CLOBBER | EX_ALL_INT64_DB | EX_ALL_INT64_API;
    const int exoid =
        ex_create(path.c_str(), mode, &cpu_word_size, &io_word_size);
    if (exoid < 0)
        throw std::runtime_error("Could not create Exodus file '" + path +
                                 "': " + ex_strerror(exoid));
    ExodusFile file(exoid);
    ex_set_int64_status(file.id(), EX_ALL_INT64_API);

    std::vector<BlockConnectivity> blocks;
    const std::vector<Quad4Element>& elements = mesh.elements();
    const std::vector<std::int64_t>& block_ids = mesh.element_block_ids();
    for (std::size_t element = 0; element < elements.size(); ++element) {
        auto block = std::find_if(blocks.begin(), blocks.end(),
                                  [block_id = block_ids[element]](
                                      const BlockConnectivity& candidate) {
                                      return candidate._id == block_id;
                                  });
        if (block == blocks.end()) {
            blocks.push_back({block_ids[element], {}});
            block = blocks.end() - 1;
        }
        block->_nodes.reserve(block->_nodes.size() + 4U);
        for (const std::size_t node : elements[element].nodes) {
            if (node >= static_cast<std::size_t>(
                            std::numeric_limits<std::int64_t>::max()))
                throw std::length_error("Exodus node ID is out of range");
            const std::int64_t exodus_node =
                static_cast<std::int64_t>(node) + 1;
            block->_nodes.push_back(exodus_node);
        }
    }

    ex_init_params parameters{};
    ex_copy_string(parameters.title, "fuelsim Quad4 mesh",
                   sizeof(parameters.title));
    parameters.num_dim = 2;
    parameters.num_nodes = checked_count(mesh.nodes().size(), "Node count");
    parameters.num_elem = checked_count(elements.size(), "Element count");
    parameters.num_elem_blk = checked_count(blocks.size(), "Block count");
    check_exodus(ex_put_init_ext(file.id(), &parameters),
                 "Could not write Exodus model parameters");

    std::vector<double> radial_coordinates;
    std::vector<double> axial_coordinates;
    radial_coordinates.reserve(mesh.nodes().size());
    axial_coordinates.reserve(mesh.nodes().size());
    for (const RzPoint& node : mesh.nodes()) {
        radial_coordinates.push_back(node.r);
        axial_coordinates.push_back(node.z);
    }
    check_exodus(ex_put_coord(file.id(), radial_coordinates.data(),
                              axial_coordinates.data(), nullptr),
                 "Could not write Exodus coordinates");
    char radial_name[] = "r";
    char axial_name[] = "z";
    char* coordinate_names[] = {radial_name, axial_name};
    check_exodus(ex_put_coord_names(file.id(), coordinate_names),
                 "Could not write Exodus coordinate names");

    for (const BlockConnectivity& block : blocks) {
        const std::int64_t block_element_count =
            checked_count(block._nodes.size() / 4U, "Block element count");
        check_exodus(ex_put_block(file.id(), EX_ELEM_BLOCK, block._id, "QUAD4",
                                  block_element_count, 4, 0, 0, 0),
                     "Could not write Exodus element block");
        check_exodus(ex_put_conn(file.id(), EX_ELEM_BLOCK, block._id,
                                 block._nodes.data(), nullptr, nullptr),
                     "Could not write Exodus Quad4 connectivity");
    }

    file.close();
}

} // namespace fuelsim
