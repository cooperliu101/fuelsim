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
    std::string _name;
    std::vector<std::int64_t> _nodes;
    std::vector<std::size_t> _source_elements;
};

std::string read_entity_name(int exoid, ex_entity_type type, std::int64_t id,
                             std::size_t maximum_length) {
    std::vector<char> name(maximum_length + 1U, '\0');
    check_exodus(ex_get_name(exoid, type, id, name.data()),
                 "Could not read Exodus entity name");
    return std::string(name.data());
}

std::vector<NodeSet> read_node_sets(int exoid, std::size_t set_count,
                                    std::size_t maximum_name_length) {
    std::vector<NodeSet> result;
    if (set_count == 0)
        return result;

    std::vector<std::int64_t> set_ids(set_count);
    check_exodus(ex_get_ids(exoid, EX_NODE_SET, set_ids.data()),
                 "Could not read Exodus node set IDs");
    result.reserve(set_count);
    for (const std::int64_t set_id : set_ids) {
        std::int64_t entry_count = 0;
        std::int64_t factor_count = 0;
        check_exodus(ex_get_set_param(exoid, EX_NODE_SET, set_id, &entry_count,
                                      &factor_count),
                     "Could not read Exodus node set parameters");
        std::vector<std::int64_t> entries(
            checked_size(entry_count, "Exodus node set size"));
        if (!entries.empty())
            check_exodus(
                ex_get_set(exoid, EX_NODE_SET, set_id, entries.data(), nullptr),
                "Could not read Exodus node set");
        NodeSet set{
            set_id,
            read_entity_name(exoid, EX_NODE_SET, set_id, maximum_name_length),
            {}};
        set.nodes.reserve(entries.size());
        for (const std::int64_t entry : entries) {
            if (entry <= 0)
                throw std::runtime_error(
                    "Exodus node set contains an invalid node ID");
            set.nodes.push_back(static_cast<std::size_t>(entry - 1));
        }
        result.push_back(std::move(set));
    }
    return result;
}

std::vector<SideSet> read_side_sets(int exoid, std::size_t set_count,
                                    std::size_t maximum_name_length) {
    std::vector<SideSet> result;
    if (set_count == 0)
        return result;

    std::vector<std::int64_t> set_ids(set_count);
    check_exodus(ex_get_ids(exoid, EX_SIDE_SET, set_ids.data()),
                 "Could not read Exodus side set IDs");
    result.reserve(set_count);
    for (const std::int64_t set_id : set_ids) {
        std::int64_t entry_count = 0;
        std::int64_t factor_count = 0;
        check_exodus(ex_get_set_param(exoid, EX_SIDE_SET, set_id, &entry_count,
                                      &factor_count),
                     "Could not read Exodus side set parameters");
        const std::size_t count =
            checked_size(entry_count, "Exodus side set size");
        std::vector<std::int64_t> elements(count);
        std::vector<std::int64_t> sides(count);
        if (count != 0)
            check_exodus(ex_get_set(exoid, EX_SIDE_SET, set_id, elements.data(),
                                    sides.data()),
                         "Could not read Exodus side set");
        SideSet set{
            set_id,
            read_entity_name(exoid, EX_SIDE_SET, set_id, maximum_name_length),
            {}};
        set.sides.reserve(count);
        for (std::size_t entry = 0; entry < count; ++entry) {
            if (elements[entry] <= 0 || sides[entry] <= 0 || sides[entry] > 4)
                throw std::runtime_error(
                    "Exodus side set contains an invalid Quad4 side");
            set.sides.push_back({static_cast<std::size_t>(elements[entry] - 1),
                                 static_cast<std::size_t>(sides[entry] - 1)});
        }
        result.push_back(std::move(set));
    }
    return result;
}

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
    const std::int64_t maximum_name_length =
        ex_inquire_int(file.id(), EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);
    if (maximum_name_length < 0)
        throw std::runtime_error("Could not read Exodus maximum name length");
    check_exodus(ex_set_max_name_length(file.id(),
                                        static_cast<int>(maximum_name_length)),
                 "Could not set Exodus maximum read name length");

    ex_init_params parameters{};
    check_exodus(ex_get_init_ext(file.id(), &parameters),
                 "Could not read Exodus model parameters");
    if (parameters.num_dim != 2)
        throw std::runtime_error("Exodus Quad4 mesh must be two-dimensional");

    const std::size_t node_count =
        checked_size(parameters.num_nodes, "Exodus node count");
    const std::size_t block_count =
        checked_size(parameters.num_elem_blk, "Exodus element block count");
    const std::size_t node_set_count =
        checked_size(parameters.num_node_sets, "Exodus node set count");
    const std::size_t side_set_count =
        checked_size(parameters.num_side_sets, "Exodus side set count");
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
    std::vector<ElementBlockInfo> element_blocks;
    element_blocks.reserve(block_count);
    elements.reserve(checked_size(parameters.num_elem, "Exodus element count"));
    element_block_ids.reserve(elements.capacity());

    for (const std::int64_t block_id : block_ids) {
        element_blocks.push_back(
            {block_id,
             read_entity_name(file.id(), EX_ELEM_BLOCK, block_id,
                              checked_size(maximum_name_length,
                                           "Exodus maximum name length"))});
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

    std::vector<NodeSet> node_sets = read_node_sets(
        file.id(), node_set_count,
        checked_size(maximum_name_length, "Exodus maximum name length"));
    std::vector<SideSet> side_sets = read_side_sets(
        file.id(), side_set_count,
        checked_size(maximum_name_length, "Exodus maximum name length"));

    file.close();
    return UnstructuredQuad4Mesh(
        std::move(nodes), std::move(elements), std::move(element_block_ids),
        std::move(element_blocks), std::move(node_sets), std::move(side_sets));
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
    blocks.reserve(mesh.element_blocks().size());
    for (const ElementBlockInfo& block : mesh.element_blocks())
        blocks.push_back({block.id, block.name, {}, {}});
    const std::vector<Quad4Element>& elements = mesh.elements();
    const std::vector<std::int64_t>& block_ids = mesh.element_block_ids();
    for (std::size_t element = 0; element < elements.size(); ++element) {
        auto block = std::find_if(blocks.begin(), blocks.end(),
                                  [block_id = block_ids[element]](
                                      const BlockConnectivity& candidate) {
                                      return candidate._id == block_id;
                                  });
        if (block == blocks.end())
            throw std::logic_error("Mesh element references an unknown block");
        block->_source_elements.push_back(element);
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
    parameters.num_node_sets =
        checked_count(mesh.node_sets().size(), "Node set count");
    parameters.num_side_sets =
        checked_count(mesh.side_sets().size(), "Side set count");
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

    std::vector<std::int64_t> written_element_ids(elements.size(), 0);
    std::int64_t next_element_id = 1;
    for (const BlockConnectivity& block : blocks) {
        const std::int64_t block_element_count =
            checked_count(block._nodes.size() / 4U, "Block element count");
        check_exodus(ex_put_block(file.id(), EX_ELEM_BLOCK, block._id, "QUAD4",
                                  block_element_count, 4, 0, 0, 0),
                     "Could not write Exodus element block");
        check_exodus(ex_put_conn(file.id(), EX_ELEM_BLOCK, block._id,
                                 block._nodes.data(), nullptr, nullptr),
                     "Could not write Exodus Quad4 connectivity");
        if (!block._name.empty())
            check_exodus(ex_put_name(file.id(), EX_ELEM_BLOCK, block._id,
                                     block._name.c_str()),
                         "Could not write Exodus element block name");
        for (const std::size_t source_element : block._source_elements)
            written_element_ids[source_element] = next_element_id++;
    }

    for (const NodeSet& set : mesh.node_sets()) {
        std::vector<std::int64_t> entries;
        entries.reserve(set.nodes.size());
        for (const std::size_t node : set.nodes) {
            if (node >= static_cast<std::size_t>(
                            std::numeric_limits<std::int64_t>::max()))
                throw std::length_error("Exodus node set ID is out of range");
            entries.push_back(static_cast<std::int64_t>(node) + 1);
        }
        check_exodus(
            ex_put_set_param(file.id(), EX_NODE_SET, set.id,
                             checked_count(entries.size(), "Node set size"), 0),
            "Could not write Exodus node set parameters");
        if (!entries.empty())
            check_exodus(ex_put_set(file.id(), EX_NODE_SET, set.id,
                                    entries.data(), nullptr),
                         "Could not write Exodus node set");
        if (!set.name.empty())
            check_exodus(
                ex_put_name(file.id(), EX_NODE_SET, set.id, set.name.c_str()),
                "Could not write Exodus node set name");
    }

    for (const SideSet& set : mesh.side_sets()) {
        std::vector<std::int64_t> set_elements;
        std::vector<std::int64_t> set_sides;
        set_elements.reserve(set.sides.size());
        set_sides.reserve(set.sides.size());
        for (const ElementSide& side : set.sides) {
            if (side.element >= static_cast<std::size_t>(
                                    std::numeric_limits<std::int64_t>::max()))
                throw std::length_error(
                    "Exodus side set element ID is out of range");
            set_elements.push_back(written_element_ids.at(side.element));
            set_sides.push_back(static_cast<std::int64_t>(side.local_side) + 1);
        }
        check_exodus(ex_put_set_param(
                         file.id(), EX_SIDE_SET, set.id,
                         checked_count(set.sides.size(), "Side set size"), 0),
                     "Could not write Exodus side set parameters");
        if (!set_elements.empty())
            check_exodus(ex_put_set(file.id(), EX_SIDE_SET, set.id,
                                    set_elements.data(), set_sides.data()),
                         "Could not write Exodus side set");
        if (!set.name.empty())
            check_exodus(
                ex_put_name(file.id(), EX_SIDE_SET, set.id, set.name.c_str()),
                "Could not write Exodus side set name");
    }

    file.close();
}

} // namespace fuelsim
