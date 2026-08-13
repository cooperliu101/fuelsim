#include "fuelsim/results_io.hpp"
#include "fuelsim/transient_problem.hpp"
#include "problem_backend_access.hpp"
#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exodusII.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
namespace fuelsim::exodus_detail {
void check_exodus(int status, const std::string& operation) {
    if (status < 0) throw std::runtime_error(operation + ": " + ex_strerror(status));
}
class ExodusFile final {
  public:
    explicit ExodusFile(int id) : _id(id) {}
    ExodusFile(const ExodusFile&) = delete;
    ExodusFile& operator=(const ExodusFile&) = delete;
    ~ExodusFile() {
        if (_id >= 0) ex_close(_id);
    }
    int id() const noexcept { return _id; }
    void close(const std::string& operation = "Could not close Exodus file") {
        const int id = _id;
        _id = -1;
        check_exodus(ex_close(id), operation);
    }

  private:
    int _id;
};
} // namespace fuelsim::exodus_detail
namespace fuelsim {
namespace {
constexpr std::uint64_t fnv_offset = 14695981039346656037ULL;
constexpr std::uint64_t fnv_prime = 1099511628211ULL;
void hash_bytes(std::uint64_t& hash, const void* data, std::size_t size) {
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; ++index) {
        hash ^= static_cast<std::uint64_t>(bytes[index]);
        hash *= fnv_prime;
    }
}
void hash_size(std::uint64_t& hash, std::size_t value) {
    const std::uint64_t encoded = static_cast<std::uint64_t>(value);
    hash_bytes(hash, &encoded, sizeof(encoded));
}
void hash_integer(std::uint64_t& hash, std::int64_t value) { hash_bytes(hash, &value, sizeof(value)); }
std::uint64_t encode_double_bits(double value) {
    std::uint64_t encoded = 0;
    static_assert(sizeof(encoded) == sizeof(value));
    std::memcpy(&encoded, &value, sizeof(value));
    return encoded;
}
double decode_double_bits(std::uint64_t encoded) {
    double result = 0.0;
    static_assert(sizeof(encoded) == sizeof(result));
    std::memcpy(&result, &encoded, sizeof(result));
    return result;
}
void hash_double(std::uint64_t& hash, double value) {
    const std::uint64_t encoded = encode_double_bits(value);
    hash_bytes(hash, &encoded, sizeof(encoded));
}
void hash_string(std::uint64_t& hash, const std::string& value) {
    hash_size(hash, value.size());
    hash_bytes(hash, value.data(), value.size());
}
void hash_thermoelastic(std::uint64_t& hash, const ThermoelasticProperties& material) {
    hash_double(hash, material.reference_young_modulus);
    const std::uint64_t signature = material.functions->signature();
    hash_bytes(hash, &signature, sizeof(signature));
}
void hash_region_definition(std::uint64_t& hash, const RegionDefinition& spatial) {
    hash_string(hash, spatial.name);
    hash_string(hash, spatial.block);
    hash_integer(hash, spatial.block_id);
    hash_thermoelastic(hash, spatial.material);
    hash_double(hash, spatial.volumetric_heat_source);
    hash_double(hash, spatial.initial_temperature);
    hash_string(hash, spatial.heat_source_function);
}
void hash_boundaries(std::uint64_t& hash, const SpatialDefinition& definition, bool include_displaced_geometry) {
    for (const BoundaryConditionDefinition& boundary : definition.boundary_conditions) {
        hash_string(hash, boundary.name);
        hash_integer(hash, static_cast<std::int64_t>(boundary.type));
        hash_string(hash, boundary.boundary);
        hash_integer(hash, static_cast<std::int64_t>(boundary.field));
        hash_double(hash, boundary.value);
        hash_integer(hash, boundary.scale_with_load ? 1 : 0);
        hash_string(hash, boundary.function);
        hash_double(hash, boundary.heat_transfer_coefficient);
        hash_double(hash, boundary.ambient_temperature);
        hash_string(hash, boundary.coefficient_function);
        hash_string(hash, boundary.ambient_temperature_function);
        if (include_displaced_geometry) hash_integer(hash, boundary.use_displaced_geometry ? 1 : 0);
    }
}
void hash_time_tables(std::uint64_t& hash, const SpatialDefinition& definition) {
    for (const PiecewiseLinearTimeTable& table : definition.time_tables) {
        hash_string(hash, table.name());
        hash_size(hash, table.times().size());
        for (std::size_t entry = 0; entry < table.times().size(); ++entry) {
            hash_double(hash, table.times()[entry]);
            hash_double(hash, table.values()[entry]);
        }
    }
}
void hash_contribution_layout(std::uint64_t& hash, const TransientProblem& problem) {
    hash_size(hash, problem.contribution_count());
    std::vector<std::size_t> dofs;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        problem.contribution_dofs(contribution, dofs);
        hash_size(hash, dofs.size());
        for (const std::size_t dof : dofs) hash_size(hash, dof);
    }
}
} // namespace
std::uint64_t transient_problem_signature(const TransientProblem& problem) {
    std::uint64_t hash = fnv_offset;
    const bool cartesian = problem.is_cartesian_3d();
    hash_string(hash, cartesian ? "cartesian_3d_hex8" : "axisymmetric_rz_quad4");
    hash_string(hash, cartesian ? "xx,yy,zz,xy,yz,xz" : "rr,zz,hoop,rz");
    hash_size(hash, cartesian ? 8 : 4);
    hash_size(hash, problem.dof_count());
    if (cartesian) {
        const cartesian::SpatialAssembly& assembly = cartesian::BackendAccess::spatial(problem);
        hash_size(hash, assembly.region_count());
        const SpatialDefinition& definition = problem.definition();
        for (std::size_t region = 0; region < assembly.region_count(); ++region) {
            const RegionDefinition& spatial = definition.regions[region];
            hash_region_definition(hash, spatial);
            const Hex8RegionMesh& mesh = assembly.region_mesh(region);
            hash_size(hash, mesh.nodes().size());
            for (const CartesianPoint3& point : mesh.nodes()) {
                hash_double(hash, point.x);
                hash_double(hash, point.y);
                hash_double(hash, point.z);
            }
            hash_size(hash, mesh.elements().size());
            for (const Hex8Element& element : mesh.elements())
                for (const std::size_t node : element.nodes) hash_size(hash, node);
            for (const std::size_t source : mesh.source_node_ids()) hash_size(hash, source);
            for (const std::size_t source : mesh.source_element_ids()) hash_size(hash, source);
        }
        hash_boundaries(hash, definition, false);
        hash_time_tables(hash, definition);
        hash_contribution_layout(hash, problem);
        return hash;
    }
    const rz::TransientBackendView backend = rz::BackendAccess::transient(problem);
    hash_size(hash, backend.spatial.region_count());
    const SpatialDefinition& definition = problem.definition();
    for (std::size_t region = 0; region < backend.spatial.region_count(); ++region) {
        const RegionDefinition& spatial = definition.regions[region];
        hash_region_definition(hash, spatial);
        hash_integer(hash, static_cast<std::int64_t>(spatial.strain_formulation));
        const RegionMesh& mesh = backend.spatial.region_mesh(region);
        hash_size(hash, mesh.nodes().size());
        for (const RzPoint& point : mesh.nodes()) {
            hash_double(hash, point.r);
            hash_double(hash, point.z);
        }
        hash_size(hash, mesh.elements().size());
        for (const Quad4Element& element : mesh.elements())
            for (const std::size_t node : element.nodes) hash_size(hash, node);
        for (const std::size_t source : mesh.source_node_ids()) hash_size(hash, source);
        for (const std::size_t source : mesh.source_element_ids()) hash_size(hash, source);
    }
    for (const ContactDefinition& contact : definition.contacts) {
        hash_string(hash, contact.name);
        hash_string(hash, contact.primary);
        hash_string(hash, contact.secondary);
        hash_integer(hash, contact.thermal ? 1 : 0);
        hash_integer(hash, contact.mechanical ? 1 : 0);
        hash_double(hash, contact.gap_conductivity);
        hash_double(hash, contact.minimum_gap);
        hash_double(hash, contact.penalty);
        hash_double(hash, contact.friction_coefficient);
        hash_integer(hash, contact.automatic_penalty ? 1 : 0);
        hash_double(hash, contact.penalty_factor);
        hash_integer(hash, static_cast<std::int64_t>(contact.mechanical_formulation));
        hash_double(hash, contact.penetration_tolerance);
        hash_size(hash, contact.maximum_augmented_iterations);
    }
    hash_boundaries(hash, definition, true);
    hash_time_tables(hash, definition);
    hash_contribution_layout(hash, problem);
    return hash;
}
namespace {
using exodus_detail::check_exodus;
using exodus_detail::ExodusFile;
std::size_t checked_size(std::int64_t value, const char* description) {
    if (value < 0 || static_cast<std::uint64_t>(value) > std::numeric_limits<std::size_t>::max())
        throw std::length_error(std::string(description) + " is out of range");
    return static_cast<std::size_t>(value);
}
std::int64_t checked_count(std::size_t value, const char* description) {
    if (value > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error(std::string(description) + " is out of range");
    return static_cast<std::int64_t>(value);
}
std::string normalized_topology(const char* topology) {
    std::string result(topology);
    std::transform(result.begin(), result.end(), result.begin(),
        [](unsigned char character) { return static_cast<char>(std::toupper(character)); });
    return result;
}
struct BlockConnectivity final {
    std::int64_t _id;
    std::string _name;
    std::vector<std::int64_t> _nodes;
    std::vector<std::size_t> _source_elements;
};
std::string read_entity_name(int exoid, ex_entity_type type, std::int64_t id, std::size_t maximum_length) {
    std::vector<char> name(maximum_length + 1U, '\0');
    check_exodus(ex_get_name(exoid, type, id, name.data()), "Could not read Exodus entity name");
    return std::string(name.data());
}
std::int64_t to_exodus_id(std::size_t index, const char* description) {
    if (index >= static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::length_error(std::string(description) + " is out of range");
    return static_cast<std::int64_t>(index) + 1;
}
void write_entity_name(
    int exoid, ex_entity_type type, std::int64_t id, const std::string& name, const std::string& operation) {
    if (!name.empty()) check_exodus(ex_put_name(exoid, type, id, name.c_str()), operation);
}
std::vector<NodeSet> read_node_sets(int exoid, std::size_t set_count, std::size_t maximum_name_length) {
    std::vector<NodeSet> result;
    if (set_count == 0) return result;
    std::vector<std::int64_t> set_ids(set_count);
    check_exodus(ex_get_ids(exoid, EX_NODE_SET, set_ids.data()), "Could not read Exodus node set IDs");
    result.reserve(set_count);
    for (const std::int64_t set_id : set_ids) {
        std::int64_t entry_count = 0, factor_count = 0;
        check_exodus(ex_get_set_param(exoid, EX_NODE_SET, set_id, &entry_count, &factor_count),
            "Could not read Exodus node set parameters");
        std::vector<std::int64_t> entries(checked_size(entry_count, "Exodus node set size"));
        if (!entries.empty())
            check_exodus(
                ex_get_set(exoid, EX_NODE_SET, set_id, entries.data(), nullptr), "Could not read Exodus node set");
        NodeSet set{set_id, read_entity_name(exoid, EX_NODE_SET, set_id, maximum_name_length), {}};
        set.nodes.reserve(entries.size());
        for (const std::int64_t entry : entries) {
            if (entry <= 0) throw std::runtime_error("Exodus node set contains an invalid node ID");
            set.nodes.push_back(static_cast<std::size_t>(entry - 1));
        }
        result.push_back(std::move(set));
    }
    return result;
}
std::vector<SideSet> read_side_sets(int exoid, std::size_t set_count, std::size_t maximum_name_length,
    std::int64_t maximum_side, const char* topology_name) {
    std::vector<SideSet> result;
    if (set_count == 0) return result;
    std::vector<std::int64_t> set_ids(set_count);
    check_exodus(ex_get_ids(exoid, EX_SIDE_SET, set_ids.data()), "Could not read Exodus side set IDs");
    result.reserve(set_count);
    for (const std::int64_t set_id : set_ids) {
        std::int64_t entry_count = 0, factor_count = 0;
        check_exodus(ex_get_set_param(exoid, EX_SIDE_SET, set_id, &entry_count, &factor_count),
            "Could not read Exodus side set parameters");
        const std::size_t count = checked_size(entry_count, "Exodus side set size");
        std::vector<std::int64_t> elements(count);
        std::vector<std::int64_t> sides(count);
        if (count != 0)
            check_exodus(ex_get_set(exoid, EX_SIDE_SET, set_id, elements.data(), sides.data()),
                "Could not read Exodus side set");
        SideSet set{set_id, read_entity_name(exoid, EX_SIDE_SET, set_id, maximum_name_length), {}};
        set.sides.reserve(count);
        for (std::size_t entry = 0; entry < count; ++entry) {
            if (elements[entry] <= 0 || sides[entry] <= 0 || sides[entry] > maximum_side)
                throw std::runtime_error(std::string("Exodus side set contains an invalid ") + topology_name + " side");
            set.sides.push_back(
                {static_cast<std::size_t>(elements[entry] - 1), static_cast<std::size_t>(sides[entry] - 1)});
        }
        result.push_back(std::move(set));
    }
    return result;
}
struct ExodusMeshData final {
    std::size_t dimension, nodes_per_element, maximum_side;
    const char* topology;
    const char* alternate_topology;
    const char* title;
    std::array<const char*, 3> coordinate_names;
    std::vector<std::array<double, 3>> nodes;
    std::vector<std::array<std::size_t, 8>> elements;
    std::vector<std::int64_t> element_block_ids;
    std::vector<ElementBlockInfo> element_blocks;
    std::vector<NodeSet> node_sets;
    std::vector<SideSet> side_sets;
};
ExodusMeshData read_exodus_mesh(const std::string& path, std::size_t dimension, std::size_t nodes_per_element,
    std::size_t maximum_side, const char* topology, const char* alternate_topology, const char* title) {
    int cpu_word_size = static_cast<int>(sizeof(double)), io_word_size = 0;
    float version = 0.0F;
    const int exoid = ex_open(path.c_str(), EX_READ, &cpu_word_size, &io_word_size, &version);
    if (exoid < 0) throw std::runtime_error("Could not open Exodus file '" + path + "': " + ex_strerror(exoid));
    ExodusFile file(exoid);
    ex_set_int64_status(file.id(), EX_ALL_INT64_API);
    const std::int64_t maximum_name_length = ex_inquire_int(file.id(), EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);
    if (maximum_name_length < 0) throw std::runtime_error("Could not read Exodus maximum name length");
    check_exodus(ex_set_max_name_length(file.id(), static_cast<int>(maximum_name_length)),
        "Could not set Exodus maximum read name length");
    const std::size_t name_length = checked_size(maximum_name_length, "Exodus maximum name length");
    ex_init_params parameters{};
    check_exodus(ex_get_init_ext(file.id(), &parameters), "Could not read Exodus model parameters");
    if (parameters.num_dim != static_cast<std::int64_t>(dimension))
        throw std::runtime_error(std::string("Exodus ") + topology + " mesh has the wrong dimension");
    const std::size_t node_count = checked_size(parameters.num_nodes, "Exodus node count"),
                      block_count = checked_size(parameters.num_elem_blk, "Exodus element block count");
    if (node_count == 0 || block_count == 0)
        throw std::runtime_error(std::string("Exodus ") + topology + " mesh must contain nodes and element blocks");
    ExodusMeshData result{
        dimension, nodes_per_element, maximum_side, topology, alternate_topology, title, {}, {}, {}, {}, {}, {}, {}};
    std::array<std::vector<double>, 3> coordinates = {
        std::vector<double>(node_count), std::vector<double>(node_count), std::vector<double>(node_count)};
    check_exodus(ex_get_coord(file.id(), coordinates[0].data(), coordinates[1].data(),
                     dimension == 3 ? coordinates[2].data() : nullptr),
        "Could not read Exodus coordinates");
    result.nodes.resize(node_count);
    for (std::size_t node = 0; node < node_count; ++node)
        result.nodes[node] = {coordinates[0][node], coordinates[1][node], dimension == 3 ? coordinates[2][node] : 0.0};
    std::vector<std::int64_t> block_ids(block_count);
    check_exodus(ex_get_ids(file.id(), EX_ELEM_BLOCK, block_ids.data()), "Could not read Exodus element block IDs");
    result.element_blocks.reserve(block_count);
    result.elements.reserve(checked_size(parameters.num_elem, "Exodus element count"));
    result.element_block_ids.reserve(result.elements.capacity());
    for (const std::int64_t block_id : block_ids) {
        result.element_blocks.push_back({block_id, read_entity_name(file.id(), EX_ELEM_BLOCK, block_id, name_length)});
        ex_block block{};
        block.id = block_id;
        block.type = EX_ELEM_BLOCK;
        check_exodus(ex_get_block_param(file.id(), &block), "Could not read Exodus element block");
        const std::string actual_topology = normalized_topology(block.topology);
        if (block.num_nodes_per_entry != static_cast<std::int64_t>(nodes_per_element) ||
            (actual_topology != topology && actual_topology != alternate_topology))
            throw std::runtime_error(std::string("Exodus element blocks must contain only ") + topology + " elements");
        const std::size_t block_elements = checked_size(block.num_entry, "Exodus block element count");
        if (block_elements > std::numeric_limits<std::size_t>::max() / nodes_per_element)
            throw std::length_error("Exodus connectivity is too large");
        std::vector<std::int64_t> connectivity(block_elements * nodes_per_element);
        check_exodus(ex_get_conn(file.id(), EX_ELEM_BLOCK, block_id, connectivity.data(), nullptr, nullptr),
            "Could not read Exodus connectivity");
        for (std::size_t element = 0; element < block_elements; ++element) {
            std::array<std::size_t, 8> nodes{};
            for (std::size_t local = 0; local < nodes_per_element; ++local) {
                const std::int64_t exodus_node = connectivity[element * nodes_per_element + local];
                if (exodus_node <= 0 || static_cast<std::uint64_t>(exodus_node) > node_count)
                    throw std::runtime_error("Exodus connectivity contains an invalid node ID");
                nodes[local] = static_cast<std::size_t>(exodus_node - 1);
            }
            result.elements.push_back(nodes);
            result.element_block_ids.push_back(block_id);
        }
    }
    if (result.elements.size() != checked_size(parameters.num_elem, "Exodus element count"))
        throw std::runtime_error("Exodus element count does not match its element blocks");
    result.node_sets =
        read_node_sets(file.id(), checked_size(parameters.num_node_sets, "Exodus node set count"), name_length);
    result.side_sets = read_side_sets(file.id(), checked_size(parameters.num_side_sets, "Exodus side set count"),
        name_length, static_cast<std::int64_t>(maximum_side), topology);
    file.close();
    return result;
}
void write_exodus_mesh(const std::string& path, const ExodusMeshData& mesh) {
    int cpu_word_size = static_cast<int>(sizeof(double)), io_word_size = static_cast<int>(sizeof(double));
    const int exoid =
        ex_create(path.c_str(), EX_CLOBBER | EX_ALL_INT64_DB | EX_ALL_INT64_API, &cpu_word_size, &io_word_size);
    if (exoid < 0) throw std::runtime_error("Could not create Exodus file '" + path + "': " + ex_strerror(exoid));
    ExodusFile file(exoid);
    ex_set_int64_status(file.id(), EX_ALL_INT64_API);
    std::vector<BlockConnectivity> blocks;
    for (const ElementBlockInfo& block : mesh.element_blocks) blocks.push_back({block.id, block.name, {}, {}});
    for (std::size_t element = 0; element < mesh.elements.size(); ++element) {
        const std::int64_t block_id = mesh.element_block_ids.at(element);
        auto block = std::find_if(blocks.begin(), blocks.end(),
            [block_id](const BlockConnectivity& candidate) { return candidate._id == block_id; });
        if (block == blocks.end()) throw std::logic_error("Mesh element references an unknown block");
        block->_source_elements.push_back(element);
        for (std::size_t local = 0; local < mesh.nodes_per_element; ++local)
            block->_nodes.push_back(to_exodus_id(mesh.elements[element][local], "Exodus node ID"));
    }
    ex_init_params parameters{};
    ex_copy_string(parameters.title, mesh.title, sizeof(parameters.title));
    parameters.num_dim = checked_count(mesh.dimension, "Dimension");
    parameters.num_nodes = checked_count(mesh.nodes.size(), "Node count");
    parameters.num_elem = checked_count(mesh.elements.size(), "Element count");
    parameters.num_elem_blk = checked_count(blocks.size(), "Block count");
    parameters.num_node_sets = checked_count(mesh.node_sets.size(), "Node set count");
    parameters.num_side_sets = checked_count(mesh.side_sets.size(), "Side set count");
    check_exodus(ex_put_init_ext(file.id(), &parameters), "Could not write Exodus model parameters");
    std::array<std::vector<double>, 3> coordinates;
    for (std::size_t component = 0; component < mesh.dimension; ++component) {
        coordinates[component].reserve(mesh.nodes.size());
        for (const auto& node : mesh.nodes) coordinates[component].push_back(node[component]);
    }
    check_exodus(ex_put_coord(file.id(), coordinates[0].data(), coordinates[1].data(),
                     mesh.dimension == 3 ? coordinates[2].data() : nullptr),
        "Could not write Exodus coordinates");
    std::array<char*, 3> coordinate_names{};
    for (std::size_t component = 0; component < mesh.dimension; ++component)
        coordinate_names[component] = const_cast<char*>(mesh.coordinate_names[component]);
    check_exodus(ex_put_coord_names(file.id(), coordinate_names.data()), "Could not write Exodus coordinate names");
    std::vector<std::int64_t> written_element_ids(mesh.elements.size(), 0);
    std::int64_t next_element_id = 1;
    for (const BlockConnectivity& block : blocks) {
        check_exodus(ex_put_block(file.id(), EX_ELEM_BLOCK, block._id, mesh.topology,
                         checked_count(block._source_elements.size(), "Block element count"),
                         checked_count(mesh.nodes_per_element, "Nodes per element"), 0, 0, 0),
            "Could not write Exodus element block");
        check_exodus(ex_put_conn(file.id(), EX_ELEM_BLOCK, block._id, block._nodes.data(), nullptr, nullptr),
            "Could not write Exodus connectivity");
        write_entity_name(
            file.id(), EX_ELEM_BLOCK, block._id, block._name, "Could not write Exodus element block name");
        for (const std::size_t source : block._source_elements) written_element_ids[source] = next_element_id++;
    }
    for (const NodeSet& set : mesh.node_sets) {
        std::vector<std::int64_t> entries;
        for (const std::size_t node : set.nodes) entries.push_back(to_exodus_id(node, "Exodus node set ID"));
        check_exodus(
            ex_put_set_param(file.id(), EX_NODE_SET, set.id, checked_count(entries.size(), "Node set size"), 0),
            "Could not write Exodus node set parameters");
        if (!entries.empty())
            check_exodus(
                ex_put_set(file.id(), EX_NODE_SET, set.id, entries.data(), nullptr), "Could not write Exodus node set");
        write_entity_name(file.id(), EX_NODE_SET, set.id, set.name, "Could not write Exodus node set name");
    }
    for (const SideSet& set : mesh.side_sets) {
        std::vector<std::int64_t> elements, sides;
        for (const ElementSide& side : set.sides) {
            elements.push_back(written_element_ids.at(side.element));
            sides.push_back(static_cast<std::int64_t>(side.local_side) + 1);
        }
        check_exodus(
            ex_put_set_param(file.id(), EX_SIDE_SET, set.id, checked_count(set.sides.size(), "Side set size"), 0),
            "Could not write Exodus side set parameters");
        if (!elements.empty())
            check_exodus(ex_put_set(file.id(), EX_SIDE_SET, set.id, elements.data(), sides.data()),
                "Could not write Exodus side set");
        write_entity_name(file.id(), EX_SIDE_SET, set.id, set.name, "Could not write Exodus side set name");
    }
    file.close();
}
} // namespace
UnstructuredQuad4Mesh read_exodus_quad4(const std::string& path) {
    ExodusMeshData data = read_exodus_mesh(path, 2, 4, 4, "QUAD4", "QUAD", "fuelsim Quad4 mesh");
    std::vector<RzPoint> nodes;
    for (const auto& node : data.nodes) nodes.push_back({node[0], node[1]});
    std::vector<Quad4Element> elements;
    for (const auto& value : data.elements) elements.push_back({{{value[0], value[1], value[2], value[3]}}});
    return UnstructuredQuad4Mesh(std::move(nodes), std::move(elements), std::move(data.element_block_ids),
        std::move(data.element_blocks), std::move(data.node_sets), std::move(data.side_sets));
}
void write_exodus_quad4(const std::string& path, const UnstructuredQuad4Mesh& mesh) {
    ExodusMeshData data{2, 4, 4, "QUAD4", "QUAD", "fuelsim Quad4 mesh", {"r", "z", nullptr}, {}, {},
        mesh.element_block_ids(), mesh.element_blocks(), mesh.node_sets(), mesh.side_sets()};
    for (const RzPoint& node : mesh.nodes()) data.nodes.push_back({node.r, node.z, 0.0});
    for (const Quad4Element& element : mesh.elements())
        data.elements.push_back({element.nodes[0], element.nodes[1], element.nodes[2], element.nodes[3]});
    write_exodus_mesh(path, data);
}
UnstructuredHex8Mesh read_exodus_hex8(const std::string& path) {
    ExodusMeshData data = read_exodus_mesh(path, 3, 8, 6, "HEX8", "HEX", "fuelsim HEX8 mesh");
    std::vector<CartesianPoint3> nodes;
    for (const auto& node : data.nodes) nodes.push_back({node[0], node[1], node[2]});
    std::vector<Hex8Element> elements;
    for (const auto& value : data.elements) elements.push_back({value});
    return UnstructuredHex8Mesh(std::move(nodes), std::move(elements), std::move(data.element_block_ids),
        std::move(data.element_blocks), std::move(data.node_sets), std::move(data.side_sets));
}
void write_exodus_hex8(const std::string& path, const UnstructuredHex8Mesh& mesh) {
    ExodusMeshData data{3, 8, 6, "HEX8", "HEX", "fuelsim HEX8 mesh", {"x", "y", "z"}, {}, {}, mesh.element_block_ids(),
        mesh.element_blocks(), mesh.node_sets(), mesh.side_sets()};
    for (const CartesianPoint3& node : mesh.nodes()) data.nodes.push_back({node.x, node.y, node.z});
    for (const Hex8Element& element : mesh.elements()) data.elements.push_back(element.nodes);
    write_exodus_mesh(path, data);
}
namespace {
using exodus_detail::check_exodus;
using exodus_detail::ExodusFile;
constexpr std::array<const char*, 4> stress_components = {"rr", "zz", "hoop", "rz"};
constexpr std::array<const char*, 6> cartesian_stress_components = {"xx", "yy", "zz", "xy", "yz", "xz"};
int checked_int(std::size_t value, const std::string& quantity) {
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error(quantity + " exceeds the Exodus int range");
    return static_cast<int>(value);
}
ExodusFile open_results(const std::string& path) {
    int cpu_word_size = static_cast<int>(sizeof(double)), io_word_size = 0;
    float version = 0.0F;
    const int exoid = ex_open(path.c_str(), EX_WRITE, &cpu_word_size, &io_word_size, &version);
    if (exoid < 0) throw std::runtime_error("Could not open Exodus results file '" + path + "': " + ex_strerror(exoid));
    ex_set_int64_status(exoid, EX_ALL_INT64_API);
    return ExodusFile(exoid);
}
std::vector<char*> variable_name_pointers(const std::vector<std::string>& names) {
    std::vector<char*> result;
    result.reserve(names.size());
    for (const std::string& name : names) result.push_back(const_cast<char*>(name.c_str()));
    return result;
}
void define_variable_names(
    int exoid, ex_entity_type type, const std::vector<std::string>& names, const std::string& category) {
    std::vector<char*> pointers = variable_name_pointers(names);
    check_exodus(ex_put_variable_param(exoid, type, static_cast<int>(pointers.size())),
        "Could not define Exodus " + category + " variables");
    check_exodus(ex_put_variable_names(exoid, type, static_cast<int>(pointers.size()), pointers.data()),
        "Could not name Exodus " + category + " variables");
}
std::vector<std::string> nodal_variable_names(const std::vector<ContactDefinition>& contacts) {
    std::vector<std::string> result = {"temperature", "displacement_r", "displacement_z"};
    for (const ContactDefinition& contact : contacts) {
        result.push_back("contact_gap_" + contact.name);
        result.push_back("contact_pressure_" + contact.name);
        result.push_back("contact_tangential_traction_" + contact.name);
        result.push_back("contact_elastic_tangential_slip_" + contact.name);
        result.push_back("contact_sliding_" + contact.name);
    }
    return result;
}
std::vector<std::string> cartesian_nodal_variable_names() {
    return {"temperature", "displacement_x", "displacement_y", "displacement_z"};
}
std::vector<std::string> global_variable_names(const std::vector<ContactDefinition>& contacts) {
    std::vector<std::string> result = {"load_factor"};
    for (const ContactDefinition& contact : contacts) {
        result.push_back("contact_heat_rate_" + contact.name);
        result.push_back("contact_force_" + contact.name);
        result.push_back("contact_tangential_force_" + contact.name);
    }
    return result;
}
void append_component_variable_names(std::vector<std::string>& result, const char* prefix, std::size_t q) {
    for (const char* component : stress_components)
        result.push_back(std::string(prefix) + std::string(component) + "_q" + std::to_string(q));
}
std::vector<std::string> stress_variable_names() {
    std::vector<std::string> result;
    result.reserve(16);
    for (std::size_t q = 0; q < 4; ++q) append_component_variable_names(result, "stress_", q);
    return result;
}
std::vector<std::string> transient_element_variable_names() {
    std::vector<std::string> result = stress_variable_names();
    result.reserve(56);
    for (std::size_t q = 0; q < 4; ++q) {
        append_component_variable_names(result, "plastic_", q);
        append_component_variable_names(result, "creep_", q);
        result.push_back("equiv_plastic_q" + std::to_string(q));
        result.push_back("equiv_creep_q" + std::to_string(q));
    }
    return result;
}
std::vector<std::string> cartesian_stress_variable_names() {
    std::vector<std::string> result;
    result.reserve(48);
    for (std::size_t q = 0; q < 8; ++q)
        for (const char* component : cartesian_stress_components)
            result.push_back("stress_" + std::string(component) + "_q" + std::to_string(q));
    return result;
}
struct ResultsMeshView {
    std::size_t node_count, element_count;
    const std::vector<ElementBlockInfo>& blocks;
    const std::vector<std::int64_t>& element_block_ids;
};
ResultsMeshView results_mesh_view(const UnstructuredQuad4Mesh& mesh) {
    return {mesh.nodes().size(), mesh.elements().size(), mesh.element_blocks(), mesh.element_block_ids()};
}
ResultsMeshView results_mesh_view(const UnstructuredHex8Mesh& mesh) {
    return {mesh.nodes().size(), mesh.elements().size(), mesh.element_blocks(), mesh.element_block_ids()};
}
void define_result_variables(const std::string& path, const ResultsMeshView& mesh,
    const std::vector<std::string>& nodal_names, const std::vector<std::string>& element_names,
    const std::vector<std::string>& global_names) {
    ExodusFile file = open_results(path);
    check_exodus(ex_set_max_name_length(file.id(), 64), "Could not set Exodus result-name length");
    define_variable_names(file.id(), EX_NODAL, nodal_names, "nodal");
    define_variable_names(file.id(), EX_GLOBAL, global_names, "global");
    define_variable_names(file.id(), EX_ELEM_BLOCK, element_names, "element");
    std::vector<int> truth(mesh.blocks.size() * element_names.size(), 1);
    check_exodus(ex_put_truth_table(file.id(), EX_ELEM_BLOCK, static_cast<int>(mesh.blocks.size()),
                     static_cast<int>(element_names.size()), truth.data()),
        "Could not define Exodus element-variable truth table");
    file.close();
}
void define_variables(const std::string& path, const UnstructuredQuad4Mesh& mesh, std::vector<std::string> nodal_names,
    std::vector<std::string> element_names, std::vector<std::string> global_names) {
    write_exodus_quad4(path, mesh);
    define_result_variables(path, results_mesh_view(mesh), nodal_names, element_names, global_names);
}
void define_variables(const std::string& path, const UnstructuredHex8Mesh& mesh, std::vector<std::string> nodal_names,
    std::vector<std::string> element_names, std::vector<std::string> global_names) {
    write_exodus_hex8(path, mesh);
    define_result_variables(path, results_mesh_view(mesh), nodal_names, element_names, global_names);
}
std::vector<std::size_t> block_elements(const ResultsMeshView& mesh, std::int64_t block_id) {
    std::vector<std::size_t> result;
    for (std::size_t element = 0; element < mesh.element_count; ++element)
        if (mesh.element_block_ids[element] == block_id) result.push_back(element);
    return result;
}
void write_result_step(const std::string& path, const ResultsMeshView& mesh, std::size_t step, double time,
    const std::vector<std::vector<double>>& nodal_values, const std::vector<std::vector<double>>& element_values,
    const std::vector<double>& global_values) {
    if (step == 0 || nodal_values.empty() || element_values.empty() || global_values.empty())
        throw std::invalid_argument("Exodus result step is incomplete");
    ExodusFile file = open_results(path);
    const int exodus_step = checked_int(step, "Exodus result step");
    check_exodus(ex_put_time(file.id(), exodus_step, &time), "Could not write Exodus result time");
    check_exodus(ex_put_var(file.id(), exodus_step, EX_GLOBAL, 1, 0, static_cast<std::int64_t>(global_values.size()),
                     global_values.data()),
        "Could not write Exodus global results");
    for (std::size_t variable = 0; variable < nodal_values.size(); ++variable) {
        if (nodal_values[variable].size() != mesh.node_count)
            throw std::invalid_argument("Exodus nodal result size does not match mesh");
        check_exodus(ex_put_var(file.id(), exodus_step, EX_NODAL, static_cast<int>(variable + 1), 1,
                         static_cast<std::int64_t>(mesh.node_count), nodal_values[variable].data()),
            "Could not write Exodus nodal results");
    }
    std::vector<std::vector<std::size_t>> block_element_lists;
    block_element_lists.reserve(mesh.blocks.size());
    for (const ElementBlockInfo& block : mesh.blocks) block_element_lists.push_back(block_elements(mesh, block.id));
    for (std::size_t variable = 0; variable < element_values.size(); ++variable) {
        if (element_values[variable].size() != mesh.element_count)
            throw std::invalid_argument("Exodus element result size does not match mesh");
        for (std::size_t block = 0; block < mesh.blocks.size(); ++block) {
            const std::vector<std::size_t>& elements = block_element_lists[block];
            std::vector<double> values;
            values.reserve(elements.size());
            for (const std::size_t element : elements) values.push_back(element_values[variable][element]);
            check_exodus(ex_put_var(file.id(), exodus_step, EX_ELEM_BLOCK, static_cast<int>(variable + 1),
                             mesh.blocks[block].id, static_cast<std::int64_t>(values.size()), values.data()),
                "Could not write Exodus element results");
        }
    }
    check_exodus(ex_update(file.id()), "Could not flush Exodus results");
    file.close();
}
void fill_region_nodal_values(const std::vector<std::size_t>& source_nodes, std::size_t region_offset,
    const DofMap& dof_map, const std::vector<double>& state, std::vector<bool>& present,
    std::vector<std::vector<double>>& values) {
    for (std::size_t local = 0; local < source_nodes.size(); ++local) {
        const std::size_t source = source_nodes[local];
        if (present.at(source)) throw std::invalid_argument("Exodus result mapping contains a shared source node");
        present[source] = true;
        const std::size_t global = region_offset + local;
        const std::vector<FieldDescriptor>& fields = dof_map.field_layout();
        for (std::size_t field = 0; field < fields.size(); ++field)
            values[field][source] = state.at(fields[field].begin + global);
    }
}
void fill_contact_nodal_values(std::size_t contact, const std::vector<std::size_t>& nodes,
    const std::vector<ContactNodeSummary>& summary, std::vector<std::vector<double>>& values) {
    if (nodes.size() != summary.size()) throw std::logic_error("Contact result mapping size mismatch");
    const double missing = std::numeric_limits<double>::quiet_NaN();
    const std::size_t base = 3 + 5 * contact;
    for (std::size_t node = 0; node < nodes.size(); ++node) {
        values[base].at(nodes[node]) = summary[node].projected ? summary[node].gap : missing;
        values[base + 1].at(nodes[node]) = summary[node].projected ? summary[node].pressure : missing;
        values[base + 2].at(nodes[node]) = summary[node].projected ? summary[node].tangential_traction : missing;
        values[base + 3].at(nodes[node]) = summary[node].projected ? summary[node].elastic_tangential_slip : missing;
        values[base + 4].at(nodes[node]) = summary[node].projected ? (summary[node].sliding ? 1.0 : 0.0) : missing;
    }
}
void fill_rz_nodal(const UnstructuredQuad4Mesh& mesh, const rz::SpatialAssembly& spatial,
    const std::vector<double>& state, std::vector<std::vector<double>>& values) {
    const double missing = std::numeric_limits<double>::quiet_NaN();
    values.assign(
        nodal_variable_names(spatial.definition().contacts).size(), std::vector<double>(mesh.nodes().size(), missing));
    std::vector<bool> present(mesh.nodes().size(), false);
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        fill_region_nodal_values(spatial.region_mesh(region).source_node_ids(), spatial.region_node_offset(region),
            spatial.dof_map(), state, present, values);
    for (std::size_t contact = 0; contact < spatial.contact_count(); ++contact) {
        const std::vector<std::size_t> nodes = spatial.contact_secondary_source_nodes(contact);
        const std::vector<ContactNodeSummary> summary = spatial.summarize_contact_nodes(contact, state);
        fill_contact_nodal_values(contact, nodes, summary, values);
    }
}
std::vector<double> rz_globals(
    const rz::SpatialAssembly& spatial, const std::vector<double>& state, double load_factor) {
    std::vector<double> result = {load_factor};
    for (std::size_t contact = 0; contact < spatial.contact_count(); ++contact) {
        const InterfaceSummary summary = spatial.summarize_interface(contact, state);
        result.push_back(summary.total_heat_rate);
        result.push_back(summary.total_contact_force);
        result.push_back(summary.total_tangential_force);
    }
    return result;
}
void store_stress_values(std::size_t source, const std::array<AxisymmetricStressValues, 4>& stresses,
    std::vector<std::vector<double>>& values) {
    for (std::size_t q = 0; q < stresses.size(); ++q) {
        const std::size_t offset = 4 * q;
        values[offset][source] = stresses[q].rr;
        values[offset + 1][source] = stresses[q].zz;
        values[offset + 2][source] = stresses[q].hoop;
        values[offset + 3][source] = stresses[q].rz;
    }
}
std::vector<std::vector<double>> steady_elements(
    const UnstructuredQuad4Mesh& mesh, const SteadyProblem& problem, const std::vector<double>& state) {
    const rz::SteadyBackendView backend = rz::BackendAccess::steady(problem);
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> result(
        stress_variable_names().size(), std::vector<double>(mesh.elements().size(), missing));
    for (std::size_t region = 0; region < backend.spatial.region_count(); ++region) {
        const RegionMesh& region_mesh = backend.spatial.region_mesh(region);
        const std::size_t contribution_offset = backend.spatial.region_element_offset(region);
        for (std::size_t element = 0; element < region_mesh.elements().size(); ++element) {
            const LocalDofs dofs = backend.spatial.contribution_dofs(contribution_offset + element);
            LocalValues local{};
            for (std::size_t index = 0; index < dofs.size(); ++index) local[index] = state.at(dofs[index]);
            const auto stresses =
                backend.kernels[region].stress_values(backend.spatial.region_element_geometry(region, element), local);
            const std::size_t source = region_mesh.source_element_ids().at(element);
            store_stress_values(source, stresses, result);
        }
    }
    return result;
}
std::vector<std::vector<double>> transient_elements(
    const UnstructuredQuad4Mesh& mesh, const TransientProblem& problem) {
    const rz::TransientBackendView backend = rz::BackendAccess::transient(problem);
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> result(
        transient_element_variable_names().size(), std::vector<double>(mesh.elements().size(), missing));
    for (std::size_t region = 0; region < backend.spatial.region_count(); ++region) {
        const RegionMesh& region_mesh = backend.spatial.region_mesh(region);
        for (std::size_t element = 0; element < region_mesh.elements().size(); ++element) {
            const std::size_t source = region_mesh.source_element_ids().at(element);
            const Quad4MaterialHistory& history = backend.histories.at(region).at(element);
            const auto& stresses = backend.stresses.at(region).at(element);
            store_stress_values(source, stresses, result);
            for (std::size_t q = 0; q < 4; ++q) {
                const std::size_t history_offset = 16 + 10 * q;
                for (std::size_t component = 0; component < 4; ++component) {
                    result[history_offset + component][source] = history[q].plastic_strain[component];
                    result[history_offset + 4 + component][source] = history[q].creep_strain[component];
                }
                result[history_offset + 8][source] = history[q].equivalent_plastic_strain;
                result[history_offset + 9][source] = history[q].equivalent_creep_strain;
            }
        }
    }
    return result;
}
void fill_cartesian_nodal(const UnstructuredHex8Mesh& mesh, const cartesian::SpatialAssembly& spatial,
    const std::vector<double>& state, std::vector<std::vector<double>>& values) {
    const double missing = std::numeric_limits<double>::quiet_NaN();
    values.assign(4, std::vector<double>(mesh.nodes().size(), missing));
    std::vector<bool> present(mesh.nodes().size(), false);
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        fill_region_nodal_values(spatial.region_mesh(region).source_node_ids(), spatial.region_node_offset(region),
            spatial.dof_map(), state, present, values);
}
void store_cartesian_stress_values(std::size_t source, const std::array<SymmetricTensor3Values, 8>& stresses,
    std::vector<std::vector<double>>& values) {
    for (std::size_t q = 0; q < 8; ++q) {
        const std::size_t offset = 6 * q;
        values[offset][source] = stresses[q].xx;
        values[offset + 1][source] = stresses[q].yy;
        values[offset + 2][source] = stresses[q].zz;
        values[offset + 3][source] = stresses[q].xy;
        values[offset + 4][source] = stresses[q].yz;
        values[offset + 5][source] = stresses[q].xz;
    }
}
std::vector<std::vector<double>> cartesian_elements(
    const UnstructuredHex8Mesh& mesh, const cartesian::SpatialAssembly& spatial, const std::vector<double>& state) {
    const double missing = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> result(48, std::vector<double>(mesh.elements().size(), missing));
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const Hex8RegionMesh& region_mesh = spatial.region_mesh(region);
        for (std::size_t element = 0; element < region_mesh.elements().size(); ++element) {
            const auto stresses = spatial.stress(region, element, state);
            store_cartesian_stress_values(region_mesh.source_element_ids()[element], stresses, result);
        }
    }
    return result;
}
} // namespace
std::string next_results_segment_path(const std::string& configured_path) {
    if (configured_path.empty()) throw std::invalid_argument("Results segment path requires a configured path");
    const std::filesystem::path configured(configured_path);
    const std::filesystem::path directory = configured.parent_path();
    const std::string stem = configured.stem().string();
    const std::string extension = configured.extension().string();
    for (std::size_t segment = 1;; ++segment) {
        const std::filesystem::path candidate = directory / (stem + ".part" + std::to_string(segment) + extension);
        if (!std::filesystem::exists(candidate)) return candidate.string();
    }
}
EngineeringHistoryWriter::EngineeringHistoryWriter(std::string path, const TransientProblem& problem)
    : _path(std::move(path)), _problem_signature(transient_problem_signature(problem)) {
    if (_path.empty()) throw std::invalid_argument("Engineering history path must not be empty");
    _stream.open(_path, std::ios::out | std::ios::trunc);
    if (!_stream) throw std::runtime_error("Could not open engineering history file '" + _path + "'");
    _stream.exceptions(std::ios::badbit | std::ios::failbit);
    _stream << "time,time_step,next_time_step,load_factor,nonlinear_iterations";
    for (const TransientConservationField& field : transient_conservation_fields) _stream << ',' << field.name;
    for (const RegionDefinition& region : problem.definition().regions) {
        const std::string prefix = ",region_" + region.name;
        _stream << prefix << "_maximum_temperature" << prefix << "_maximum_equivalent_plastic_strain" << prefix
                << "_maximum_equivalent_creep_strain";
    }
    for (const ContactDefinition& contact : problem.definition().contacts) {
        const std::string prefix = ",contact_" + contact.name;
        _stream << prefix << "_minimum_gap" << prefix << "_maximum_pressure" << prefix << "_total_heat_rate" << prefix
                << "_total_force" << prefix << "_total_tangential_force";
    }
    _stream << '\n' << std::scientific << std::setprecision(12);
}
void EngineeringHistoryWriter::append(
    const TransientProblem& problem, double time_step, double next_time_step, int nonlinear_iterations) {
    if (problem.time_step_active())
        throw std::logic_error("Engineering history cannot be written during an active time step");
    if (transient_problem_signature(problem) != _problem_signature)
        throw std::invalid_argument("Engineering history problem does not match writer model");
    _stream << problem.committed_time() << ',' << time_step << ',' << next_time_step << ','
            << problem.committed_load_factor() << ',' << nonlinear_iterations;
    const TransientConservationSummary& conservation = problem.last_conservation_summary();
    for (const TransientConservationField& field : transient_conservation_fields)
        _stream << ',' << conservation.*field.member;
    const std::vector<double>& state = problem.committed_solution();
    for (std::size_t region = 0; region < problem.definition().regions.size(); ++region) {
        const RegionStateSummary summary = problem.summarize_region(region);
        _stream << ',' << summary.maximum_temperature << ',' << summary.maximum_equivalent_plastic_strain << ','
                << summary.maximum_equivalent_creep_strain;
    }
    if (!problem.is_cartesian_3d()) {
        const rz::TransientBackendView backend = rz::BackendAccess::transient(problem);
        for (std::size_t contact = 0; contact < problem.definition().contacts.size(); ++contact) {
            const InterfaceSummary summary = backend.spatial.summarize_interface(contact, state);
            _stream << ',' << summary.minimum_gap << ',' << summary.maximum_contact_pressure << ','
                    << summary.total_heat_rate << ',' << summary.total_contact_force << ','
                    << summary.total_tangential_force;
        }
    }
    _stream << '\n';
    _stream.flush();
}
void write_steady_results(const std::string& path, const UnstructuredQuad4Mesh& mesh, const SteadyProblem& problem,
    const std::vector<double>& state) {
    if (path.empty()) throw std::invalid_argument("Exodus result path must not be empty");
    const SpatialDefinition& definition = rz::BackendAccess::steady(problem).spatial.definition();
    const std::vector<std::string> nodal = nodal_variable_names(definition.contacts);
    const std::vector<std::string> element = stress_variable_names();
    const std::vector<std::string> global = global_variable_names(definition.contacts);
    define_variables(path, mesh, nodal, element, global);
    std::vector<std::vector<double>> nodal_values;
    fill_rz_nodal(mesh, rz::BackendAccess::steady(problem).spatial, state, nodal_values);
    write_result_step(path, results_mesh_view(mesh), 1, 1.0, nodal_values, steady_elements(mesh, problem, state),
        rz_globals(rz::BackendAccess::steady(problem).spatial, state, problem.load_factor()));
}
void write_steady_results(const std::string& path, const UnstructuredHex8Mesh& mesh, const SteadyProblem& problem,
    const std::vector<double>& state) {
    if (path.empty()) throw std::invalid_argument("Exodus result path must not be empty");
    define_variables(path, mesh, cartesian_nodal_variable_names(), cartesian_stress_variable_names(), {"load_factor"});
    std::vector<std::vector<double>> nodal_values;
    const cartesian::SpatialAssembly& spatial = cartesian::BackendAccess::spatial(problem);
    fill_cartesian_nodal(mesh, spatial, state, nodal_values);
    write_result_step(path, results_mesh_view(mesh), 1, 1.0, nodal_values, cartesian_elements(mesh, spatial, state),
        {problem.load_factor()});
}
ExodusTransientResultsWriter::ExodusTransientResultsWriter(
    std::string path, UnstructuredQuad4Mesh mesh, const TransientProblem& problem)
    : _path(std::move(path)), _rz_mesh(std::make_unique<UnstructuredQuad4Mesh>(std::move(mesh))),
      _problem_signature(transient_problem_signature(problem)), _step_count(0) {
    if (_path.empty()) throw std::invalid_argument("Exodus result path must not be empty");
    const std::vector<ContactDefinition>& contacts = problem.definition().contacts;
    define_variables(_path, *_rz_mesh, nodal_variable_names(contacts), transient_element_variable_names(),
        global_variable_names(contacts));
}
ExodusTransientResultsWriter::ExodusTransientResultsWriter(
    std::string path, UnstructuredHex8Mesh mesh, const TransientProblem& problem)
    : _path(std::move(path)), _hex_mesh(std::make_unique<UnstructuredHex8Mesh>(std::move(mesh))),
      _problem_signature(transient_problem_signature(problem)), _step_count(0) {
    if (_path.empty()) throw std::invalid_argument("Exodus result path must not be empty");
    define_variables(
        _path, *_hex_mesh, cartesian_nodal_variable_names(), cartesian_stress_variable_names(), {"load_factor"});
}
void ExodusTransientResultsWriter::append(const TransientProblem& problem) {
    if (problem.time_step_active())
        throw std::logic_error("Exodus results cannot be written during an active time step");
    if (transient_problem_signature(problem) != _problem_signature)
        throw std::invalid_argument("Exodus result problem does not match writer model");
    ++_step_count;
    std::vector<std::vector<double>> nodal_values;
    if (_hex_mesh) {
        const cartesian::SpatialAssembly& spatial = cartesian::BackendAccess::spatial(problem);
        fill_cartesian_nodal(*_hex_mesh, spatial, problem.committed_solution(), nodal_values);
        write_result_step(_path, results_mesh_view(*_hex_mesh), _step_count, problem.committed_time(), nodal_values,
            cartesian_elements(*_hex_mesh, spatial, problem.committed_solution()), {problem.committed_load_factor()});
    } else {
        const rz::SpatialAssembly& spatial = rz::BackendAccess::transient(problem).spatial;
        fill_rz_nodal(*_rz_mesh, spatial, problem.committed_solution(), nodal_values);
        write_result_step(_path, results_mesh_view(*_rz_mesh), _step_count, problem.committed_time(), nodal_values,
            transient_elements(*_rz_mesh, problem),
            rz_globals(spatial, problem.committed_solution(), problem.committed_load_factor()));
    }
}
std::size_t ExodusTransientResultsWriter::step_count() const noexcept { return _step_count; }
namespace {
constexpr std::array<unsigned char, 16> checkpoint_magic = {
    'F', 'U', 'E', 'L', 'S', 'I', 'M', '_', 'C', 'H', 'E', 'C', 'K', 'P', 'T', '\0'};
constexpr std::uint32_t checkpoint_version = 8U;
constexpr std::uint32_t rz_geometry_tag = 1U;
constexpr std::uint32_t cartesian_geometry_tag = 2U;
constexpr std::uint32_t endian_marker = 0x01020304U;
constexpr std::uint64_t maximum_checkpoint_bytes = 16ULL * 1024ULL * 1024ULL * 1024ULL;
std::uint64_t checksum(const std::vector<unsigned char>& bytes) {
    std::uint64_t value = fnv_offset;
    hash_bytes(value, bytes.data(), bytes.size());
    return value;
}
class BinaryBuffer final {
  public:
    void append_u32(std::uint32_t value) {
        for (std::size_t byte = 0; byte < 4; ++byte) _bytes.push_back(static_cast<unsigned char>(value >> (8U * byte)));
    }
    void append_u64(std::uint64_t value) {
        for (std::size_t byte = 0; byte < 8; ++byte) _bytes.push_back(static_cast<unsigned char>(value >> (8U * byte)));
    }
    void append_double(double value) { append_u64(encode_double_bits(value)); }
    void append_bytes(const unsigned char* data, std::size_t size) { _bytes.insert(_bytes.end(), data, data + size); }
    const std::vector<unsigned char>& bytes() const noexcept { return _bytes; }

  private:
    std::vector<unsigned char> _bytes;
};
class BinaryCursor final {
  public:
    explicit BinaryCursor(const std::vector<unsigned char>& bytes) : _bytes(bytes), _position(0) {}
    std::uint32_t read_u32() {
        require(4);
        std::uint32_t result = 0;
        for (std::size_t byte = 0; byte < 4; ++byte)
            result |= static_cast<std::uint32_t>(_bytes[_position++]) << (8U * byte);
        return result;
    }
    std::uint64_t read_u64() {
        require(8);
        std::uint64_t result = 0;
        for (std::size_t byte = 0; byte < 8; ++byte)
            result |= static_cast<std::uint64_t>(_bytes[_position++]) << (8U * byte);
        return result;
    }
    double read_double() { return decode_double_bits(read_u64()); }
    void read_bytes(unsigned char* destination, std::size_t size) {
        require(size);
        std::memcpy(destination, _bytes.data() + _position, size);
        _position += size;
    }
    bool at_end() const noexcept { return _position == _bytes.size(); }

  private:
    void require(std::size_t size) const {
        if (size > _bytes.size() - _position) throw std::runtime_error("Checkpoint payload is truncated");
    }
    const std::vector<unsigned char>& _bytes;
    std::size_t _position;
};
void append_strains(BinaryBuffer& payload, const std::array<double, 4>& values) {
    for (const double value : values) payload.append_double(value);
}
void read_strains(BinaryCursor& payload, std::array<double, 4>& values) {
    for (double& value : values) value = payload.read_double();
}
void append_material_point(
    BinaryBuffer& payload, const MaterialPointState& state, const AxisymmetricStressValues& stress) {
    append_strains(payload, state.elastic_strain);
    append_strains(payload, state.plastic_strain);
    append_strains(payload, state.creep_strain);
    payload.append_double(state.equivalent_plastic_strain);
    payload.append_double(state.equivalent_creep_strain);
    payload.append_double(stress.rr);
    payload.append_double(stress.zz);
    payload.append_double(stress.hoop);
    payload.append_double(stress.rz);
}
void read_material_point(BinaryCursor& payload, MaterialPointState& state, AxisymmetricStressValues& stress) {
    read_strains(payload, state.elastic_strain);
    read_strains(payload, state.plastic_strain);
    read_strains(payload, state.creep_strain);
    state.equivalent_plastic_strain = payload.read_double();
    state.equivalent_creep_strain = payload.read_double();
    stress.rr = payload.read_double();
    stress.zz = payload.read_double();
    stress.hoop = payload.read_double();
    stress.rz = payload.read_double();
}
void append_conservation(BinaryBuffer& payload, const TransientConservationSummary& summary) {
    for (const TransientConservationField& field : transient_conservation_fields)
        payload.append_double(summary.*field.member);
}
TransientConservationSummary read_conservation(BinaryCursor& payload) {
    TransientConservationSummary result;
    for (const TransientConservationField& field : transient_conservation_fields) {
        result.*field.member = payload.read_double();
        if (!std::isfinite(result.*field.member))
            throw std::runtime_error("Checkpoint conservation summary is invalid");
    }
    return result;
}
BinaryBuffer state_payload(const TransientProblem& problem, double next_time_step) {
    BinaryBuffer payload;
    payload.append_u64(transient_problem_signature(problem));
    const bool cartesian = problem.is_cartesian_3d();
    payload.append_u32(cartesian ? cartesian_geometry_tag : rz_geometry_tag);
    const TransientCommittedState state =
        cartesian ? cartesian::BackendAccess::committed_state(problem) : rz::BackendAccess::committed_state(problem);
    payload.append_double(state.time);
    payload.append_double(state.load_factor);
    payload.append_double(next_time_step);
    append_conservation(payload, state.conservation);
    payload.append_u64(static_cast<std::uint64_t>(state.solution.size()));
    for (const double value : state.solution) payload.append_double(value);
    if (cartesian) return payload;
    payload.append_u64(static_cast<std::uint64_t>(state.contact_histories.size()));
    for (const auto& contact : state.contact_histories) {
        payload.append_u64(static_cast<std::uint64_t>(contact.size()));
        for (const ContactPointHistory& history : contact) {
            payload.append_double(history.elastic_tangential_slip);
            payload.append_u32(history.sliding ? 1U : 0U);
            payload.append_double(history.normal_multiplier);
        }
    }
    payload.append_u64(static_cast<std::uint64_t>(state.material_histories.size()));
    for (std::size_t region = 0; region < state.material_histories.size(); ++region) {
        payload.append_u64(static_cast<std::uint64_t>(state.material_histories[region].size()));
        for (std::size_t element = 0; element < state.material_histories[region].size(); ++element)
            for (std::size_t q = 0; q < 4; ++q)
                append_material_point(
                    payload, state.material_histories[region][element][q], state.material_stresses[region][element][q]);
    }
    return payload;
}
void append_checkpoint_header(BinaryBuffer& file, const std::vector<unsigned char>& payload) {
    file.append_bytes(checkpoint_magic.data(), checkpoint_magic.size());
    file.append_u32(checkpoint_version);
    file.append_u32(endian_marker);
    file.append_u64(static_cast<std::uint64_t>(payload.size()));
    file.append_u64(checksum(payload));
}
std::vector<unsigned char> read_file(const std::string& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) throw std::runtime_error("Could not open checkpoint '" + path + "'");
    const std::streamoff length = input.tellg();
    if (length < 0 || static_cast<std::uint64_t>(length) > maximum_checkpoint_bytes)
        throw std::runtime_error("Checkpoint size is invalid: " + path);
    input.seekg(0, std::ios::beg);
    std::vector<unsigned char> result(static_cast<std::size_t>(length));
    if (!result.empty()) input.read(reinterpret_cast<char*>(result.data()), length);
    if (!input) throw std::runtime_error("Could not read checkpoint '" + path + "'");
    return result;
}
void write_atomic(const std::string& path, const std::vector<unsigned char>& bytes) {
    const std::string temporary = path + ".tmp";
    try {
        std::ofstream output(temporary, std::ios::binary | std::ios::out | std::ios::trunc);
        if (!output) throw std::runtime_error("Could not create checkpoint '" + temporary + "'");
        output.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
        output.close();
        if (!output) throw std::runtime_error("Could not finish checkpoint '" + temporary + "'");
        if (std::rename(temporary.c_str(), path.c_str()) != 0)
            throw std::runtime_error("Could not install checkpoint '" + path + "': " + std::strerror(errno));
    } catch (...) {
        (void)std::remove(temporary.c_str());
        throw;
    }
}
} // namespace
void write_transient_checkpoint(const std::string& path, const TransientProblem& problem, double next_time_step) {
    if (path.empty()) throw std::invalid_argument("Checkpoint path must not be empty");
    if (problem.time_step_active()) throw std::logic_error("Checkpoint cannot be written during an active time step");
    if (!std::isfinite(next_time_step) || !(next_time_step > 0.0))
        throw std::invalid_argument("Checkpoint next time step must be finite and positive");
    const BinaryBuffer payload = state_payload(problem, next_time_step);
    BinaryBuffer file;
    append_checkpoint_header(file, payload.bytes());
    file.append_bytes(payload.bytes().data(), payload.bytes().size());
    write_atomic(path, file.bytes());
}
double restore_transient_checkpoint(const std::string& path, TransientProblem& problem) {
    const std::vector<unsigned char> file = read_file(path);
    BinaryCursor header(file);
    std::array<unsigned char, checkpoint_magic.size()> magic{};
    header.read_bytes(magic.data(), magic.size());
    if (magic != checkpoint_magic) throw std::runtime_error("Checkpoint magic does not match fuelsim");
    if (header.read_u32() != checkpoint_version) throw std::runtime_error("Checkpoint version is not supported");
    if (header.read_u32() != endian_marker) throw std::runtime_error("Checkpoint endian marker is invalid");
    const std::uint64_t payload_size = header.read_u64();
    const std::uint64_t expected_checksum = header.read_u64();
    constexpr std::size_t header_size = 16U + 4U + 4U + 8U + 8U;
    if (payload_size != file.size() - header_size) throw std::runtime_error("Checkpoint payload length is invalid");
    std::vector<unsigned char> payload_bytes(file.begin() + static_cast<std::ptrdiff_t>(header_size), file.end());
    if (checksum(payload_bytes) != expected_checksum)
        throw std::runtime_error("Checkpoint payload checksum does not match");
    BinaryCursor payload(payload_bytes);
    if (payload.read_u64() != transient_problem_signature(problem))
        throw std::runtime_error("Checkpoint model signature does not match the current problem");
    const std::uint32_t geometry_tag = payload.read_u32();
    const bool cartesian = problem.is_cartesian_3d();
    if (geometry_tag != (cartesian ? cartesian_geometry_tag : rz_geometry_tag))
        throw std::runtime_error("Checkpoint geometry does not match the current problem");
    TransientCommittedState state;
    state.time = payload.read_double();
    state.load_factor = payload.read_double();
    const double next_time_step = payload.read_double();
    state.conservation = read_conservation(payload);
    if (!std::isfinite(next_time_step) || !(next_time_step > 0.0))
        throw std::runtime_error("Checkpoint next time step is invalid");
    if (payload.read_u64() != problem.dof_count())
        throw std::runtime_error("Checkpoint solution size does not match the current problem");
    state.solution.resize(problem.dof_count());
    for (double& value : state.solution) value = payload.read_double();
    if (cartesian) {
        if (!payload.at_end()) throw std::runtime_error("Checkpoint payload contains trailing data");
        cartesian::BackendAccess::restore_committed_state(problem, std::move(state));
        return next_time_step;
    }
    const std::uint64_t contact_count = payload.read_u64();
    const TransientCommittedState expected_state = rz::BackendAccess::committed_state(problem);
    const auto& expected_histories = expected_state.contact_histories;
    if (contact_count != expected_histories.size())
        throw std::runtime_error("Checkpoint contact count does not match the current problem");
    state.contact_histories.resize(expected_histories.size());
    for (std::size_t contact = 0; contact < expected_histories.size(); ++contact) {
        const std::uint64_t node_count = payload.read_u64();
        if (node_count != expected_histories[contact].size())
            throw std::runtime_error("Checkpoint contact-node count does not match the current problem");
        state.contact_histories[contact].resize(expected_histories[contact].size());
        for (ContactPointHistory& history : state.contact_histories[contact]) {
            history.elastic_tangential_slip = payload.read_double();
            const std::uint32_t sliding = payload.read_u32();
            if (sliding > 1U) throw std::runtime_error("Checkpoint friction state is invalid");
            history.sliding = sliding == 1U;
            history.normal_multiplier = payload.read_double();
            if (!std::isfinite(history.normal_multiplier) || history.normal_multiplier < 0.0)
                throw std::runtime_error("Checkpoint normal contact multiplier is invalid");
        }
    }
    const std::uint64_t region_count = payload.read_u64();
    const rz::TransientBackendView backend = rz::BackendAccess::transient(problem);
    if (region_count != backend.spatial.region_count())
        throw std::runtime_error("Checkpoint region count does not match the current problem");
    state.material_histories.resize(backend.spatial.region_count());
    state.material_stresses.resize(backend.spatial.region_count());
    for (std::size_t region = 0; region < backend.spatial.region_count(); ++region) {
        const std::size_t elements = backend.spatial.region_mesh(region).elements().size();
        if (payload.read_u64() != elements)
            throw std::runtime_error("Checkpoint element count does not match the current problem");
        state.material_histories[region].resize(elements);
        state.material_stresses[region].resize(elements);
        for (std::size_t element = 0; element < elements; ++element)
            for (std::size_t q = 0; q < 4; ++q)
                read_material_point(
                    payload, state.material_histories[region][element][q], state.material_stresses[region][element][q]);
    }
    if (!payload.at_end()) throw std::runtime_error("Checkpoint payload contains trailing data");
    rz::BackendAccess::restore_committed_state(problem, std::move(state));
    return next_time_step;
}
} // namespace fuelsim
