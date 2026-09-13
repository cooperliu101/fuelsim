#include "support/exodus_result_reader.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exodusII.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using fuelsim::test::ExodusResults;

void require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

void checked(int status, const char* message) {
    require(status >= 0, message);
}

std::size_t size(std::int64_t count) {
    require(count >= 0 && static_cast<std::uint64_t>(count) <= std::numeric_limits<std::size_t>::max(),
        "Invalid Exodus geometry count");
    return static_cast<std::size_t>(count);
}

std::size_t product(std::int64_t first, std::int64_t second) {
    const auto a = size(first), b = size(second);
    require(b == 0 || a <= std::numeric_limits<std::size_t>::max() / b, "Exodus geometry count overflow");
    return a * b;
}

class GeometryFile final {
  public:
    explicit GeometryFile(const std::string& path) {
        int cpu = static_cast<int>(sizeof(double)), disk = 0;
        float version = 0.0F;
        _id = ex_open(path.c_str(), EX_READ, &cpu, &disk, &version);
        require(_id >= 0, "Cannot open geometry in " + path);
        checked(ex_set_int64_status(_id, EX_ALL_INT64_API), "Cannot configure Exodus integer reads");
    }

    ~GeometryFile() { ex_close(_id); }

    GeometryFile(const GeometryFile&) = delete;
    GeometryFile& operator=(const GeometryFile&) = delete;

    int id() const noexcept { return _id; }

  private:
    int _id = -1;
};

struct GeometryRecords final {
    std::vector<std::int64_t> integers;
    std::vector<double> values;
    std::vector<std::string> names;
};

// The shared history reader supplies coordinates and all result fields. These
// additional records cover connectivity, axial-control attributes and sets.
GeometryRecords geometry_records(const std::string& path) {
    const GeometryFile file(path);
    const auto id = file.id();
    GeometryRecords result;
    const auto maximum_name = size(ex_inquire_int(id, EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH));
    require(maximum_name < static_cast<std::size_t>(std::numeric_limits<int>::max()),
        "Exodus name length is too large");
    checked(ex_set_max_name_length(id, static_cast<int>(maximum_name)), "Cannot configure full geometry names");
    const auto read_name = [&](ex_entity_type kind, std::int64_t entity) {
        std::vector<char> name(maximum_name + 1, '\0');
        checked(ex_get_name(id, kind, entity, name.data()), "Cannot read an Exodus geometry name");
        return std::string(name.data());
    };
    for (const auto inquiry : {EX_INQ_DIM, EX_INQ_NODES, EX_INQ_ELEM}) {
        const auto count = ex_inquire_int(id, inquiry);
        (void)size(count);
        result.integers.push_back(count);
    }
    for (const auto kind : {EX_NODE_MAP, EX_ELEM_MAP}) {
        const auto count = kind == EX_NODE_MAP ? result.integers[1] : result.integers[2];
        std::vector<std::int64_t> map(size(count));
        if (!map.empty())
            checked(ex_get_id_map(id, kind, map.data()), "Cannot read Exodus node or element numbering");
        result.integers.insert(result.integers.end(), map.begin(), map.end());
    }
    const auto block_count = ex_inquire_int(id, EX_INQ_ELEM_BLK);
    result.integers.push_back(block_count);
    std::vector<std::int64_t> blocks(size(block_count));
    if (!blocks.empty())
        checked(ex_get_ids(id, EX_ELEM_BLOCK, blocks.data()), "Cannot read Exodus block IDs");
    for (const auto block_id : blocks) {
        ex_block block{};
        block.type = EX_ELEM_BLOCK;
        block.id = block_id;
        checked(ex_get_block_param(id, &block), "Cannot read Exodus block topology");
        result.names.push_back(read_name(EX_ELEM_BLOCK, block_id));
        result.names.emplace_back(block.topology);
        for (const auto count : {block_id,
                 block.num_entry,
                 block.num_nodes_per_entry,
                 block.num_edges_per_entry,
                 block.num_faces_per_entry,
                 block.num_attribute})
            result.integers.push_back(count);
        std::vector<std::int64_t> nodes(product(block.num_entry, block.num_nodes_per_entry));
        std::vector<std::int64_t> edges(product(block.num_entry, block.num_edges_per_entry));
        std::vector<std::int64_t> faces(product(block.num_entry, block.num_faces_per_entry));
        checked(ex_get_conn(id,
                    EX_ELEM_BLOCK,
                    block_id,
                    nodes.data(),
                    edges.empty() ? nullptr : edges.data(),
                    faces.empty() ? nullptr : faces.data()),
            "Cannot read complete Exodus element connectivity");
        for (const auto* entries : {&nodes, &edges, &faces})
            result.integers.insert(result.integers.end(), entries->begin(), entries->end());
        if (block.num_attribute > 0) {
            std::vector<std::vector<char>> names(size(block.num_attribute), std::vector<char>(maximum_name + 1, '\0'));
            std::vector<char*> pointers;
            for (auto& name : names)
                pointers.push_back(name.data());
            checked(ex_get_attr_names(id, EX_ELEM_BLOCK, block_id, pointers.data()),
                "Cannot read block attribute names");
            for (const auto& name : names)
                result.names.emplace_back(name.data());
            std::vector<double> attributes(product(block.num_entry, block.num_attribute));
            checked(ex_get_attr(id, EX_ELEM_BLOCK, block_id, attributes.data()),
                "Cannot read axial-control attributes");
            result.values.insert(result.values.end(), attributes.begin(), attributes.end());
        }
    }
    for (const auto kind : {EX_NODE_SET, EX_SIDE_SET}) {
        const auto count = ex_inquire_int(id, kind == EX_NODE_SET ? EX_INQ_NODE_SETS : EX_INQ_SIDE_SETS);
        result.integers.push_back(count);
        std::vector<std::int64_t> sets(size(count));
        if (!sets.empty())
            checked(ex_get_ids(id, kind, sets.data()), "Cannot read Exodus set IDs");
        for (const auto set : sets) {
            std::int64_t entries = 0, factors = 0;
            checked(ex_get_set_param(id, kind, set, &entries, &factors), "Cannot read Exodus set dimensions");
            result.integers.insert(result.integers.end(), {set, entries, factors});
            result.names.push_back(read_name(kind, set));
            std::vector<std::int64_t> members(size(entries)), sides(kind == EX_SIDE_SET ? size(entries) : 0);
            if (entries > 0)
                checked(ex_get_set(id, kind, set, members.data(), sides.empty() ? nullptr : sides.data()),
                    "Cannot read Exodus set membership");
            result.integers.insert(result.integers.end(), members.begin(), members.end());
            result.integers.insert(result.integers.end(), sides.begin(), sides.end());
            std::vector<double> weights(size(factors));
            if (!weights.empty())
                checked(ex_get_set_dist_fact(id, kind, set, weights.data()), "Cannot read Exodus set weights");
            result.values.insert(result.values.end(), weights.begin(), weights.end());
        }
    }
    for (const auto value : result.values)
        require(std::isfinite(value), "Exodus geometry attributes and weights must be finite");
    return result;
}

enum class Quantity : std::size_t {
    exact,
    temperature,
    length,
    strain,
    stress,
    force,
    heat_rate,
    heat_flux,
    energy,
    area,
    count
};

struct Metric final {
    const char* name;
    const char* unit;
    double absolute_tolerance;
    std::size_t finite_count = 0, nan_count = 0;
    double maximum_absolute_error = 0.0, maximum_relative_error = 0.0;
};

constexpr double relative_tolerance = 1e-10;
using Metrics = std::array<Metric, static_cast<std::size_t>(Quantity::count)>;

bool starts(const std::string& value, const char* prefix) {
    return value.rfind(prefix, 0) == 0;
}

Quantity quantity(const std::string& name) {
    if (name == "node_role" || name == "material_point_count" || name == "load_factor" || starts(name, "reference_")
        || starts(name, "contact_sliding_"))
        return Quantity::exact;
    if (name == "temperature" || starts(name, "temperature_q"))
        return Quantity::temperature;
    if (starts(name, "displacement_") || starts(name, "contact_gap_")
        || starts(name, "contact_elastic_tangential_slip_") || starts(name, "contact_total_tangential_slip_"))
        return Quantity::length;
    if (name == "axial_strain" || starts(name, "elastic_") || starts(name, "plastic_") || starts(name, "creep_")
        || starts(name, "equiv_") || starts(name, "conservation_relative_"))
        return Quantity::strain;
    if (starts(name, "stress_") || starts(name, "contact_pressure_") || starts(name, "contact_tangential_traction_"))
        return Quantity::stress;
    if (name == "axial_force" || starts(name, "reaction_force_") || starts(name, "contact_force_")
        || starts(name, "contact_tangential_force_") || name == "conservation_unconstrained_mechanical_residual_l2")
        return Quantity::force;
    if (starts(name, "contact_area_"))
        return Quantity::area;
    if (starts(name, "contact_heat_flux_"))
        return Quantity::heat_flux;
    if (name == "reaction_heat_flux" || starts(name, "contact_heat_rate_") || name == "conservation_generated_heat_rate"
        || name == "conservation_stored_heat_rate" || name == "conservation_convection_heat_rate"
        || name == "conservation_surface_heat_input_rate" || name == "conservation_interface_heat_imbalance"
        || name == "conservation_dirichlet_heat_input_rate" || name == "conservation_global_thermal_balance"
        || name == "conservation_unconstrained_thermal_residual_l2")
        return Quantity::heat_rate;
    if (name == "conservation_internal_mechanical_work_increment"
        || name == "conservation_pressure_traction_work_increment" || name == "conservation_body_force_work_increment"
        || name == "conservation_dirichlet_reaction_work_increment" || name == "conservation_contact_work_increment"
        || name == "conservation_mechanical_work_balance" || name == "conservation_elastic_energy_change"
        || name == "conservation_plastic_dissipation_increment" || name == "conservation_creep_dissipation_increment"
        || name == "conservation_friction_dissipation_increment"
        || name == "conservation_trapezoidal_pressure_traction_work_increment"
        || name == "conservation_trapezoidal_dirichlet_reaction_work_increment"
        || name == "conservation_mechanical_hourglass_energy"
        || name == "conservation_mechanical_hourglass_energy_change")
        return Quantity::energy;
    throw std::runtime_error("No physical tolerance is defined for result field '" + name + "'");
}

void compare_value(double expected,
    double actual,
    const std::string& field,
    double time,
    std::size_t entry,
    Metrics& metrics) {
    const auto kind = quantity(field);
    auto& metric = metrics[static_cast<std::size_t>(kind)];
    const auto location = field + " at time " + std::to_string(time) + ", entry " + std::to_string(entry);
    require(std::isnan(expected) == std::isnan(actual), "NaN mask mismatch in " + location);
    if (std::isnan(expected)) {
        ++metric.nan_count;
        return;
    }
    require(std::isfinite(expected) && std::isfinite(actual), "Nonfinite result in " + location);
    const double difference = std::abs(actual - expected);
    // Zero references use the stated SI absolute tolerance. No reference
    // denominator floor or field-dependent data normalization is introduced.
    const double bound =
        kind == Quantity::exact ? 0.0 : metric.absolute_tolerance + relative_tolerance * std::abs(expected);
    if (difference > bound) {
        std::cerr << std::setprecision(17) << location << ": reference=" << expected << " other=" << actual
                  << " absolute_error=" << difference << " allowed=" << bound << ' ' << metric.unit << '\n';
        throw std::runtime_error("Production radial results differ beyond their physical tolerance");
    }
    ++metric.finite_count;
    metric.maximum_absolute_error = std::max(metric.maximum_absolute_error, difference);
    if (expected != 0.0)
        metric.maximum_relative_error = std::max(metric.maximum_relative_error, difference / std::abs(expected));
}

void unique_names(const std::vector<std::string>& names) {
    require(!names.empty(), "Production result variable lists must not be empty");
    const std::set<std::string> unique(names.begin(), names.end());
    require(unique.size() == names.size() && unique.count("") == 0,
        "Production result variable names must be unique and nonempty");
}

void compare_frame(const ExodusResults& serial, const ExodusResults& other, Metrics& metrics) {
    for (const auto* nodes : {&serial.nodes, &other.nodes})
        for (const auto& node : *nodes)
            for (const auto coordinate : node)
                require(std::isfinite(coordinate), "Production geometry coordinates must be finite");
    require(serial.nodes == other.nodes && serial.block_names == other.block_names
                && serial.block_element_counts == other.block_element_counts
                && serial.side_set_names == other.side_set_names && serial.side_set_sizes == other.side_set_sizes
                && serial.side_set_face_nodes == other.side_set_face_nodes,
        "Production result geometry or named block/side-set layout differs");
    require(serial.nodal_variable_names == other.nodal_variable_names
                && serial.element_variable_names == other.element_variable_names
                && serial.global_variable_names == other.global_variable_names,
        "Production result variable names or their order differs");
    std::size_t elements = 0;
    for (const auto count : serial.block_element_counts)
        elements += count;
    for (std::size_t category = 0; category < 2; ++category) {
        const auto& names = category == 0 ? serial.nodal_variable_names : serial.element_variable_names;
        const auto& expected = category == 0 ? serial.nodal_variables : serial.element_variables;
        const auto& actual = category == 0 ? other.nodal_variables : other.element_variables;
        const auto count = category == 0 ? serial.nodes.size() : elements;
        unique_names(names);
        require(expected.size() == names.size() && actual.size() == names.size(),
            "Production variable array count differs");
        for (std::size_t field = 0; field < names.size(); ++field) {
            require(expected[field].size() == count && actual[field].size() == count,
                "Production node/element value count differs in " + names[field]);
            for (std::size_t entry = 0; entry < count; ++entry)
                compare_value(expected[field][entry], actual[field][entry], names[field], serial.time, entry, metrics);
        }
    }
    unique_names(serial.global_variable_names);
    require(serial.global_variables.size() == serial.global_variable_names.size()
                && other.global_variables.size() == serial.global_variable_names.size(),
        "Production global value count differs");
    for (std::size_t field = 0; field < serial.global_variable_names.size(); ++field)
        compare_value(serial.global_variables[field],
            other.global_variables[field],
            serial.global_variable_names[field],
            serial.time,
            field,
            metrics);
}

void run(const std::string& serial_path, const std::string& other_path) {
    const auto serial_geometry = geometry_records(serial_path), other_geometry = geometry_records(other_path);
    require(serial_geometry.integers == other_geometry.integers && serial_geometry.values == other_geometry.values
                && serial_geometry.names == other_geometry.names,
        "Exodus topology, connectivity, axial controls, numbering or named sets differ");
    const auto serial = fuelsim::test::read_exodus_history(serial_path);
    const auto other = fuelsim::test::read_exodus_history(other_path);
    require(!serial.empty() && !other.empty(), "Both production outputs must contain frames");
    for (const auto* history : {&serial, &other})
        for (std::size_t frame = 0; frame < history->size(); ++frame) {
            require(std::isfinite((*history)[frame].time), "Production output time must be finite");
            if (frame > 0)
                require((*history)[frame].time > (*history)[frame - 1].time,
                    "Production output times must increase strictly");
        }
    constexpr double time_tolerance = 1e-12; // Seconds, independent of field-value tolerances.
    std::size_t first = serial.size();
    for (std::size_t frame = 0; frame < serial.size(); ++frame)
        if (std::abs(serial[frame].time - other.front().time) <= time_tolerance) {
            require(first == serial.size(), "The restart beginning matches multiple serial frames");
            first = frame;
        }
    require(first < serial.size() && other.size() == serial.size() - first,
        "Other output must contain a complete serial time suffix without omitted frames");
    Metrics metrics{{{"exact", "exact", 0.0},
        {"temperature", "K", 1e-10},
        {"length", "m", 1e-12},
        {"strain_or_relative_balance", "1", 1e-12},
        {"stress", "Pa", 1e-4},
        {"force", "N", 1e-8},
        {"heat_rate", "W", 1e-10},
        {"heat_flux", "W/m2", 1e-6},
        {"energy", "J", 1e-10},
        {"area", "m2", 1e-16}}};
    for (std::size_t frame = 0; frame < other.size(); ++frame) {
        require(std::abs(serial[first + frame].time - other[frame].time) <= time_tolerance,
            "Production output time differs from the corresponding serial frame");
        compare_frame(serial[first + frame], other[frame], metrics);
    }
    std::cout << std::scientific << std::setprecision(12) << "compared_frames=" << other.size()
              << " serial_start_frame=" << first << " compared_nodes=" << serial.front().nodes.size()
              << " relative_tolerance=" << relative_tolerance << '\n';
    for (const auto& metric : metrics)
        std::cout << metric.name << " finite_values=" << metric.finite_count
                  << " matched_nan_values=" << metric.nan_count
                  << " maximum_absolute_error=" << metric.maximum_absolute_error
                  << " maximum_nonzero_relative_error=" << metric.maximum_relative_error
                  << " absolute_tolerance=" << metric.absolute_tolerance << ' ' << metric.unit << '\n';
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 3)
        return 2;
    try {
        run(argv[1], argv[2]);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
