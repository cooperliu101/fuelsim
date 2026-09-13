#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exodusII.h>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using fuelsim::test::ExodusResults;
constexpr double relative_tolerance = 1e-9;
constexpr double zero_absolute_tolerance = 1e-10; // W, only for an exactly zero reference.

void require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

enum class Topology { hex8, quad4, quad8 };

struct BlockMaterial final {
    std::string name;
    std::size_t element_count;
    double initial_temperature;
    double initial_density;
    double specific_heat;
    double specific_heat_temperature_coefficient;
};

struct CaseSettings final {
    Topology topology;
    std::size_t material_point_count;
    std::size_t quadrature_order;
    std::size_t node_count;
    std::size_t increments;
    double time_step;
    std::vector<BlockMaterial> blocks;
};

// These constants are the material definitions in the corresponding complete,
// tracked verification/fuelsim/transient_*.fsi inputs. In particular, initial
// density uses the region's initial temperature, not a committed/output state.
CaseSettings settings(const std::string& name) {
    if (name == "b534")
        return {Topology::hex8, 1, 2, 12, 1, 1.0, {{"solid", 2, 300.0, 2000.0, 3000.0, 0.0}}};
    if (name == "b539")
        return {Topology::hex8,
            1,
            2,
            24,
            20,
            0.02,
            {{"primary", 2, 300.0, 100.0, 1.0, 0.001}, {"secondary", 2, 300.0, 100.0, 1.0, 0.001}}};
    if (name == "b544")
        return {Topology::hex8, 1, 2, 30, 10, 100000.0, {{"solid", 8, 300.0, 2000.0, 500.0, 1.25}}};
    if (name == "b61")
        return {Topology::hex8,
            1,
            2,
            1785,
            5,
            2.0,
            {{"clad", 1120, 600.0, 6500.0, 330.0, 0.0}, {"meat", 80, 600.0, 10970.0, 300.0, 0.0}}};
    const std::array<std::string, 4> models = {"cax4t", "cax4rt", "cax8t", "cax8rt"};
    for (std::size_t model = 0; model < models.size(); ++model) {
        const bool finite = name == "b13_finite_" + models[model];
        const bool integrated = name == "b13_integrated_" + models[model];
        if (!finite && !integrated)
            continue;
        const bool quadratic = model >= 2;
        const std::size_t points = model == 0 ? 4 : model == 1 ? 1 : model == 2 ? 9 : 4;
        return {quadratic ? Topology::quad8 : Topology::quad4,
            points,
            model == 2 ? 3U : 2U,
            finite ? (quadratic ? 138U : 53U) : (quadratic ? 402U : 150U),
            finite ? 20U : 96U,
            finite ? 1.0 : 0.0625,
            {{"fuel", finite ? 24U : 48U, 600.0, 10970.0, 300.0, 0.0},
                {"clad", finite ? 10U : 56U, 600.0, 6500.0, 330.0, 0.0}}};
    }
    throw std::invalid_argument("Unknown initial-mass verification case: " + name);
}

class ExodusFile final {
  public:
    explicit ExodusFile(const std::string& path) {
        int cpu_word_size = static_cast<int>(sizeof(double));
        int io_word_size = 0;
        float version = 0.0F;
        _id = ex_open(path.c_str(), EX_READ, &cpu_word_size, &io_word_size, &version);
        require(_id >= 0, "Could not open Exodus connectivity: " + path);
        if (ex_set_int64_status(_id, EX_ALL_INT64_API) < 0) {
            ex_close(_id);
            throw std::runtime_error("Could not select 64-bit Exodus connectivity");
        }
    }

    ~ExodusFile() { ex_close(_id); }

    ExodusFile(const ExodusFile&) = delete;
    ExodusFile& operator=(const ExodusFile&) = delete;

    int id() const { return _id; }

  private:
    int _id = -1;
};

struct Shape final {
    std::array<double, 8> value{};
    std::array<std::array<double, 3>, 8> derivative{};
};

constexpr std::array<std::array<double, 3>, 8> signs = {
    {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}}};

Shape quad4_shape(double xi, double eta) {
    Shape result;
    for (std::size_t n = 0; n < 4; ++n) {
        const double a = signs[n][0], b = signs[n][1];
        result.value[n] = 0.25 * (1 + a * xi) * (1 + b * eta);
        result.derivative[n][0] = 0.25 * a * (1 + b * eta);
        result.derivative[n][1] = 0.25 * b * (1 + a * xi);
    }
    return result;
}

Shape geometry_shape(Topology topology, double xi, double eta, double zeta) {
    if (topology == Topology::quad4)
        return quad4_shape(xi, eta);
    Shape result;
    if (topology == Topology::hex8) {
        for (std::size_t n = 0; n < 8; ++n) {
            const double a = signs[n][0], b = signs[n][1], c = signs[n][2];
            result.value[n] = 0.125 * (1 + a * xi) * (1 + b * eta) * (1 + c * zeta);
            result.derivative[n][0] = 0.125 * a * (1 + b * eta) * (1 + c * zeta);
            result.derivative[n][1] = 0.125 * b * (1 + a * xi) * (1 + c * zeta);
            result.derivative[n][2] = 0.125 * c * (1 + a * xi) * (1 + b * eta);
        }
        return result;
    }
    for (std::size_t n = 0; n < 4; ++n) {
        const double a = signs[n][0], b = signs[n][1];
        result.value[n] = 0.25 * (1 + a * xi) * (1 + b * eta) * (a * xi + b * eta - 1);
        result.derivative[n][0] = 0.25 * a * (1 + b * eta) * (2 * a * xi + b * eta);
        result.derivative[n][1] = 0.25 * b * (1 + a * xi) * (a * xi + 2 * b * eta);
    }
    result.value[4] = 0.5 * (1 - xi * xi) * (1 - eta);
    result.value[5] = 0.5 * (1 + xi) * (1 - eta * eta);
    result.value[6] = 0.5 * (1 - xi * xi) * (1 + eta);
    result.value[7] = 0.5 * (1 - xi) * (1 - eta * eta);
    result.derivative[4] = {-xi * (1 - eta), -0.5 * (1 - xi * xi), 0};
    result.derivative[5] = {0.5 * (1 - eta * eta), -(1 + xi) * eta, 0};
    result.derivative[6] = {-xi * (1 + eta), 0.5 * (1 - xi * xi), 0};
    result.derivative[7] = {-0.5 * (1 - eta * eta), -(1 - xi) * eta, 0};
    return result;
}

struct CapacityPoint final {
    double reference_measure = 0.0;
    std::array<double, 8> temperature_shape{};
};

struct Element final {
    std::array<std::size_t, 8> nodes{};
    std::size_t temperature_nodes = 0;
    BlockMaterial material;
    std::vector<CapacityPoint> points;
    std::array<double, 8> nodal_reference_volume{};
};

void integrate_reference_geometry(Element& element,
    const CaseSettings& config,
    const std::vector<std::array<double, 3>>& coordinates) {
    const double g2 = 1.0 / std::sqrt(3.0), g3 = std::sqrt(3.0 / 5.0);
    const std::vector<double> locations =
        config.quadrature_order == 3 ? std::vector<double>{-g3, 0, g3} : std::vector<double>{-g2, g2};
    const std::vector<double> weights =
        config.quadrature_order == 3 ? std::vector<double>{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0} : std::vector<double>{1, 1};
    const bool cartesian = config.topology == Topology::hex8;
    const std::size_t geometry_nodes = config.topology == Topology::quad4 ? 4 : 8;
    const std::size_t dimension = cartesian ? 3 : 2;
    for (std::size_t i = 0; i < locations.size(); ++i)
        for (std::size_t j = 0; j < locations.size(); ++j)
            for (std::size_t k = 0; k < (cartesian ? locations.size() : 1); ++k) {
                const Shape shape =
                    geometry_shape(config.topology, locations[i], locations[j], cartesian ? locations[k] : 0.0);
                std::array<std::array<double, 3>, 3> jacobian{};
                double radius = 0.0;
                for (std::size_t n = 0; n < geometry_nodes; ++n) {
                    const auto& position = coordinates.at(element.nodes[n]);
                    radius += shape.value[n] * position[0];
                    for (std::size_t a = 0; a < dimension; ++a)
                        for (std::size_t b = 0; b < dimension; ++b)
                            jacobian[a][b] += position[a] * shape.derivative[n][b];
                }
                const auto& a = jacobian;
                const double determinant = cartesian ? a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                                                           - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                                                           + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
                                                     : a[0][0] * a[1][1] - a[0][1] * a[1][0];
                require(std::isfinite(determinant) && determinant > 0.0, "Invalid independent reference Jacobian");
                require(cartesian || (std::isfinite(radius) && radius > 0.0), "Invalid independent reference radius");
                CapacityPoint point;
                point.reference_measure =
                    determinant * weights[i] * weights[j] * (cartesian ? weights[k] : 2.0 * std::acos(-1.0) * radius);
                point.temperature_shape = cartesian ? shape.value : quad4_shape(locations[i], locations[j]).value;
                require(std::isfinite(point.reference_measure) && point.reference_measure > 0,
                    "Invalid reference capacity measure");
                element.points.push_back(point);
                for (std::size_t n = 0; n < element.temperature_nodes; ++n)
                    element.nodal_reference_volume[n] += point.reference_measure * point.temperature_shape[n];
            }
}

std::vector<Element> read_elements(const std::string& path, const CaseSettings& config, const ExodusResults& frame) {
    const ExodusFile file(path);
    const auto block_count = ex_inquire_int(file.id(), EX_INQ_ELEM_BLK);
    require(block_count == static_cast<std::int64_t>(config.blocks.size()), "Incorrect number of material blocks");
    std::vector<std::int64_t> ids(config.blocks.size());
    require(ex_get_ids(file.id(), EX_ELEM_BLOCK, ids.data()) >= 0, "Could not read material block IDs");
    const auto maximum_name_length = ex_inquire_int(file.id(), EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH);
    require(maximum_name_length > 0 && maximum_name_length < std::numeric_limits<int>::max(),
        "Invalid Exodus maximum name length");
    require(ex_set_max_name_length(file.id(), static_cast<int>(maximum_name_length)) >= 0,
        "Could not configure full block names");
    std::vector<char> name(static_cast<std::size_t>(maximum_name_length) + 1, '\0');
    const std::string expected_topology = config.topology == Topology::hex8    ? "HEX8"
                                          : config.topology == Topology::quad4 ? "QUAD4"
                                                                               : "QUAD8";
    const std::size_t nodes_per_element = config.topology == Topology::quad4 ? 4 : 8;
    std::vector<Element> result;
    std::vector<bool> found(config.blocks.size(), false);
    for (std::size_t b = 0; b < ids.size(); ++b) {
        std::fill(name.begin(), name.end(), '\0');
        require(ex_get_name(file.id(), EX_ELEM_BLOCK, ids[b], name.data()) >= 0, "Could not read material block name");
        const auto material = std::find_if(config.blocks.begin(), config.blocks.end(), [&](const BlockMaterial& item) {
            return item.name == name.data();
        });
        require(material != config.blocks.end(), "Unexpected material block: " + std::string(name.data()));
        const auto material_index = static_cast<std::size_t>(material - config.blocks.begin());
        require(!found[material_index], "Duplicate material block");
        found[material_index] = true;
        ex_block block{};
        block.type = EX_ELEM_BLOCK;
        block.id = ids[b];
        require(ex_get_block_param(file.id(), &block) >= 0, "Could not read block topology");
        require(block.topology == expected_topology
                    && block.num_nodes_per_entry == static_cast<std::int64_t>(nodes_per_element),
            "Unexpected element topology");
        require(block.num_entry == static_cast<std::int64_t>(material->element_count), "Incomplete element coverage");
        require(frame.block_names.at(b) == material->name
                    && frame.block_element_counts.at(b) == material->element_count,
            "Connectivity and result block metadata disagree");
        std::vector<std::int64_t> connectivity(material->element_count * nodes_per_element);
        require(ex_get_conn(file.id(), EX_ELEM_BLOCK, ids[b], connectivity.data(), nullptr, nullptr) >= 0,
            "Could not read reference element connectivity");
        for (std::size_t e = 0; e < material->element_count; ++e) {
            Element element;
            element.material = *material;
            element.temperature_nodes = config.topology == Topology::hex8 ? 8 : 4;
            for (std::size_t n = 0; n < nodes_per_element; ++n) {
                const auto node = connectivity[e * nodes_per_element + n];
                require(node > 0 && node <= static_cast<std::int64_t>(frame.nodes.size()), "Invalid connectivity node");
                element.nodes[n] = static_cast<std::size_t>(node - 1);
                for (std::size_t earlier = 0; earlier < n; ++earlier)
                    require(element.nodes[earlier] != element.nodes[n], "Repeated node within an element");
            }
            integrate_reference_geometry(element, config, frame.nodes);
            result.push_back(std::move(element));
        }
    }
    return result;
}

double specific_heat(const BlockMaterial& material, double temperature) {
    const double result =
        material.specific_heat
        + material.specific_heat_temperature_coefficient * (temperature - material.initial_temperature);
    require(std::isfinite(result) && result > 0.0, "Invalid independently evaluated specific heat");
    return result;
}

double reference_storage(const std::vector<Element>& elements,
    Topology topology,
    const std::vector<double>& current,
    const std::vector<double>& previous,
    double dt) {
    double result = 0.0;
    for (const Element& element : elements) {
        if (topology == Topology::quad8) {
            for (const CapacityPoint& point : element.points) {
                double temperature = 0.0, old_temperature = 0.0;
                for (std::size_t n = 0; n < element.temperature_nodes; ++n) {
                    temperature += point.temperature_shape[n] * current.at(element.nodes[n]);
                    old_temperature += point.temperature_shape[n] * previous.at(element.nodes[n]);
                }
                result += point.reference_measure * element.material.initial_density
                          * specific_heat(element.material, temperature) * (temperature - old_temperature) / dt;
            }
        } else {
            for (std::size_t n = 0; n < element.temperature_nodes; ++n) {
                const std::size_t node = element.nodes[n];
                result += element.nodal_reference_volume[n] * element.material.initial_density
                          * specific_heat(element.material, current.at(node)) * (current.at(node) - previous.at(node))
                          / dt;
            }
        }
    }
    require(std::isfinite(result), "Nonfinite independently reconstructed storage rate");
    return result;
}

void validate_frame(const ExodusResults& frame, const CaseSettings& config, std::size_t element_count) {
    require(frame.nodes.size() == config.node_count, "Incomplete nodal coverage");
    const std::vector<std::string> fields = config.topology == Topology::hex8 ? std::vector<std::string>{"temperature",
                                                                                    "displacement_x",
                                                                                    "displacement_y",
                                                                                    "displacement_z",
                                                                                    "reaction_heat_flux",
                                                                                    "reaction_force_x",
                                                                                    "reaction_force_y",
                                                                                    "reaction_force_z"}
                                                                              : std::vector<std::string>{"temperature",
                                                                                    "displacement_r",
                                                                                    "displacement_z",
                                                                                    "reaction_heat_flux",
                                                                                    "reaction_force_r",
                                                                                    "reaction_force_z"};
    for (const auto& field : fields) {
        const auto& values = frame.nodal(field);
        require(values.size() == config.node_count, "Incomplete nodal field: " + field);
        for (double value : values)
            require(std::isfinite(value), "Nonfinite nodal field: " + field);
    }
    for (double value : frame.global_variables)
        require(std::isfinite(value), "Nonfinite global output");
    if (config.topology != Topology::hex8) {
        const auto& active = frame.element("material_point_count");
        require(active.size() == element_count, "Incomplete material-point coverage");
        for (double value : active)
            require(value == static_cast<double>(config.material_point_count),
                "Incorrect number of active material points");
    }
    // Current HEX8 output repeats the single C3D8RT history in eight slots and
    // has no material_point_count field. These are output slots, not eight
    // independent material samples. RZ output explicitly marks inactive slots.
    const std::size_t finite_output_slots = config.topology == Topology::hex8 ? 8 : config.material_point_count;
    for (std::size_t v = 0; v < frame.element_variable_names.size(); ++v) {
        const std::string& field = frame.element_variable_names[v];
        const auto suffix = field.rfind("_q");
        if (suffix == std::string::npos)
            continue;
        const auto point = static_cast<std::size_t>(std::stoul(field.substr(suffix + 2)));
        const auto& values = frame.element_variables.at(v);
        require(values.size() == element_count, "Incomplete material field: " + field);
        for (double value : values)
            require(point < finite_output_slots ? std::isfinite(value) : std::isnan(value),
                "Invalid active/inactive material field: " + field);
    }
}

double check_thermal_balance(const ExodusResults& frame) {
    const double stored = frame.global("conservation_stored_heat_rate");
    const double convection = frame.global("conservation_convection_heat_rate");
    const double interface = frame.global("conservation_interface_heat_imbalance");
    const double generated = frame.global("conservation_generated_heat_rate");
    const double surface = frame.global("conservation_surface_heat_input_rate");
    const double dirichlet = frame.global("conservation_dirichlet_heat_input_rate");
    const double reconstructed = stored + convection + interface - generated - surface - dirichlet;
    const double scale = std::abs(stored) + std::abs(convection) + std::abs(interface) + std::abs(generated)
                         + std::abs(surface) + std::abs(dirichlet);
    const double balance = frame.global("conservation_global_thermal_balance");
    const double relative = scale == 0 ? 0 : std::abs(reconstructed) / scale;
    require(scale == 0 ? std::abs(reconstructed) <= zero_absolute_tolerance : relative < relative_tolerance,
        "Global thermal balance exceeds 1e-9 relative tolerance");
    require(std::abs(balance - reconstructed) <= 32.0 * std::numeric_limits<double>::epsilon() * scale,
        "Reported global thermal balance does not equal its heat-rate components");
    require(std::abs(frame.global("conservation_relative_thermal_balance") - relative)
                <= 32.0 * std::numeric_limits<double>::epsilon(),
        "Reported relative thermal balance is inconsistent");
    return relative;
}

void check(const std::string& name, const std::string& path) {
    const CaseSettings config = settings(name);
    const auto history = fuelsim::test::read_exodus_history(path);
    require(history.size() == config.increments + 1, "Expected the initial frame and every accepted increment");
    require(history.front().time == 0.0, "Initial temperature frame is missing");
    for (const auto& coordinate : history.front().nodes)
        for (double component : coordinate)
            require(std::isfinite(component), "Nonfinite reference coordinate");
    const auto elements = read_elements(path, config, history.front());
    const auto& initial_temperature = history.front().nodal("temperature");
    std::vector<bool> active_temperature(config.node_count, false);
    for (const Element& element : elements)
        for (std::size_t n = 0; n < element.temperature_nodes; ++n) {
            const std::size_t node = element.nodes[n];
            active_temperature.at(node) = true;
            require(initial_temperature.at(node) == element.material.initial_temperature,
                "Initial output temperature differs from this case's prescribed initial state");
        }
    fuelsim::test::FieldErrorMetrics storage;
    double maximum_relative_balance = 0.0;
    for (std::size_t step = 0; step < history.size(); ++step) {
        const auto& frame = history[step];
        const double expected_time = static_cast<double>(step) * config.time_step;
        const double time_tolerance =
            64.0 * std::numeric_limits<double>::epsilon() * static_cast<double>(config.increments) * config.time_step;
        require(std::isfinite(frame.time) && std::abs(frame.time - expected_time) <= time_tolerance,
            "Incomplete or unexpected accepted time sequence");
        validate_frame(frame, config, elements.size());
        if (config.topology == Topology::quad8) {
            const auto& mask = frame.nodal("temperature_active");
            for (std::size_t node = 0; node < active_temperature.size(); ++node)
                require(mask.at(node) == (active_temperature[node] ? 1.0 : 0.0),
                    "Incorrect active temperature-node mask");
        }
        if (step == 0)
            continue;
        const double dt = frame.time - history[step - 1].time;
        require(dt > 0.0, "Nonpositive accepted time increment");
        const double reference = reference_storage(elements,
            config.topology,
            frame.nodal("temperature"),
            history[step - 1].nodal("temperature"),
            dt);
        storage.add(frame.global("conservation_stored_heat_rate"), reference);
        maximum_relative_balance = std::max(maximum_relative_balance, check_thermal_balance(frame));
    }
    fuelsim::test::print_relative_metrics(name + "_independent_initial_mass_storage", storage);
    std::cout << name << " accepted_increments=" << config.increments << " elements=" << elements.size()
              << " active_temperature_nodes=" << std::count(active_temperature.begin(), active_temperature.end(), true)
              << " maximum_relative_thermal_balance=" << maximum_relative_balance
              << " zero_reference_absolute_tolerance_W=" << zero_absolute_tolerance << '\n';
    require(storage.value_count == config.increments, "Incomplete stored-heat-rate history comparison");
    require(fuelsim::test::relative_metrics_below(storage, relative_tolerance),
        "Independent reference-mass storage exceeds 1e-9 relative tolerance");
    require(storage.maximum_zero_reference_difference <= zero_absolute_tolerance,
        "Zero-reference storage exceeds its separate absolute tolerance");
}
} // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3, "Usage: fuelsim_production_initial_mass <case_name> <results.e>");
        check(argv[1], argv[2]);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Independent initial-mass production verification failed: " << error.what() << '\n';
        return 1;
    }
}
