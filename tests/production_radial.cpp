#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exodusII.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
using Row = std::map<std::string, double>;
using fuelsim::test::ExodusResults;
using fuelsim::test::FieldErrorMetrics;
using fuelsim::test::GroupedFieldErrorMetrics;
constexpr double relative_tolerance = 0.001;

std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) {
        if (!value.empty() && value.back() == '\r')
            value.pop_back();
        fields.push_back(value);
    }
    return fields;
}

double number(const std::string& value) {
    std::size_t count = 0;
    const double parsed = std::stod(value, &count);
    if (count != value.size() || !std::isfinite(parsed))
        throw std::runtime_error("Invalid finite numerical reference value");
    return parsed;
}

std::vector<Row> read_rows(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("Cannot read reference " + path.string());
    std::string line;
    std::getline(stream, line);
    const auto header = split(line);
    std::vector<Row> result;
    while (std::getline(stream, line)) {
        const auto fields = split(line);
        if (fields.size() != header.size())
            throw std::runtime_error("Incomplete reference row in " + path.string());
        Row row;
        for (std::size_t i = 0; i < header.size(); ++i)
            if (!row.emplace(header[i], number(fields[i])).second)
                throw std::runtime_error("Duplicate reference column");
        result.push_back(std::move(row));
    }
    if (result.empty())
        throw std::runtime_error("Reference contains no accepted results");
    return result;
}

std::map<std::string, std::string> read_summary(const std::string& path) {
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("Missing production summary");
    std::string line;
    std::getline(stream, line);
    if (line != "metric,value")
        throw std::runtime_error("Unexpected production summary header");
    std::map<std::string, std::string> values;
    while (std::getline(stream, line)) {
        const auto fields = split(line);
        if (fields.size() != 2 || !values.emplace(fields[0], fields[1]).second)
            throw std::runtime_error("Invalid production summary row");
    }
    if (values.at("completed") != "true" && values.at("completed") != "1")
        throw std::runtime_error("Production solve did not complete");
    return values;
}

void require_summary(const std::string& path, double final_time) {
    const auto values = read_summary(path);
    if (number(values.at("accepted_steps")) != 10 || number(values.at("rejected_steps")) != 0
        || std::abs(number(values.at("committed_time")) - final_time) > 1e-12)
        throw std::runtime_error("Production accepted history differs from the native reference");
}

void require_time(double actual, double reference) {
    if (!std::isfinite(actual) || !std::isfinite(reference)
        || std::abs(actual - reference) > 8e-8 * std::max(std::abs(actual), std::abs(reference)) + 1e-16)
        throw std::runtime_error("Production and native reference times do not match");
}

double constrained_zero(double native, double tolerance) {
    if (!std::isfinite(native) || std::abs(native) > tolerance)
        throw std::runtime_error("Native result violates an exact prescribed-zero constraint");
    return 0.0;
}

bool check(const std::string& name, const FieldErrorMetrics& metric, double zero_tolerance) {
    if (metric.value_count == 0)
        throw std::runtime_error("Empty scalar acceptance metric");
    if (metric.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, metric);
    else
        fuelsim::test::print_absolute_metrics(name, metric);
    std::cout << name << "_zero_reference_count=" << metric.zero_reference_count << '\n'
              << name << "_maximum_zero_reference_absolute_difference=" << metric.maximum_zero_reference_difference
              << '\n';
    return (!metric.has_relative_norm() || fuelsim::test::relative_metrics_below(metric, relative_tolerance))
           && metric.maximum_zero_reference_difference <= zero_tolerance;
}

bool check(const std::string& name, const GroupedFieldErrorMetrics& metric, double zero_tolerance) {
    if (metric.group_count == 0)
        throw std::runtime_error("Empty tensor acceptance metric");
    if (metric.has_relative_norm())
        fuelsim::test::print_grouped_relative_metrics(name, metric);
    std::cout << name << "_zero_reference_count=" << metric.zero_reference_count << '\n'
              << name << "_maximum_zero_reference_absolute_difference=" << metric.maximum_zero_reference_difference
              << '\n';
    return (!metric.has_relative_norm() || fuelsim::test::grouped_relative_metrics_below(metric, relative_tolerance))
           && metric.maximum_zero_reference_difference <= zero_tolerance;
}

void add_tensor(GroupedFieldErrorMetrics& metric,
    const ExodusResults& frame,
    const Row& reference,
    const std::string& prefix,
    std::size_t q,
    std::size_t element = 0) {
    const std::array<std::string, 4> components{"rr", "zz", "hoop", "rz"};
    std::array<double, 4> actual{}, expected{};
    for (std::size_t i = 0; i < 4; ++i) {
        const std::string field = prefix + "_" + components[i];
        actual[i] = frame.element(field + "_q" + std::to_string(q)).at(element);
        expected[i] = reference.at(field);
        if (i == 3) {
            actual[i] *= std::sqrt(2.0);
            expected[i] *= std::sqrt(2.0);
        }
    }
    metric.add(actual.data(), expected.data(), actual.size());
}

void require_uniform_control_attributes(const std::string& output) {
    int cpu_word_size = 8, disk_word_size = 0;
    float version = 0.0F;
    int file = ex_open(output.c_str(), EX_READ | EX_ALL_INT64_API, &cpu_word_size, &disk_word_size, &version);
    if (file < 0)
        throw std::runtime_error("Cannot open steady result control-node attributes");
    try {
        ex_init_params parameters{};
        if (ex_get_init_ext(file, &parameters) < 0 || parameters.num_dim != 2 || parameters.num_nodes != 4
            || parameters.num_elem != 1 || parameters.num_elem_blk != 1)
            throw std::runtime_error("Steady radial result must retain the complete source mesh");
        std::int64_t block_id = 0;
        if (ex_get_ids(file, EX_ELEM_BLOCK, &block_id) < 0)
            throw std::runtime_error("Cannot read steady radial element block");
        ex_block block{};
        block.type = EX_ELEM_BLOCK;
        block.id = block_id;
        if (ex_get_block_param(file, &block) < 0 || std::string(block.topology) != "BAR2" || block.num_entry != 1
            || block.num_nodes_per_entry != 2 || block.num_attribute != 2)
            throw std::runtime_error("Steady radial output lost its BAR2 topology or axial attributes");
        std::array<std::int64_t, 2> connectivity{};
        std::array<double, 2> attributes{};
        std::array<std::array<char, MAX_STR_LENGTH + 1>, 2> names{};
        std::array<char*, 2> name_pointers{names[0].data(), names[1].data()};
        if (ex_get_conn(file, EX_ELEM_BLOCK, block_id, connectivity.data(), nullptr, nullptr) < 0
            || ex_get_attr_names(file, EX_ELEM_BLOCK, block_id, name_pointers.data()) < 0
            || ex_get_attr(file, EX_ELEM_BLOCK, block_id, attributes.data()) < 0)
            throw std::runtime_error("Cannot read steady radial connectivity and axial attributes");
        const std::map<std::string, double> references{{names[0].data(), attributes[0]},
            {names[1].data(), attributes[1]}};
        if (connectivity != std::array<std::int64_t, 2>{1, 2}
            || references != std::map<std::string, double>{{"axial_lower_node", 3.0}, {"axial_upper_node", 4.0}})
            throw std::runtime_error("Steady radial output changed its radial and axial source-node mapping");
        const int status = ex_close(file);
        file = -1;
        if (status < 0)
            throw std::runtime_error("Cannot close steady radial result after attribute checks");
    } catch (...) {
        if (file >= 0)
            (void)ex_close(file);
        throw;
    }
}

bool steady_finite(const std::string& output, const std::string& summary) {
    const auto values = read_summary(summary);
    if (values.at("problem") != "steady" || number(values.at("load_steps_completed")) != 5.0
        || number(values.at("rejected_load_steps")) != 0.0 || number(values.at("load_cutbacks")) != 0.0
        || number(values.at("regions")) != 1.0 || number(values.at("contacts")) != 0.0)
        throw std::runtime_error("Steady radial production solve did not complete the prescribed five load steps");
    require_uniform_control_attributes(output);
    const auto frames = fuelsim::test::read_exodus_history(output);
    if (frames.size() != 1 || frames.front().time != 0.0 || frames.front().step_count != 1)
        throw std::runtime_error("Steady radial output must contain one complete final frame");
    const auto& frame = frames.front();
    if (frame.nodal_variable_names
            != std::vector<std::string>{"temperature", "displacement_r", "displacement_z", "node_role"}
        || frame.global_variable_names != std::vector<std::string>{"load_factor"}
        || frame.element("material_point_count") != std::vector<double>{2.0})
        throw std::runtime_error("Steady radial result fields or active material-point count changed");
    constexpr std::array<std::array<double, 3>, 4> coordinates{
        {{{0.004, 0.005, 0.0}}, {{0.005, 0.005, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, 0.01, 0.0}}}};
    if (frame.nodes.size() != coordinates.size())
        throw std::runtime_error("Steady radial result has the wrong source-node count");
    FieldErrorMetrics temperature, radial_displacement, axial_displacement, axial_strain, load_factor;
    FieldErrorMetrics zero_stress, zero_strain, zero_axial_force, geometry, constrained_displacement;
    const double imposed_strain = 1e-4 * (1100.0 - 600.0);
    const double stretch = (2.0 + imposed_strain) / (2.0 - imposed_strain);
    const double extension = stretch - 1.0;
    for (std::size_t n = 0; n < frame.nodes.size(); ++n) {
        for (std::size_t component = 0; component < 3; ++component)
            geometry.add(frame.nodes[n][component], coordinates[n][component]);
        if (frame.nodal("node_role").at(n) != (n < 2 ? 1.0 : 2.0))
            throw std::runtime_error("Steady radial output lost the radial and axial node roles");
        if (n < 2) {
            temperature.add(frame.nodal("temperature").at(n), 1100.0);
            radial_displacement.add(frame.nodal("displacement_r").at(n), coordinates[n][0] * extension);
        } else if (!std::isnan(frame.nodal("temperature").at(n)) || !std::isnan(frame.nodal("displacement_r").at(n)))
            throw std::runtime_error(
                "Steady axial control nodes must keep inactive radial and temperature fields undefined");
        if (n == 2)
            constrained_displacement.add(frame.nodal("displacement_z").at(n), 0.0);
        else
            axial_displacement.add(frame.nodal("displacement_z").at(n), coordinates[n][1] * extension);
    }
    geometry.add(frame.element("reference_height").at(0), 0.01);
    axial_strain.add(frame.element("axial_strain").at(0), extension);
    zero_axial_force.add(frame.element("axial_force").at(0), 0.0);
    load_factor.add(frame.global("load_factor"), 1.0);
    const double gauss = 1.0 / std::sqrt(3.0);
    for (std::size_t q = 0; q < 2; ++q) {
        const std::string suffix = "_q" + std::to_string(q);
        const double radius = 0.0045 + (q == 0 ? -gauss : gauss) * 0.0005;
        geometry.add(frame.element("reference_r" + suffix).at(0), radius);
        geometry.add(frame.element("reference_z" + suffix).at(0), 0.005);
        geometry.add(frame.element("reference_measure" + suffix).at(0), std::acos(-1.0) * radius * 0.001 * 0.01);
        temperature.add(frame.element("temperature" + suffix).at(0), 1100.0);
        for (const char* component : {"rr", "zz", "hoop", "rz"}) {
            zero_stress.add(frame.element(std::string("stress_") + component + suffix).at(0), 0.0);
            for (const char* prefix : {"elastic_", "plastic_", "creep_"})
                zero_strain.add(frame.element(std::string(prefix) + component + suffix).at(0), 0.0);
        }
        zero_strain.add(frame.element("equiv_plastic" + suffix).at(0), 0.0);
        zero_strain.add(frame.element("equiv_creep" + suffix).at(0), 0.0);
    }
    bool passed = true;
    passed = check("steady_temperature", temperature, 1e-11) && passed;
    passed = check("steady_radial_displacement", radial_displacement, 1e-13) && passed;
    passed = check("steady_axial_displacement", axial_displacement, 1e-13) && passed;
    passed = check("steady_axial_strain", axial_strain, 1e-13) && passed;
    passed = check("steady_load_factor", load_factor, 1e-13) && passed;
    passed = check("steady_reference_geometry", geometry, 1e-14) && passed;
    passed = check("steady_constrained_axial_displacement", constrained_displacement, 1e-13) && passed;
    passed = check("steady_stress", zero_stress, 1e-4) && passed;
    passed = check("steady_elastic_and_inelastic_strain", zero_strain, 1e-13) && passed;
    passed = check("steady_axial_force", zero_axial_force, 1e-8) && passed;
    std::cout << "steady_analytic_stretch=" << stretch
              << "\nradial_gps_steady_finite_qualification=" << (passed ? "passed" : "failed") << '\n';
    return passed;
}

bool uniform(const std::string& output,
    const std::string& summary,
    const std::filesystem::path& references,
    bool finite) {
    require_summary(summary, 1.0);
    const auto frames = fuelsim::test::read_exodus_history(output);
    const std::string case_name = finite ? "gps_uniform_axial_finite" : "gps_uniform_axial";
    const auto native_nodes = read_rows(references / (case_name + "_nodes.csv"));
    const auto native_points = read_rows(references / (case_name + "_points.csv"));
    if (frames.size() != 11 || frames.front().time != 0.0 || native_nodes.size() != 40 || native_points.size() != 40)
        throw std::runtime_error("Uniform radial qualification requires all ten increments and every native sample");
    FieldErrorMetrics temperature, radial_reaction, axial_reaction, heat_reaction, axial_force, axial_strain;
    FieldErrorMetrics midpoint_displacement, zero_mechanism, prescribed_zero_displacement, zero_stress, zero_shear;
    GroupedFieldErrorMetrics displacement, stress, elastic;
    constexpr std::array<std::size_t, 4> radial_node{0, 1, 1, 0};
    constexpr std::array<std::size_t, 4> axial_node{2, 2, 3, 3};
    for (std::size_t step = 1; step < frames.size(); ++step) {
        const auto& frame = frames[step];
        require_time(frame.time, static_cast<double>(step) / 10.0);
        if (frame.nodes.size() != 4 || frame.element("material_point_count").size() != 1
            || frame.element("material_point_count").at(0) != 2.0)
            throw std::runtime_error("Uniform radial geometry or active quadrature count changed");
        const auto& t = frame.nodal("temperature");
        const auto& ur = frame.nodal("displacement_r");
        const auto& uz = frame.nodal("displacement_z");
        const auto& role = frame.nodal("node_role");
        const auto& rr = frame.nodal("reaction_force_r");
        const auto& rz = frame.nodal("reaction_force_z");
        const auto& heat = frame.nodal("reaction_heat_flux");
        for (std::size_t n = 0; n < 4; ++n) {
            if (role.at(n) != (n < 2 ? 1.0 : 2.0))
                throw std::runtime_error("Radial and axial control node roles changed");
            if ((n < 2 && !std::isnan(rz.at(n)))
                || (n >= 2
                    && (!std::isnan(t.at(n)) || !std::isnan(ur.at(n)) || !std::isnan(rr.at(n))
                        || !std::isnan(heat.at(n)))))
                throw std::runtime_error("Inactive control-node fields must remain explicitly undefined");
            const auto& reference = native_nodes[(step - 1) * 4 + n];
            require_time(frame.time, reference.at("time"));
            if (reference.at("node") != static_cast<double>(n + 1))
                throw std::runtime_error("Native annulus node correspondence changed");
            temperature.add(t.at(radial_node[n]), reference.at("temperature"));
            const std::array<double, 2> actual{ur.at(radial_node[n]), uz.at(axial_node[n])};
            const std::array<double, 2> expected{reference.at("ur"),
                n < 2 ? constrained_zero(reference.at("uz"), 1e-13) : reference.at("uz")};
            displacement.add(actual.data(), expected.data(), actual.size());
            if (n < 2)
                prescribed_zero_displacement.add(uz.at(axial_node[n]), expected[1]);
        }
        const auto& bottom_inner = native_nodes[(step - 1) * 4];
        const auto& bottom_outer = native_nodes[(step - 1) * 4 + 1];
        const auto& top_outer = native_nodes[(step - 1) * 4 + 2];
        const auto& top_inner = native_nodes[(step - 1) * 4 + 3];
        for (std::size_t n = 0; n < 2; ++n) {
            midpoint_displacement.add(uz.at(n), 0.5 * (uz.at(2) + uz.at(3)));
            radial_reaction.add(rr.at(n), constrained_zero(native_nodes[(step - 1) * 4 + n].at("rf_r"), 1e-8));
        }
        heat_reaction.add(heat.at(0), bottom_inner.at("reaction_heat") + top_inner.at("reaction_heat"));
        heat_reaction.add(heat.at(1), bottom_outer.at("reaction_heat") + top_outer.at("reaction_heat"));
        axial_reaction.add(rz.at(2), bottom_inner.at("rf_z") + bottom_outer.at("rf_z"));
        axial_reaction.add(rz.at(3), constrained_zero(top_inner.at("rf_z") + top_outer.at("rf_z"), 1e-8));
        axial_strain.add(frame.element("axial_strain").at(0), (top_inner.at("uz") - bottom_inner.at("uz")) / 0.01);
        double reference_force = 0.0;
        for (std::size_t point = 0; point < 4; ++point) {
            const auto& reference = native_points[(step - 1) * 4 + point];
            require_time(frame.time, reference.at("time"));
            if (reference.at("element") != 1.0 || reference.at("point") != static_cast<double>(point + 1))
                throw std::runtime_error("Native annulus material-point correspondence changed");
            const std::size_t q = point % 2;
            add_tensor(stress, frame, reference, "stress", q);
            add_tensor(elastic, frame, reference, "elastic", q);
            for (const char* component : {"rr", "hoop", "rz"}) {
                const std::string field = std::string("stress_") + component;
                zero_stress.add(frame.element(field + "_q" + std::to_string(q)).at(0),
                    constrained_zero(reference.at(field), 1e-4));
            }
            zero_shear.add(frame.element("elastic_rz_q" + std::to_string(q)).at(0),
                constrained_zero(reference.at("elastic_rz"), 1e-13));
            reference_force += reference.at("stress_zz") * reference.at("volume")
                               / (finite ? 0.01 + top_inner.at("uz") - bottom_inner.at("uz") : 0.01);
        }
        axial_force.add(frame.element("axial_force").at(0), reference_force);
        for (std::size_t q = 0; q < 2; ++q) {
            const std::string suffix = "_q" + std::to_string(q);
            for (const char* prefix : {"plastic_", "creep_"})
                for (const char* component : {"rr", "zz", "hoop", "rz"})
                    zero_mechanism.add(frame.element(std::string(prefix) + component + suffix).at(0), 0.0);
            zero_mechanism.add(frame.element("equiv_plastic" + suffix).at(0), 0.0);
            zero_mechanism.add(frame.element("equiv_creep" + suffix).at(0), 0.0);
        }
    }
    bool passed = true;
    passed = check("temperature", temperature, 1e-11) && passed;
    passed = check("displacement", displacement, 1e-13) && passed;
    passed = check("stress", stress, 1e-4) && passed;
    passed = check("elastic_strain", elastic, 1e-13) && passed;
    passed = check("radial_reaction", radial_reaction, 1e-8) && passed;
    passed = check("axial_reaction", axial_reaction, 1e-8) && passed;
    passed = check("heat_reaction", heat_reaction, 1e-11) && passed;
    passed = check("axial_force", axial_force, 1e-8) && passed;
    passed = check("axial_strain", axial_strain, 1e-13) && passed;
    passed = check("radial_node_axial_interpolation", midpoint_displacement, 1e-13) && passed;
    passed = check("prescribed_zero_displacement", prescribed_zero_displacement, 1e-13) && passed;
    passed = check("analytical_zero_stress", zero_stress, 1e-4) && passed;
    passed = check("analytical_zero_elastic_shear", zero_shear, 1e-13) && passed;
    passed = check("inactive_inelastic_history", zero_mechanism, 1e-13) && passed;
    std::cout << "radial_gps_abaqus_qualification=" << (passed ? "passed" : "failed") << '\n';
    return passed;
}

using ContactKey = std::tuple<double, std::string, std::string, std::size_t, std::size_t>;

std::map<ContactKey, double> read_contact(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream)
        throw std::runtime_error("Cannot read native contact reference");
    std::string line;
    std::getline(stream, line);
    if (line != "time,field,node,component,value")
        throw std::runtime_error("Unknown native contact reference format");
    std::map<ContactKey, double> result;
    while (std::getline(stream, line)) {
        const auto fields = split(line);
        if (fields.size() != 5)
            throw std::runtime_error("Invalid native contact row");
        std::istringstream label(fields[1]);
        std::string name, pair, extra;
        if (!(label >> name >> pair) || (label >> extra))
            throw std::runtime_error("Invalid native contact field name");
        const double node = number(fields[2]), component = number(fields[3]);
        if (node < 1 || std::floor(node) != node || component < 0 || component > 1
            || std::floor(component) != component)
            throw std::runtime_error("Invalid native contact node or component");
        const ContactKey key{number(fields[0]),
            pair,
            name,
            static_cast<std::size_t>(node),
            static_cast<std::size_t>(component)};
        if (!result.emplace(key, number(fields[4])).second)
            throw std::runtime_error("Duplicate native contact sample");
    }
    return result;
}

bool contact_small(const std::string& output, const std::string& summary, const std::filesystem::path& references) {
    require_summary(summary, 10.0);
    const auto frames = fuelsim::test::read_exodus_history(output);
    const auto native_nodes = read_rows(references / "gps_two_slice_contact_nodes.csv");
    const auto native_points = read_rows(references / "gps_two_slice_contact_points.csv");
    const auto native_contacts = read_contact(references / "gps_two_slice_contact_contact.csv");
    if (frames.size() != 11 || frames.front().time != 0 || native_nodes.size() != 160 || native_points.size() != 160)
        throw std::runtime_error("Contact operator requires every one of the ten accepted increments");
    constexpr std::array<std::size_t, 16> radial_node{0, 1, 1, 0, 2, 3, 3, 2, 4, 5, 5, 4, 6, 7, 7, 6};
    constexpr std::array<std::size_t, 16> axial_node{8, 8, 9, 9, 11, 11, 12, 12, 9, 9, 10, 10, 12, 12, 13, 13};
    constexpr std::array<std::size_t, 4> source_element{0, 2, 1, 3};
    const std::array<std::string, 2> pair{"SECONDARY_ONE/PRIMARY_ONE", "SECONDARY_TWO/PRIMARY_TWO"};
    const std::array<std::array<std::size_t, 2>, 2> secondary{{{2, 3}, {10, 11}}};
    const std::array<std::array<std::size_t, 2>, 2> primary{{{5, 8}, {13, 16}}};
    FieldErrorMetrics temperature, pressure, gap, traction, slip, elastic_slip, area, heat_flux, heat_rate;
    FieldErrorMetrics normal_force, tangent_force, heat_reaction, zero_mechanism, native_force_balance;
    FieldErrorMetrics body_axial_reaction, production_axial_force_balance;
    FieldErrorMetrics prescribed_zero_displacement, midpoint_displacement;
    GroupedFieldErrorMetrics displacement, diagnostic_stress, diagnostic_elastic;
    std::size_t positive_sliding = 0, negative_sliding = 0, sticking = 0, opening = 0;
    for (std::size_t step = 1; step < frames.size(); ++step) {
        const auto& frame = frames[step];
        require_time(frame.time, static_cast<double>(step));
        if (frame.nodes.size() != 14 || frame.element("material_point_count") != std::vector<double>(4, 2.0))
            throw std::runtime_error("Two-slice radial source layout or material-point counts changed");
        const auto& t = frame.nodal("temperature");
        const auto& ur = frame.nodal("displacement_r");
        const auto& uz = frame.nodal("displacement_z");
        const auto& node_role = frame.nodal("node_role");
        for (std::size_t n = 0; n < 14; ++n) {
            if (node_role.at(n) != (n < 8 ? 1.0 : 2.0))
                throw std::runtime_error("Two-slice source node roles changed");
            if ((n < 8 && !std::isnan(frame.nodal("reaction_force_z").at(n)))
                || (n >= 8
                    && (!std::isnan(t.at(n)) || !std::isnan(ur.at(n))
                        || !std::isnan(frame.nodal("reaction_force_r").at(n))
                        || !std::isnan(frame.nodal("reaction_heat_flux").at(n)))))
                throw std::runtime_error("Inactive contact control-node fields must remain explicitly undefined");
        }
        for (std::size_t n = 0; n < 16; ++n) {
            const auto& reference = native_nodes[(step - 1) * 16 + n];
            require_time(frame.time, reference.at("time"));
            if (reference.at("node") != static_cast<double>(n + 1))
                throw std::runtime_error("Native contact node mapping changed");
            temperature.add(t.at(radial_node[n]), reference.at("temperature"));
            const bool outer = (n / 4) % 2 == 1;
            const std::array<double, 2> actual{ur.at(radial_node[n]), uz.at(axial_node[n])};
            const std::array<double, 2> expected{
                outer || step == 6 || step == 7 ? constrained_zero(reference.at("ur"), 1e-13) : reference.at("ur"),
                outer || step == 1 ? constrained_zero(reference.at("uz"), 1e-13) : reference.at("uz")};
            displacement.add(actual.data(), expected.data(), actual.size());
            if (outer || step == 6 || step == 7)
                prescribed_zero_displacement.add(actual[0], expected[0]);
            if (outer || step == 1)
                prescribed_zero_displacement.add(actual[1], expected[1]);
        }
        for (std::size_t n = 0; n < 8; ++n) {
            const std::size_t lower = 8 + n / 4 + ((n % 4) / 2) * 3;
            midpoint_displacement.add(uz.at(n), 0.5 * (uz.at(lower) + uz.at(lower + 1)));
        }
        for (std::size_t point = 0; point < 16; ++point) {
            const auto& reference = native_points[(step - 1) * 16 + point];
            require_time(frame.time, reference.at("time"));
            if (reference.at("element") != static_cast<double>(point / 4 + 1)
                || reference.at("point") != static_cast<double>(point % 4 + 1))
                throw std::runtime_error("Native contact material-point mapping changed");
            const std::size_t element = source_element[point / 4];
            add_tensor(diagnostic_stress, frame, reference, "stress", point % 2, element);
            add_tensor(diagnostic_elastic, frame, reference, "elastic", point % 2, element);
        }
        double total_reference_tangent = 0.0;
        for (std::size_t layer = 0; layer < 2; ++layer) {
            const std::string suffix = "_layer" + std::to_string(layer + 1);
            const double height = layer == 0 ? 0.01 : 0.02;
            const double point_area = std::acos(-1.0) * 0.005 * height;
            const double allowable = height * 1e-5;
            const double native_time = native_nodes[(step - 1) * 16].at("time");
            const auto native = [&](const std::string& field, std::size_t node, std::size_t component = 0) {
                return native_contacts.at({native_time, pair[layer], field, node, component});
            };
            double reference_normal = 0.0, reference_tangent = 0.0, reference_heat = 0.0;
            for (std::size_t q = 0; q < 2; ++q) {
                const std::size_t node = secondary[layer][q];
                const std::string quadrature = suffix + "_q" + std::to_string(q);
                const auto actual = [&](const std::string& field) {
                    return frame.element("contact_" + field + quadrature).at(layer);
                };
                const double p = native("CPRESS", node);
                const bool open = step == 6 || step == 7;
                const bool zero_traction = open || step == 1 || step == 8;
                const double tau =
                    zero_traction ? constrained_zero(native("CSHEAR1", node), 1e-4) : native("CSHEAR1", node);
                const double total_slip =
                    step == 1 ? constrained_zero(native("CSLIP1", node), 1e-13) : native("CSLIP1", node);
                const double spring = zero_traction ? 0.0 : tau * allowable / (0.3 * p);
                pressure.add(actual("pressure"), open ? constrained_zero(p, 1e-4) : p);
                gap.add(actual("gap"), native("COPEN", node));
                traction.add(actual("tangential_traction"), tau);
                slip.add(actual("total_tangential_slip"), total_slip);
                elastic_slip.add(actual("elastic_tangential_slip"), spring);
                area.add(actual("area"), point_area);
                heat_flux.add(actual("heat_flux"), native("HFL", node));
                const bool sliding = !open && std::abs(tau) / (0.3 * p) > 1.0 - 1e-10;
                if (actual("sliding") != (sliding ? 1.0 : 0.0))
                    throw std::runtime_error("Production and native Coulomb branch classification differ");
                if (open)
                    ++opening;
                else if (!sliding)
                    ++sticking;
                else if (tau > 0)
                    ++positive_sliding;
                else
                    ++negative_sliding;
                reference_normal -= native("CNORMF", node, 0);
                reference_tangent -= native("CSHEARF", node, 1);
                reference_heat += native("HFL", node) * point_area;
            }
            // Native surface contact can distribute an opposite radial force couple
            // to the primary endpoints. Conservation is an integrated pair property.
            for (const char* field : {"CNORMF", "CSHEARF"})
                for (std::size_t component = 0; component < 2; ++component) {
                    double balance = 0.0;
                    for (std::size_t q = 0; q < 2; ++q)
                        balance +=
                            native(field, secondary[layer][q], component) + native(field, primary[layer][q], component);
                    native_force_balance.add(balance, 0.0);
                }
            normal_force.add(frame.global("contact_force" + suffix), reference_normal);
            const double expected_tangent = step == 1 || step == 6 || step == 7 || step == 8
                                                ? constrained_zero(reference_tangent, 1e-8)
                                                : reference_tangent;
            tangent_force.add(frame.global("contact_tangential_force" + suffix), expected_tangent);
            total_reference_tangent += expected_tangent;
            heat_rate.add(frame.global("contact_heat_rate" + suffix), reference_heat);
            const std::size_t native_offset = (step - 1) * 16 + layer * 8;
            heat_reaction.add(frame.nodal("reaction_heat_flux").at(layer * 4),
                native_nodes[native_offset].at("reaction_heat") + native_nodes[native_offset + 3].at("reaction_heat"));
            heat_reaction.add(frame.nodal("reaction_heat_flux").at(layer * 4 + 3),
                native_nodes[native_offset + 5].at("reaction_heat")
                    + native_nodes[native_offset + 6].at("reaction_heat"));
        }
        // Every axial control is prescribed, so its residual does not constrain
        // the solve. Check actual reactions on both bodies independently: each
        // body's internal axial forces cancel when its three controls are summed.
        // This remains valid despite the different coarse body stress operators.
        double inner_axial_reaction = 0.0, outer_axial_reaction = 0.0;
        for (std::size_t station = 0; station < 3; ++station) {
            inner_axial_reaction += frame.nodal("reaction_force_z").at(8 + station);
            outer_axial_reaction += frame.nodal("reaction_force_z").at(11 + station);
        }
        body_axial_reaction.add(inner_axial_reaction, total_reference_tangent);
        body_axial_reaction.add(outer_axial_reaction, -total_reference_tangent);
        production_axial_force_balance.add(inner_axial_reaction + outer_axial_reaction, 0.0);
        for (std::size_t element = 0; element < 4; ++element)
            for (std::size_t q = 0; q < 2; ++q) {
                const std::string suffix = "_q" + std::to_string(q);
                for (const char* prefix : {"plastic_", "creep_"})
                    for (const char* component : {"rr", "zz", "hoop", "rz"})
                        zero_mechanism.add(frame.element(std::string(prefix) + component + suffix).at(element), 0.0);
                zero_mechanism.add(frame.element("equiv_plastic" + suffix).at(element), 0.0);
                zero_mechanism.add(frame.element("equiv_creep" + suffix).at(element), 0.0);
            }
    }
    if (positive_sliding != 8 || negative_sliding != 4 || sticking != 20 || opening != 8)
        throw std::runtime_error("Required sticking, reverse sliding, opening and recontact samples are missing");
    bool passed = true;
    passed = check("temperature", temperature, 1e-11) && passed;
    passed = check("prescribed_displacement", displacement, 1e-13) && passed;
    passed = check("prescribed_zero_displacement", prescribed_zero_displacement, 1e-13) && passed;
    passed = check("radial_node_axial_interpolation", midpoint_displacement, 1e-13) && passed;
    passed = check("contact_pressure", pressure, 1e-4) && passed;
    passed = check("contact_gap", gap, 1e-13) && passed;
    passed = check("contact_tangential_traction", traction, 1e-4) && passed;
    passed = check("contact_total_slip", slip, 1e-13) && passed;
    passed = check("contact_elastic_slip", elastic_slip, 1e-13) && passed;
    passed = check("contact_area", area, 1e-13) && passed;
    passed = check("contact_heat_flux", heat_flux, 1e-7) && passed;
    passed = check("contact_heat_rate", heat_rate, 1e-8) && passed;
    passed = check("contact_normal_force", normal_force, 1e-8) && passed;
    passed = check("contact_tangential_force", tangent_force, 1e-8) && passed;
    passed = check("thermal_boundary_reaction", heat_reaction, 1e-8) && passed;
    passed = check("native_contact_force_conservation", native_force_balance, 1e-8) && passed;
    passed = check("production_body_axial_reaction", body_axial_reaction, 1e-8) && passed;
    passed = check("production_axial_force_conservation", production_axial_force_balance, 1e-8) && passed;
    passed = check("inactive_inelastic_history", zero_mechanism, 1e-13) && passed;
    fuelsim::test::print_grouped_relative_metrics("diagnostic_body_stress", diagnostic_stress);
    fuelsim::test::print_grouped_relative_metrics("diagnostic_body_elastic_strain", diagnostic_elastic);
    std::cout << "positive_sliding_samples=" << positive_sliding << "\nnegative_sliding_samples=" << negative_sliding
              << "\nsticking_samples=" << sticking << "\nopen_samples=" << opening
              << "\nqualification_scope=prescribed_motion_contact_and_thermal_operators\n"
              << "radial_gps_contact_abaqus_qualification=" << (passed ? "passed" : "failed") << '\n';
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: production_radial <uniform-small|uniform-finite|contact-small|steady-finite> <results.e> "
                     "<summary.csv> "
                     "<reference_dir>\n";
        return 2;
    }
    try {
        const std::string mode = argv[1];
        if (mode == "uniform-small" || mode == "uniform-finite")
            return uniform(argv[2], argv[3], argv[4], mode == "uniform-finite") ? 0 : 1;
        if (mode == "contact-small")
            return contact_small(argv[2], argv[3], argv[4]) ? 0 : 1;
        if (mode == "steady-finite")
            return steady_finite(argv[2], argv[3]) ? 0 : 1;
        throw std::invalid_argument("Unsupported radial qualification mode: " + mode);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
