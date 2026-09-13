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

enum class NonuniformReference { finite_constant, finite_elastic_temperature, small_thermal, finite_thermal };

bool nonuniform(const std::string& output,
    const std::string& summary,
    const std::filesystem::path& references,
    NonuniformReference reference_case) {
    const bool variable_thermal =
        reference_case == NonuniformReference::small_thermal || reference_case == NonuniformReference::finite_thermal;
    const bool finite = reference_case != NonuniformReference::small_thermal;
    const bool temperature_dependent_elasticity = reference_case == NonuniformReference::finite_elastic_temperature;
    const std::string native_name =
        variable_thermal ? (finite ? "gps_thermal_finite" : "gps_thermal_small") : "gps_nonuniform_finite";
    const auto reference_directory =
        temperature_dependent_elasticity ? references / "nonuniform_elastic_temperature_probe" : references;
    require_summary(summary, 1.0);
    const auto frames = fuelsim::test::read_exodus_history(output);
    const auto native_nodes = read_rows(reference_directory / (native_name + "_nodes.csv"));
    const auto native_points = read_rows(reference_directory / (native_name + "_points.csv"));
    if (frames.size() != 11 || frames.front().time != 0.0 || native_nodes.size() != 40 || native_points.size() != 40)
        throw std::runtime_error("Nonuniform finite qualification requires all ten increments and every native sample");
    constexpr std::array<std::size_t, 4> radial_node{0, 1, 1, 0};
    constexpr std::array<std::size_t, 4> axial_node{2, 2, 3, 3};
    FieldErrorMetrics temperature, radial_reaction, axial_reaction, axial_force, axial_strain;
    FieldErrorMetrics midpoint_displacement, prescribed_zero_displacement, zero_stress_shear, zero_elastic_shear;
    FieldErrorMetrics zero_mechanism, elastic_energy_change;
    FieldErrorMetrics heat_reaction, stored_heat_rate, generated_heat_rate, dirichlet_heat_rate, thermal_balance;
    FieldErrorMetrics source_node_reaction_change;
    GroupedFieldErrorMetrics displacement, stress, elastic;
    GroupedFieldErrorMetrics reconstructed_heat_flux_diagnostic;
    double previous_native_elastic_energy = 0.0;
    for (std::size_t step = 1; step < frames.size(); ++step) {
        const auto& frame = frames[step];
        require_time(frame.time, static_cast<double>(step) / 10.0);
        if (frame.nodes.size() != 4 || frame.element("material_point_count") != std::vector<double>{2.0})
            throw std::runtime_error("Nonuniform finite geometry or active quadrature count changed");
        const auto& t = frame.nodal("temperature");
        const auto& ur = frame.nodal("displacement_r");
        const auto& uz = frame.nodal("displacement_z");
        const auto& role = frame.nodal("node_role");
        const auto& rr = frame.nodal("reaction_force_r");
        const auto& rz = frame.nodal("reaction_force_z");
        const auto& heat = frame.nodal("reaction_heat_flux");
        for (std::size_t n = 0; n < 4; ++n) {
            if (role.at(n) != (n < 2 ? 1.0 : 2.0))
                throw std::runtime_error("Nonuniform finite radial and axial control node roles changed");
            if ((n < 2 && !std::isnan(rz.at(n)))
                || (n >= 2
                    && (!std::isnan(t.at(n)) || !std::isnan(ur.at(n)) || !std::isnan(rr.at(n))
                        || !std::isnan(frame.nodal("reaction_heat_flux").at(n)))))
                throw std::runtime_error("Inactive nonuniform finite control-node fields must remain undefined");
            const auto& reference = native_nodes[(step - 1) * 4 + n];
            require_time(frame.time, reference.at("time"));
            if (reference.at("node") != static_cast<double>(n + 1))
                throw std::runtime_error("Native nonuniform finite node correspondence changed");
            temperature.add(t.at(radial_node[n]), reference.at("temperature"));
            const std::array<double, 2> actual{ur.at(radial_node[n]), uz.at(axial_node[n])};
            const std::array<double, 2> expected{reference.at("ur"),
                n < 2 ? constrained_zero(reference.at("uz"), 1e-13) : reference.at("uz")};
            displacement.add(actual.data(), expected.data(), actual.size());
            if (n < 2)
                prescribed_zero_displacement.add(actual[1], expected[1]);
        }
        const auto& bottom_inner = native_nodes[(step - 1) * 4];
        const auto& bottom_outer = native_nodes[(step - 1) * 4 + 1];
        const auto& top_outer = native_nodes[(step - 1) * 4 + 2];
        const auto& top_inner = native_nodes[(step - 1) * 4 + 3];
        const std::array<double, 2> native_heat = {bottom_inner.at("reaction_heat") + top_inner.at("reaction_heat"),
            bottom_outer.at("reaction_heat") + top_outer.at("reaction_heat")};
        for (std::size_t n = 0; n < 2; ++n)
            heat_reaction.add(heat.at(n), native_heat[n]);
        if (variable_thermal && step == 8) {
            const auto& old_heat = frames[step - 1].nodal("reaction_heat_flux");
            const std::size_t offset = (step - 2) * 4;
            source_node_reaction_change.add(old_heat.at(0) - heat.at(0),
                native_nodes[offset].at("reaction_heat") + native_nodes[offset + 3].at("reaction_heat")
                    - native_heat[0]);
            source_node_reaction_change.add(old_heat.at(1) - heat.at(1),
                native_nodes[offset + 1].at("reaction_heat") + native_nodes[offset + 2].at("reaction_heat")
                    - native_heat[1]);
        }
        radial_reaction.add(rr.at(0), bottom_inner.at("rf_r") + top_inner.at("rf_r"));
        radial_reaction.add(rr.at(1), bottom_outer.at("rf_r") + top_outer.at("rf_r"));
        const double lower_force = bottom_inner.at("rf_z") + bottom_outer.at("rf_z");
        const double upper_force = top_inner.at("rf_z") + top_outer.at("rf_z");
        axial_reaction.add(rz.at(2), lower_force);
        axial_reaction.add(rz.at(3), upper_force);
        // Compare the actual native end force directly, without inferring
        // material-point integration weights from deformed coordinates.
        axial_force.add(frame.element("axial_force").at(0), upper_force);
        const double reference_height = top_inner.at("z") - bottom_inner.at("z");
        if (!(reference_height > 0.0))
            throw std::runtime_error("Nonuniform finite reference slice height must be positive");
        axial_strain.add(frame.element("axial_strain").at(0),
            variable_thermal && step == 5
                ? constrained_zero((top_inner.at("uz") - bottom_inner.at("uz")) / reference_height, 1e-13)
                : (top_inner.at("uz") - bottom_inner.at("uz")) / reference_height);
        for (std::size_t n = 0; n < 2; ++n)
            midpoint_displacement.add(uz.at(n), 0.5 * (uz.at(2) + uz.at(3)));
        const double reference_inner_radius = bottom_inner.at("r");
        const double reference_outer_radius = bottom_outer.at("r");
        const double current_inner_radius = reference_inner_radius + bottom_inner.at("ur");
        const double current_outer_radius = reference_outer_radius + bottom_outer.at("ur");
        const double current_height = reference_height + top_inner.at("uz") - bottom_inner.at("uz");
        if (!(reference_outer_radius > reference_inner_radius && reference_inner_radius > 0.0
                && current_outer_radius > current_inner_radius && current_inner_radius > 0.0 && current_height > 0.0))
            throw std::runtime_error("Nonuniform finite native annulus geometry is invalid");
        const double current_volume =
            std::acos(-1.0)
            * (current_outer_radius * current_outer_radius - current_inner_radius * current_inner_radius)
            * current_height;
        const double reference_volume =
            std::acos(-1.0)
            * (reference_outer_radius * reference_outer_radius - reference_inner_radius * reference_inner_radius)
            * reference_height;
        const double active_volume = finite ? current_volume : reference_volume;
        // The native BF amplitude is averaged over each accepted interval.
        const double source = variable_thermal && step >= 8 ? (step == 8 ? 5e6 : 1e7) : 0.0;
        const double native_generated = source * active_volume;
        const double native_total_heat = native_heat[0] + native_heat[1];
        const double native_stored = native_total_heat + native_generated;
        generated_heat_rate.add(frame.global("conservation_generated_heat_rate"), native_generated);
        stored_heat_rate.add(frame.global("conservation_stored_heat_rate"),
            variable_thermal && step >= 7 ? constrained_zero(native_stored, 1e-10) : native_stored);
        dirichlet_heat_rate.add(frame.global("conservation_dirichlet_heat_input_rate"),
            variable_thermal && step == 7 ? constrained_zero(native_total_heat, 1e-10) : native_total_heat);
        thermal_balance.add(frame.global("conservation_global_thermal_balance"), 0.0);
        const auto& previous_ur = frames[step - 1].nodal("displacement_r");
        const double conduction_length =
            reference_outer_radius - reference_inner_radius
            + (finite ? 0.5 * (ur.at(1) - ur.at(0) + previous_ur.at(1) - previous_ur.at(0)) : 0.0);
        double native_elastic_energy = 0.0;
        for (std::size_t point = 0; point < 4; ++point) {
            const auto& reference = native_points[(step - 1) * 4 + point];
            require_time(frame.time, reference.at("time"));
            if (reference.at("element") != 1.0 || reference.at("point") != static_cast<double>(point + 1))
                throw std::runtime_error("Native nonuniform finite material-point correspondence changed");
            const std::size_t q = point % 2;
            // Reconstruct exported nodal values for diagnosis only: the actual
            // GPS body heat flux is not a production output field.
            const double conductivity = variable_thermal ? 10.0 + 0.1 * (t.at(q) - 600.0) : 10.0;
            const std::array<double, 2> diagnostic_flux = {-conductivity * (t.at(1) - t.at(0)) / conduction_length,
                0.0};
            const std::array<double, 2> native_flux = {reference.at("heat_flux_r"), reference.at("heat_flux_z")};
            reconstructed_heat_flux_diagnostic.add(diagnostic_flux.data(), native_flux.data(), 2);
            add_tensor(stress, frame, reference, "stress", q);
            add_tensor(elastic, frame, reference, "elastic", q);
            zero_stress_shear.add(frame.element("stress_rz_q" + std::to_string(q)).at(0),
                constrained_zero(reference.at("stress_rz"), 1e-4));
            zero_elastic_shear.add(frame.element("elastic_rz_q" + std::to_string(q)).at(0),
                constrained_zero(reference.at("elastic_rz"), 1e-13));
            // Reconstruct energy from native fields using the same mechanical
            // weights: each reference Gauss volume fraction times total current
            // volume. This is not a comparison with native SENER or point IVOL.
            const double eta = 0.5 * (1.0 + (q == 0 ? -1.0 : 1.0) / std::sqrt(3.0));
            const double reference_radius = (1.0 - eta) * reference_inner_radius + eta * reference_outer_radius;
            const double mechanical_measure =
                active_volume * reference_radius / (2.0 * (reference_inner_radius + reference_outer_radius));
            double stress_elastic_product = 0.0;
            for (const char* component : {"rr", "zz", "hoop", "rz"})
                stress_elastic_product += reference.at(std::string("stress_") + component)
                                          * reference.at(std::string("elastic_") + component)
                                          * (std::string(component) == "rz" ? 2.0 : 1.0);
            native_elastic_energy += 0.5 * stress_elastic_product * mechanical_measure;
        }
        elastic_energy_change.add(frame.global("conservation_elastic_energy_change"),
            variable_thermal && step >= 7
                ? constrained_zero(native_elastic_energy - previous_native_elastic_energy, 1e-12)
                : native_elastic_energy - previous_native_elastic_energy);
        previous_native_elastic_energy = native_elastic_energy;
        for (std::size_t q = 0; q < 2; ++q) {
            const std::string suffix = "_q" + std::to_string(q);
            for (const char* prefix : {"plastic_", "creep_"})
                for (const char* component : {"rr", "zz", "hoop", "rz"})
                    zero_mechanism.add(frame.element(std::string(prefix) + component + suffix).at(0), 0.0);
            zero_mechanism.add(frame.element("equiv_plastic" + suffix).at(0), 0.0);
            zero_mechanism.add(frame.element("equiv_creep" + suffix).at(0), 0.0);
        }
    }
    if (temperature.value_count != 40 || displacement.group_count != 40 || stress.group_count != 40
        || elastic.group_count != 40 || zero_stress_shear.value_count != 40 || zero_elastic_shear.value_count != 40
        || radial_reaction.value_count != 20 || axial_reaction.value_count != 20 || axial_force.value_count != 10
        || axial_strain.value_count != 10 || elastic_energy_change.value_count != 10
        || zero_mechanism.value_count != 200 || heat_reaction.value_count != 20 || stored_heat_rate.value_count != 10
        || generated_heat_rate.value_count != 10 || dirichlet_heat_rate.value_count != 10
        || thermal_balance.value_count != 10 || reconstructed_heat_flux_diagnostic.group_count != 40
        || (variable_thermal && source_node_reaction_change.value_count != 2))
        throw std::runtime_error("Nonuniform finite comparison lost native field or inactive-history samples");
    bool passed = true;
    passed = check("temperature", temperature, 1e-11) && passed;
    passed = check("displacement", displacement, 1e-13) && passed;
    passed = check("stress", stress, 1e-4) && passed;
    passed = check("elastic_strain", elastic, 1e-13) && passed;
    passed = check("radial_reaction", radial_reaction, 1e-8) && passed;
    passed = check("axial_reaction", axial_reaction, 1e-8) && passed;
    passed = check("axial_force", axial_force, 1e-8) && passed;
    passed = check("axial_strain", axial_strain, 1e-13) && passed;
    passed = check("mechanical_measure_elastic_energy_change", elastic_energy_change, 1e-12) && passed;
    passed = check("radial_node_axial_interpolation", midpoint_displacement, 1e-13) && passed;
    passed = check("prescribed_zero_displacement", prescribed_zero_displacement, 1e-13) && passed;
    passed = check("analytical_zero_stress_shear", zero_stress_shear, 1e-4) && passed;
    passed = check("analytical_zero_elastic_shear", zero_elastic_shear, 1e-13) && passed;
    passed = check("inactive_inelastic_history", zero_mechanism, 1e-13) && passed;
    passed = check("body_boundary_heat_reaction", heat_reaction, 1e-10) && passed;
    passed = check("stored_heat_rate", stored_heat_rate, 1e-10) && passed;
    passed = check("generated_heat_rate", generated_heat_rate, 1e-10) && passed;
    passed = check("dirichlet_heat_input_rate", dirichlet_heat_rate, 1e-10) && passed;
    passed = check("thermal_balance", thermal_balance, 1e-10) && passed;
    if (variable_thermal)
        passed = check("source_node_reaction_change", source_node_reaction_change, 1e-10) && passed;
    fuelsim::test::print_grouped_relative_metrics("reconstructed_heat_flux_diagnostic",
        reconstructed_heat_flux_diagnostic);
    std::cout << "native_node_samples=" << displacement.group_count
              << "\nnative_material_point_samples=" << stress.group_count
              << "\nqualification_scope=prescribed_nonuniform_" << (finite ? "finite" : "small")
              << "_thermomechanics_with_"
              << (temperature_dependent_elasticity ? "paired_endpoint_material_temperature" : "constant_elasticity")
              << "\nthermal_operator_cross_comparison=nodal_heat_reactions_and_actual_conservation_outputs\n"
              << "body_heat_flux_comparison=diagnostic_reconstruction_not_production_field_output\n"
              << "radial_gps_nonuniform_abaqus_qualification=" << (passed ? "passed" : "failed") << '\n';
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

// Native CAX4T chain has shared middle nodes and one axial displacement per section.
bool chain_small(const std::string& output, const std::string& summary, const std::filesystem::path& references) {
    require_summary(summary, 1.0);
    const auto frames = fuelsim::test::read_exodus_history(output);
    const auto nodes = read_rows(references / "gps_two_slice_chain_nodes.csv");
    const auto points = read_rows(references / "gps_two_slice_chain_points.csv");
    const auto contacts = read_contact(references / "gps_two_slice_chain_contact.csv");
    if (frames.size() != 11 || nodes.size() != 120 || points.size() != 160)
        throw std::runtime_error("Continuous chain requires all ten frames, twelve nodes and sixteen material points");
    constexpr std::array<std::size_t, 12> labels{1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 15, 16};
    constexpr std::array<std::size_t, 12> radial{0, 1, 1, 0, 2, 3, 3, 2, 5, 4, 7, 6};
    constexpr std::array<std::size_t, 12> axial{8, 8, 9, 9, 11, 11, 12, 12, 10, 10, 13, 13};
    constexpr std::array<std::size_t, 4> elements{0, 2, 1, 3};
    constexpr std::array<std::array<std::size_t, 2>, 2> secondary{{{2, 3}, {3, 11}}};
    constexpr std::array<std::array<std::size_t, 2>, 2> primary{{{5, 8}, {8, 16}}};
    const std::array<std::string, 2> pairs{"SECONDARY_ONE/PRIMARY_ONE", "SECONDARY_TWO/PRIMARY_TWO"};
    FieldErrorMetrics temperature, displacement, middle, reaction, axial_strain, axial_force;
    FieldErrorMetrics pressure, traction, gap, normal_force, tangent_force, balance, slip;
    GroupedFieldErrorMetrics stress, elastic;
    for (std::size_t step = 1; step < frames.size(); ++step) {
        const auto& f = frames[step];
        const double time = nodes[(step - 1) * 12].at("time");
        require_time(f.time, time);
        if (f.nodes.size() != 14 || f.element("material_point_count") != std::vector<double>(4, 2.0))
            throw std::runtime_error("Continuous chain production layout changed");
        std::map<std::size_t, Row> by_label;
        std::array<double, 6> native_reactions{};
        for (std::size_t n = 0; n < 12; ++n) {
            const auto& r = nodes[(step - 1) * 12 + n];
            require_time(f.time, r.at("time"));
            if (r.at("node") != static_cast<double>(labels[n]))
                throw std::runtime_error("Continuous chain native node mapping changed");
            by_label.emplace(labels[n], r);
            temperature.add(f.nodal("temperature").at(radial[n]), r.at("temperature"));
            const bool outer = axial[n] >= 11;
            displacement.add(f.nodal("displacement_r").at(radial[n]),
                outer ? constrained_zero(r.at("ur"), 1e-13) : r.at("ur"));
            const bool zero = axial[n] == 8 || axial[n] == 11 || axial[n] == 13 || (outer && step <= 5);
            const double expected = zero ? constrained_zero(r.at("uz"), 1e-13) : r.at("uz");
            displacement.add(f.nodal("displacement_z").at(axial[n]), expected);
            if (axial[n] == 9 || axial[n] == 12)
                middle.add(f.nodal("displacement_z").at(axial[n]), expected);
            native_reactions[axial[n] - 8] += r.at("rf_z");
        }
        double total_reaction = 0;
        for (std::size_t n = 0; n < 6; ++n) {
            const double expected = n == 1 || n == 4 || (n >= 3 && step <= 5)
                                        ? constrained_zero(native_reactions[n], 1e-8)
                                        : native_reactions[n];
            reaction.add(f.nodal("reaction_force_z").at(8 + n), expected);
            total_reaction += f.nodal("reaction_force_z").at(8 + n);
        }
        balance.add(total_reaction, 0.0);
        for (std::size_t e = 0; e < 4; ++e) {
            double force = 0;
            for (std::size_t q = 0; q < 4; ++q) {
                Row r = points[(step - 1) * 16 + e * 4 + q];
                require_time(f.time, r.at("time"));
                if (r.at("element") != static_cast<double>(e + 1) || r.at("point") != static_cast<double>(q + 1))
                    throw std::runtime_error("Continuous chain material point mapping changed");
                r["stress_rz"] = constrained_zero(r.at("stress_rz"), 1e-4);
                r["elastic_rz"] = constrained_zero(r.at("elastic_rz"), 1e-13);
                if (e % 2 == 1 && step <= 5)
                    for (const char* c : {"rr", "zz", "hoop"}) {
                        r[std::string("stress_") + c] = constrained_zero(r.at(std::string("stress_") + c), 1e-4);
                        r[std::string("elastic_") + c] = constrained_zero(r.at(std::string("elastic_") + c), 1e-13);
                    }
                add_tensor(stress, f, r, "stress", q % 2, elements[e]);
                add_tensor(elastic, f, r, "elastic", q % 2, elements[e]);
                force += r.at("stress_zz") * r.at("volume") / (e < 2 ? 0.01 : 0.02);
                axial_strain.add(f.element("axial_strain").at(elements[e]), r.at("elastic_zz"));
            }
            axial_force.add(f.element("axial_force").at(elements[e]), force);
        }
        double normal = 0, tangent = 0;
        for (std::size_t layer = 0; layer < 2; ++layer) {
            const auto native = [&](const std::string& field, std::size_t node, std::size_t component = 0) {
                return contacts.at({time, pairs[layer], field, node, component});
            };
            for (std::size_t q = 0; q < 2; ++q) {
                const auto node = secondary[layer][q];
                const std::string suffix = "_interface_q" + std::to_string(q);
                const auto zero_contact = [&](double v) {
                    return step <= 5 ? constrained_zero(v, 1e-4) : v;
                };
                pressure.add(f.element("contact_pressure" + suffix).at(layer), zero_contact(native("CPRESS", node)));
                traction.add(f.element("contact_tangential_traction" + suffix).at(layer),
                    zero_contact(native("CSHEAR1", node)));
                const double g = native("COPEN", node);
                gap.add(f.element("contact_gap" + suffix).at(layer), step == 5 ? constrained_zero(g, 1e-13) : g);
                normal -= native("CNORMF", node, 0);
                tangent -= native("CSHEARF", node, 1);
                // Compare total relative motion at the actual Gauss locations from native U.
                // Native CSLIP1 is a recovered contact-history field, not this kinematic quantity.
                const double eta = 0.5 * (1.0 + (q == 0 ? -1.0 : 1.0) / std::sqrt(3.0));
                double expected = 0;
                for (std::size_t end = 0; end < 2; ++end)
                    expected +=
                        (end == 0 ? 1.0 - eta : eta)
                        * (by_label.at(secondary[layer][end]).at("uz") - by_label.at(primary[layer][end]).at("uz"));
                const auto& uz = f.nodal("displacement_z");
                const double actual =
                    (1.0 - eta) * (uz.at(8 + layer) - uz.at(11 + layer)) + eta * (uz.at(9 + layer) - uz.at(12 + layer));
                slip.add(actual, expected);
            }
        }
        normal_force.add(f.global("contact_force_interface"), step <= 5 ? constrained_zero(normal, 1e-8) : normal);
        tangent_force.add(f.global("contact_tangential_force_interface"),
            step <= 5 ? constrained_zero(tangent, 1e-8) : tangent);
    }
    bool passed = true;
    for (const auto& entry : std::vector<std::tuple<std::string, const FieldErrorMetrics*, double>>{
             {"chain_temperature", &temperature, 1e-11},
             {"chain_displacement", &displacement, 1e-13},
             {"chain_free_middle_displacement", &middle, 1e-13},
             {"chain_axial_reaction", &reaction, 1e-8},
             {"chain_axial_strain", &axial_strain, 1e-13},
             {"chain_axial_force", &axial_force, 1e-8},
             {"chain_pressure", &pressure, 1e-4},
             {"chain_traction", &traction, 1e-4},
             {"chain_gap", &gap, 1e-13},
             {"chain_normal_force", &normal_force, 1e-8},
             {"chain_tangent_force", &tangent_force, 1e-8},
             {"chain_axial_balance", &balance, 1e-8},
             {"chain_relative_axial_motion", &slip, 1e-13}})
        passed = check(std::get<0>(entry), *std::get<1>(entry), std::get<2>(entry)) && passed;
    passed = check("chain_stress", stress, 1e-4) && passed;
    passed = check("chain_elastic_strain", elastic, 1e-13) && passed;
    return passed;
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
    FieldErrorMetrics radial_reaction, axial_control_reaction, axial_force;
    FieldErrorMetrics body_axial_reaction, production_axial_force_balance;
    FieldErrorMetrics prescribed_zero_displacement, midpoint_displacement;
    FieldErrorMetrics zero_stress_shear, zero_elastic_shear;
    GroupedFieldErrorMetrics displacement, stress, elastic;
    GroupedFieldErrorMetrics outer_stress, outer_elastic, inner_stress, inner_elastic;
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
        std::array<double, 8> reference_radial_reaction{};
        std::array<double, 6> reference_axial_reaction{};
        for (std::size_t n = 0; n < 16; ++n) {
            const auto& reference = native_nodes[(step - 1) * 16 + n];
            require_time(frame.time, reference.at("time"));
            if (reference.at("node") != static_cast<double>(n + 1))
                throw std::runtime_error("Native contact node mapping changed");
            temperature.add(t.at(radial_node[n]), reference.at("temperature"));
            reference_radial_reaction.at(radial_node[n]) += reference.at("rf_r");
            reference_axial_reaction.at(axial_node[n] - 8) += reference.at("rf_z");
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
        // A radial BAR2 node represents both native axial endpoints at the same
        // radius. Shared axial controls collect both adjoining native end faces.
        for (std::size_t n = 0; n < reference_radial_reaction.size(); ++n)
            radial_reaction.add(frame.nodal("reaction_force_r").at(n), reference_radial_reaction[n]);
        for (std::size_t n = 0; n < reference_axial_reaction.size(); ++n) {
            const bool unloaded_middle = (n == 1 || n == 4) && (step == 1 || step == 6 || step == 7 || step == 8);
            axial_control_reaction.add(frame.nodal("reaction_force_z").at(n + 8),
                unloaded_middle ? constrained_zero(reference_axial_reaction[n], 1e-8) : reference_axial_reaction[n]);
        }
        for (std::size_t n = 0; n < 8; ++n) {
            const std::size_t lower = 8 + n / 4 + ((n % 4) / 2) * 3;
            midpoint_displacement.add(uz.at(n), 0.5 * (uz.at(lower) + uz.at(lower + 1)));
        }
        std::array<double, 4> reference_axial_force{};
        for (std::size_t point = 0; point < 16; ++point) {
            const auto& reference = native_points[(step - 1) * 16 + point];
            require_time(frame.time, reference.at("time"));
            if (reference.at("element") != static_cast<double>(point / 4 + 1)
                || reference.at("point") != static_cast<double>(point % 4 + 1))
                throw std::runtime_error("Native contact material-point mapping changed");
            const std::size_t element = source_element[point / 4];
            add_tensor(stress, frame, reference, "stress", point % 2, element);
            add_tensor(elastic, frame, reference, "elastic", point % 2, element);
            if ((point / 4) % 2 == 1) {
                add_tensor(outer_stress, frame, reference, "stress", point % 2, element);
                add_tensor(outer_elastic, frame, reference, "elastic", point % 2, element);
            } else {
                add_tensor(inner_stress, frame, reference, "stress", point % 2, element);
                add_tensor(inner_elastic, frame, reference, "elastic", point % 2, element);
            }
            const std::string suffix = "_q" + std::to_string(point % 2);
            zero_stress_shear.add(frame.element("stress_rz" + suffix).at(element),
                constrained_zero(reference.at("stress_rz"), 1e-4));
            zero_elastic_shear.add(frame.element("elastic_rz" + suffix).at(element),
                constrained_zero(reference.at("elastic_rz"), 1e-13));
            const double height = point < 8 ? 0.01 : 0.02;
            reference_axial_force[element] += reference.at("stress_zz") * reference.at("volume") / height;
        }
        for (std::size_t element = 0; element < reference_axial_force.size(); ++element)
            axial_force.add(frame.element("axial_force").at(element), reference_axial_force[element]);
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
    if (outer_stress.group_count != 80 || outer_elastic.group_count != 80 || inner_stress.group_count != 80
        || inner_elastic.group_count != 80 || stress.group_count != 160 || elastic.group_count != 160
        || zero_stress_shear.value_count != 160 || zero_elastic_shear.value_count != 160)
        throw std::runtime_error("Contact body comparison must retain every native material point in all ten steps");
    if (radial_reaction.value_count != 80 || axial_control_reaction.value_count != 60 || axial_force.value_count != 40)
        throw std::runtime_error("Contact force comparison must retain every active node and section in all ten steps");
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
    passed = check("body_radial_reaction", radial_reaction, 1e-8) && passed;
    passed = check("body_axial_control_reaction", axial_control_reaction, 1e-8) && passed;
    passed = check("body_axial_section_force", axial_force, 1e-8) && passed;
    passed = check("production_body_axial_reaction", body_axial_reaction, 1e-8) && passed;
    passed = check("production_axial_force_conservation", production_axial_force_balance, 1e-8) && passed;
    passed = check("inactive_inelastic_history", zero_mechanism, 1e-13) && passed;
    passed = check("body_stress", stress, 1e-4) && passed;
    passed = check("body_elastic_strain", elastic, 1e-13) && passed;
    passed = check("inner_body_stress", inner_stress, 1e-4) && passed;
    passed = check("inner_body_elastic_strain", inner_elastic, 1e-13) && passed;
    passed = check("constrained_outer_body_stress", outer_stress, 1e-4) && passed;
    passed = check("constrained_outer_body_elastic_strain", outer_elastic, 1e-13) && passed;
    passed = check("body_zero_stress_shear", zero_stress_shear, 1e-4) && passed;
    passed = check("body_zero_elastic_shear", zero_elastic_shear, 1e-13) && passed;
    std::cout
        << "positive_sliding_samples=" << positive_sliding << "\nnegative_sliding_samples=" << negative_sliding
        << "\nsticking_samples=" << sticking << "\nopen_samples=" << opening
        << "\nqualification_scope=prescribed_motion_contact_and_thermal_operators_and_all_body_fields_and_reactions\n"
        << "native_body_material_point_samples=" << stress.group_count
        << "\nradial_reaction_samples=" << radial_reaction.value_count
        << "\naxial_control_reaction_samples=" << axial_control_reaction.value_count
        << "\naxial_section_force_samples=" << axial_force.value_count << '\n'
        << "radial_gps_contact_abaqus_qualification=" << (passed ? "passed" : "failed") << '\n';
    return passed;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: production_radial "
                     "<uniform-small|uniform-finite|nonuniform-finite|material-temperature-finite|thermal-small|"
                     "thermal-finite|contact-small|chain-"
                     "small|steady-finite> "
                     "<results.e> "
                     "<summary.csv> "
                     "<reference_dir>\n";
        return 2;
    }
    try {
        const std::string mode = argv[1];
        if (mode == "uniform-small" || mode == "uniform-finite")
            return uniform(argv[2], argv[3], argv[4], mode == "uniform-finite") ? 0 : 1;
        if (mode == "nonuniform-finite")
            return nonuniform(argv[2], argv[3], argv[4], NonuniformReference::finite_constant) ? 0 : 1;
        if (mode == "material-temperature-finite")
            return nonuniform(argv[2], argv[3], argv[4], NonuniformReference::finite_elastic_temperature) ? 0 : 1;
        if (mode == "thermal-small" || mode == "thermal-finite")
            return nonuniform(argv[2],
                       argv[3],
                       argv[4],
                       mode == "thermal-small" ? NonuniformReference::small_thermal
                                               : NonuniformReference::finite_thermal)
                       ? 0
                       : 1;
        if (mode == "chain-small")
            return chain_small(argv[2], argv[3], argv[4]) ? 0 : 1;
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
