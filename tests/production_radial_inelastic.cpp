#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using fuelsim::test::ExodusResults;
using fuelsim::test::FieldErrorMetrics;
using Tensor = std::array<double, 4>;
constexpr std::array<const char*, 4> components{"rr", "zz", "hoop", "rz"};
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double young_modulus = 1e9, poisson_ratio = 0.3;
constexpr double shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio));
constexpr double bulk_modulus = young_modulus / (3.0 * (1.0 - 2.0 * poisson_ratio));
constexpr double yield_stress = 1e6, hardening_modulus = 1e7;
constexpr double creep_coefficient = 0.005 / 1e6;
constexpr double time_step = 0.1, final_extension = 0.0004, height = 0.01;
constexpr double inner_radius = 0.004, outer_radius = 0.005;
constexpr double relative_tolerance = 0.001;
constexpr double strain_zero_tolerance = 1e-12, stress_zero_tolerance = 1e-4;

struct Reference final {
    double total_axial_strain = 0.0;
    double equivalent_stress = 0.0;
    double equivalent_plastic_strain = 0.0;
    double equivalent_creep_strain = 0.0;
    double plastic_increment = 0.0;
    double creep_increment = 0.0;
    Tensor stress{}, elastic{}, plastic{}, creep{};
};

// This proportional, diagonal deformation has no spin. With d denoting the
// axial Hughes-Winget increment, q_trial = q_old + 2 G d. Norton exponent one
// gives dc = dt b q; simultaneous plastic consistency gives q = Y_old + H dp.
// Substituting both into q_trial - q = 3 G (dp + dc) yields the rational
// expression below. It requires no production material functions or root solve.
Reference advance(const Reference& old, std::size_t step) {
    Reference next = old;
    const double stretch_old = 1.0 + final_extension / height * time_step * static_cast<double>(step - 1);
    const double stretch_new = 1.0 + final_extension / height * time_step * static_cast<double>(step);
    const double increment = 2.0 * (stretch_new - stretch_old) / (stretch_new + stretch_old);
    next.total_axial_strain += increment;
    const double trial = old.equivalent_stress + 2.0 * shear_modulus * increment;
    const double old_yield = yield_stress + hardening_modulus * old.equivalent_plastic_strain;
    const double relaxation = 1.0 + 3.0 * shear_modulus * time_step * creep_coefficient;
    next.plastic_increment = (trial - relaxation * old_yield) / (3.0 * shear_modulus + relaxation * hardening_modulus);
    next.equivalent_stress = old_yield + hardening_modulus * next.plastic_increment;
    next.creep_increment = time_step * creep_coefficient * next.equivalent_stress;
    if (!(next.plastic_increment > 0.0) || !(next.creep_increment > 0.0))
        throw std::runtime_error("The independent reference must activate plasticity and creep at every increment");
    next.equivalent_plastic_strain += next.plastic_increment;
    next.equivalent_creep_strain += next.creep_increment;
    constexpr Tensor direction{-0.5, 1.0, -0.5, 0.0};
    const double mean_stress = bulk_modulus * next.total_axial_strain;
    for (std::size_t component = 0; component < 4; ++component) {
        next.plastic[component] = direction[component] * next.equivalent_plastic_strain;
        next.creep[component] = direction[component] * next.equivalent_creep_strain;
        next.elastic[component] =
            (component == 1 ? next.total_axial_strain : 0.0) - next.plastic[component] - next.creep[component];
        next.stress[component] =
            component == 3 ? 0.0 : mean_stress + 2.0 / 3.0 * direction[component] * next.equivalent_stress;
    }
    return next;
}

double number(const std::string& value) {
    std::size_t count = 0;
    const double parsed = std::stod(value, &count);
    if (count != value.size() || !std::isfinite(parsed))
        throw std::runtime_error("Summary must contain finite numerical values");
    return parsed;
}

void require_summary(const std::string& path) {
    std::ifstream stream(path);
    std::string line;
    if (!stream || !std::getline(stream, line) || line != "metric,value")
        throw std::runtime_error("Missing or malformed production summary");
    std::map<std::string, std::string> values;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        const auto comma = line.find(',');
        if (comma == std::string::npos || line.find(',', comma + 1) != std::string::npos
            || !values.emplace(line.substr(0, comma), line.substr(comma + 1)).second)
            throw std::runtime_error("Malformed production summary row");
    }
    if ((values.at("completed") != "true" && values.at("completed") != "1")
        || number(values.at("accepted_steps")) != 10.0 || number(values.at("rejected_steps")) != 0.0
        || std::abs(number(values.at("committed_time")) - 1.0) > 1e-12)
        throw std::runtime_error("Finite inelastic qualification requires ten completed, unrejected time steps");
}

bool check(const std::string& name, const FieldErrorMetrics& metric, double zero_tolerance) {
    if (metric.value_count == 0)
        throw std::runtime_error("Empty material-history acceptance metric");
    if (metric.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, metric);
    else
        fuelsim::test::print_absolute_metrics(name, metric);
    return (!metric.has_relative_norm() || fuelsim::test::relative_metrics_below(metric, relative_tolerance))
           && metric.maximum_zero_reference_difference <= zero_tolerance;
}

Tensor point_tensor(const ExodusResults& frame, const std::string& prefix, std::size_t point) {
    Tensor tensor{};
    for (std::size_t component = 0; component < 4; ++component)
        tensor[component] = frame.element(prefix + "_" + components[component] + "_q" + std::to_string(point)).at(0);
    return tensor;
}

void add_tensor(std::map<std::string, FieldErrorMetrics>& metrics,
    const std::string& prefix,
    const Tensor& actual,
    const Tensor& expected) {
    for (std::size_t component = 0; component < 4; ++component)
        metrics[prefix + "_" + components[component]].add(actual[component], expected[component]);
}

bool qualify(const std::string& output, const std::string& summary) {
    require_summary(summary);
    const auto frames = fuelsim::test::read_exodus_history(output);
    if (frames.size() != 11)
        throw std::runtime_error("Finite inelastic qualification requires the initial frame and all ten increments");
    std::map<std::string, FieldErrorMetrics> metrics;
    Reference reference;
    std::array<double, 2> old_plastic{}, old_creep{};
    std::size_t simultaneous_point_increments = 0;
    double maximum_hw_absolute_difference = 0.0;
    constexpr std::array<std::array<double, 3>, 4> coordinates{
        {{{inner_radius, 0.005, 0.0}}, {{outer_radius, 0.005, 0.0}}, {{0.0, 0.0, 0.0}}, {{0.0, height, 0.0}}}};
    for (std::size_t step = 0; step < frames.size(); ++step) {
        const auto& frame = frames[step];
        const double time = time_step * static_cast<double>(step);
        if (!std::isfinite(frame.time) || std::abs(frame.time - time) > 1e-12
            || frame.nodes.size() != coordinates.size() || frame.block_names != std::vector<std::string>{"solid"}
            || frame.block_element_counts != std::vector<std::size_t>{1}
            || frame.element("material_point_count") != std::vector<double>{2.0})
            throw std::runtime_error("Finite inelastic history has unexpected times, mesh or material-point count");
        for (std::size_t node = 0; node < coordinates.size(); ++node)
            if (frame.nodes[node] != coordinates[node])
                throw std::runtime_error("Finite inelastic qualification requires the explicit gps_uniform.e geometry");
        if (step != 0)
            reference = advance(reference, step);
        const double extension = final_extension * time;
        for (std::size_t node = 0; node < coordinates.size(); ++node) {
            if (frame.nodal("node_role").at(node) != (node < 2 ? 1.0 : 2.0))
                throw std::runtime_error("Finite inelastic output changed the radial and axial node roles");
            const double expected_axial = node < 2 ? 0.5 * extension : (node == 2 ? 0.0 : extension);
            metrics["displacement_z"].add(frame.nodal("displacement_z").at(node), expected_axial);
            if (node < 2) {
                metrics["temperature"].add(frame.nodal("temperature").at(node), 600.0);
                metrics["displacement_r"].add(frame.nodal("displacement_r").at(node), 0.0);
                metrics["reaction_heat_flux"].add(frame.nodal("reaction_heat_flux").at(node), 0.0);
                if (!std::isnan(frame.nodal("reaction_force_z").at(node)))
                    throw std::runtime_error("Axial reactions are undefined on radial temperature nodes");
            } else if (!std::isnan(frame.nodal("temperature").at(node))
                       || !std::isnan(frame.nodal("displacement_r").at(node))
                       || !std::isnan(frame.nodal("reaction_force_r").at(node))
                       || !std::isnan(frame.nodal("reaction_heat_flux").at(node)))
                throw std::runtime_error("Unused axial control-node fields must remain NaN");
        }
        metrics["axial_strain"].add(frame.element("axial_strain").at(0), extension / height);
        const double area = pi * (outer_radius * outer_radius - inner_radius * inner_radius);
        const double force = area * reference.stress[1];
        metrics["axial_force"].add(frame.element("axial_force").at(0), force);
        metrics["reaction_force_z"].add(frame.nodal("reaction_force_z").at(2), -force);
        metrics["reaction_force_z"].add(frame.nodal("reaction_force_z").at(3), force);
        metrics["reaction_force_r"].add(frame.nodal("reaction_force_r").at(0),
            -2.0 * pi * inner_radius * (height + extension) * reference.stress[0]);
        metrics["reaction_force_r"].add(frame.nodal("reaction_force_r").at(1),
            2.0 * pi * outer_radius * (height + extension) * reference.stress[0]);
        for (std::size_t point = 0; point < 2; ++point) {
            const std::string suffix = "_q" + std::to_string(point);
            const Tensor stress = point_tensor(frame, "stress", point);
            const Tensor elastic = point_tensor(frame, "elastic", point);
            const Tensor plastic = point_tensor(frame, "plastic", point);
            const Tensor creep = point_tensor(frame, "creep", point);
            add_tensor(metrics, "stress", stress, reference.stress);
            add_tensor(metrics, "elastic", elastic, reference.elastic);
            add_tensor(metrics, "plastic", plastic, reference.plastic);
            add_tensor(metrics, "creep", creep, reference.creep);
            metrics["temperature"].add(frame.element("temperature" + suffix).at(0), 600.0);
            const double equivalent_plastic = frame.element("equiv_plastic" + suffix).at(0);
            const double equivalent_creep = frame.element("equiv_creep" + suffix).at(0);
            metrics["equiv_plastic"].add(equivalent_plastic, reference.equivalent_plastic_strain);
            metrics["equiv_creep"].add(equivalent_creep, reference.equivalent_creep_strain);
            if (step != 0) {
                const double dp = equivalent_plastic - old_plastic[point];
                const double dc = equivalent_creep - old_creep[point];
                if (!(dp > 0.0) || !(dc > 0.0))
                    throw std::runtime_error(
                        "Both inelastic mechanisms must grow at each identical point and increment");
                metrics["equiv_plastic_increment"].add(dp, reference.plastic_increment);
                metrics["equiv_creep_increment"].add(dc, reference.creep_increment);
                ++simultaneous_point_increments;
            }
            old_plastic[point] = equivalent_plastic;
            old_creep[point] = equivalent_creep;
            for (std::size_t component = 0; component < 4; ++component) {
                const double actual = elastic[component] + plastic[component] + creep[component];
                const double expected = component == 1 ? reference.total_axial_strain : 0.0;
                metrics[std::string("accumulated_hw_") + components[component]].add(actual, expected);
                maximum_hw_absolute_difference = std::max(maximum_hw_absolute_difference, std::abs(actual - expected));
            }
        }
    }
    // A 0.1% material-field gate alone cannot distinguish the discrete HW sum
    // from log(lambda_z): they differ by 5.03e-8 here. Check the exact discrete
    // accumulated strain separately, including its prescribed-zero components.
    bool accepted = simultaneous_point_increments == 20 && maximum_hw_absolute_difference <= strain_zero_tolerance;
    std::cout << "simultaneous_plastic_creep_point_increments=" << simultaneous_point_increments << '\n'
              << "accumulated_hw_maximum_absolute_difference=" << maximum_hw_absolute_difference << '\n';
    for (const auto& item : metrics) {
        double zero_tolerance = strain_zero_tolerance;
        if (item.first.rfind("stress_", 0) == 0)
            zero_tolerance = stress_zero_tolerance;
        else if (item.first == "axial_force" || item.first.rfind("reaction_force_", 0) == 0)
            zero_tolerance = 1e-8;
        else if (item.first == "reaction_heat_flux")
            zero_tolerance = 1e-10;
        accepted = check("gps_inelastic_finite_" + item.first, item.second, zero_tolerance) && accepted;
    }
    return accepted;
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3)
            throw std::invalid_argument("Usage: production_radial_inelastic <results.e> <summary.csv>");
        std::cout << std::scientific << std::setprecision(12);
        return qualify(argv[1], argv[2]) ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
