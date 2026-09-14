#include "support/production_hex8_full_field.hpp"
#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace fuelsim::test {
namespace {
struct CartesianPoint3 {
    double x = 0, y = 0, z = 0;
};

struct SymmetricTensor3Values {
    double xx = 0, yy = 0, zz = 0, xy = 0, yz = 0, xz = 0;
};

struct NodeReference final {
    std::size_t increment = 0, node = 0;
    double time = 0.0;
    std::array<double, 8> fields{};
};

struct IntegrationReference final {
    std::size_t increment = 0, element = 0, point = 0;
    double time = 0.0;
    CartesianPoint3 position{};
    double temperature = 0.0;
    std::array<double, 3> heat_flux{};
    SymmetricTensor3Values stress{}, logarithmic_strain{}, elastic_strain{}, plastic_strain{}, creep_strain{};
    double equivalent_plastic_strain = 0.0, equivalent_creep_strain = 0.0, integration_volume = 0.0;
};

struct ContactReference final {
    std::size_t increment = 0, node = 0, state = 0;
    double time = 0.0;
    CartesianPoint3 position{};
    double opening = 0.0, pressure = 0.0, slip_first = 0.0, slip_second = 0.0;
    std::array<double, 3> normal_force{}, shear_force{}, tangent_first{}, tangent_second{};
    double heat_flux = 0.0, shear_traction_first = 0.0, shear_traction_second = 0.0;
};

struct EnergyReference final {
    std::size_t increment = 0;
    double time = 0.0, internal = 0.0, elastic = 0.0, plastic = 0.0, creep = 0.0, friction = 0.0, external_work = 0.0,
           boundary_heat_rate = 0.0, artificial = 0.0;
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ','))
        result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size())
        throw std::invalid_argument("Incomplete Abaqus HEX8 full-field row in " + path);
    std::size_t parsed = 0;
    const double result = std::stod(values[index], &parsed);
    if (parsed != values[index].size() || !std::isfinite(result))
        throw std::invalid_argument("Invalid Abaqus HEX8 full-field number in " + path);
    return std::abs(result) < 1.0e-20 ? 0.0 : result;
}

std::size_t positive_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0)
        throw std::invalid_argument("Abaqus HEX8 full-field index is not a positive integer in " + path);
    return static_cast<std::size_t>(rounded);
}

SymmetricTensor3Values
tensor(const std::vector<std::string>& values, std::size_t start, const std::string& path, bool engineering_shear) {
    const double scale = engineering_shear ? 0.5 : 1.0;
    return {number(values, start, path),
        number(values, start + 1, path),
        number(values, start + 2, path),
        scale * number(values, start + 3, path),
        scale * number(values, start + 5, path),
        scale * number(values, start + 4, path)};
}

std::vector<NodeReference> read_nodes(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus HEX8 nodal reference: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "increment,time_s,node,temperature_k,u1_m,u2_m,u3_m,reaction_heat_flux_w,rf1_n,rf2_n,rf3_n")
        throw std::invalid_argument("Unexpected Abaqus HEX8 nodal header in " + path);
    std::vector<NodeReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 11)
            throw std::invalid_argument("Unexpected Abaqus HEX8 nodal columns in " + path);
        NodeReference reference;
        reference.increment = positive_integer(number(values, 0, path), path);
        reference.time = number(values, 1, path);
        reference.node = positive_integer(number(values, 2, path), path);
        for (std::size_t field = 0; field < reference.fields.size(); ++field)
            reference.fields[field] = number(values, field + 3, path);
        result.push_back(reference);
    }
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus HEX8 integration reference: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected =
        "increment,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,hfl1_w_m2,hfl2_w_m2,hfl3_w_m2,"
        "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa,le11,le22,le33,le12_engineering,le13_engineering,"
        "le23_engineering,ee11,ee22,ee33,ee12_engineering,ee13_engineering,ee23_engineering,pe11,pe22,pe33,"
        "pe12_engineering,pe13_engineering,pe23_engineering,peeq,ce11,ce22,ce33,ce12_engineering,"
        "ce13_engineering,ce23_engineering,ceeq,ivol_m3";
    if (line != expected)
        throw std::invalid_argument("Unexpected Abaqus HEX8 integration header in " + path);
    std::vector<IntegrationReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 44)
            throw std::invalid_argument("Unexpected Abaqus HEX8 integration columns in " + path);
        IntegrationReference reference;
        reference.increment = positive_integer(number(values, 0, path), path);
        reference.time = number(values, 1, path);
        reference.element = positive_integer(number(values, 2, path), path);
        reference.point = positive_integer(number(values, 3, path), path);
        reference.position = {number(values, 4, path), number(values, 5, path), number(values, 6, path)};
        reference.temperature = number(values, 7, path);
        reference.heat_flux = {number(values, 8, path), number(values, 9, path), number(values, 10, path)};
        reference.stress = tensor(values, 11, path, false);
        reference.logarithmic_strain = tensor(values, 17, path, true);
        reference.elastic_strain = tensor(values, 23, path, true);
        reference.plastic_strain = tensor(values, 29, path, true);
        reference.equivalent_plastic_strain = number(values, 35, path);
        reference.creep_strain = tensor(values, 36, path, true);
        reference.equivalent_creep_strain = number(values, 42, path);
        reference.integration_volume = number(values, 43, path);
        result.push_back(reference);
    }
    return result;
}

std::size_t nonnegative_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 0.0)
        throw std::invalid_argument("Abaqus HEX8 full-field index is not a nonnegative integer in " + path);
    return static_cast<std::size_t>(rounded);
}

std::vector<ContactReference> read_contact(const std::string& path, bool integrated_baseline = false) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus HEX8 contact reference: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected =
        integrated_baseline
            ? "increment,time_s,node,x_m,y_m,z_m,opening_m,pressure_pa,slip1_m,slip2_m,normal_force1_n,"
              "normal_force2_n,normal_force3_n,shear_force1_n,shear_force2_n,shear_force3_n,contact_heat_flux_w,state"
            : "increment,time_s,node,x_m,y_m,z_m,opening_m,pressure_pa,slip1_m,slip2_m,normal_force1_n,"
              "normal_force2_n,normal_force3_n,shear_force1_n,shear_force2_n,shear_force3_n,contact_heat_flux_w,"
              "shear_traction1_pa,shear_traction2_pa,tangent1_x,tangent1_y,tangent1_z,tangent2_x,tangent2_y,"
              "tangent2_z,state";
    if (line != expected)
        throw std::invalid_argument("Unexpected Abaqus HEX8 contact header in " + path);
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != (integrated_baseline ? 18U : 26U))
            throw std::invalid_argument("Unexpected Abaqus HEX8 contact columns in " + path);
        ContactReference reference;
        reference.increment = positive_integer(number(values, 0, path), path);
        reference.time = number(values, 1, path);
        reference.node = positive_integer(number(values, 2, path), path);
        reference.position = {number(values, 3, path), number(values, 4, path), number(values, 5, path)};
        reference.opening = number(values, 6, path);
        reference.pressure = number(values, 7, path);
        reference.slip_first = number(values, 8, path);
        reference.slip_second = number(values, 9, path);
        for (std::size_t component = 0; component < 3; ++component) {
            reference.normal_force[component] = number(values, 10 + component, path);
            reference.shear_force[component] = number(values, 13 + component, path);
            if (!integrated_baseline) {
                reference.tangent_first[component] = number(values, 19 + component, path);
                reference.tangent_second[component] = number(values, 22 + component, path);
            }
        }
        reference.heat_flux = number(values, 16, path);
        if (integrated_baseline) {
            reference.tangent_first = {0, 0, 1};
            reference.tangent_second = {0, -1, 0};
        } else {
            reference.shear_traction_first = number(values, 17, path);
            reference.shear_traction_second = number(values, 18, path);
        }
        reference.state = nonnegative_integer(number(values, integrated_baseline ? 17 : 25, path), path);
        result.push_back(reference);
    }
    return result;
}

void require_contact_free_reference(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus HEX8 contact reference: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected =
        "increment,time_s,node,x_m,y_m,z_m,opening_m,pressure_pa,slip1_m,slip2_m,normal_force1_n,"
        "normal_force2_n,normal_force3_n,shear_force1_n,shear_force2_n,shear_force3_n,contact_heat_flux_w,"
        "shear_traction1_pa,shear_traction2_pa,tangent1_x,tangent1_y,tangent1_z,tangent2_x,tangent2_y,"
        "tangent2_z,state";
    if (line != expected)
        throw std::invalid_argument("Unexpected Abaqus HEX8 contact header in " + path);
    while (std::getline(input, line))
        if (!line.empty())
            throw std::invalid_argument("Bulk result comparison cannot omit contact reference rows");
}

std::vector<EnergyReference> read_energy(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read Abaqus HEX8 energy reference: " + path);
    std::string line;
    std::getline(input, line);
    const std::string legacy_header =
        "increment,time_s,allie_j,allse_j,allpd_j,allcd_j,allfd_j,allwk_j,boundary_heat_rate_w";
    const bool has_artificial_energy = line == legacy_header + ",allae_j";
    if (line != legacy_header && !has_artificial_energy)
        throw std::invalid_argument("Unexpected Abaqus HEX8 energy header in " + path);
    std::vector<EnergyReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != (has_artificial_energy ? 10 : 9))
            throw std::invalid_argument("Unexpected Abaqus HEX8 energy columns in " + path);
        result.push_back({positive_integer(number(values, 0, path), path),
            number(values, 1, path),
            number(values, 2, path),
            number(values, 3, path),
            number(values, 4, path),
            number(values, 5, path),
            number(values, 6, path),
            number(values, 7, path),
            number(values, 8, path),
            has_artificial_energy ? number(values, 9, path) : 0.0});
    }
    return result;
}

std::array<double, 6> components(const SymmetricTensor3Values& value) {
    return {value.xx, value.yy, value.zz, value.xy, value.yz, value.xz};
}

using Matrix3 = std::array<std::array<double, 3>, 3>;

bool metrics_pass(const FieldErrorMetrics& metrics,
    double aggregate_tolerance,
    double pointwise_tolerance,
    double zero_tolerance,
    double qualified_pointwise_absolute_tolerance = 0.0,
    double peak_relative_tolerance = 0.0) {
    if (metrics.has_relative_norm()) {
        const bool aggregate_passed =
            metrics.relative_l2() < aggregate_tolerance
            && metrics.relative_absolute_peak()
                   < (peak_relative_tolerance > 0.0 ? peak_relative_tolerance : aggregate_tolerance);
        const double maximum_pointwise_absolute_difference =
            std::abs(metrics.maximum_pointwise_relative_actual - metrics.maximum_pointwise_relative_reference);
        const bool pointwise_passed =
            metrics.maximum_pointwise_relative < pointwise_tolerance
            || (qualified_pointwise_absolute_tolerance > 0.0
                && maximum_pointwise_absolute_difference < qualified_pointwise_absolute_tolerance);
        if (!aggregate_passed || !pointwise_passed)
            return false;
    }
    return metrics.maximum_zero_reference_difference < zero_tolerance;
}

bool grouped_metrics_pass(const GroupedFieldErrorMetrics& metrics,
    double aggregate_tolerance,
    double pointwise_tolerance,
    double zero_tolerance,
    double qualified_pointwise_absolute_tolerance = 0.0) {
    if (metrics.has_relative_norm()) {
        const double pointwise_absolute_difference =
            metrics.maximum_pointwise_relative * metrics.maximum_pointwise_reference_norm;
        const bool pointwise_passed = metrics.maximum_pointwise_relative < pointwise_tolerance
                                      || (qualified_pointwise_absolute_tolerance > 0.0
                                          && pointwise_absolute_difference < qualified_pointwise_absolute_tolerance);
        if (!(metrics.relative_l2() < aggregate_tolerance && metrics.relative_absolute_peak() < aggregate_tolerance
                && pointwise_passed))
            return false;
    }
    return metrics.maximum_zero_reference_difference < zero_tolerance;
}

bool report_metric(const std::string& name,
    const FieldErrorMetrics& metrics,
    double aggregate_tolerance,
    double pointwise_tolerance,
    double zero_tolerance,
    bool gate,
    double qualified_pointwise_absolute_tolerance = 0.0,
    double peak_relative_tolerance = 0.0) {
    if (metrics.has_relative_norm())
        print_relative_metrics(name, metrics);
    else
        print_absolute_metrics(name, metrics);
    if (!gate)
        return true;
    const bool passed = metrics_pass(metrics,
        aggregate_tolerance,
        pointwise_tolerance,
        zero_tolerance,
        qualified_pointwise_absolute_tolerance,
        peak_relative_tolerance);
    if (!passed)
        std::cerr << "[FAIL] " << name << " exceeds its full-field gate\n";
    return passed;
}

bool report_grouped(const std::string& name,
    const GroupedFieldErrorMetrics& metrics,
    double aggregate_tolerance,
    double pointwise_tolerance,
    double zero_tolerance,
    bool gate,
    double qualified_pointwise_absolute_tolerance = 0.0) {
    if (metrics.has_relative_norm())
        print_grouped_relative_metrics(name, metrics);
    else {
        std::cout << name << "_zero_reference_count=" << metrics.zero_reference_count << '\n'
                  << name << "_maximum_zero_reference_absolute_difference=" << metrics.maximum_zero_reference_difference
                  << '\n'
                  << name << "_maximum_difference_index=" << metrics.maximum_difference_index << '\n';
    }
    if (!gate)
        return true;
    const bool passed = grouped_metrics_pass(metrics,
        aggregate_tolerance,
        pointwise_tolerance,
        zero_tolerance,
        qualified_pointwise_absolute_tolerance);
    if (!passed)
        std::cerr << "[FAIL] " << name << " exceeds its complete-vector or complete-tensor gate\n";
    return passed;
}

bool compare_contact_history(const std::vector<ExodusResults>& frames, const ProductionHex8FullFieldOptions& options) {
    const bool integrated_baseline = options.case_name == "b523";
    const auto refs = read_contact(options.reference_prefix + "_contact.csv", integrated_baseline);
    std::array<FieldErrorMetrics, 3> baseline_resultants;
    const auto prefix = options.case_name + "_";
    const auto suffix = "_" + options.contact_name;
    if (refs.empty())
        throw std::invalid_argument("Contact reference must not be empty");
    std::array<GroupedFieldErrorMetrics, 8> vectors;
    FieldErrorMetrics pressure, normal_magnitude;
    std::set<std::pair<std::size_t, std::size_t>> seen;
    std::map<std::size_t, double> initial_faces;
    std::size_t state_count = 0, matched_states = 0, sticking_seen = 0, sliding_seen = 0;
    bool active_seen = false, reopened = false, recontact = false, crossed = false, final_active = false,
         final_sticking = false;
    double minimum_y = 0, maximum_y = 0, maximum_basis_error = 0;
    const std::array<std::string, 3> axes = {"x", "y", "z"};
    for (std::size_t step = 1; step < frames.size(); ++step) {
        const auto& frame = frames[step];
        std::array<double, 3> force{}, expected_force{}, moment{}, expected_moment{}, center{}, expected_center{};
        double weight = 0, expected_weight = 0, tangent_y = 0;
        std::size_t active = 0, sliding = 0, count = 0, output_count = 0;
        for (double flag : frame.nodal("contact_projected" + suffix))
            if (!std::isnan(flag))
                ++output_count;
        for (const auto& ref : refs) {
            if (ref.increment != step)
                continue;
            if (ref.node == 0 || ref.node > frame.nodes.size() || std::abs(ref.time - frame.time) > 1e-7
                || !seen.emplace(step, ref.node).second)
                throw std::invalid_argument("Contact reference association is invalid");
            ++count;
            const auto n = ref.node - 1;
            std::array<double, 3> p{}, ep = {ref.position.x, ref.position.y, ref.position.z}, nf{},
                                       enf = ref.normal_force, sf{}, esf = ref.shear_force, f{}, ef{}, slip{}, eslip{};
            for (std::size_t c = 0; c < 3; ++c) {
                p[c] = frame.nodes[n][c] + frame.nodal("displacement_" + axes[c]).at(n);
                nf[c] = -frame.nodal("contact_normal_force_" + axes[c] + suffix).at(n);
                sf[c] = -frame.nodal("contact_tangential_force_" + axes[c] + suffix).at(n);
                f[c] = nf[c] + sf[c];
                ef[c] = enf[c] + esf[c];
                slip[c] = frame.nodal("contact_total_slip_" + axes[c] + suffix).at(n);
                eslip[c] = ref.slip_first * ref.tangent_first[c] + ref.slip_second * ref.tangent_second[c];
                force[c] += f[c];
                expected_force[c] += ef[c];
            }
            vectors[0].add(p.data(), ep.data(), 3);
            if (ref.state != 0)
                vectors[1].add(slip.data(), eslip.data(), 3);
            vectors[2].add(nf.data(), enf.data(), 3);
            vectors[3].add(f.data(), ef.data(), 3);
            vectors[7].add(sf.data(), esf.data(), 3);
            const double w = std::hypot(nf[0], nf[1], nf[2]), ew = std::hypot(enf[0], enf[1], enf[2]);
            normal_magnitude.add(w, ew);
            pressure.add(frame.nodal("contact_pressure" + suffix).at(n), ref.pressure);
            weight += w;
            expected_weight += ew;
            for (std::size_t c = 0; c < 3; ++c) {
                const auto next = (c + 1) % 3, last = (c + 2) % 3;
                moment[c] += p[next] * f[last] - p[last] * f[next];
                expected_moment[c] += ep[next] * ef[last] - ep[last] * ef[next];
                center[c] += w * p[c];
                expected_center[c] += ew * ep[c];
            }
            double first_norm = 0, second_norm = 0, orthogonal = 0;
            for (std::size_t c = 0; c < 3; ++c) {
                first_norm += ref.tangent_first[c] * ref.tangent_first[c];
                second_norm += ref.tangent_second[c] * ref.tangent_second[c];
                orthogonal += ref.tangent_first[c] * ref.tangent_second[c];
            }
            maximum_basis_error = std::max(
                {maximum_basis_error, std::abs(first_norm - 1), std::abs(second_norm - 1), std::abs(orthogonal)});
            const std::size_t state = frame.nodal("contact_pressure" + suffix).at(n) <= 0  ? 0
                                      : frame.nodal("contact_sliding" + suffix).at(n) == 1 ? 2
                                                                                           : 1;
            ++state_count;
            if (state == ref.state)
                ++matched_states;
            if (state != 0) {
                ++active;
                if (state == 2)
                    ++sliding;
                tangent_y -= sf[1];
                const double face = frame.nodal("contact_primary_face" + suffix).at(n);
                if (!std::isfinite(face))
                    throw std::invalid_argument("Active contact has no primary face");
                if (!initial_faces.emplace(n, face).second && initial_faces.at(n) != face)
                    crossed = true;
            }
        }
        if (count == 0 || count != output_count)
            throw std::invalid_argument("Contact reference does not cover every output constraint");
        vectors[4].add(force.data(), expected_force.data(), 3);
        vectors[5].add(moment.data(), expected_moment.data(), 3);
        for (std::size_t c = 0; c < 3; ++c) {
            baseline_resultants[0].add(force[c], expected_force[c]);
            baseline_resultants[1].add(moment[c], expected_moment[c]);
        }
        if (weight > 0 && expected_weight > 0) {
            for (std::size_t c = 0; c < 3; ++c) {
                center[c] /= weight;
                expected_center[c] /= expected_weight;
            }
            vectors[6].add(center.data(), expected_center.data(), 3);
            for (std::size_t c = 0; c < 3; ++c)
                baseline_resultants[2].add(center[c], expected_center[c]);
        }
        if (active == 0) {
            if (active_seen)
                reopened = true;
        } else {
            if (reopened)
                recontact = true;
            active_seen = true;
        }
        sticking_seen += active - sliding;
        sliding_seen += sliding;
        minimum_y = std::min(minimum_y, tangent_y);
        maximum_y = std::max(maximum_y, tangent_y);
        final_active = active > 0;
        final_sticking = final_active && sliding == 0;
    }
    const double pointwise = options.contact_pointwise_relative_tolerance > 0
                                 ? options.contact_pointwise_relative_tolerance
                                 : options.contact_relative_tolerance;
    const std::array<std::string, 8> names = {"contact_position",
        "contact_slip_vector",
        "contact_normal_force_vector",
        "contact_complete_force_vector",
        "contact_resultant",
        "contact_moment",
        "contact_normal_force_center",
        "contact_shear_force_vector"};
    bool passed = true;
    for (std::size_t i = 0; i < vectors.size(); ++i)
        passed = report_grouped(prefix + names[i],
                     vectors[i],
                     options.contact_relative_tolerance,
                     pointwise,
                     i == 0 || i == 6 ? options.coordinate_tolerance
                     : i == 1         ? 1e-10
                                      : 1e-2,
                     !integrated_baseline && i != 7 && (i != 1 || options.gate_contact_slip),
                     i == 1 ? options.contact_slip_pointwise_absolute_tolerance : 0)
                 && passed;
    passed = report_metric(prefix + "contact_pressure",
                 pressure,
                 options.contact_relative_tolerance,
                 pointwise,
                 1,
                 options.gate_contact_pressure)
             && passed;
    passed = report_metric(prefix + "contact_normal_force_magnitude",
                 normal_magnitude,
                 options.contact_relative_tolerance,
                 pointwise,
                 1e-2,
                 !integrated_baseline)
             && passed;
    if (integrated_baseline)
        for (std::size_t i = 0; i < 3; ++i)
            passed = report_metric(prefix + names[4 + i] + "_components",
                         baseline_resultants[i],
                         5e-3,
                         5e-3,
                         i == 2 ? 1e-8 : 1,
                         true)
                     && passed;
    const double fraction = static_cast<double>(matched_states) / static_cast<double>(state_count);
    bool transition = true;
    if (options.contact_transition == "cycle")
        transition = active_seen && reopened && recontact && final_active;
    else if (options.contact_transition == "nonmatching")
        transition = active_seen && reopened && recontact && crossed && final_active;
    else if (options.contact_transition == "reversal")
        transition = sticking_seen > 0 && sliding_seen > 0 && minimum_y < 0 && maximum_y > 0 && final_sticking;
    std::cout << prefix << "contact_state_match_fraction=" << fraction << '\n'
              << prefix << "contact_transition_verified=" << transition << '\n'
              << prefix << "maximum_tangent_basis_error=" << maximum_basis_error << '\n';
    if (seen.size() != refs.size() || !final_active || !transition || maximum_basis_error >= 1e-9
        || (options.gate_contact_state && fraction < options.minimum_contact_state_match_fraction)) {
        std::cerr << "[FAIL] " << prefix << "contact coverage, state, transition, or reference basis check failed\n";
        passed = false;
    }
    if (options.case_name == "b527_medium" && !(std::abs(frames.back().global("contact_heat_rate" + suffix)) > 0)) {
        std::cerr << "[FAIL] B5.27 requires active thermal contact\n";
        passed = false;
    }
    return passed;
}

std::array<double, 6> tensor_output(const ExodusResults& frame, const std::string& name, std::size_t q, std::size_t e) {
    std::array<double, 6> result{};
    const std::array<std::string, 6> names = {"xx", "yy", "zz", "xy", "yz", "xz"};
    for (std::size_t c = 0; c < 6; ++c)
        result[c] = frame.element(name + names[c] + "_q" + std::to_string(q)).at(e);
    return result;
}
} // namespace

bool compare_production_hex8_full_field(const std::string& output_path, const ProductionHex8FullFieldOptions& options) {
    if (options.case_name.empty() || options.reference_prefix.empty() || options.expected_steps == 0
        || !(options.time_step > 0.0))
        throw std::invalid_argument("Abaqus HEX8 full-field options are incomplete");
    const double bulk_pointwise_tolerance = options.bulk_pointwise_relative_tolerance > 0.0
                                                ? options.bulk_pointwise_relative_tolerance
                                                : options.bulk_relative_tolerance,
                 energy_pointwise_tolerance = options.energy_pointwise_relative_tolerance > 0.0
                                                  ? options.energy_pointwise_relative_tolerance
                                                  : options.energy_relative_tolerance;
    const auto pointwise_tolerance = [](double configured, double fallback) {
        return configured > 0.0 ? configured : fallback;
    };
    const double displacement_pointwise_tolerance =
                     pointwise_tolerance(options.displacement_pointwise_relative_tolerance, bulk_pointwise_tolerance),
                 reaction_heat_flux_pointwise_tolerance =
                     pointwise_tolerance(options.reaction_heat_flux_pointwise_relative_tolerance,
                         bulk_pointwise_tolerance),
                 reaction_pointwise_tolerance =
                     pointwise_tolerance(options.reaction_pointwise_relative_tolerance, bulk_pointwise_tolerance),
                 stress_pointwise_tolerance =
                     pointwise_tolerance(options.stress_pointwise_relative_tolerance, bulk_pointwise_tolerance),
                 logarithmic_strain_pointwise_tolerance =
                     pointwise_tolerance(options.logarithmic_strain_pointwise_relative_tolerance,
                         bulk_pointwise_tolerance),
                 elastic_strain_pointwise_tolerance =
                     pointwise_tolerance(options.elastic_strain_pointwise_relative_tolerance, bulk_pointwise_tolerance),
                 inelastic_pointwise_tolerance =
                     pointwise_tolerance(options.inelastic_pointwise_relative_tolerance, bulk_pointwise_tolerance);
    if (!(options.bulk_relative_tolerance > 0.0) || !(options.energy_relative_tolerance > 0.0)
        || !(options.reaction_zero_absolute_tolerance > 0.0))
        throw std::invalid_argument("Abaqus HEX8 full-field tolerances are invalid");
    const std::vector<NodeReference> nodes = read_nodes(options.reference_prefix + "_nodal.csv");
    const std::vector<IntegrationReference> integration =
        read_integration(options.reference_prefix + "_integration.csv");
    if (options.contact_name.empty())
        require_contact_free_reference(options.reference_prefix + "_contact.csv");
    const std::vector<EnergyReference> energy = read_energy(options.reference_prefix + "_energy.csv");
    const std::size_t integration_points_per_element = options.reduced_integration ? 1 : 8;
    const auto frames = read_exodus_history(output_path);
    if (frames.size() != options.expected_steps + 1)
        throw std::invalid_argument("Unexpected production history frame count");
    const auto& final = frames.back();
    const std::size_t node_count = final.nodes.size();
    std::size_t element_count = 0;
    for (auto count : final.block_element_counts)
        element_count += count;
    if (nodes.size() != options.expected_steps * node_count
        || integration.size() != options.expected_steps * element_count * integration_points_per_element
        || energy.size() != options.expected_steps)
        throw std::invalid_argument("Reference row count does not match production history");
    for (std::size_t step = 1; step < frames.size(); ++step)
        if (std::abs(frames[step].time - options.time_step * static_cast<double>(step)) > 1e-7)
            throw std::invalid_argument("Unexpected production time");
    const std::string prefix = options.case_name + "_";
    bool passed = true;
    std::array<FieldErrorMetrics, 8> nodal_metrics;
    GroupedFieldErrorMetrics displacement_vector, reaction_force_vector;
    const std::array<std::string, 8> output_names = {"temperature",
        "displacement_x",
        "displacement_y",
        "displacement_z",
        "reaction_heat_flux",
        "reaction_force_x",
        "reaction_force_y",
        "reaction_force_z"};
    const std::array<std::string, 4> constraints = {"dirichlet_temperature",
        "dirichlet_displacement_x",
        "dirichlet_displacement_y",
        "dirichlet_displacement_z"};
    std::set<std::pair<std::size_t, std::size_t>> mapped_nodes;
    for (const auto& reference : nodes) {
        if (reference.increment < 1 || reference.increment >= frames.size() || reference.node < 1
            || reference.node > node_count || !mapped_nodes.emplace(reference.increment, reference.node).second)
            throw std::invalid_argument("Invalid or repeated reference node");
        const auto& frame = frames.at(reference.increment);
        if (std::abs(frame.time - reference.time) > 1e-7)
            throw std::invalid_argument("Reference nodal time differs");
        std::array<double, 3> actual_displacement{}, expected_displacement{}, actual_reaction{}, expected_reaction{};
        for (std::size_t f = 0; f < 8; ++f) {
            double value = frame.nodal(output_names[f]).at(reference.node - 1);
            if (f >= 4 && frame.nodal(constraints[f - 4]).at(reference.node - 1) == 0)
                value = 0;
            nodal_metrics[f].add(value, reference.fields[f]);
            if (f == 4 && options.gate_reaction_heat_flux
                && options.reaction_heat_flux_pointwise_absolute_tolerance > 0.0 && reference.fields[f] != 0.0) {
                // Apply the original relative-or-absolute qualification to every nonzero sample.
                const double absolute_difference = std::abs(value - reference.fields[f]);
                const double relative_difference = absolute_difference / std::abs(reference.fields[f]);
                if (!(relative_difference < reaction_heat_flux_pointwise_tolerance
                        || absolute_difference < options.reaction_heat_flux_pointwise_absolute_tolerance)) {
                    std::cerr << "[FAIL] " << prefix << "reaction_heat_flux sample time_s=" << reference.time
                              << " node=" << reference.node << " relative_error=" << relative_difference
                              << " absolute_error_W=" << absolute_difference
                              << " relative_tolerance=" << reaction_heat_flux_pointwise_tolerance
                              << " absolute_tolerance_W=" << options.reaction_heat_flux_pointwise_absolute_tolerance
                              << '\n';
                    passed = false;
                }
            }
            if (f >= 1 && f < 4) {
                actual_displacement[f - 1] = value;
                expected_displacement[f - 1] = reference.fields[f];
            }
            if (f >= 5) {
                actual_reaction[f - 5] = value;
                expected_reaction[f - 5] = reference.fields[f];
            }
        }
        displacement_vector.add(actual_displacement.data(), expected_displacement.data(), 3);
        reaction_force_vector.add(actual_reaction.data(), expected_reaction.data(), 3);
    }
    const std::array<std::string, 8> nodal_names = {"temperature",
        "displacement_x",
        "displacement_y",
        "displacement_z",
        "reaction_heat_flux",
        "reaction_force_x",
        "reaction_force_y",
        "reaction_force_z"};
    for (std::size_t field = 0; field < nodal_metrics.size(); ++field)
        passed = report_metric(prefix + nodal_names[field],
                     nodal_metrics[field],
                     options.bulk_relative_tolerance,
                     field == 4 ? reaction_heat_flux_pointwise_tolerance : bulk_pointwise_tolerance,
                     field == 0   ? 1.0e-8
                     : field < 4  ? 1.0e-10
                     : field == 4 ? 1.0e-2
                                  : 1.0,
                     field == 0 || (field == 4 && options.gate_reaction_heat_flux),
                     field == 4 ? options.reaction_heat_flux_pointwise_absolute_tolerance : 0.0,
                     field == 4 ? options.reaction_heat_flux_peak_relative_tolerance : 0.0)
                 && passed;
    passed = report_grouped(prefix + "displacement_vector",
                 displacement_vector,
                 options.bulk_relative_tolerance,
                 displacement_pointwise_tolerance,
                 1.0e-10,
                 true,
                 options.displacement_pointwise_absolute_tolerance)
             && passed;
    passed = report_grouped(prefix + "reaction_force_vector",
                 reaction_force_vector,
                 options.bulk_relative_tolerance,
                 reaction_pointwise_tolerance,
                 options.reaction_zero_absolute_tolerance,
                 true,
                 options.reaction_pointwise_absolute_tolerance)
             && passed;
    std::array<FieldErrorMetrics, 37> integration_metrics;
    GroupedFieldErrorMetrics integration_position, heat_flux_vector, stress_tensor, logarithmic_strain_tensor,
        elastic_strain_tensor, plastic_strain_tensor, creep_strain_tensor;
    double maximum_integration_coordinate_difference = 0;
    std::set<std::tuple<std::size_t, std::size_t, std::size_t>> mapped_integration_points;
    for (const auto& reference : integration) {
        if (reference.increment < 1 || reference.increment >= frames.size() || reference.element < 1
            || reference.element > element_count)
            throw std::invalid_argument("Invalid integration reference");
        const auto& frame = frames.at(reference.increment);
        if (std::abs(frame.time - reference.time) > 1e-7)
            throw std::invalid_argument("Reference integration time differs");
        const auto e = reference.element - 1;
        std::size_t closest = 0;
        double distance = std::numeric_limits<double>::infinity();
        std::array<double, 3> closest_position{};
        for (std::size_t q = 0; q < integration_points_per_element; ++q) {
            const auto suffix = "_q" + std::to_string(q);
            const std::array<double, 3> p = {frame.element("current_x" + suffix).at(e),
                frame.element("current_y" + suffix).at(e),
                frame.element("current_z" + suffix).at(e)};
            const double d =
                std::hypot(p[0] - reference.position.x, p[1] - reference.position.y, p[2] - reference.position.z);
            if (d < distance) {
                distance = d;
                closest = q;
                closest_position = p;
            }
        }
        if (!std::isfinite(distance) || !mapped_integration_points.emplace(reference.increment, e, closest).second)
            throw std::invalid_argument("Integration association is not unique");
        maximum_integration_coordinate_difference = std::max(maximum_integration_coordinate_difference, distance);
        const auto suffix = "_q" + std::to_string(closest);
        const std::array<double, 3> actual_heat_flux = {frame.element("heat_flux_x" + suffix).at(e),
            frame.element("heat_flux_y" + suffix).at(e),
            frame.element("heat_flux_z" + suffix).at(e)};
        const auto actual_stress = tensor_output(frame, "stress_", closest, e);
        const auto actual_logarithmic = tensor_output(frame, "logarithmic_strain_", closest, e);
        const auto elastic = tensor_output(frame, "elastic_", closest, e),
                   plastic = tensor_output(frame, "plastic_", closest, e),
                   creep = tensor_output(frame, "creep_", closest, e);
        const auto expected_stress = components(reference.stress),
                   expected_logarithmic = components(reference.logarithmic_strain),
                   expected_elastic = components(reference.elastic_strain),
                   expected_plastic = components(reference.plastic_strain),
                   expected_creep = components(reference.creep_strain);
        const std::array<double, 3> expected_position = {reference.position.x,
            reference.position.y,
            reference.position.z};
        integration_position.add(closest_position.data(), expected_position.data(), 3);
        for (std::size_t c = 0; c < 3; ++c)
            integration_metrics[c].add(actual_heat_flux[c], reference.heat_flux[c]);
        for (std::size_t c = 0; c < 6; ++c) {
            integration_metrics[3 + c].add(actual_stress[c], expected_stress[c]);
            integration_metrics[9 + c].add(actual_logarithmic[c], expected_logarithmic[c]);
            integration_metrics[15 + c].add(elastic[c], expected_elastic[c]);
            integration_metrics[21 + c].add(plastic[c], expected_plastic[c]);
            integration_metrics[28 + c].add(creep[c], expected_creep[c]);
        }
        heat_flux_vector.add(actual_heat_flux.data(), reference.heat_flux.data(), 3);
        stress_tensor.add(actual_stress.data(), expected_stress.data(), 6);
        logarithmic_strain_tensor.add(actual_logarithmic.data(), expected_logarithmic.data(), 6);
        elastic_strain_tensor.add(elastic.data(), expected_elastic.data(), 6);
        plastic_strain_tensor.add(plastic.data(), expected_plastic.data(), 6);
        creep_strain_tensor.add(creep.data(), expected_creep.data(), 6);
        integration_metrics[27].add(frame.element("equiv_plastic" + suffix).at(e), reference.equivalent_plastic_strain);
        integration_metrics[34].add(frame.element("equiv_creep" + suffix).at(e), reference.equivalent_creep_strain);
        integration_metrics[35].add(frame.element("material_temperature" + suffix).at(e), reference.temperature);
        integration_metrics[36].add(frame.element("integration_measure" + suffix).at(e), reference.integration_volume);
    }
    const std::array<std::string, 37> integration_names = {"heat_flux_x",
        "heat_flux_y",
        "heat_flux_z",
        "stress_xx",
        "stress_yy",
        "stress_zz",
        "stress_xy",
        "stress_yz",
        "stress_xz",
        "logarithmic_strain_xx",
        "logarithmic_strain_yy",
        "logarithmic_strain_zz",
        "logarithmic_strain_xy",
        "logarithmic_strain_yz",
        "logarithmic_strain_xz",
        "elastic_strain_xx",
        "elastic_strain_yy",
        "elastic_strain_zz",
        "elastic_strain_xy",
        "elastic_strain_yz",
        "elastic_strain_xz",
        "plastic_strain_xx",
        "plastic_strain_yy",
        "plastic_strain_zz",
        "plastic_strain_xy",
        "plastic_strain_yz",
        "plastic_strain_xz",
        "equivalent_plastic_strain",
        "creep_strain_xx",
        "creep_strain_yy",
        "creep_strain_zz",
        "creep_strain_xy",
        "creep_strain_yz",
        "creep_strain_xz",
        "equivalent_creep_strain",
        "material_temperature",
        "integration_volume"};
    for (std::size_t field = 0; field < integration_metrics.size(); ++field)
        report_metric(prefix + integration_names[field],
            integration_metrics[field],
            options.bulk_relative_tolerance,
            bulk_pointwise_tolerance,
            field < 3    ? 1.0e-6
            : field < 9  ? 1.0
            : field < 35 ? 1.0e-12
                         : 1.0e-10,
            false);
    passed = report_grouped(prefix + "integration_position",
                 integration_position,
                 options.bulk_relative_tolerance,
                 bulk_pointwise_tolerance,
                 options.coordinate_tolerance,
                 true)
             && passed;
    passed = report_grouped(prefix + "stress_tensor",
                 stress_tensor,
                 options.bulk_relative_tolerance,
                 stress_pointwise_tolerance,
                 options.case_name == "b523" ? 1e-12 : 1.0,
                 true,
                 options.stress_pointwise_absolute_tolerance)
             && passed;
    passed = report_grouped(prefix + "logarithmic_strain_tensor",
                 logarithmic_strain_tensor,
                 options.bulk_relative_tolerance,
                 logarithmic_strain_pointwise_tolerance,
                 1.0e-12,
                 true,
                 options.logarithmic_strain_pointwise_absolute_tolerance)
             && passed;
    passed = report_grouped(prefix + "elastic_strain_tensor",
                 elastic_strain_tensor,
                 options.bulk_relative_tolerance,
                 elastic_strain_pointwise_tolerance,
                 1.0e-12,
                 true,
                 options.elastic_strain_pointwise_absolute_tolerance)
             && passed;
    passed = report_grouped(prefix + "plastic_strain_tensor",
                 plastic_strain_tensor,
                 options.bulk_relative_tolerance,
                 inelastic_pointwise_tolerance,
                 1.0e-12,
                 true)
             && passed;
    passed = report_grouped(prefix + "creep_strain_tensor",
                 creep_strain_tensor,
                 options.bulk_relative_tolerance,
                 inelastic_pointwise_tolerance,
                 1.0e-12,
                 true)
             && passed;
    passed = report_metric(prefix + "equivalent_plastic_strain",
                 integration_metrics[27],
                 options.bulk_relative_tolerance,
                 inelastic_pointwise_tolerance,
                 1.0e-12,
                 true)
             && passed;
    passed = report_metric(prefix + "equivalent_creep_strain",
                 integration_metrics[34],
                 options.bulk_relative_tolerance,
                 inelastic_pointwise_tolerance,
                 1.0e-12,
                 true)
             && passed;
    passed = report_metric(prefix + "material_temperature",
                 integration_metrics[35],
                 options.bulk_relative_tolerance,
                 bulk_pointwise_tolerance,
                 options.case_name == "b523" ? 1e-12 : 1.0e-8,
                 true)
             && passed;
    passed = report_metric(prefix + "integration_volume",
                 integration_metrics[36],
                 options.bulk_relative_tolerance,
                 bulk_pointwise_tolerance,
                 1.0e-15,
                 true)
             && passed;
    report_grouped(prefix + "heat_flux_vector",
        heat_flux_vector,
        options.bulk_relative_tolerance,
        bulk_pointwise_tolerance,
        1.0e-6,
        false);
    std::cout << prefix << "maximum_integration_coordinate_difference=" << maximum_integration_coordinate_difference
              << '\n';
    if (maximum_integration_coordinate_difference >= (options.case_name == "b523" ? 1e-4 : 1e-3)) {
        std::cerr << "[FAIL] " << options.case_name << " integration-point coordinates do not map uniquely\n";
        passed = false;
    }

    std::array<FieldErrorMetrics, 8> energy_metrics;
    FieldErrorMetrics stored_heat;
    double cumulative_stored_heat = 0, cumulative_reference_heat = 0;
    double maximum_abaqus_artificial_energy = 0.0, maximum_abaqus_artificial_energy_fraction = 0.0;
    double cumulative_elastic = 0.0, cumulative_plastic = 0.0, cumulative_creep = 0.0, cumulative_friction = 0.0,
           cumulative_external_work = 0.0;
    for (std::size_t increment = 1; increment <= options.expected_steps; ++increment) {
        const auto& snapshot = frames.at(increment);
        const EnergyReference& expected = energy.at(increment - 1);
        if (expected.increment != increment || std::abs(expected.time - snapshot.time) > 1.0e-7)
            throw std::invalid_argument(options.case_name + " Abaqus energy time is invalid");
        cumulative_elastic += snapshot.global("conservation_elastic_energy_change");
        cumulative_plastic += snapshot.global("conservation_plastic_dissipation_increment");
        cumulative_creep += snapshot.global("conservation_creep_dissipation_increment");
        cumulative_friction += snapshot.global("conservation_friction_dissipation_increment");
        cumulative_external_work += snapshot.global("conservation_trapezoidal_pressure_traction_work_increment")
                                    + snapshot.global("conservation_trapezoidal_dirichlet_reaction_work_increment")
                                    + snapshot.global("conservation_body_force_work_increment");
        maximum_abaqus_artificial_energy = std::max(maximum_abaqus_artificial_energy, std::abs(expected.artificial));
        if (expected.internal != 0.0)
            maximum_abaqus_artificial_energy_fraction =
                std::max(maximum_abaqus_artificial_energy_fraction, std::abs(expected.artificial / expected.internal));
        energy_metrics[0].add(cumulative_elastic + cumulative_plastic + cumulative_creep,
            expected.internal - expected.artificial);
        energy_metrics[1].add(cumulative_elastic, expected.elastic);
        energy_metrics[2].add(cumulative_plastic, expected.plastic);
        energy_metrics[3].add(cumulative_creep, expected.creep);
        energy_metrics[4].add(cumulative_friction, expected.friction);
        energy_metrics[5].add(cumulative_external_work, expected.external_work);
        energy_metrics[6].add(snapshot.global("conservation_dirichlet_heat_input_rate"), expected.boundary_heat_rate);
        energy_metrics[7].add(snapshot.global("conservation_mechanical_hourglass_energy"), expected.artificial);
        if (options.case_name == "b523") {
            cumulative_stored_heat += snapshot.global("conservation_stored_heat_rate") * options.time_step;
            cumulative_reference_heat += expected.boundary_heat_rate * options.time_step;
            stored_heat.add(cumulative_stored_heat, cumulative_reference_heat);
        }
    }
    const std::array<std::string, 8> energy_names = {"internal_energy",
        "elastic_energy",
        "plastic_dissipation",
        "creep_dissipation",
        "friction_dissipation",
        "external_work",
        "boundary_heat_rate",
        "mechanical_hourglass_energy"};
    for (std::size_t field = 0; field < energy_metrics.size(); ++field) {
        const bool comparable = field != 4;
        passed = report_metric(prefix + energy_names[field],
                     energy_metrics[field],
                     options.energy_relative_tolerance,
                     energy_pointwise_tolerance,
                     field == 6 && options.case_name != "b523" ? 1.0e-2 : 1.0e-8,
                     comparable,
                     field == 5   ? options.external_work_pointwise_absolute_tolerance
                     : field == 7 ? options.hourglass_energy_pointwise_absolute_tolerance
                                  : 0.0)
                 && passed;
    }
    std::cout << prefix << "abaqus_artificial_energy_maximum_absolute=" << maximum_abaqus_artificial_energy << '\n'
              << prefix
              << "abaqus_artificial_energy_maximum_internal_fraction=" << maximum_abaqus_artificial_energy_fraction
              << '\n';
    std::cout << prefix << "compared_nodal_rows=" << nodes.size() << '\n'
              << prefix << "compared_integration_rows=" << integration.size() << '\n'
              << prefix << "compared_energy_rows=" << energy.size() << '\n';
    if (!options.contact_name.empty())
        passed = compare_contact_history(frames, options) && passed;
    if (options.case_name == "b523") {
        passed = report_metric(prefix + "stored_heat", stored_heat, 1e-2, 1e-2, 1e-8, true) && passed;
        const auto& final = frames.back();
        const auto active = std::count_if(final.nodal("contact_pressure_coupled_contact").begin(),
            final.nodal("contact_pressure_coupled_contact").end(),
            [](double p) { return p > 0; });
        const auto& plastic = final.element("equiv_plastic_q0");
        const auto& creep = final.element("equiv_creep_q0");
        if (active != 4 || !(*std::max_element(plastic.begin(), plastic.end()) > 0)
            || !(*std::max_element(creep.begin(), creep.end()) > 0)
            || !(cumulative_friction > 0 && energy.back().friction > 0)) {
            std::cerr << "[FAIL] B5.23 must activate both inelastic mechanisms, friction, and all four contact nodes\n";
            passed = false;
        }
    }
    if (passed)
        std::cout << "[PASS] " << options.case_name << " Abaqus full-field comparison\n";
    return passed;
}
} // namespace fuelsim::test
