#include "support/exodus_result_reader.hpp"
#include "support/field_error_metrics.hpp"
#include "support/production_checks.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr double comparison_tolerance = 5.0e-3;
// Individual tensor components cross zero during reversal.  Keep their
// componentwise pointwise errors as diagnostics, while acceptance uses the
// complete axisymmetric tensor field with the rz component counted twice in
// the Frobenius norm.
constexpr double time_tolerance = 1.0e-12;
enum class ExpectedBehavior {
    elastic,
    plastic,
    creep,
    coupled,
};

struct VariantConfig final {
    std::string name;
    ExpectedBehavior behavior = ExpectedBehavior::elastic;
    std::size_t element_count = 0;
    std::size_t step_count = 0;
    double comparison_tolerance = 0.0;
    double plastic_trace_tolerance = 0.0;
    double creep_trace_tolerance = 0.0;
};

VariantConfig variant_config(const std::string& name) {
    if (name == "production")
        return {name, ExpectedBehavior::coupled, 4, 100, comparison_tolerance, 7.0e-6, 2.0e-7};
    if (name == "elastic_displacement")
        return {name, ExpectedBehavior::elastic, 4, 100, comparison_tolerance, 1.0e-14, 1.0e-14};
    if (name == "plastic_displacement")
        return {name, ExpectedBehavior::plastic, 4, 100, comparison_tolerance, 7.0e-6, 1.0e-14};
    if (name == "creep_displacement")
        return {name, ExpectedBehavior::creep, 4, 100, comparison_tolerance, 1.0e-14, 1.3e-5};
    if (name == "coupled_displacement")
        return {name, ExpectedBehavior::coupled, 4, 100, comparison_tolerance, 7.0e-6, 2.0e-7};
    if (name == "coupled_pressure")
        return {name, ExpectedBehavior::coupled, 4, 100, comparison_tolerance, 7.0e-6, 2.0e-7};
    if (name == "material_oracle")
        return {name, ExpectedBehavior::coupled, 1, 100, 5.0e-6, 2.0e-6, 1.0e-7};
    throw std::invalid_argument("Unknown M4.3 comparison variant: " + name);
}

bool plastic_active(ExpectedBehavior behavior) noexcept {
    return behavior == ExpectedBehavior::plastic || behavior == ExpectedBehavior::coupled;
}

bool creep_active(ExpectedBehavior behavior) noexcept {
    return behavior == ExpectedBehavior::creep || behavior == ExpectedBehavior::coupled;
}

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    for (;;) {
        const std::size_t end = line.find(',', begin);
        fields.push_back(line.substr(begin, end - begin));
        if (end == std::string::npos)
            return fields;
        begin = end + 1;
    }
}

std::size_t column_index(const std::vector<std::string>& header, const std::string& name) {
    const auto found = std::find(header.begin(), header.end(), name);
    if (found == header.end())
        throw std::invalid_argument("M4.3 CSV is missing column '" + name + "'");
    return static_cast<std::size_t>(found - header.begin());
}

double csv_value(const std::vector<std::string>& fields, std::size_t column, const std::string& path) {
    if (column >= fields.size())
        throw std::invalid_argument("M4.3 CSV row is incomplete: " + path);
    std::size_t parsed = 0;
    const double value = std::stod(fields[column], &parsed);
    if (parsed != fields[column].size() || !std::isfinite(value))
        throw std::invalid_argument("M4.3 CSV value is invalid: " + path);
    return value;
}

std::size_t csv_id(const std::vector<std::string>& fields, std::size_t column, const std::string& path) {
    const double value = csv_value(fields, column, path);
    if (value < 0.0 || value > static_cast<double>(std::numeric_limits<std::size_t>::max())
        || std::floor(value) != value)
        throw std::invalid_argument("M4.3 CSV element ID is invalid: " + path);
    return static_cast<std::size_t>(value);
}

struct StressValues final {
    double rr = 0.0, zz = 0.0, hoop = 0.0, rz = 0.0;
};

struct MaterialValues final {
    std::array<double, 4> elastic_strain{}, plastic_strain{}, creep_strain{};
    double equivalent_plastic_strain = 0.0, equivalent_creep_strain = 0.0;
};

struct ElementSnapshot final {
    std::size_t element_id = 0;
    double radius = 0.0;
    double axial_coordinate = 0.0;
    StressValues stress{};
    MaterialValues state{};
};

struct HistorySnapshot final {
    double time = 0.0;
    double reference_height = 0.0;
    double top_radial_displacement = 0.0;
    double top_axial_displacement = 0.0;
    std::vector<ElementSnapshot> elements;
};

class OutputHistory final {
  public:
    OutputHistory(const std::string& path, std::size_t expected_elements) {
        const auto final = fuelsim::test::read_final_exodus_results(path);
        if (final.block_element_counts.size() != 1 || final.block_element_counts[0] != expected_elements)
            throw std::runtime_error("M4.3 production result has incorrect region or element count");
        const std::array<std::string, 4> components = {"rr", "zz", "hoop", "rz"};
        for (std::size_t step = 1; step <= final.step_count; ++step) {
            const auto result = fuelsim::test::read_exodus_results(path, step);
            if (result.time <= 0.0)
                continue;
            HistorySnapshot snapshot;
            snapshot.time = result.time;
            double minimum_z = result.nodes.front()[1], maximum_z = minimum_z;
            for (const auto& point : result.nodes) {
                minimum_z = std::min(minimum_z, point[1]);
                maximum_z = std::max(maximum_z, point[1]);
            }
            snapshot.reference_height = maximum_z - minimum_z;
            std::size_t count = 0;
            for (std::size_t node = 0; node < result.nodes.size(); ++node)
                if (std::abs(result.nodes[node][1] - maximum_z) <= time_tolerance) {
                    ++count;
                    snapshot.top_radial_displacement += result.nodal("displacement_r")[node];
                    snapshot.top_axial_displacement += result.nodal("displacement_z")[node];
                }
            if (count == 0)
                throw std::runtime_error("M4.3 result has no top nodes");
            snapshot.top_radial_displacement /= static_cast<double>(count);
            snapshot.top_axial_displacement /= static_cast<double>(count);
            for (std::size_t element = 0; element < expected_elements; ++element) {
                ElementSnapshot value;
                value.element_id = element;
                double total_weight = 0.0;
                for (std::size_t q = 0; q < 4; ++q) {
                    const auto suffix = "_q" + std::to_string(q);
                    // The mean of four symmetric Gauss coordinates is the Quad4 nodal centroid.
                    value.radius += result.element("reference_r" + suffix)[element] / 4.0;
                    value.axial_coordinate += result.element("reference_z" + suffix)[element] / 4.0;
                    const double weight = result.element("reference_measure" + suffix)[element];
                    if (!(weight > 0.0))
                        throw std::runtime_error("M4.3 result has nonpositive reference measure");
                    total_weight += weight;
                    value.stress.rr += weight * result.element("stress_rr" + suffix)[element];
                    value.stress.zz += weight * result.element("stress_zz" + suffix)[element];
                    value.stress.hoop += weight * result.element("stress_hoop" + suffix)[element];
                    value.stress.rz += weight * result.element("stress_rz" + suffix)[element];
                    double plastic_trace = 0.0, creep_trace = 0.0;
                    for (std::size_t c = 0; c < 4; ++c) {
                        const double elastic = result.element("elastic_" + components[c] + suffix)[element];
                        const double plastic = result.element("plastic_" + components[c] + suffix)[element];
                        const double creep = result.element("creep_" + components[c] + suffix)[element];
                        value.state.elastic_strain[c] += weight * elastic;
                        value.state.plastic_strain[c] += weight * plastic;
                        value.state.creep_strain[c] += weight * creep;
                        if (c < 3) {
                            plastic_trace += plastic;
                            creep_trace += creep;
                        }
                    }
                    _maximum_plastic_trace = std::max(_maximum_plastic_trace, std::abs(plastic_trace));
                    _maximum_creep_trace = std::max(_maximum_creep_trace, std::abs(creep_trace));
                    value.state.equivalent_plastic_strain += weight * result.element("equiv_plastic" + suffix)[element];
                    value.state.equivalent_creep_strain += weight * result.element("equiv_creep" + suffix)[element];
                }
                value.stress.rr /= total_weight;
                value.stress.zz /= total_weight;
                value.stress.hoop /= total_weight;
                value.stress.rz /= total_weight;
                for (std::size_t c = 0; c < 4; ++c) {
                    value.state.elastic_strain[c] /= total_weight;
                    value.state.plastic_strain[c] /= total_weight;
                    value.state.creep_strain[c] /= total_weight;
                }
                value.state.equivalent_plastic_strain /= total_weight;
                value.state.equivalent_creep_strain /= total_weight;
                snapshot.elements.push_back(value);
            }
            _snapshots.push_back(snapshot);
        }
    }

    const std::vector<HistorySnapshot>& snapshots() const noexcept { return _snapshots; }

    double maximum_plastic_trace() const noexcept { return _maximum_plastic_trace; }

    double maximum_creep_trace() const noexcept { return _maximum_creep_trace; }

  private:
    std::vector<HistorySnapshot> _snapshots;
    double _maximum_plastic_trace = 0.0, _maximum_creep_trace = 0.0;
};

struct ReferenceSnapshot final {
    double time = 0.0;
    std::size_t element_id = 0;
    double radius = 0.0;
    double axial_coordinate = 0.0;
    std::array<double, 4> stress{};
    std::array<double, 4> elastic{};
    std::array<double, 4> combined_inelastic{};
    double equivalent_plastic = 0.0;
    double equivalent_creep = 0.0;
};

std::vector<ReferenceSnapshot> read_reference_history(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Could not read M4.3 MOOSE history: " + path);
    std::string line;
    if (!std::getline(input, line))
        throw std::invalid_argument("M4.3 MOOSE history is empty: " + path);
    const std::vector<std::string> header = split_csv(line);
    const auto column = [&header](const std::string& name) {
        return column_index(header, name);
    };
    const std::size_t time = column("sample_time");
    const std::size_t id = column("id");
    const std::size_t radius = column("x");
    const std::size_t axial_coordinate = column("y");
    const std::array<std::size_t, 4> stress = {column("stress_rr"),
        column("stress_zz"),
        column("stress_hoop"),
        column("stress_rz")};
    const std::array<std::size_t, 4> elastic = {column("elastic_rr"),
        column("elastic_zz"),
        column("elastic_hoop"),
        column("elastic_rz")};
    const std::array<std::size_t, 4> combined_inelastic = {column("combined_rr"),
        column("combined_zz"),
        column("combined_hoop"),
        column("combined_rz")};
    const std::size_t equivalent_plastic = column("effective_plastic");
    const std::size_t equivalent_creep = column("effective_creep");
    std::vector<ReferenceSnapshot> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> fields = split_csv(line);
        ReferenceSnapshot value;
        value.time = csv_value(fields, time, path);
        value.element_id = csv_id(fields, id, path);
        value.radius = csv_value(fields, radius, path);
        value.axial_coordinate = csv_value(fields, axial_coordinate, path);
        for (std::size_t component = 0; component < 4; ++component) {
            value.stress[component] = csv_value(fields, stress[component], path);
            value.elastic[component] = csv_value(fields, elastic[component], path);
            value.combined_inelastic[component] = csv_value(fields, combined_inelastic[component], path);
        }
        value.equivalent_plastic = csv_value(fields, equivalent_plastic, path);
        value.equivalent_creep = csv_value(fields, equivalent_creep, path);
        result.push_back(value);
    }
    if (result.empty())
        throw std::invalid_argument("M4.3 MOOSE history has no transient rows: " + path);
    std::sort(result.begin(), result.end(), [](const ReferenceSnapshot& left, const ReferenceSnapshot& right) {
        if (left.time != right.time)
            return left.time < right.time;
        return left.element_id < right.element_id;
    });
    return result;
}

std::array<double, 4> stress_components(const StressValues& stress) {
    return {stress.rr, stress.zz, stress.hoop, stress.rz};
}

bool check_metrics(const std::string& name,
    const fuelsim::test::FieldErrorMetrics& metrics,
    double zero_reference_tolerance,
    double aggregate_tolerance,
    double pointwise_tolerance) {
    fuelsim::test::print_relative_metrics(name, metrics);
    std::cout << name << "_maximum_absolute_difference=" << metrics.maximum_absolute_difference << '\n';
    std::cout << name << "_relative_l2_tolerance=" << aggregate_tolerance << '\n';
    std::cout << name << "_relative_absolute_peak_tolerance=" << aggregate_tolerance << '\n';
    std::cout << name << "_maximum_pointwise_relative_tolerance=" << pointwise_tolerance << '\n';
    return check(metrics.relative_l2() < aggregate_tolerance && metrics.relative_absolute_peak() < aggregate_tolerance
                     && metrics.maximum_pointwise_relative_error() < pointwise_tolerance
                     && metrics.maximum_zero_reference_difference < zero_reference_tolerance,
        name + " three MOOSE metrics and zero-reference error pass");
}

void print_tensor_metric_locations(const std::string& name,
    const fuelsim::test::FieldErrorMetrics& metrics,
    const std::vector<ReferenceSnapshot>& reference) {
    constexpr std::array<const char*, 4> components = {"rr", "zz", "hoop", "rz"};
    const auto print_location = [&name, &reference, &components](const std::string& metric, std::size_t flat_index) {
        const std::size_t row = flat_index / components.size();
        const std::size_t component = flat_index % components.size();
        if (row >= reference.size())
            throw std::logic_error("M4.3 metric index exceeds history size");
        std::cout << name << '_' << metric << "_time=" << reference[row].time << '\n';
        std::cout << name << '_' << metric << "_element_id=" << reference[row].element_id << '\n';
        std::cout << name << '_' << metric << "_component=" << components[component] << '\n';
    };
    print_location("maximum_pointwise_relative_location", metrics.maximum_pointwise_relative_index);
    print_location("maximum_absolute_difference_location", metrics.maximum_absolute_difference_index);
}

const HistorySnapshot& snapshot_at(const std::vector<HistorySnapshot>& values, double time) {
    const auto found = std::find_if(values.begin(), values.end(), [time](const HistorySnapshot& value) {
        return std::abs(value.time - time) < time_tolerance;
    });
    if (found == values.end())
        throw std::invalid_argument("M4.3 accepted history misses event time");
    return *found;
}

bool check_load_path(const std::vector<HistorySnapshot>& snapshots, ExpectedBehavior behavior) {
    const HistorySnapshot& first_stretch = snapshot_at(snapshots, 1.0);
    const HistorySnapshot& positive_shear = snapshot_at(snapshots, 2.0);
    const HistorySnapshot& axial_reversal = snapshot_at(snapshots, 3.0);
    const HistorySnapshot& shear_reversal = snapshot_at(snapshots, 4.0);
    const HistorySnapshot& final_stretch = snapshot_at(snapshots, 5.0);
    const auto close = [](double actual, double expected) {
        return std::abs(actual - expected) < 1.0e-14;
    };
    bool passed =
        check(close(first_stretch.top_radial_displacement, 0.0) && close(first_stretch.top_axial_displacement, 2.0e-4)
                  && close(positive_shear.top_radial_displacement, 2.0e-3)
                  && close(positive_shear.top_axial_displacement, 2.0e-4)
                  && close(axial_reversal.top_radial_displacement, 2.0e-3)
                  && close(axial_reversal.top_axial_displacement, -1.0e-4)
                  && close(shear_reversal.top_radial_displacement, -1.5e-3)
                  && close(shear_reversal.top_axial_displacement, -1.0e-4)
                  && close(final_stretch.top_radial_displacement, -1.5e-3)
                  && close(final_stretch.top_axial_displacement, 6.0e-5),
            "M4.3 hits every noncoaxial load-path event");
    const double shear = positive_shear.top_radial_displacement / positive_shear.reference_height;
    const double axial_stretch = 1.0 + positive_shear.top_axial_displacement / positive_shear.reference_height;
    const double positive_polar_rotation = std::abs(std::atan2(-shear, 1.0 + axial_stretch));
    std::cout << "m43_positive_polar_rotation=" << positive_polar_rotation << '\n';
    passed = check(positive_polar_rotation > 0.44, "M4.3 positive-shear stage exceeds 25 degrees rotation") && passed;
    double maximum_plastic_shear = 0.0;
    double maximum_creep_shear = 0.0;
    for (const ElementSnapshot& element : positive_shear.elements) {
        maximum_plastic_shear = std::max(maximum_plastic_shear, std::abs(element.state.plastic_strain[3]));
        maximum_creep_shear = std::max(maximum_creep_shear, std::abs(element.state.creep_strain[3]));
    }
    if (plastic_active(behavior))
        passed = check(maximum_plastic_shear > 0.1, "M4.3 activates rotated plastic shear history") && passed;
    else
        passed = check(maximum_plastic_shear == 0.0, "M4.3 inactive plastic history remains zero") && passed;
    if (creep_active(behavior))
        passed = check(maximum_creep_shear > 1.0e-5, "M4.3 activates rotated creep shear history") && passed;
    else
        passed = check(maximum_creep_shear == 0.0, "M4.3 inactive creep history remains zero") && passed;
    return passed;
}

} // namespace

bool fuelsim::test::check_rz_noncoaxial(const std::string& output_path,
    const std::string& variant_name,
    const std::string& history_reference_path) {
    const auto variant = variant_config(variant_name);
    OutputHistory observer(output_path, variant.element_count);
    bool passed = check(observer.snapshots().size() == variant.step_count, "M4.3 compares every fixed step");
    passed = check_load_path(observer.snapshots(), variant.behavior) && passed;
    const std::vector<ReferenceSnapshot> reference = read_reference_history(history_reference_path);
    const std::vector<HistorySnapshot>& actual = observer.snapshots();
    passed = check(reference.size() == actual.size() * variant.element_count,
                 "M4.3 compares every element at every accepted step")
             && passed;
    if (reference.size() != actual.size() * variant.element_count)
        return false;
    fuelsim::test::FieldErrorMetrics stress;
    fuelsim::test::FieldErrorMetrics elastic;
    fuelsim::test::FieldErrorMetrics combined_inelastic;
    fuelsim::test::GroupedFieldErrorMetrics stress_tensor;
    fuelsim::test::GroupedFieldErrorMetrics elastic_tensor;
    fuelsim::test::GroupedFieldErrorMetrics combined_inelastic_tensor;
    fuelsim::test::FieldErrorMetrics equivalent_plastic;
    fuelsim::test::FieldErrorMetrics equivalent_creep;
    double maximum_time_difference = 0.0;
    double maximum_coordinate_difference = 0.0;
    std::size_t reference_row = 0;
    for (std::size_t step = 0; step < actual.size(); ++step) {
        passed =
            check(actual[step].elements.size() == variant.element_count, "M4.3 observer records every element per step")
            && passed;
        if (actual[step].elements.size() != variant.element_count)
            return false;
        for (const ElementSnapshot& element : actual[step].elements) {
            const ReferenceSnapshot& expected = reference[reference_row];
            maximum_time_difference = std::max(maximum_time_difference, std::abs(actual[step].time - expected.time));
            maximum_coordinate_difference = std::max({maximum_coordinate_difference,
                std::abs(element.radius - expected.radius),
                std::abs(element.axial_coordinate - expected.axial_coordinate)});
            passed =
                check(element.element_id == expected.element_id, "M4.3 element IDs match at every time step") && passed;
            const std::array<double, 4> actual_stress = stress_components(element.stress);
            std::array<double, 5> actual_stress_tensor{};
            std::array<double, 5> expected_stress_tensor{};
            std::array<double, 5> actual_elastic_tensor{};
            std::array<double, 5> expected_elastic_tensor{};
            std::array<double, 5> actual_combined_inelastic_tensor{};
            std::array<double, 5> expected_combined_inelastic_tensor{};
            for (std::size_t component = 0; component < 4; ++component) {
                stress.add(actual_stress[component], expected.stress[component]);
                elastic.add(element.state.elastic_strain[component], expected.elastic[component]);
                combined_inelastic.add(element.state.plastic_strain[component] + element.state.creep_strain[component],
                    expected.combined_inelastic[component]);
                actual_stress_tensor[component] = actual_stress[component];
                expected_stress_tensor[component] = expected.stress[component];
                actual_elastic_tensor[component] = element.state.elastic_strain[component];
                expected_elastic_tensor[component] = expected.elastic[component];
                actual_combined_inelastic_tensor[component] =
                    element.state.plastic_strain[component] + element.state.creep_strain[component];
                expected_combined_inelastic_tensor[component] = expected.combined_inelastic[component];
            }
            actual_stress_tensor[4] = actual_stress[3];
            expected_stress_tensor[4] = expected.stress[3];
            actual_elastic_tensor[4] = element.state.elastic_strain[3];
            expected_elastic_tensor[4] = expected.elastic[3];
            actual_combined_inelastic_tensor[4] = element.state.plastic_strain[3] + element.state.creep_strain[3];
            expected_combined_inelastic_tensor[4] = expected.combined_inelastic[3];
            stress_tensor.add(actual_stress_tensor.data(), expected_stress_tensor.data(), actual_stress_tensor.size());
            elastic_tensor.add(actual_elastic_tensor.data(),
                expected_elastic_tensor.data(),
                actual_elastic_tensor.size());
            combined_inelastic_tensor.add(actual_combined_inelastic_tensor.data(),
                expected_combined_inelastic_tensor.data(),
                actual_combined_inelastic_tensor.size());
            equivalent_plastic.add(element.state.equivalent_plastic_strain, expected.equivalent_plastic);
            equivalent_creep.add(element.state.equivalent_creep_strain, expected.equivalent_creep);
            ++reference_row;
        }
    }
    passed = check(maximum_time_difference < time_tolerance, "M4.3 MOOSE and fuelsim history times match") && passed;
    passed = check(maximum_coordinate_difference < 1.0e-12, "M4.3 MOOSE and fuelsim element centroids match") && passed;
    std::cout << "m43_compared_element_time_rows=" << reference_row << '\n';
    std::cout << "m43_maximum_element_coordinate_difference=" << maximum_coordinate_difference << '\n';
    std::cout << "m43_maximum_plastic_trace=" << observer.maximum_plastic_trace() << '\n';
    std::cout << "m43_maximum_creep_trace=" << observer.maximum_creep_trace() << '\n';
    passed = check(observer.maximum_plastic_trace() < variant.plastic_trace_tolerance
                       && observer.maximum_creep_trace() < variant.creep_trace_tolerance,
                 "M4.3 default-Rashid accumulated trace drift stays below "
                 "its acceptance limits")
             && passed;
    const std::string prefix = "m43_" + variant.name + "_element_qp_average_";
    print_tensor_metric_locations(prefix + "stress", stress, reference);
    print_tensor_metric_locations(prefix + "elastic_strain", elastic, reference);
    print_tensor_metric_locations(prefix + "combined_inelastic_strain", combined_inelastic, reference);
    fuelsim::test::print_relative_metrics(prefix + "stress_component_diagnostic", stress);
    fuelsim::test::print_relative_metrics(prefix + "elastic_strain_component_diagnostic", elastic);
    fuelsim::test::print_grouped_relative_metrics(prefix + "stress_tensor", stress_tensor);
    fuelsim::test::print_grouped_relative_metrics(prefix + "elastic_strain_tensor", elastic_tensor);
    passed = check(fuelsim::test::grouped_relative_metrics_below(stress_tensor, variant.comparison_tolerance)
                       && stress.maximum_zero_reference_difference < 1.0e-3,
                 prefix + "stress complete-tensor three metrics and component zero-reference error pass")
             && passed;
    passed = check(fuelsim::test::grouped_relative_metrics_below(elastic_tensor, variant.comparison_tolerance)
                       && elastic.maximum_zero_reference_difference < 1.0e-12,
                 prefix + "elastic strain complete-tensor three metrics and component zero-reference error pass")
             && passed;
    if (variant.behavior == ExpectedBehavior::elastic) {
        fuelsim::test::print_absolute_metrics(prefix + "combined_inelastic_strain", combined_inelastic);
        passed =
            check(combined_inelastic.maximum_absolute_difference < 1.0e-14
                      && combined_inelastic.maximum_actual < 1.0e-14 && combined_inelastic.maximum_reference < 1.0e-14,
                "M4.3 elastic combined inelastic history stays zero")
            && passed;
    } else {
        fuelsim::test::print_relative_metrics(prefix + "combined_inelastic_strain_component_diagnostic",
            combined_inelastic);
        fuelsim::test::print_grouped_relative_metrics(prefix + "combined_inelastic_strain_tensor",
            combined_inelastic_tensor);
        passed =
            check(fuelsim::test::grouped_relative_metrics_below(combined_inelastic_tensor, variant.comparison_tolerance)
                      && combined_inelastic.maximum_zero_reference_difference < 1.0e-12,
                prefix
                    + "combined inelastic strain complete-tensor three metrics and component zero-reference "
                      "error pass")
            && passed;
    }
    if (plastic_active(variant.behavior))
        passed = check_metrics(prefix + "equivalent_plastic",
                     equivalent_plastic,
                     1.0e-12,
                     variant.comparison_tolerance,
                     variant.comparison_tolerance)
                 && passed;
    else {
        fuelsim::test::print_absolute_metrics(prefix + "equivalent_plastic", equivalent_plastic);
        passed =
            check(equivalent_plastic.maximum_absolute_difference < 1.0e-14
                      && equivalent_plastic.maximum_actual < 1.0e-14 && equivalent_plastic.maximum_reference < 1.0e-14,
                "M4.3 inactive equivalent plastic history stays zero")
            && passed;
    }
    if (creep_active(variant.behavior))
        passed = check_metrics(prefix + "equivalent_creep",
                     equivalent_creep,
                     1.0e-12,
                     variant.comparison_tolerance,
                     variant.comparison_tolerance)
                 && passed;
    else {
        fuelsim::test::print_absolute_metrics(prefix + "equivalent_creep", equivalent_creep);
        passed = check(equivalent_creep.maximum_absolute_difference < 1.0e-14
                           && equivalent_creep.maximum_actual < 1.0e-14 && equivalent_creep.maximum_reference < 1.0e-14,
                     "M4.3 inactive equivalent creep history stays zero")
                 && passed;
    }
    return passed;
}
