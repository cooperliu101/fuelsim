#include "support/field_error_metrics.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace fuelsim::test {
void FieldErrorMetrics::add(double actual, double reference) {
    if (!std::isfinite(actual) || !std::isfinite(reference))
        throw std::invalid_argument("Full-field values must be finite");
    const double difference = actual - reference;
    difference_squared += difference * difference;
    reference_squared += reference * reference;
    maximum_actual = std::max(maximum_actual, std::abs(actual));
    maximum_reference = std::max(maximum_reference, std::abs(reference));
    if (std::abs(difference) > maximum_absolute_difference) {
        maximum_absolute_difference = std::abs(difference);
        maximum_absolute_difference_index = value_count;
        maximum_absolute_difference_actual = actual;
        maximum_absolute_difference_reference = reference;
    }
    if (reference != 0.0) {
        const double relative = std::abs(difference) / std::abs(reference);
        if (relative > maximum_pointwise_relative) {
            maximum_pointwise_relative = relative;
            maximum_pointwise_relative_index = value_count;
            maximum_pointwise_relative_actual = actual;
            maximum_pointwise_relative_reference = reference;
        }
        ++nonzero_reference_count;
    } else {
        maximum_zero_reference_difference = std::max(maximum_zero_reference_difference, std::abs(difference));
        ++zero_reference_count;
    }
    ++value_count;
}

bool FieldErrorMetrics::has_relative_norm() const noexcept {
    return reference_squared > 0.0 && maximum_reference > 0.0 && nonzero_reference_count > 0;
}

double FieldErrorMetrics::relative_l2() const {
    if (!has_relative_norm()) throw std::domain_error("Relative full-field norm is undefined");
    return std::sqrt(difference_squared / reference_squared);
}

double FieldErrorMetrics::relative_absolute_peak() const {
    if (!has_relative_norm()) throw std::domain_error("Relative full-field peak is undefined");
    return std::abs(maximum_actual - maximum_reference) / maximum_reference;
}

double FieldErrorMetrics::maximum_pointwise_relative_error() const {
    if (!has_relative_norm()) throw std::domain_error("Pointwise relative full-field error is undefined");
    return maximum_pointwise_relative;
}

double FieldErrorMetrics::absolute_l2() const noexcept { return std::sqrt(difference_squared); }

double FieldErrorMetrics::absolute_peak() const noexcept { return std::abs(maximum_actual - maximum_reference); }

void GroupedFieldErrorMetrics::add(const double* actual, const double* reference, std::size_t component_count) {
    if (actual == nullptr || reference == nullptr || component_count == 0)
        throw std::invalid_argument("Grouped full-field values require nonempty component arrays");
    double group_difference_squared = 0.0, group_reference_squared = 0.0, group_actual_squared = 0.0;
    for (std::size_t component = 0; component < component_count; ++component) {
        if (!std::isfinite(actual[component]) || !std::isfinite(reference[component]))
            throw std::invalid_argument("Grouped full-field values must be finite");
        const double difference = actual[component] - reference[component];
        group_difference_squared += difference * difference;
        group_reference_squared += reference[component] * reference[component];
        group_actual_squared += actual[component] * actual[component];
    }
    const double difference_norm = std::sqrt(group_difference_squared),
                 reference_norm = std::sqrt(group_reference_squared), actual_norm = std::sqrt(group_actual_squared);
    difference_squared += group_difference_squared;
    reference_squared += group_reference_squared;
    maximum_reference = std::max(maximum_reference, reference_norm);
    if (difference_norm > maximum_difference) {
        maximum_difference = difference_norm;
        maximum_difference_index = group_count;
    }
    if (reference_norm != 0.0) {
        const double relative = difference_norm / reference_norm;
        if (relative > maximum_pointwise_relative) {
            maximum_pointwise_relative = relative;
            maximum_pointwise_relative_index = group_count;
            maximum_pointwise_actual_norm = actual_norm;
            maximum_pointwise_reference_norm = reference_norm;
        }
        ++nonzero_reference_count;
    } else {
        maximum_zero_reference_difference = std::max(maximum_zero_reference_difference, difference_norm);
        ++zero_reference_count;
    }
    ++group_count;
}

bool GroupedFieldErrorMetrics::has_relative_norm() const noexcept {
    return reference_squared > 0.0 && maximum_reference > 0.0 && nonzero_reference_count > 0;
}

double GroupedFieldErrorMetrics::relative_l2() const {
    if (!has_relative_norm()) throw std::domain_error("Grouped relative full-field norm is undefined");
    return std::sqrt(difference_squared / reference_squared);
}

double GroupedFieldErrorMetrics::relative_absolute_peak() const {
    if (!has_relative_norm()) throw std::domain_error("Grouped relative full-field peak is undefined");
    return maximum_difference / maximum_reference;
}

bool relative_metrics_below(const FieldErrorMetrics& metrics, double tolerance) {
    return metrics.relative_l2() < tolerance && metrics.relative_absolute_peak() < tolerance &&
           metrics.maximum_pointwise_relative_error() < tolerance;
}

bool relative_metrics_below_with_pointwise_tolerance(
    const FieldErrorMetrics& metrics, double aggregate_tolerance, double pointwise_tolerance) {
    return metrics.relative_l2() < aggregate_tolerance && metrics.relative_absolute_peak() < aggregate_tolerance &&
           metrics.maximum_pointwise_relative_error() < pointwise_tolerance;
}

bool absolute_metrics_below(const FieldErrorMetrics& metrics, double tolerance) {
    return metrics.absolute_l2() < tolerance && metrics.absolute_peak() < tolerance &&
           metrics.maximum_absolute_difference < tolerance;
}

bool grouped_relative_metrics_below(const GroupedFieldErrorMetrics& metrics, double tolerance) {
    return metrics.relative_l2() < tolerance && metrics.relative_absolute_peak() < tolerance &&
           metrics.maximum_pointwise_relative < tolerance;
}

void print_relative_metrics(const std::string& name, const FieldErrorMetrics& metrics) {
    std::cout << name << "_relative_l2=" << metrics.relative_l2() << '\n';
    std::cout << name << "_relative_absolute_peak=" << metrics.relative_absolute_peak() << '\n';
    std::cout << name << "_maximum_pointwise_relative=" << metrics.maximum_pointwise_relative_error() << '\n';
    std::cout << name << "_maximum_pointwise_relative_index=" << metrics.maximum_pointwise_relative_index << '\n';
    std::cout << name << "_maximum_pointwise_relative_actual=" << metrics.maximum_pointwise_relative_actual << '\n';
    std::cout << name << "_maximum_pointwise_relative_reference=" << metrics.maximum_pointwise_relative_reference
              << '\n';
    std::cout << name << "_maximum_absolute_difference_index=" << metrics.maximum_absolute_difference_index << '\n';
    std::cout << name << "_maximum_absolute_difference_actual=" << metrics.maximum_absolute_difference_actual << '\n';
    std::cout << name << "_maximum_absolute_difference_reference=" << metrics.maximum_absolute_difference_reference
              << '\n';
    std::cout << name << "_zero_reference_count=" << metrics.zero_reference_count << '\n';
    std::cout << name << "_maximum_zero_reference_absolute_difference=" << metrics.maximum_zero_reference_difference
              << '\n';
}

void print_absolute_metrics(const std::string& name, const FieldErrorMetrics& metrics) {
    std::cout << name << "_absolute_l2=" << metrics.absolute_l2() << '\n';
    std::cout << name << "_absolute_peak=" << metrics.absolute_peak() << '\n';
    std::cout << name << "_maximum_pointwise_absolute=" << metrics.maximum_absolute_difference << '\n';
}

void print_grouped_relative_metrics(const std::string& name, const GroupedFieldErrorMetrics& metrics) {
    std::cout << name << "_relative_l2=" << metrics.relative_l2() << '\n';
    std::cout << name << "_relative_absolute_peak=" << metrics.relative_absolute_peak() << '\n';
    std::cout << name << "_maximum_pointwise_relative=" << metrics.maximum_pointwise_relative << '\n';
    std::cout << name << "_maximum_pointwise_relative_index=" << metrics.maximum_pointwise_relative_index << '\n';
    std::cout << name << "_maximum_pointwise_actual_norm=" << metrics.maximum_pointwise_actual_norm << '\n';
    std::cout << name << "_maximum_pointwise_reference_norm=" << metrics.maximum_pointwise_reference_norm << '\n';
    std::cout << name << "_maximum_difference_index=" << metrics.maximum_difference_index << '\n';
    std::cout << name << "_maximum_absolute_group_difference=" << metrics.maximum_difference << '\n';
    std::cout << name << "_zero_reference_count=" << metrics.zero_reference_count << '\n';
    std::cout << name << "_maximum_zero_reference_absolute_difference=" << metrics.maximum_zero_reference_difference
              << '\n';
}
} // namespace fuelsim::test
