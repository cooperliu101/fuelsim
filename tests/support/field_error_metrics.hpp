#pragma once
#include <cstddef>
#include <string>

namespace fuelsim::test {
struct FieldErrorMetrics final {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    double maximum_actual = 0.0;
    double maximum_reference = 0.0;
    double maximum_absolute_difference = 0.0;
    double maximum_pointwise_relative = 0.0;
    double maximum_zero_reference_difference = 0.0;
    std::size_t value_count = 0;
    std::size_t nonzero_reference_count = 0;
    std::size_t zero_reference_count = 0;
    std::size_t maximum_absolute_difference_index = 0;
    std::size_t maximum_pointwise_relative_index = 0;
    double maximum_pointwise_relative_actual = 0.0;
    double maximum_pointwise_relative_reference = 0.0;
    double maximum_absolute_difference_actual = 0.0;
    double maximum_absolute_difference_reference = 0.0;
    void add(double actual, double reference);
    bool has_relative_norm() const noexcept;
    double relative_l2() const;
    double relative_absolute_peak() const;
    double maximum_pointwise_relative_error() const;
    double absolute_l2() const noexcept;
    double absolute_peak() const noexcept;
};

struct GroupedFieldErrorMetrics final {
    double difference_squared = 0.0;
    double reference_squared = 0.0;
    double maximum_difference = 0.0;
    double maximum_reference = 0.0;
    double maximum_pointwise_relative = 0.0;
    double maximum_zero_reference_difference = 0.0;
    std::size_t group_count = 0;
    std::size_t nonzero_reference_count = 0;
    std::size_t zero_reference_count = 0;
    std::size_t maximum_pointwise_relative_index = 0;
    std::size_t maximum_difference_index = 0;
    double maximum_pointwise_actual_norm = 0.0;
    double maximum_pointwise_reference_norm = 0.0;
    void add(const double* actual, const double* reference, std::size_t component_count);
    bool has_relative_norm() const noexcept;
    double relative_l2() const;
    double relative_absolute_peak() const;
};

bool relative_metrics_below(const FieldErrorMetrics& metrics, double tolerance);
bool relative_metrics_below_with_pointwise_tolerance(
    const FieldErrorMetrics& metrics, double aggregate_tolerance, double pointwise_tolerance);
bool absolute_metrics_below(const FieldErrorMetrics& metrics, double tolerance);
bool grouped_relative_metrics_below(const GroupedFieldErrorMetrics& metrics, double tolerance);
void print_relative_metrics(const std::string& name, const FieldErrorMetrics& metrics);
void print_absolute_metrics(const std::string& name, const FieldErrorMetrics& metrics);
void print_grouped_relative_metrics(const std::string& name, const GroupedFieldErrorMetrics& metrics);
} // namespace fuelsim::test
