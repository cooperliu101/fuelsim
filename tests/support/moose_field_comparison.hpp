#ifndef FUELSIM_TEST_MOOSE_FIELD_COMPARISON_HPP
#define FUELSIM_TEST_MOOSE_FIELD_COMPARISON_HPP
#include "fuelsim/spatial_definition.hpp"
#include "fuelsim/steady_problem.hpp"
#include "fuelsim/transient_problem.hpp"
#include <cstddef>
#include <string>
#include <vector>
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
    void add(double actual, double reference);
    bool has_relative_norm() const noexcept;
    double relative_l2() const;
    double relative_absolute_peak() const;
    double maximum_pointwise_relative_error() const;
    double absolute_l2() const noexcept;
    double absolute_peak() const noexcept;
};
struct NodalFieldReference final {
    double radius;
    double axial_coordinate;
    double temperature;
    double radial_displacement;
    double axial_displacement;
};
struct NodalFieldComparison final {
    FieldErrorMetrics temperature;
    FieldErrorMetrics radial_displacement;
    FieldErrorMetrics axial_displacement;
    double maximum_coordinate_difference = 0.0;
    std::size_t node_count = 0;
};
std::vector<NodalFieldReference> read_moose_nodal_reference(const std::string& path);
NodalFieldComparison compare_moose_nodal_fields(
    const SteadyProblem& problem, const std::vector<double>& state, const std::vector<NodalFieldReference>& reference);
NodalFieldComparison compare_moose_nodal_fields(const TransientProblem& problem, const std::vector<double>& state,
    const std::vector<NodalFieldReference>& reference);
std::vector<double> read_moose_contact_pressure_reference(const std::string& path, std::vector<double>& coordinates);
FieldErrorMetrics compare_moose_contact_pressure(const std::vector<ContactNodeSummary>& actual,
    const std::vector<double>& reference, const std::vector<double>& reference_coordinates,
    double coordinate_tolerance);
bool relative_metrics_below(const FieldErrorMetrics& metrics, double tolerance);
bool relative_metrics_below_with_pointwise_tolerance(
    const FieldErrorMetrics& metrics, double aggregate_tolerance, double pointwise_tolerance);
bool absolute_metrics_below(const FieldErrorMetrics& metrics, double tolerance);
void print_relative_metrics(const std::string& name, const FieldErrorMetrics& metrics);
void print_absolute_metrics(const std::string& name, const FieldErrorMetrics& metrics);
} // namespace fuelsim::test
#endif
