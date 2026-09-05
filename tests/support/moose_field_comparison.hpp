#ifndef FUELSIM_TEST_MOOSE_FIELD_COMPARISON_HPP
#define FUELSIM_TEST_MOOSE_FIELD_COMPARISON_HPP
#include "fuelsim/core/spatial_definition.hpp"
#include "fuelsim/core/steady_problem.hpp"
#include "fuelsim/core/transient_problem.hpp"
#include "support/field_error_metrics.hpp"
#include <cstddef>
#include <string>
#include <vector>

namespace fuelsim::test {

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
} // namespace fuelsim::test
#endif
