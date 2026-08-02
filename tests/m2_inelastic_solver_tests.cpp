#include "fuelsim/case_input.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/problem_solver.hpp"
#include "support/moose_field_comparison.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

constexpr double moose_relative_tolerance = 1.0e-3;

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

double relative_error(double actual, double expected) {
    if (!std::isfinite(actual) || !std::isfinite(expected) || expected == 0.0)
        return std::numeric_limits<double>::infinity();
    return std::abs(actual - expected) / std::abs(expected);
}

bool check_scalar_metrics(const std::string& name, double actual,
                          double reference, double tolerance) {
    const double relative_l2 = relative_error(actual, reference);
    const double relative_absolute_peak =
        std::abs(std::abs(actual) - std::abs(reference)) / std::abs(reference);
    const double maximum_pointwise_relative = relative_l2;
    std::cout << name << "_relative_l2=" << relative_l2 << '\n';
    std::cout << name << "_relative_absolute_peak=" << relative_absolute_peak
              << '\n';
    std::cout << name
              << "_maximum_pointwise_relative=" << maximum_pointwise_relative
              << '\n';
    return check(relative_l2 < tolerance &&
                     relative_absolute_peak < tolerance &&
                     maximum_pointwise_relative < tolerance,
                 name + " three MOOSE error metrics pass");
}

fuelsim::SolverOptions
solver_options(const fuelsim::FuelSimCaseDefinition& definition) {
    return {definition.solver.absolute_tolerance,
            definition.solver.relative_tolerance,
            definition.solver.step_tolerance,
            definition.solver.maximum_iterations};
}

fuelsim::TransientTimeOptions
time_options(const fuelsim::FuelSimCaseDefinition& definition) {
    return {definition.transient_execution.end_time,
            definition.transient_execution.initial_time_step,
            definition.transient_execution.minimum_time_step,
            definition.transient_execution.maximum_time_step,
            definition.transient_execution.growth_factor,
            definition.transient_execution.cutback_factor,
            definition.transient_execution.maximum_cutbacks,
            definition.transient_execution.load_ramp_time};
}

class TransientCaseRun final {
  public:
    explicit TransientCaseRun(const std::string& input_path)
        : _definition(fuelsim::CaseInputReader::read(input_path)),
          _source(fuelsim::ExodusMeshIo::read_quad4(_definition.mesh_file)),
          _problem(_definition.transient_definition(), _source),
          _result(fuelsim::solve_transient(_problem, time_options(_definition),
                                           solver_options(_definition))) {
        if (_definition.problem != fuelsim::CaseProblem::transient)
            throw std::invalid_argument(
                "M2.2 comparison requires a transient input card");
        if (_problem.region_count() != 1 ||
            _problem.region_mesh(0).elements().size() != 1)
            throw std::invalid_argument(
                "M2.2 comparison input requires one region and one Quad4");
    }

    const fuelsim::FuelSimCaseDefinition& definition() const noexcept {
        return _definition;
    }

    const fuelsim::UnstructuredQuad4Mesh& source() const noexcept {
        return _source;
    }

    const fuelsim::TransientProblem& problem() const noexcept {
        return _problem;
    }

    const fuelsim::TransientResult& result() const noexcept {
        return _result;
    }

  private:
    fuelsim::FuelSimCaseDefinition _definition;
    fuelsim::UnstructuredQuad4Mesh _source;
    fuelsim::TransientProblem _problem;
    fuelsim::TransientResult _result;
};

double average_axial_stress(const TransientCaseRun& run) {
    double value = 0.0;
    for (const fuelsim::AxisymmetricStressValues& stress :
         run.problem().material_stress(0, 0))
        value += stress.zz / 4.0;
    return value;
}

double average_equivalent_plastic(const TransientCaseRun& run) {
    double value = 0.0;
    for (const fuelsim::MaterialPointState& point :
         run.problem().material_history(0, 0))
        value += point.equivalent_plastic_strain / 4.0;
    return value;
}

double average_equivalent_creep(const TransientCaseRun& run) {
    double value = 0.0;
    for (const fuelsim::MaterialPointState& point :
         run.problem().material_history(0, 0))
        value += point.equivalent_creep_strain / 4.0;
    return value;
}

double maximum_inelastic_trace(const TransientCaseRun& run,
                               bool plastic_strain) {
    double maximum = 0.0;
    for (const fuelsim::MaterialPointState& point :
         run.problem().material_history(0, 0)) {
        const std::array<double, 4>& strain =
            plastic_strain ? point.plastic_strain : point.creep_strain;
        maximum =
            std::max(maximum, std::abs(strain[0] + strain[1] + strain[2]));
    }
    return maximum;
}

double average_top_displacement(const TransientCaseRun& run) {
    const fuelsim::RegionBoundary top =
        run.problem().region_mesh(0).map_side_set(run.source(), "top");
    double value = 0.0;
    for (std::size_t node : top.nodes) {
        value += run.result().committed_state.at(
            run.problem().dof_map().axial_displacement(node));
    }
    return value / static_cast<double>(top.nodes.size());
}

bool temperatures_are_600(const TransientCaseRun& run) {
    for (std::size_t node = 0; node < run.problem().dof_map().node_count();
         ++node) {
        if (std::abs(run.result().committed_state.at(
                         run.problem().dof_map().temperature(node)) -
                     600.0) > 1.0e-12)
            return false;
    }
    return true;
}

bool common_run_checks(const std::string& name, const TransientCaseRun& run,
                       std::size_t expected_steps,
                       const std::string& nodal_reference_path) {
    bool passed =
        check(run.result().completed, name + " input-card load path converged");
    passed = check(run.result().accepted_steps.size() == expected_steps,
                   name + " commits every configured time step") &&
             passed;
    passed = check(run.result().aggregate_timing.workspace_setups == 1,
                   name + " reuses one PETSc workspace") &&
             passed;
    passed = check(temperatures_are_600(run),
                   name + " committed temperature remains 600 K") &&
             passed;
    const std::vector<fuelsim::test::NodalFieldReference> reference =
        fuelsim::test::read_moose_nodal_reference(nodal_reference_path);
    const fuelsim::test::NodalFieldComparison fields =
        fuelsim::test::compare_moose_nodal_fields(
            run.problem(), run.result().committed_state, reference);
    passed =
        check(fields.node_count == run.source().nodes().size() &&
                  fields.maximum_coordinate_difference < 1.0e-12,
              name + " compares every MOOSE node at matching coordinates") &&
        passed;
    passed = check(fuelsim::test::relative_metrics_below(
                       fields.temperature, moose_relative_tolerance),
                   name + " full-field temperature three errors pass") &&
             passed;
    passed =
        check(fuelsim::test::relative_metrics_below(fields.radial_displacement,
                                                    moose_relative_tolerance),
              name + " full-field radial displacement three errors pass") &&
        passed;
    passed = check(fuelsim::test::relative_metrics_below(
                       fields.axial_displacement, moose_relative_tolerance),
                   name + " full-field axial displacement three errors pass") &&
             passed;
    fuelsim::test::print_relative_metrics("m22_" + name + "_temperature",
                                          fields.temperature);
    fuelsim::test::print_relative_metrics(
        "m22_" + name + "_radial_displacement", fields.radial_displacement);
    fuelsim::test::print_relative_metrics("m22_" + name + "_axial_displacement",
                                          fields.axial_displacement);
    return passed;
}

bool test_j2_moose_comparison(const std::string& input_path,
                              const std::string& nodal_reference_path) {
    constexpr double moose_axial_stress = 201980198.0198;
    constexpr double moose_equivalent_plastic = 0.000990099009901;
    constexpr double moose_top_displacement = 2.0e-6;
    const TransientCaseRun run(input_path);
    const double axial_stress = average_axial_stress(run);
    const double equivalent_plastic = average_equivalent_plastic(run);
    const double top_displacement = average_top_displacement(run);

    bool passed = common_run_checks("j2", run, 10, nodal_reference_path);
    passed =
        check_scalar_metrics("m22_j2_axial_stress", axial_stress,
                             moose_axial_stress, moose_relative_tolerance) &&
        passed;
    passed = check_scalar_metrics("m22_j2_equivalent_plastic",
                                  equivalent_plastic, moose_equivalent_plastic,
                                  moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_j2_top_displacement", top_displacement,
                                  moose_top_displacement,
                                  moose_relative_tolerance) &&
             passed;
    passed = check(maximum_inelastic_trace(run, true) < 1.0e-12,
                   "J2 plastic strain is trace-free") &&
             passed;
    return passed;
}

bool test_norton_moose_comparison(const std::string& input_path,
                                  const std::string& nodal_reference_path) {
    constexpr double moose_axial_stress = 99998007.620195;
    constexpr double moose_equivalent_creep = 9.9991036503199e-5;
    constexpr double moose_top_displacement = 5.9997808603447e-7;
    const TransientCaseRun run(input_path);
    const double axial_stress = average_axial_stress(run);
    const double equivalent_creep = average_equivalent_creep(run);
    const double top_displacement = average_top_displacement(run);

    bool passed = common_run_checks("norton", run, 10, nodal_reference_path);
    passed =
        check_scalar_metrics("m22_norton_axial_stress", axial_stress,
                             moose_axial_stress, moose_relative_tolerance) &&
        passed;
    passed = check_scalar_metrics("m22_norton_equivalent_creep",
                                  equivalent_creep, moose_equivalent_creep,
                                  moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_norton_top_displacement",
                                  top_displacement, moose_top_displacement,
                                  moose_relative_tolerance) &&
             passed;
    passed = check(maximum_inelastic_trace(run, false) < 1.0e-12,
                   "Norton creep strain is trace-free") &&
             passed;
    return passed;
}

bool test_coupled_moose_comparison(const std::string& displacement_input,
                                   const std::string& traction_input,
                                   const std::string& displacement_reference,
                                   const std::string& traction_reference) {
    constexpr double moose_axial_stress = 200999992.08159;
    constexpr double moose_equivalent_plastic = 4.9999406119027e-4;
    constexpr double moose_equivalent_creep = 2.4564701918764e-4;
    constexpr double moose_top_displacement = 1.7506410289082e-6;
    constexpr double analytic_equivalent_plastic = 5.0e-4;
    constexpr double analytic_equivalent_creep = 2.4564818025e-4;
    constexpr double analytic_top_displacement = 1.75064818025e-6;
    const TransientCaseRun traction(traction_input);
    const double axial_stress = average_axial_stress(traction);
    const double equivalent_plastic = average_equivalent_plastic(traction);
    const double equivalent_creep = average_equivalent_creep(traction);
    const double top_displacement = average_top_displacement(traction);

    bool passed =
        common_run_checks("coupled_traction", traction, 10, traction_reference);
    passed =
        check_scalar_metrics("m22_coupled_traction_axial_stress", axial_stress,
                             moose_axial_stress, moose_relative_tolerance) &&
        passed;
    passed = check_scalar_metrics("m22_coupled_traction_equivalent_plastic",
                                  equivalent_plastic, moose_equivalent_plastic,
                                  moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_traction_equivalent_creep",
                                  equivalent_creep, moose_equivalent_creep,
                                  moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_traction_top_displacement",
                                  top_displacement, moose_top_displacement,
                                  moose_relative_tolerance) &&
             passed;
    passed = check(relative_error(equivalent_plastic,
                                  analytic_equivalent_plastic) < 1.0e-8 &&
                       relative_error(equivalent_creep,
                                      analytic_equivalent_creep) < 1.0e-8 &&
                       relative_error(top_displacement,
                                      analytic_top_displacement) < 1.0e-8,
                   "coupled traction matches independent uniaxial history") &&
             passed;
    passed = check(maximum_inelastic_trace(traction, true) < 1.0e-12 &&
                       maximum_inelastic_trace(traction, false) < 1.0e-12,
                   "coupled traction inelastic strains are trace-free") &&
             passed;

    constexpr double displacement_moose_stress = 200963368.63564;
    constexpr double displacement_moose_plastic = 4.8168428343307e-4;
    constexpr double displacement_moose_creep = 5.1349887359505e-4;
    constexpr double displacement_moose_displacement = 2.0e-6;
    const TransientCaseRun displacement(displacement_input);
    passed = common_run_checks("coupled_displacement", displacement, 10,
                               displacement_reference) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_displacement_axial_stress",
                                  average_axial_stress(displacement),
                                  displacement_moose_stress,
                                  moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_displacement_equivalent_plastic",
                                  average_equivalent_plastic(displacement),
                                  displacement_moose_plastic,
                                  moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_displacement_equivalent_creep",
                                  average_equivalent_creep(displacement),
                                  displacement_moose_creep,
                                  moose_relative_tolerance) &&
             passed;
    passed = check_scalar_metrics("m22_coupled_displacement_top_displacement",
                                  average_top_displacement(displacement),
                                  displacement_moose_displacement,
                                  moose_relative_tolerance) &&
             passed;
    return passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 9) {
        std::cerr << "Usage: fuelsim_m2_inelastic_solver_tests "
                     "<j2.fsi> <norton.fsi> <coupled_displacement.fsi> "
                     "<coupled_traction.fsi> <j2-nodes.csv> "
                     "<norton-nodes.csv> <coupled-displacement-nodes.csv> "
                     "<coupled-traction-nodes.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(
            argc, argv, "fuelsim input-card M2.2 MOOSE comparison tests\n");
        bool passed = test_j2_moose_comparison(argv[1], argv[5]);
        passed = test_norton_moose_comparison(argv[2], argv[6]) && passed;
        passed =
            test_coupled_moose_comparison(argv[3], argv[4], argv[7], argv[8]) &&
            passed;
        if (!passed)
            return 1;
        std::cout << "[PASS] input-card M2.2 MOOSE comparison tests\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M2.2 tests raised: " << error.what() << '\n';
        return 1;
    }
}
