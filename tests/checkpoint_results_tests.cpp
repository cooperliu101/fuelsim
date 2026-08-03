#include "fuelsim/case_input.hpp"
#include "fuelsim/checkpoint_io.hpp"
#include "fuelsim/exodus_mesh_io.hpp"
#include "fuelsim/petsc_solver.hpp"
#include "fuelsim/problem_solver.hpp"
#include "fuelsim/results_io.hpp"

#include <exodusII.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

bool check(bool condition, const std::string& message) {
    if (condition)
        return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

bool expect_failure(const std::function<void()>& function,
                    const std::string& expected, const std::string& message) {
    try {
        function();
    } catch (const std::exception& error) {
        return check(std::string(error.what()).find(expected) !=
                         std::string::npos,
                     message);
    }
    return check(false, message);
}

bool nearly_equal(double left, double right) {
    const double scale = std::max({1.0, std::abs(left), std::abs(right)});
    return std::abs(left - right) <= 2.0e-13 * scale;
}

bool compare_committed_states(const fuelsim::TransientCommittedState& left,
                              const fuelsim::TransientCommittedState& right) {
    bool passed =
        check(nearly_equal(left.time, right.time) &&
                  nearly_equal(left.load_factor, right.load_factor),
              "restart preserves committed time and load") &&
        check(left.solution.size() == right.solution.size(),
              "restart preserves nodal-state layout") &&
        check(left.material_histories.size() == right.material_histories.size(),
              "restart preserves material-region layout");
    if (!passed)
        return false;
    for (std::size_t dof = 0; dof < left.solution.size(); ++dof)
        passed = check(nearly_equal(left.solution[dof], right.solution[dof]),
                       "restart reproduces the uninterrupted nodal state") &&
                 passed;
    for (std::size_t region = 0; region < left.material_histories.size();
         ++region) {
        if (!check(left.material_histories[region].size() ==
                       right.material_histories[region].size(),
                   "restart preserves material-element layout"))
            return false;
        for (std::size_t element = 0;
             element < left.material_histories[region].size(); ++element) {
            for (std::size_t q = 0; q < 4; ++q) {
                const fuelsim::MaterialPointState& a =
                    left.material_histories[region][element][q];
                const fuelsim::MaterialPointState& b =
                    right.material_histories[region][element][q];
                for (std::size_t component = 0; component < 4; ++component) {
                    passed =
                        check(nearly_equal(a.elastic_strain[component],
                                           b.elastic_strain[component]) &&
                                  nearly_equal(a.plastic_strain[component],
                                           b.plastic_strain[component]) &&
                                  nearly_equal(a.creep_strain[component],
                                               b.creep_strain[component]),
                              "restart reproduces committed tensor history") &&
                        passed;
                }
                passed = check(nearly_equal(a.equivalent_plastic_strain,
                                            b.equivalent_plastic_strain) &&
                                   nearly_equal(a.equivalent_creep_strain,
                                                b.equivalent_creep_strain),
                               "restart reproduces committed scalar history") &&
                         passed;
                const fuelsim::AxisymmetricStressValues& stress_a =
                    left.material_stresses[region][element][q];
                const fuelsim::AxisymmetricStressValues& stress_b =
                    right.material_stresses[region][element][q];
                passed = check(nearly_equal(stress_a.rr, stress_b.rr) &&
                                   nearly_equal(stress_a.zz, stress_b.zz) &&
                                   nearly_equal(stress_a.hoop, stress_b.hoop) &&
                                   nearly_equal(stress_a.rz, stress_b.rz),
                               "restart reproduces committed stresses") &&
                         passed;
            }
        }
    }
    return passed;
}

fuelsim::TransientTimeOptions time_options(double end_time) {
    return {end_time, 1.0, 0.125, 1.0, 1.0, 0.5, 3, 20.0};
}

class ResultsObserver final : public fuelsim::TransientStepObserver {
  public:
    explicit ResultsObserver(fuelsim::ExodusTransientResultsWriter& writer)
        : _writer(writer), _steps(0) {}

    void accepted_step(const fuelsim::TransientProblem& problem,
                       const fuelsim::TransientAcceptedStep&) override {
        _writer.append(problem);
        ++_steps;
    }

    std::size_t steps() const noexcept {
        return _steps;
    }

  private:
    fuelsim::ExodusTransientResultsWriter& _writer;
    std::size_t _steps;
};

bool verify_exodus(const std::string& path,
                   const fuelsim::UnstructuredQuad4Mesh& mesh,
                   const fuelsim::TransientProblem& problem,
                   std::size_t expected_steps) {
    int cpu_word_size = static_cast<int>(sizeof(double));
    int io_word_size = 0;
    float version = 0.0F;
    const int exoid =
        ex_open(path.c_str(), EX_READ, &cpu_word_size, &io_word_size, &version);
    if (!check(exoid >= 0, "transient Exodus result can be reopened"))
        return false;
    ex_set_int64_status(exoid, EX_ALL_INT64_API);
    int nodal_variables = 0;
    int element_variables = 0;
    int global_variables = 0;
    bool passed =
        check(ex_inquire_int(exoid, EX_INQ_TIME) ==
                  static_cast<std::int64_t>(expected_steps),
              "Exodus stores the initial and every accepted committed step") &&
        check(ex_get_variable_param(exoid, EX_NODAL, &nodal_variables) == 0 &&
                  nodal_variables == 5,
              "Exodus defines temperature, displacement, gap and pressure") &&
        check(ex_get_variable_param(exoid, EX_ELEM_BLOCK, &element_variables) ==
                      0 &&
                  element_variables == 56,
              "Exodus defines stress and inelastic integration-point fields") &&
        check(ex_get_variable_param(exoid, EX_GLOBAL, &global_variables) == 0 &&
                  global_variables == 3,
              "Exodus defines load and conservative interface totals");

    const int last_step = static_cast<int>(expected_steps);
    double time = 0.0;
    std::vector<double> temperatures(mesh.nodes().size(), 0.0);
    passed = check(ex_get_time(exoid, last_step, &time) == 0 &&
                       nearly_equal(time, problem.committed_time()),
                   "Exodus last time equals the committed physical time") &&
             check(ex_get_var(exoid, last_step, EX_NODAL, 1, 1,
                              static_cast<std::int64_t>(temperatures.size()),
                              temperatures.data()) == 0,
                   "Exodus temperature field is readable") &&
             passed;
    for (std::size_t region = 0; region < problem.region_count(); ++region) {
        const fuelsim::RegionMesh& region_mesh = problem.region_mesh(region);
        const std::size_t offset = problem.region_node_offset(region);
        for (std::size_t node = 0; node < region_mesh.nodes().size(); ++node) {
            const std::size_t source = region_mesh.source_node_ids()[node];
            const double expected = problem.committed_solution().at(
                problem.dof_map().temperature(offset + node));
            passed =
                check(nearly_equal(temperatures.at(source), expected),
                      "Exodus nodal temperature uses source-mesh mapping") &&
                passed;
        }
    }
    passed =
        check(ex_close(exoid) == 0, "Exodus result closes cleanly") && passed;
    return passed;
}

bool verify_steady_results(const fuelsim::FuelSimCaseDefinition& input,
                           const fuelsim::UnstructuredQuad4Mesh& mesh,
                           const std::string& results_path) {
    fuelsim::SteadyProblem problem(input.steady_definition(), mesh);
    const fuelsim::SolverOptions solver{
        input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};
    const fuelsim::SteadyResult result = fuelsim::solve_steady(
        problem,
        {input.steady_execution.load_steps,
         input.steady_execution.cutback_factor,
         input.steady_execution.maximum_cutbacks,
         input.steady_execution.minimum_load_increment},
        solver);
    if (!check(result.completed && result.solve.converged,
               "steady result fixture converges"))
        return false;
    fuelsim::ExodusResultsIo::write_steady(results_path, mesh, problem,
                                           result.solve.state);

    int cpu_word_size = static_cast<int>(sizeof(double));
    int io_word_size = 0;
    float version = 0.0F;
    const int exoid = ex_open(results_path.c_str(), EX_READ, &cpu_word_size,
                              &io_word_size, &version);
    if (!check(exoid >= 0, "steady Exodus result can be reopened"))
        return false;
    int nodal_variables = 0;
    int element_variables = 0;
    int global_variables = 0;
    const bool passed =
        check(ex_inquire_int(exoid, EX_INQ_TIME) == 1,
              "steady Exodus result contains one final state") &&
        check(ex_get_variable_param(exoid, EX_NODAL, &nodal_variables) == 0 &&
                  nodal_variables == 5,
              "steady Exodus result contains nodal contact fields") &&
        check(ex_get_variable_param(exoid, EX_ELEM_BLOCK, &element_variables) ==
                      0 &&
                  element_variables == 16,
              "steady Exodus result contains four-point stresses") &&
        check(ex_get_variable_param(exoid, EX_GLOBAL, &global_variables) == 0 &&
                  global_variables == 3,
              "steady Exodus result contains interface totals") &&
        check(ex_close(exoid) == 0, "steady Exodus result closes cleanly");
    return passed;
}

bool run_tests(const std::string& steady_input_path,
               const std::string& transient_input_path,
               const std::string& checkpoint_path,
               const std::string& results_path) {
    const fuelsim::FuelSimCaseDefinition steady_input =
        fuelsim::CaseInputReader::read(steady_input_path);
    const std::filesystem::path configured_results(results_path);
    const std::filesystem::path first_segment =
        configured_results.parent_path() /
        (configured_results.stem().string() + ".part1" +
         configured_results.extension().string());
    {
        std::ofstream occupied(first_segment);
        occupied << "previous segment";
    }
    const std::filesystem::path second_segment =
        configured_results.parent_path() /
        (configured_results.stem().string() + ".part2" +
         configured_results.extension().string());
    bool passed =
        check(fuelsim::next_results_segment_path(results_path) ==
                  second_segment.string(),
              "restart result segmentation preserves occupied earlier files");
    std::filesystem::remove(first_segment);
    const fuelsim::UnstructuredQuad4Mesh steady_mesh =
        fuelsim::ExodusMeshIo::read_quad4(steady_input.mesh_file);
    passed = verify_steady_results(steady_input, steady_mesh, results_path) &&
             passed;

    const fuelsim::FuelSimCaseDefinition input =
        fuelsim::CaseInputReader::read(transient_input_path);
    const fuelsim::UnstructuredQuad4Mesh mesh =
        fuelsim::ExodusMeshIo::read_quad4(input.mesh_file);
    const fuelsim::SolverOptions solver{
        input.solver.absolute_tolerance, input.solver.relative_tolerance,
        input.solver.step_tolerance, input.solver.maximum_iterations};

    fuelsim::TransientProblem uninterrupted(input.transient_definition(), mesh);
    const fuelsim::TransientResult full =
        fuelsim::solve_transient(uninterrupted, time_options(20.0), solver);
    passed =
        check(full.completed, "uninterrupted PCMI solve completes") && passed;

    fuelsim::TransientProblem split(input.transient_definition(), mesh);
    fuelsim::ExodusTransientResultsWriter writer(results_path, mesh, split);
    writer.append(split);
    ResultsObserver observer(writer);
    const fuelsim::TransientResult first =
        fuelsim::solve_transient(split, time_options(10.0), solver, &observer);
    passed = check(first.completed &&
                       observer.steps() == first.accepted_steps.size(),
                   "first restart segment observes every accepted step") &&
             passed;
    const std::string history_path = results_path + ".history.csv";
    {
        fuelsim::EngineeringHistoryWriter history(history_path, split);
        history.append(split, 1.0, 1.0, first.last_attempt.nonlinear_iterations);
    }
    {
        std::ifstream history(history_path);
        std::string header;
        std::string values;
        std::getline(history, header);
        std::getline(history, values);
        passed =
            check(header.find("region_cladding_maximum_temperature") !=
                          std::string::npos &&
                      header.find("contact_fuel_cladding_minimum_gap") !=
                          std::string::npos &&
                      !values.empty(),
                  "engineering history uses named region/contact columns") &&
            passed;
    }
    std::filesystem::remove(history_path);
    fuelsim::TransientCheckpointIo::write(checkpoint_path, split,
                                          first.next_time_step);
    const fuelsim::TransientCommittedState split_state =
        split.committed_state();

    fuelsim::TransientProblem restarted(input.transient_definition(), mesh);
    const double restored_time_step =
        fuelsim::TransientCheckpointIo::restore(checkpoint_path, restarted);
    passed = check(restored_time_step == first.next_time_step,
                   "restart preserves the committed controller step") &&
             passed;
    passed =
        compare_committed_states(split_state, restarted.committed_state()) &&
        passed;
    fuelsim::TransientTimeOptions restart_options = time_options(20.0);
    restart_options.initial_time_step = restored_time_step;
    const fuelsim::TransientResult second = fuelsim::solve_transient(
        restarted, restart_options, solver);
    passed = check(second.completed, "restarted PCMI solve reaches end time") &&
             compare_committed_states(uninterrupted.committed_state(),
                                      restarted.committed_state()) &&
             passed;
    passed =
        verify_exodus(results_path, mesh, split, writer.step_count()) && passed;

    fuelsim::TransientProblem mismatch(input.transient_definition(), mesh);
    fuelsim::TransientProblemDefinition changed = input.transient_definition();
    changed.regions[1].material.plasticity.yield_stress *= 1.01;
    fuelsim::TransientProblem changed_problem(std::move(changed), mesh);
    passed = expect_failure(
                 [&]() {
                     (void)fuelsim::TransientCheckpointIo::restore(
                         checkpoint_path, changed_problem);
                 },
                 "signature", "checkpoint rejects a changed material model") &&
             passed;

    mismatch.begin_time_step({1.0, 0.05});
    passed = expect_failure(
                 [&]() {
                     fuelsim::TransientCheckpointIo::write(checkpoint_path,
                                                           mismatch, 1.0);
                 },
                 "active time step",
                 "checkpoint cannot capture uncommitted trial state") &&
             passed;
    mismatch.rollback_time_step();
    passed =
        expect_failure(
            [&]() {
                fuelsim::TransientCheckpointIo::write(checkpoint_path, mismatch,
                                                      0.0);
            },
            "next time step", "checkpoint rejects invalid controller state") &&
        passed;

    {
        std::fstream file(checkpoint_path,
                          std::ios::binary | std::ios::in | std::ios::out);
        if (!file)
            return check(false, "checkpoint can be opened for corruption test");
        file.seekg(48, std::ios::beg);
        char byte = 0;
        file.read(&byte, 1);
        file.clear();
        file.seekp(48, std::ios::beg);
        byte = static_cast<char>(byte ^ 0x5A);
        file.write(&byte, 1);
    }
    passed = expect_failure(
                 [&]() {
                     (void)fuelsim::TransientCheckpointIo::restore(
                         checkpoint_path, mismatch);
                 },
                 "checksum", "checkpoint detects payload corruption") &&
             passed;

    const int checkpoint_remove = std::remove(checkpoint_path.c_str());
    const int results_remove = std::remove(results_path.c_str());
    return check(checkpoint_remove == 0 && results_remove == 0,
                 "M3.0 test artifacts are removed") &&
           passed;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: checkpoint_results_tests <steady.fsi> "
                     "<transient.fsi> <checkpoint> <results.e>\n";
        return 2;
    }
    try {
        fuelsim::PetscSession session(
            argc, argv, "fuelsim M3.0 checkpoint/results tests\n");
        if (!run_tests(argv[1], argv[2], argv[3], argv[4]))
            return 1;
        std::cout << "[PASS] fuelsim M3.0 checkpoint and Exodus results\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] M3.0 tests raised: " << error.what() << '\n';
        return 1;
    }
}
