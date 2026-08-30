#include "fuelsim/solver/solve_workflows.hpp"
#include "support/cartesian3d_problem_access.hpp"
#include "support/material_factory.hpp"
#include "support/moose_field_comparison.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr double standard_time_step = 0.1;
constexpr std::array<double, 10> plastic_displacement_path = {
    0.0005, 0.002, 0.004, 0.002, 0.0, -0.002, -0.004, -0.001, 0.003, 0.0};
constexpr std::array<double, 10> coupled_displacement_path = {
    0.0015, 0.0020, 0.0025, 0.0030, 0.0035, 0.0040, 0.0045, 0.0050, 0.0055, 0.0060};
constexpr std::array<double, 20> noncoaxial_axial_path = {0.0012, 0.00165, 0.0021, 0.00255, 0.003, 0.003, 0.003, 0.003,
    0.003, 0.003, 0.002, 0.001, 0.0, -0.001, -0.002, -0.002, -0.002, -0.002, -0.002, -0.002};
constexpr std::array<double, 20> noncoaxial_shear_path = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0008, 0.0016, 0.0024, 0.0032,
    0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.0024, 0.0008, -0.0008, -0.0024, -0.004};

enum class Branch {
    plastic,
    coupled,
    noncoaxial,
    finite_plastic,
    finite_coupled,
    finite_noncoaxial,
    finite_reduced_plastic,
    finite_reduced_coupled
};

struct NodeReference final {
    std::size_t stage, node;
    double time;
    std::array<double, 8> fields;
};

struct IntegrationReference final {
    std::size_t stage;
    double time;
    fuelsim::CartesianPoint3 position;
    double temperature;
    fuelsim::SymmetricTensor3Values stress, strain, elastic_strain, plastic_strain, creep_strain;
    double equivalent_plastic_strain, equivalent_creep_strain, integration_volume;
};

struct EnergyReference final {
    std::size_t stage;
    double time, internal_energy, plastic_dissipation, creep_dissipation, external_work;
};

const char* case_id(Branch branch) {
    switch (branch) {
    case Branch::plastic: return "B5.10";
    case Branch::coupled: return "B5.12";
    case Branch::noncoaxial: return "B5.13";
    case Branch::finite_plastic: return "B5.15";
    case Branch::finite_coupled: return "B5.17";
    case Branch::finite_noncoaxial: return "B5.18";
    case Branch::finite_reduced_plastic: return "B5.35";
    case Branch::finite_reduced_coupled: return "B5.37";
    }
    throw std::logic_error("Unknown Abaqus inelastic branch");
}

const char* metric_prefix(Branch branch) {
    switch (branch) {
    case Branch::plastic: return "b510_";
    case Branch::coupled: return "b512_";
    case Branch::noncoaxial: return "b513_";
    case Branch::finite_plastic: return "b515_";
    case Branch::finite_coupled: return "b517_";
    case Branch::finite_noncoaxial: return "b518_";
    case Branch::finite_reduced_plastic: return "b535_";
    case Branch::finite_reduced_coupled: return "b537_";
    }
    throw std::logic_error("Unknown Abaqus inelastic branch");
}

bool finite_strain(Branch branch) {
    return branch == Branch::finite_plastic || branch == Branch::finite_coupled ||
           branch == Branch::finite_noncoaxial || branch == Branch::finite_reduced_plastic ||
           branch == Branch::finite_reduced_coupled;
}

bool reduced_integration(Branch branch) {
    return branch == Branch::finite_reduced_plastic || branch == Branch::finite_reduced_coupled;
}

bool plastic_only(Branch branch) {
    return branch == Branch::plastic || branch == Branch::finite_plastic || branch == Branch::finite_reduced_plastic;
}

bool noncoaxial(Branch branch) { return branch == Branch::noncoaxial || branch == Branch::finite_noncoaxial; }

bool monotonic_coupled(Branch branch) {
    return branch == Branch::coupled || branch == Branch::finite_coupled || branch == Branch::finite_reduced_coupled;
}

std::size_t stage_count(Branch branch) { return noncoaxial(branch) ? 20 : 10; }

double time_step(Branch branch) { return noncoaxial(branch) ? 0.001 : standard_time_step; }

Branch parse_branch(const std::string& value) {
    if (value == "plastic") return Branch::plastic;
    if (value == "coupled") return Branch::coupled;
    if (value == "noncoaxial") return Branch::noncoaxial;
    if (value == "finite_plastic") return Branch::finite_plastic;
    if (value == "finite_coupled") return Branch::finite_coupled;
    if (value == "finite_noncoaxial") return Branch::finite_noncoaxial;
    if (value == "finite_reduced_plastic") return Branch::finite_reduced_plastic;
    if (value == "finite_reduced_coupled") return Branch::finite_reduced_coupled;
    throw std::invalid_argument("Abaqus inelastic branch is not recognized");
}

struct StepSnapshot final {
    double time;
    std::vector<double> state;
    std::vector<fuelsim::CartesianMaterialPointState> material;
    fuelsim::TransientConservationSummary conservation;
};

class SnapshotObserver final : public fuelsim::TransientStepObserver {
  public:
    void accepted_step(const fuelsim::TransientProblem& problem, const fuelsim::TransientAcceptedStep& step) override {
        StepSnapshot snapshot{step.time, problem.committed_solution(), {}, step.conservation};
        const fuelsim::CartesianMaterialHistory& history =
            fuelsim::cartesian::ProblemAccess::material_history(problem, 0, 0);
        if (history.size() != 1 && history.size() != 8)
            throw std::invalid_argument("Abaqus inelastic material history does not contain one or eight points");
        snapshot.material.assign(history.begin(), history.end());
        _snapshots.push_back(std::move(snapshot));
    }

    const std::vector<StepSnapshot>& snapshots() const noexcept { return _snapshots; }

  private:
    std::vector<StepSnapshot> _snapshots;
};

bool check(bool condition, const std::string& message) {
    if (condition) return true;
    std::cerr << "[FAIL] " << message << '\n';
    return false;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string value;
    while (std::getline(stream, value, ',')) result.push_back(value);
    return result;
}

double number(const std::vector<std::string>& values, std::size_t index, const std::string& path) {
    if (index >= values.size()) throw std::invalid_argument("Incomplete Abaqus B5.10 row in " + path);
    return std::stod(values[index]);
}

std::size_t positive_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0)
        throw std::invalid_argument("Abaqus B5.10 index is not a positive integer in " + path);
    return static_cast<std::size_t>(rounded);
}

double zero_noise(double value, double threshold) { return std::abs(value) < threshold ? 0.0 : value; }

fuelsim::SymmetricTensor3Values tensor(
    const std::vector<std::string>& values, std::size_t first, const std::string& path, double zero_threshold) {
    fuelsim::SymmetricTensor3Values result = {number(values, first, path), number(values, first + 1, path),
        number(values, first + 2, path), 0.5 * number(values, first + 3, path), 0.5 * number(values, first + 5, path),
        0.5 * number(values, first + 4, path)};
    for (double* value : {&result.xx, &result.yy, &result.zz, &result.xy, &result.yz, &result.xz})
        *value = zero_noise(*value, zero_threshold);
    return result;
}

std::vector<NodeReference> read_nodes(const std::string& path, Branch branch) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.10 nodes: " + path);
    std::string line;
    std::getline(input, line);
    if (line != "stage,time_s,node,temperature_k,ux_m,uy_m,uz_m,reaction_heat_flux_w,rfx_n,rfy_n,rfz_n")
        throw std::invalid_argument("Unexpected Abaqus B5.10 nodal header in " + path);
    std::vector<NodeReference> result;
    double maximum_active_reaction = 0.0;
    double maximum_raw_theoretical_zero_reaction = 0.0;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 11) throw std::invalid_argument("Unexpected Abaqus B5.10 nodal column count in " + path);
        std::array<double, 8> fields{};
        for (std::size_t field = 0; field < fields.size(); ++field) {
            fields[field] = number(values, field + 3, path);
            if (field == 5 || field == 6)
                maximum_active_reaction = std::max(maximum_active_reaction, std::abs(fields[field]));
            const bool theoretical_zero_reaction =
                (noncoaxial(branch) && field == 7) || (!noncoaxial(branch) && (field == 6 || field == 7));
            if (theoretical_zero_reaction) {
                maximum_raw_theoretical_zero_reaction =
                    std::max(maximum_raw_theoretical_zero_reaction, std::abs(fields[field]));
                fields[field] = 0.0;
            } else if (field != 0)
                fields[field] = zero_noise(fields[field], field >= 5 ? 1.0e-3 : (field == 4 ? 1.0e-10 : 1.0e-16));
        }
        result.push_back({positive_integer(number(values, 0, path), path),
            positive_integer(number(values, 2, path), path), number(values, 1, path), fields});
    }
    if (result.size() != 8 * stage_count(branch))
        throw std::invalid_argument("Abaqus small-strain nodal reference has an unexpected row count");
    {
        const double relative_residual = maximum_raw_theoretical_zero_reaction / maximum_active_reaction;
        std::cout << metric_prefix(branch)
                  << "abaqus_raw_theoretical_zero_reaction_maximum_absolute=" << maximum_raw_theoretical_zero_reaction
                  << '\n';
        std::cout << metric_prefix(branch)
                  << "abaqus_raw_theoretical_zero_reaction_relative_active_scale=" << relative_residual << '\n';
        if (!(maximum_active_reaction > 0.0 && relative_residual < 1.0e-7))
            throw std::invalid_argument(
                "Abaqus theoretical-zero reaction residual is too large to canonicalize as zero");
    }
    return result;
}

std::vector<IntegrationReference> read_integration(const std::string& path, Branch branch) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.10 integration points: " + path);
    std::string line;
    std::getline(input, line);
    const std::string strain_header = finite_strain(branch)
                                          ? "le11,le22,le33,le12_engineering,le13_engineering,le23_engineering,"
                                          : "e11,e22,e33,e12_engineering,e13_engineering,e23_engineering,";
    const std::string common = "stage,time_s,element,integration_point,x_m,y_m,z_m,temperature_k,"
                               "s11_pa,s22_pa,s33_pa,s12_pa,s13_pa,s23_pa," +
                               strain_header +
                               "ee11,ee22,ee33,ee12_engineering,ee13_engineering,ee23_engineering,"
                               "pe11,pe22,pe33,pe12_engineering,pe13_engineering,pe23_engineering,peeq,";
    const std::string expected_header =
        common + (plastic_only(branch)
                         ? "ivol_m3"
                         : "ce11,ce22,ce33,ce12_engineering,ce13_engineering,ce23_engineering,ceeq,ivol_m3");
    if (line != expected_header)
        throw std::invalid_argument("Unexpected Abaqus B5.10 integration-point header in " + path);
    std::vector<IntegrationReference> result;
    double maximum_active_stress = 0.0;
    double maximum_raw_free_z_stress = 0.0;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != (plastic_only(branch) ? 34 : 41))
            throw std::invalid_argument("Unexpected Abaqus B5.10 integration-point column count in " + path);
        if (positive_integer(number(values, 2, path), path) != 1)
            throw std::invalid_argument("Abaqus B5.10 reference must contain one element");
        fuelsim::SymmetricTensor3Values stress = {number(values, 8, path), number(values, 9, path),
            number(values, 10, path), number(values, 11, path), number(values, 13, path), number(values, 12, path)};
        maximum_active_stress =
            std::max({maximum_active_stress, std::abs(stress.xx), std::abs(stress.yy), std::abs(stress.xy)});
        if (noncoaxial(branch)) maximum_raw_free_z_stress = std::max(maximum_raw_free_z_stress, std::abs(stress.zz));
        for (double* value : {&stress.xx, &stress.yy, &stress.zz, &stress.xy, &stress.yz, &stress.xz})
            *value = zero_noise(*value, 1.0);
        if (noncoaxial(branch)) stress.zz = 0.0;
        const fuelsim::SymmetricTensor3Values creep =
            plastic_only(branch) ? fuelsim::SymmetricTensor3Values{} : tensor(values, 33, path, 1.0e-15);
        result.push_back({positive_integer(number(values, 0, path), path), number(values, 1, path),
            {number(values, 4, path), number(values, 5, path), number(values, 6, path)}, number(values, 7, path),
            stress, tensor(values, 14, path, 1.0e-15), tensor(values, 20, path, 1.0e-15),
            tensor(values, 26, path, 1.0e-15), creep, number(values, 32, path),
            plastic_only(branch) ? 0.0 : number(values, 39, path),
            number(values, plastic_only(branch) ? 33 : 40, path)});
    }
    if (result.size() != (reduced_integration(branch) ? 1 : 8) * stage_count(branch))
        throw std::invalid_argument("Abaqus small-strain integration-point reference has an unexpected row count");
    if (noncoaxial(branch)) {
        const double relative_residual = maximum_raw_free_z_stress / maximum_active_stress;
        std::cout << metric_prefix(branch) << "abaqus_raw_free_z_stress_maximum_absolute=" << maximum_raw_free_z_stress
                  << '\n';
        std::cout << metric_prefix(branch) << "abaqus_raw_free_z_stress_relative_active_scale=" << relative_residual
                  << '\n';
        if (!(maximum_active_stress > 0.0 && relative_residual < 1.0e-7))
            throw std::invalid_argument("Abaqus B5.13 free-z stress residual is too large to canonicalize as zero");
    }
    return result;
}

std::vector<EnergyReference> read_energy(const std::string& path, Branch branch) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Could not read Abaqus B5.10 energy: " + path);
    std::string line;
    std::getline(input, line);
    const std::string expected_header =
        plastic_only(branch)
            ? "stage,time_s,internal_energy_j,plastic_dissipation_j,external_work_j"
            : "stage,time_s,internal_energy_j,plastic_dissipation_j,creep_dissipation_j,external_work_j";
    if (line != expected_header) throw std::invalid_argument("Unexpected Abaqus B5.10 energy header in " + path);
    std::vector<EnergyReference> result;
    while (std::getline(input, line)) {
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != (plastic_only(branch) ? 5 : 6))
            throw std::invalid_argument("Unexpected Abaqus B5.10 energy column count in " + path);
        result.push_back({positive_integer(number(values, 0, path), path), number(values, 1, path),
            number(values, 2, path), number(values, 3, path), plastic_only(branch) ? 0.0 : number(values, 4, path),
            number(values, plastic_only(branch) ? 4 : 5, path)});
    }
    if (result.size() != stage_count(branch))
        throw std::invalid_argument("Abaqus small-strain energy reference has an unexpected row count");
    return result;
}

fuelsim::UnstructuredHex8Mesh mesh() {
    return fuelsim::UnstructuredHex8Mesh({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
                                             {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}},
        {{{{0, 1, 2, 3, 4, 5, 6, 7}}}}, {1}, {{1, "solid"}}, {},
        {{11, "left", {{0, 3}}}, {12, "right", {{0, 1}}}, {13, "y0", {{0, 0}}}, {14, "z0", {{0, 4}}}});
}

fuelsim::ThermoelasticProperties material(Branch branch) {
    fuelsim::ThermoelasticProperties result = fuelsim::test::with_plasticity(
        fuelsim::test::thermoelastic(0.0, 1.0, 2.0e11, 0.3, 0.0, 600.0, 0.0, 0.0, 0.0, 1.0, 1.0), 2.0e8, 2.0e9, 600.0);
    if (!plastic_only(branch))
        result = fuelsim::test::with_norton(
            std::move(result), monotonic_coupled(branch) ? 1.0e-4 : 1.0e-5, 1.0e8, 3.0, 600.0);
    return result;
}

fuelsim::SpatialDefinition definition(Branch branch) {
    fuelsim::SpatialDefinition result;
    result.regions.push_back({"solid", "solid", material(branch), 0.0, 600.0, -1, "",
        finite_strain(branch) ? fuelsim::StrainFormulation::finite : fuelsim::StrainFormulation::small});
    if (reduced_integration(branch))
        result.regions.back().hex8_element_formulation = fuelsim::Hex8ElementFormulation::c3d8rt;
    if (noncoaxial(branch)) {
        std::vector<double> times{0.0}, axial{0.0}, shear{0.0};
        for (std::size_t stage = 0; stage < stage_count(branch); ++stage) {
            times.push_back(time_step(branch) * static_cast<double>(stage + 1));
            axial.push_back(noncoaxial_axial_path[stage]);
            shear.push_back(noncoaxial_shear_path[stage]);
        }
        result.time_tables.emplace_back("right_axial", times, std::move(axial));
        result.time_tables.emplace_back("right_shear", std::move(times), std::move(shear));
        result.boundary_conditions = {
            {"temperature_left", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::temperature, 600.0},
            {"temperature_right", fuelsim::BoundaryConditionType::dirichlet, "right", fuelsim::Field::temperature,
                600.0},
            {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::displacement_x, 0.0},
            {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::displacement_y, 0.0},
            {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0},
            {"right_x", fuelsim::BoundaryConditionType::dirichlet, "right", fuelsim::Field::displacement_x, 1.0, false,
                "right_axial"},
            {"right_y", fuelsim::BoundaryConditionType::dirichlet, "right", fuelsim::Field::displacement_y, 1.0, false,
                "right_shear"},
        };
        return result;
    }
    const std::array<double, 10>& displacement_path =
        plastic_only(branch) ? plastic_displacement_path : coupled_displacement_path;
    std::vector<double> times{0.0}, displacements{0.0};
    for (std::size_t stage = 0; stage < displacement_path.size(); ++stage) {
        times.push_back(time_step(branch) * static_cast<double>(stage + 1));
        displacements.push_back(displacement_path[stage]);
    }
    result.time_tables.emplace_back("right_displacement", std::move(times), std::move(displacements));
    result.boundary_conditions = {
        {"temperature", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::temperature, 600.0},
        {"fix_x", fuelsim::BoundaryConditionType::dirichlet, "left", fuelsim::Field::displacement_x, 0.0},
        {"fix_y", fuelsim::BoundaryConditionType::dirichlet, "y0", fuelsim::Field::displacement_y, 0.0},
        {"fix_z", fuelsim::BoundaryConditionType::dirichlet, "z0", fuelsim::Field::displacement_z, 0.0},
        {"pull_x", fuelsim::BoundaryConditionType::dirichlet, "right", fuelsim::Field::displacement_x, 1.0, false,
            "right_displacement"},
    };
    return result;
}

fuelsim::SolverOptions solver_options() {
    fuelsim::SolverOptions options;
    options.absolute_tolerance = 1.0e-10;
    options.relative_tolerance = 1.0e-12;
    options.maximum_iterations = 40;
    options.linear_solver = fuelsim::SolverOptions::LinearSolver::direct;
    options.preconditioner = fuelsim::SolverOptions::Preconditioner::lu;
    options.field_residual_scaling = true;
    options.residual_reduction_tolerance = 1.0e-11;
    options.temperature_residual_absolute_tolerance = 1.0e-10;
    options.mechanical_residual_absolute_tolerance = 1.0e-2;
    return options;
}

std::size_t stage_from_time(double time, Branch branch) {
    const double value = time / time_step(branch);
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 1.0 || rounded > static_cast<double>(stage_count(branch)))
        throw std::invalid_argument("Fuelsim small-strain time does not identify a prescribed stage");
    return static_cast<std::size_t>(rounded);
}

std::vector<double> raw_residual(fuelsim::TransientProblem& problem, const std::vector<double>& state) {
    std::vector<double> result(problem.dof_count(), 0.0);
    fuelsim::ContributionWorkspace workspace;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        problem.evaluate_contribution(contribution, state, workspace, false);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            result[workspace.dofs[local]] += workspace.residual[local];
    }
    return result;
}

bool metrics_pass(const fuelsim::test::FieldErrorMetrics& metrics, double relative_tolerance, double zero_tolerance) {
    if (metrics.has_relative_norm() && !fuelsim::test::relative_metrics_below(metrics, relative_tolerance))
        return false;
    return metrics.maximum_zero_reference_difference < zero_tolerance;
}

void print_metrics(const std::string& name, const fuelsim::test::FieldErrorMetrics& metrics) {
    if (metrics.has_relative_norm())
        fuelsim::test::print_relative_metrics(name, metrics);
    else
        fuelsim::test::print_absolute_metrics(name, metrics);
}

std::array<double, 6> components(const fuelsim::SymmetricTensor3Values& value) {
    return {value.xx, value.yy, value.zz, value.xy, value.yz, value.xz};
}

fuelsim::SymmetricTensor3Values logarithmic_strain(
    const fuelsim::Hex8QuadraturePoint& point, const fuelsim::Hex8LocalValues& state) {
    double f00 = 1.0, f01 = 0.0, f10 = 0.0, f11 = 1.0, f22 = 1.0;
    for (std::size_t node = 0; node < 8; ++node) {
        f00 += state[8 + node] * point.gradient[node][0];
        f01 += state[8 + node] * point.gradient[node][1];
        f10 += state[16 + node] * point.gradient[node][0];
        f11 += state[16 + node] * point.gradient[node][1];
        f22 += state[24 + node] * point.gradient[node][2];
    }
    const double b00 = f00 * f00 + f01 * f01;
    const double b01 = f00 * f10 + f01 * f11;
    const double b11 = f10 * f10 + f11 * f11;
    const double mean = 0.5 * (b00 + b11);
    const double radius = std::hypot(0.5 * (b00 - b11), b01);
    const double first = mean + radius, second = mean - radius;
    if (!(first > 0.0) || !(second > 0.0) || !(f22 > 0.0))
        throw std::domain_error("Finite-strain Abaqus comparison requires a positive left stretch tensor");
    if (radius <= 1.0e-14 * mean) return {0.5 * std::log(mean), 0.5 * std::log(mean), std::log(f22), 0.0, 0.0, 0.0};
    const double beta = (std::log(first) - std::log(second)) / (first - second);
    const double alpha = (first * std::log(second) - second * std::log(first)) / (first - second);
    return {0.5 * (alpha + beta * b00), 0.5 * (alpha + beta * b11), std::log(f22), 0.5 * beta * b01, 0.0, 0.0};
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: fuelsim_b510_hex8_c3d8t_small_j2_abaqus_tests "
                     "<plastic|coupled|noncoaxial|finite_plastic|finite_coupled|finite_noncoaxial|"
                     "finite_reduced_plastic|finite_reduced_coupled> "
                     "<nodes.csv> <integration.csv> <energy.csv>\n";
        return 2;
    }
    try {
        std::cout << std::scientific << std::setprecision(12);
        fuelsim::PetscSession session(argc, argv, "fuelsim Abaqus inelastic comparison\n");
        const Branch branch = parse_branch(argv[1]);
        const std::vector<NodeReference> node_reference = read_nodes(argv[2], branch);
        const std::vector<IntegrationReference> integration_reference = read_integration(argv[3], branch);
        const std::vector<EnergyReference> energy_reference = read_energy(argv[4], branch);
        const fuelsim::UnstructuredHex8Mesh input_mesh = mesh();
        fuelsim::TransientProblem problem(definition(branch), input_mesh);
        SnapshotObserver observer;
        const double step = time_step(branch);
        const std::size_t stages = stage_count(branch);
        const fuelsim::TransientResult solve = fuelsim::solve_transient(problem,
            {step * static_cast<double>(stages), step, step, step, 1.0, 0.5, 0, 0.0}, solver_options(), &observer);
        bool passed = check(solve.completed && observer.snapshots().size() == stages &&
                                solve.accepted_steps.size() == stages && solve.rejected_steps.empty(),
            std::string(case_id(branch)) + " accepts every prescribed inelastic increment without cutback");
        std::map<std::size_t, const StepSnapshot*> snapshots;
        for (const StepSnapshot& snapshot : observer.snapshots())
            snapshots.emplace(stage_from_time(snapshot.time, branch), &snapshot);

        std::array<fuelsim::test::FieldErrorMetrics, 8> nodal_metrics;
        fuelsim::TransientProblem reaction_problem(definition(branch), input_mesh);
        const auto& dofs = fuelsim::cartesian::ProblemAccess::dof_map(reaction_problem);
        const std::array<fuelsim::Field, 4> fields = {fuelsim::Field::temperature, fuelsim::Field::displacement_x,
            fuelsim::Field::displacement_y, fuelsim::Field::displacement_z};
        for (std::size_t stage = 1; stage <= stages; ++stage) {
            const StepSnapshot& snapshot = *snapshots.at(stage);
            reaction_problem.begin_time_step({snapshot.time, 1.0, true});
            const std::vector<double> reaction = raw_residual(reaction_problem, snapshot.state);
            for (const NodeReference& reference : node_reference) {
                if (reference.stage != stage) continue;
                if (reference.node < 1 || reference.node > 8 || std::abs(reference.time - snapshot.time) > 1.0e-12)
                    throw std::invalid_argument("Abaqus B5.10 nodal label or time lies outside the path");
                const std::size_t node = reference.node - 1;
                for (std::size_t field = 0; field < fields.size(); ++field) {
                    const std::size_t dof = dofs.dof(fields[field], node);
                    nodal_metrics[field].add(snapshot.state[dof], reference.fields[field]);
                    nodal_metrics[4 + field].add(reaction[dof], reference.fields[4 + field]);
                }
            }
            reaction_problem.commit_time_step(snapshot.state);
            passed = check(snapshot.conservation.plastic_dissipation_increment >= -1.0e-8 &&
                               snapshot.conservation.creep_dissipation_increment >= -1.0e-8,
                         std::string(case_id(branch)) + " inelastic dissipation is nonnegative at stage " +
                             std::to_string(stage)) &&
                     passed;
        }
        const std::array<std::string, 8> nodal_names = {"temperature", "displacement_x", "displacement_y",
            "displacement_z", "reaction_heat_flux", "reaction_force_x", "reaction_force_y", "reaction_force_z"};
        const std::array<double, 8> nodal_zero_tolerances = {
            1.0e-12, 1.0e-15, 1.0e-15, 1.0e-15, 1.0e-10, 1.0e-2, 1.0e-2, 1.0e-2};
        for (std::size_t field = 0; field < nodal_metrics.size(); ++field) {
            print_metrics(metric_prefix(branch) + nodal_names[field], nodal_metrics[field]);
            passed = check(metrics_pass(nodal_metrics[field], 1.0e-3, nodal_zero_tolerances[field]),
                         std::string(case_id(branch)) + " " + nodal_names[field] +
                             " metrics are below the acceptance limits") &&
                     passed;
        }

        std::array<fuelsim::test::FieldErrorMetrics, 34> integration_metrics;
        double maximum_coordinate_error = 0.0;
        const fuelsim::Hex8Geometry geometry = fuelsim::make_hex8_geometry({{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
            {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}}});
        for (const IntegrationReference& reference : integration_reference) {
            const StepSnapshot& snapshot = *snapshots.at(reference.stage);
            fuelsim::Hex8LocalAdValues passive{};
            fuelsim::Hex8LocalValues passive_values{};
            for (std::size_t node = 0; node < 8; ++node) {
                passive[node] = snapshot.state[dofs.dof(fuelsim::Field::temperature, node)];
                passive[8 + node] = snapshot.state[dofs.dof(fuelsim::Field::displacement_x, node)];
                passive[16 + node] = snapshot.state[dofs.dof(fuelsim::Field::displacement_y, node)];
                passive[24 + node] = snapshot.state[dofs.dof(fuelsim::Field::displacement_z, node)];
            }
            for (std::size_t local = 0; local < passive.size(); ++local) passive_values[local] = passive[local].value();
            std::size_t closest = 0;
            double closest_squared = std::numeric_limits<double>::max();
            const std::size_t point_count = reduced_integration(branch) ? 1 : 8;
            for (std::size_t q = 0; q < point_count; ++q) {
                const fuelsim::Hex8QuadraturePoint& quadrature_point =
                    reduced_integration(branch) ? geometry.reduced_point : geometry.points[q];
                fuelsim::CartesianPoint3 point = quadrature_point.position;
                if (finite_strain(branch))
                    for (std::size_t node = 0; node < 8; ++node) {
                        point.x += quadrature_point.shape[node] * passive[8 + node].value();
                        point.y += quadrature_point.shape[node] * passive[16 + node].value();
                        point.z += quadrature_point.shape[node] * passive[24 + node].value();
                    }
                const double distance_squared = std::pow(point.x - reference.position.x, 2) +
                                                std::pow(point.y - reference.position.y, 2) +
                                                std::pow(point.z - reference.position.z, 2);
                if (distance_squared < closest_squared) {
                    closest = q;
                    closest_squared = distance_squared;
                }
            }
            maximum_coordinate_error = std::max(maximum_coordinate_error, std::sqrt(closest_squared));
            const fuelsim::CartesianMaterialPointState& actual = snapshot.material[closest];
            const fuelsim::Hex8QuadraturePoint& material_point =
                reduced_integration(branch) ? geometry.reduced_point : geometry.points[closest];
            const std::array<double, 6> actual_stress = components(actual.stress);
            const std::array<double, 6> expected_stress = components(reference.stress);
            std::array<double, 6> actual_total{};
            if (finite_strain(branch))
                actual_total = components(logarithmic_strain(material_point, passive_values));
            else {
                double volume = 0.0, average_trace = 0.0;
                for (const fuelsim::Hex8QuadraturePoint& point : geometry.points) {
                    const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
                        point, passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::small);
                    volume += point.weighted_measure;
                    average_trace += point.weighted_measure *
                                     (kinematics.strain_increment.xx.value() + kinematics.strain_increment.yy.value() +
                                         kinematics.strain_increment.zz.value());
                }
                average_trace /= volume;
                const fuelsim::CartesianKinematics kinematics = fuelsim::evaluate_cartesian_incremental_kinematics(
                    geometry.points[closest], passive, fuelsim::Hex8LocalValues{}, fuelsim::StrainFormulation::small);
                fuelsim::SymmetricTensor3 actual_strain = kinematics.strain_increment;
                const adlite::Scalar correction =
                    (average_trace - actual_strain.xx - actual_strain.yy - actual_strain.zz) / 3.0;
                actual_strain.xx += correction;
                actual_strain.yy += correction;
                actual_strain.zz += correction;
                actual_total = {actual_strain.xx.value(), actual_strain.yy.value(), actual_strain.zz.value(),
                    actual_strain.xy.value(), actual_strain.yz.value(), actual_strain.xz.value()};
            }
            const std::array<double, 6> expected_total = components(reference.strain);
            const std::array<double, 6> expected_elastic = components(reference.elastic_strain);
            const std::array<double, 6> expected_plastic = components(reference.plastic_strain);
            const std::array<double, 6> expected_creep = components(reference.creep_strain);
            for (std::size_t component = 0; component < 6; ++component) {
                integration_metrics[component].add(actual_stress[component], expected_stress[component]);
                integration_metrics[6 + component].add(actual_total[component], expected_total[component]);
                integration_metrics[12 + component].add(actual.elastic_strain[component], expected_elastic[component]);
                integration_metrics[18 + component].add(actual.plastic_strain[component], expected_plastic[component]);
                integration_metrics[24 + component].add(actual.creep_strain[component], expected_creep[component]);
            }
            integration_metrics[30].add(actual.equivalent_plastic_strain, reference.equivalent_plastic_strain);
            integration_metrics[31].add(actual.equivalent_creep_strain, reference.equivalent_creep_strain);
            integration_metrics[32].add(passive[0].value(), reference.temperature);
            const double actual_measure = finite_strain(branch)
                                              ? fuelsim::evaluate_cartesian_incremental_kinematics(material_point,
                                                    passive, passive_values, fuelsim::StrainFormulation::finite)
                                                    .current_weighted_measure.value()
                                              : material_point.weighted_measure;
            integration_metrics[33].add(actual_measure, reference.integration_volume);
        }
        const std::array<std::string, 34> integration_names = {"stress_xx", "stress_yy", "stress_zz", "stress_xy",
            "stress_yz", "stress_xz", "strain_xx", "strain_yy", "strain_zz", "strain_xy", "strain_yz", "strain_xz",
            "elastic_strain_xx", "elastic_strain_yy", "elastic_strain_zz", "elastic_strain_xy", "elastic_strain_yz",
            "elastic_strain_xz", "plastic_strain_xx", "plastic_strain_yy", "plastic_strain_zz", "plastic_strain_xy",
            "plastic_strain_yz", "plastic_strain_xz", "creep_strain_xx", "creep_strain_yy", "creep_strain_zz",
            "creep_strain_xy", "creep_strain_yz", "creep_strain_xz", "equivalent_plastic_strain",
            "equivalent_creep_strain", "material_temperature", "integration_volume"};
        for (std::size_t field = 0; field < integration_metrics.size(); ++field) {
            print_metrics(metric_prefix(branch) + integration_names[field], integration_metrics[field]);
            const double zero_tolerance = field < 6 ? 1.0e-2 : 1.0e-14;
            passed =
                check(metrics_pass(integration_metrics[field], 1.0e-3, zero_tolerance),
                    std::string(case_id(branch)) + " " + integration_names[field] + " metrics are below 0.1 percent") &&
                passed;
        }
        passed = check(maximum_coordinate_error < 1.0e-11,
                     std::string(case_id(branch)) + " maps every integration point at every accepted time") &&
                 passed;
        std::cout << metric_prefix(branch)
                  << "maximum_integration_coordinate_absolute_difference=" << maximum_coordinate_error << '\n';

        std::array<fuelsim::test::FieldErrorMetrics, 4> energy_metrics;
        double cumulative_elastic = 0.0, cumulative_plastic = 0.0, cumulative_creep = 0.0;
        for (const EnergyReference& reference : energy_reference) {
            const StepSnapshot& snapshot = *snapshots.at(reference.stage);
            cumulative_elastic += snapshot.conservation.elastic_energy_change;
            cumulative_plastic += snapshot.conservation.plastic_dissipation_increment;
            cumulative_creep += snapshot.conservation.creep_dissipation_increment;
            energy_metrics[0].add(
                cumulative_elastic + cumulative_plastic + cumulative_creep, reference.internal_energy);
            energy_metrics[1].add(cumulative_plastic, reference.plastic_dissipation);
            energy_metrics[2].add(cumulative_creep, reference.creep_dissipation);
            energy_metrics[3].add(cumulative_elastic + cumulative_plastic + cumulative_creep, reference.external_work);
        }
        const std::array<std::string, 4> energy_names = {
            "internal_energy", "plastic_dissipation", "creep_dissipation", "external_work"};
        for (std::size_t field = 0; field < energy_metrics.size(); ++field) {
            print_metrics(metric_prefix(branch) + energy_names[field], energy_metrics[field]);
            passed = check(metrics_pass(energy_metrics[field], 1.0e-3, 1.0e-8),
                         std::string(case_id(branch)) + " " + energy_names[field] + " metrics are below 0.1 percent") &&
                     passed;
        }
        if (plastic_only(branch))
            passed =
                check(snapshots.at(3)->material[0].equivalent_plastic_strain > 0.0 &&
                          snapshots.at(4)->material[0].equivalent_plastic_strain ==
                              snapshots.at(3)->material[0].equivalent_plastic_strain &&
                          snapshots.at(7)->material[0].equivalent_plastic_strain >
                              snapshots.at(4)->material[0].equivalent_plastic_strain &&
                          snapshots.at(9)->material[0].equivalent_plastic_strain >
                              snapshots.at(7)->material[0].equivalent_plastic_strain,
                    std::string(case_id(branch)) +
                        " activates plasticity, preserves it during elastic unload, and accumulates it on reversal") &&
                passed;
        else if (monotonic_coupled(branch))
            passed = check(snapshots.at(1)->material[0].equivalent_plastic_strain > 0.0 &&
                               snapshots.at(1)->material[0].equivalent_creep_strain > 0.0 &&
                               snapshots.at(10)->material[0].equivalent_plastic_strain >
                                   snapshots.at(1)->material[0].equivalent_plastic_strain &&
                               snapshots.at(10)->material[0].equivalent_creep_strain >
                                   snapshots.at(1)->material[0].equivalent_creep_strain,
                         std::string(case_id(branch)) +
                             " keeps plasticity and creep active and accumulating throughout the path") &&
                     passed;
        else
            passed = check(snapshots.at(1)->material[0].equivalent_plastic_strain > 0.0 &&
                               snapshots.at(1)->material[0].equivalent_creep_strain > 0.0 &&
                               snapshots.at(stages)->material[0].equivalent_plastic_strain >
                                   snapshots.at(1)->material[0].equivalent_plastic_strain &&
                               snapshots.at(stages)->material[0].equivalent_creep_strain >
                                   snapshots.at(1)->material[0].equivalent_creep_strain &&
                               snapshots.at(10)->material[0].plastic_strain[3] *
                                       snapshots.at(stages)->material[0].plastic_strain[3] <
                                   0.0,
                         std::string(case_id(branch)) +
                             " accumulates both mechanisms and reverses the plastic shear direction") &&
                     passed;
        if (passed && session.rank() == 0) {
            std::cout << "[PASS] " << case_id(branch) << " Abaqus "
                      << (reduced_integration(branch) ? "C3D8RT " : "C3D8T ")
                      << (finite_strain(branch) ? "finite-strain" : "small-strain") << ' '
                      << (plastic_only(branch)
                                 ? "J2 load/unload/reversal"
                                 : (monotonic_coupled(branch) ? "coupled plastic-creep" : "noncoaxial reversal"))
                      << " comparison\n";
        }
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "[FAIL] Abaqus inelastic comparison raised: " << error.what() << '\n';
        return 1;
    }
}
