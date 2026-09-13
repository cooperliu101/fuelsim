#include "support/abaqus_hex8_full_field.hpp"
#include "c3d8rt.hpp"
#include "c3d8t.hpp"
#include "quad4_face.hpp"
#include "support/c3d_recovery.hpp"
#include "support/cartesian3d_problem_access.hpp"
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
#include <utility>

namespace fuelsim::test {
namespace {
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

std::size_t nonnegative_integer(double value, const std::string& path) {
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1.0e-12 || rounded < 0.0)
        throw std::invalid_argument("Abaqus HEX8 full-field index is not a nonnegative integer in " + path);
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

std::vector<ContactReference> read_contact(const std::string& path) {
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
    std::vector<ContactReference> result;
    while (std::getline(input, line)) {
        if (line.empty())
            continue;
        const std::vector<std::string> values = split_csv(line);
        if (values.size() != 26)
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
            reference.tangent_first[component] = number(values, 19 + component, path);
            reference.tangent_second[component] = number(values, 22 + component, path);
        }
        reference.heat_flux = number(values, 16, path);
        reference.shear_traction_first = number(values, 17, path);
        reference.shear_traction_second = number(values, 18, path);
        reference.state = nonnegative_integer(number(values, 25, path), path);
        result.push_back(reference);
    }
    return result;
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

SymmetricTensor3Values logarithmic_strain(const Hex8QuadraturePoint& point, const Hex8LocalValues& state) {
    Matrix3 deformation = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    for (std::size_t component = 0; component < 3; ++component)
        for (std::size_t direction = 0; direction < 3; ++direction)
            for (std::size_t node = 0; node < 8; ++node)
                deformation[component][direction] +=
                    state[8 * (component + 1) + node] * point.gradient[node][direction];
    Matrix3 left{};
    for (std::size_t row = 0; row < 3; ++row)
        for (std::size_t column = 0; column < 3; ++column)
            for (std::size_t inner = 0; inner < 3; ++inner)
                left[row][column] += deformation[row][inner] * deformation[column][inner];
    Matrix3 vectors = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    for (std::size_t sweep = 0; sweep < 32; ++sweep) {
        std::size_t p = 0, q = 1;
        double largest = std::abs(left[p][q]);
        for (std::size_t row = 0; row < 3; ++row)
            for (std::size_t column = row + 1; column < 3; ++column)
                if (std::abs(left[row][column]) > largest) {
                    p = row;
                    q = column;
                    largest = std::abs(left[row][column]);
                }
        const double scale = std::max({1.0, std::abs(left[0][0]), std::abs(left[1][1]), std::abs(left[2][2])});
        if (largest <= 1.0e-15 * scale)
            break;
        const double tau = (left[q][q] - left[p][p]) / (2.0 * left[p][q]);
        const double tangent = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent), sine = tangent * cosine;
        const double app = left[p][p], aqq = left[q][q], apq = left[p][q];
        left[p][p] = app - tangent * apq;
        left[q][q] = aqq + tangent * apq;
        left[p][q] = left[q][p] = 0.0;
        for (std::size_t row = 0; row < 3; ++row) {
            if (row == p || row == q)
                continue;
            const double arp = left[row][p], arq = left[row][q];
            left[row][p] = left[p][row] = cosine * arp - sine * arq;
            left[row][q] = left[q][row] = sine * arp + cosine * arq;
        }
        for (std::size_t row = 0; row < 3; ++row) {
            const double vrp = vectors[row][p], vrq = vectors[row][q];
            vectors[row][p] = cosine * vrp - sine * vrq;
            vectors[row][q] = sine * vrp + cosine * vrq;
        }
    }
    Matrix3 result{};
    for (std::size_t mode = 0; mode < 3; ++mode) {
        if (!(left[mode][mode] > 0.0))
            throw std::domain_error("Abaqus HEX8 logarithmic strain requires a positive left stretch tensor");
        const double value = 0.5 * std::log(left[mode][mode]);
        for (std::size_t row = 0; row < 3; ++row)
            for (std::size_t column = 0; column < 3; ++column)
                result[row][column] += value * vectors[row][mode] * vectors[column][mode];
    }
    return {result[0][0], result[1][1], result[2][2], result[0][1], result[1][2], result[0][2]};
}

std::vector<double> raw_residual(TransientProblem& problem, const std::vector<double>& state) {
    std::vector<double> result(problem.dof_count(), 0.0);
    ContributionWorkspace workspace;
    for (std::size_t contribution = 0; contribution < problem.contribution_count(); ++contribution) {
        problem.evaluate_contribution(contribution, state, workspace, false);
        for (std::size_t local = 0; local < workspace.dofs.size(); ++local)
            result[workspace.dofs[local]] += workspace.residual[local];
    }
    return result;
}

std::vector<double> thermal_contact_residual(const cartesian::SpatialAssembly& spatial,
    const std::vector<double>& state,
    double& conservation_error) {
    std::vector<double> result(state.size(), 0.0);
    for (std::size_t contribution = 0; contribution < spatial.contribution_count(); ++contribution) {
        if (spatial.contribution_type(contribution) != SpatialContributionType::thermal_contact)
            continue;
        std::vector<std::size_t> dofs;
        spatial.contribution_dofs(contribution, dofs);
        std::vector<double> local(dofs.size());
        for (std::size_t index = 0; index < dofs.size(); ++index)
            local[index] = state[dofs[index]];
        std::vector<double> residual;
        spatial.compute_contribution(contribution, local, nullptr, nullptr, 0.0, residual, nullptr);
        double balance = 0.0;
        for (std::size_t row = 0; row < dofs.size(); ++row) {
            result[dofs[row]] += residual[row];
            if (row < 8)
                balance += residual[row];
        }
        conservation_error = std::max(conservation_error, std::abs(balance));
    }
    return result;
}

bool metrics_pass(const FieldErrorMetrics& metrics,
    double aggregate_tolerance,
    double pointwise_tolerance,
    double zero_tolerance,
    double qualified_pointwise_absolute_tolerance = 0.0) {
    if (metrics.has_relative_norm()) {
        const bool aggregate_passed =
            metrics.relative_l2() < aggregate_tolerance && metrics.relative_absolute_peak() < aggregate_tolerance;
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
    double qualified_pointwise_absolute_tolerance = 0.0) {
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
        qualified_pointwise_absolute_tolerance);
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

std::pair<std::size_t, std::size_t> locate_source_element(const cartesian::SpatialAssembly& spatial,
    std::size_t source) {
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const Hex8RegionMesh& region_mesh = spatial.region_mesh(region);
        const auto found =
            std::find(region_mesh.source_element_ids().begin(), region_mesh.source_element_ids().end(), source);
        if (found != region_mesh.source_element_ids().end())
            return {region, static_cast<std::size_t>(found - region_mesh.source_element_ids().begin())};
    }
    throw std::invalid_argument("Abaqus HEX8 source element is absent from the Fuelsim problem");
}

std::size_t matching_primary_source(const UnstructuredHex8Mesh& mesh,
    const std::vector<std::size_t>& contact_sources,
    std::size_t secondary_source) {
    const CartesianPoint3& secondary = mesh.nodes().at(secondary_source);
    std::size_t result = std::numeric_limits<std::size_t>::max();
    double minimum_distance = std::numeric_limits<double>::infinity();
    for (std::size_t source = 0; source < mesh.nodes().size(); ++source) {
        if (std::find(contact_sources.begin(), contact_sources.end(), source) != contact_sources.end())
            continue;
        const CartesianPoint3& candidate = mesh.nodes().at(source);
        const double distance =
            std::hypot(candidate.x - secondary.x, candidate.y - secondary.y, candidate.z - secondary.z);
        if (distance < minimum_distance) {
            minimum_distance = distance;
            result = source;
        }
    }
    if (result == std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("Abaqus HEX8 matching primary contact node was not found");
    return result;
}

std::vector<double> secondary_nodal_areas(const UnstructuredHex8Mesh& mesh,
    const std::string& secondary_boundary,
    const cartesian::SpatialAssembly& spatial,
    const std::map<std::size_t, std::size_t>& source_to_global,
    const std::vector<double>& state) {
    static constexpr std::array<std::array<std::size_t, 4>, 6> face_nodes = {
        {{{0, 1, 5, 4}}, {{1, 2, 6, 5}}, {{2, 3, 7, 6}}, {{0, 4, 7, 3}}, {{0, 3, 2, 1}}, {{4, 5, 6, 7}}}};
    std::vector<double> result(mesh.nodes().size(), 0.0);
    const std::array<Field, 3> displacement_fields = {Field::displacement_x,
        Field::displacement_y,
        Field::displacement_z};
    for (const ElementSide& side : mesh.side_set(secondary_boundary).sides) {
        Quad4FaceCoordinates current{};
        std::array<std::size_t, 4> sources{};
        for (std::size_t local = 0; local < 4; ++local) {
            const std::size_t source =
                mesh.elements().at(side.element).nodes.at(face_nodes.at(side.local_side).at(local));
            sources[local] = source;
            current[local] = mesh.nodes().at(source);
            const std::size_t global = source_to_global.at(source);
            current[local].x += state.at(spatial.dof(displacement_fields[0], global));
            current[local].y += state.at(spatial.dof(displacement_fields[1], global));
            current[local].z += state.at(spatial.dof(displacement_fields[2], global));
        }
        const Quad4FaceGeometry geometry = make_quad4_face_geometry(current);
        for (const Quad4FaceQuadraturePoint& point : geometry.points)
            for (std::size_t local = 0; local < 4; ++local)
                result[sources[local]] += point.shape[local] * point.weighted_measure;
    }
    return result;
}
} // namespace

void AbaqusHex8SnapshotObserver::accepted_step(const TransientProblem& problem, const TransientAcceptedStep& step) {
    AbaqusHex8StepSnapshot snapshot;
    snapshot.time = step.time;
    snapshot.load_factor = step.load_factor;
    snapshot.state = problem.committed_solution();
    snapshot.conservation = step.conservation;
    const cartesian::SpatialAssembly& spatial = cartesian::ProblemAccess::view(problem);
    std::size_t source_element_count = 0;
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        for (const std::size_t source : spatial.region_mesh(region).source_element_ids())
            source_element_count = std::max(source_element_count, source + 1);
    snapshot.material_by_source_element.resize(source_element_count);
    for (std::size_t region = 0; region < spatial.region_count(); ++region) {
        const Hex8RegionMesh& region_mesh = spatial.region_mesh(region);
        for (std::size_t element = 0; element < region_mesh.elements().size(); ++element) {
            const std::size_t source = region_mesh.source_element_ids().at(element);
            const CartesianMaterialHistory& history =
                cartesian::ProblemAccess::material_history(problem, region, element);
            snapshot.material_by_source_element[source].assign(history.begin(), history.end());
        }
    }
    const auto& contact_histories = cartesian::ProblemAccess::committed_contact_histories(problem);
    if (!contact_histories.empty()) {
        snapshot.contact = cartesian::ProblemAccess::summarize_contact_nodes(problem, 0, snapshot.state);
        snapshot.contact_history = contact_histories.front();
    }
    _snapshots.push_back(std::move(snapshot));
}

bool compare_abaqus_hex8_full_field(const TransientProblem& solved_problem,
    const SpatialDefinition& definition,
    const UnstructuredHex8Mesh& mesh,
    const std::vector<AbaqusHex8StepSnapshot>& snapshots,
    const AbaqusHex8FullFieldOptions& options) {
    if (options.case_name.empty() || options.reference_prefix.empty() || options.expected_steps == 0
        || !(options.time_step > 0.0))
        throw std::invalid_argument("Abaqus HEX8 full-field options are incomplete");
    const double bulk_pointwise_tolerance = options.bulk_pointwise_relative_tolerance > 0.0
                                                ? options.bulk_pointwise_relative_tolerance
                                                : options.bulk_relative_tolerance,
                 contact_pointwise_tolerance = options.contact_pointwise_relative_tolerance > 0.0
                                                   ? options.contact_pointwise_relative_tolerance
                                                   : options.contact_relative_tolerance,
                 energy_pointwise_tolerance = options.energy_pointwise_relative_tolerance > 0.0
                                                  ? options.energy_pointwise_relative_tolerance
                                                  : options.energy_relative_tolerance;
    const auto pointwise_tolerance = [](double configured, double fallback) {
        return configured > 0.0 ? configured : fallback;
    };
    const double
        displacement_pointwise_tolerance =
            pointwise_tolerance(options.displacement_pointwise_relative_tolerance, bulk_pointwise_tolerance),
        reaction_heat_flux_pointwise_tolerance =
            pointwise_tolerance(options.reaction_heat_flux_pointwise_relative_tolerance, bulk_pointwise_tolerance),
        reaction_pointwise_tolerance =
            pointwise_tolerance(options.reaction_pointwise_relative_tolerance, bulk_pointwise_tolerance),
        stress_pointwise_tolerance =
            pointwise_tolerance(options.stress_pointwise_relative_tolerance, bulk_pointwise_tolerance),
        logarithmic_strain_pointwise_tolerance =
            pointwise_tolerance(options.logarithmic_strain_pointwise_relative_tolerance, bulk_pointwise_tolerance),
        elastic_strain_pointwise_tolerance =
            pointwise_tolerance(options.elastic_strain_pointwise_relative_tolerance, bulk_pointwise_tolerance),
        inelastic_pointwise_tolerance =
            pointwise_tolerance(options.inelastic_pointwise_relative_tolerance, bulk_pointwise_tolerance),
        contact_slip_pointwise_tolerance =
            pointwise_tolerance(options.contact_slip_pointwise_relative_tolerance, contact_pointwise_tolerance),
        contact_replayed_heat_rate_relative_tolerance =
            pointwise_tolerance(options.contact_replayed_heat_rate_relative_tolerance,
                options.contact_relative_tolerance),
        contact_replayed_heat_rate_pointwise_tolerance =
            pointwise_tolerance(options.contact_replayed_heat_rate_pointwise_relative_tolerance,
                contact_pointwise_tolerance),
        contact_total_heat_rate_relative_tolerance =
            pointwise_tolerance(options.contact_total_heat_rate_relative_tolerance, options.contact_relative_tolerance),
        contact_total_heat_rate_pointwise_tolerance =
            pointwise_tolerance(options.contact_total_heat_rate_pointwise_relative_tolerance,
                contact_total_heat_rate_relative_tolerance);
    if (!(options.bulk_relative_tolerance > 0.0) || !(options.contact_relative_tolerance > 0.0)
        || !(options.energy_relative_tolerance > 0.0) || !(options.reaction_zero_absolute_tolerance > 0.0)
        || !(options.minimum_contact_state_match_fraction >= 0.0)
        || !(options.minimum_contact_state_match_fraction <= 1.0))
        throw std::invalid_argument("Abaqus HEX8 full-field tolerances are invalid");
    const std::vector<NodeReference> nodes = read_nodes(options.reference_prefix + "_nodal.csv");
    const std::vector<IntegrationReference> integration =
        read_integration(options.reference_prefix + "_integration.csv");
    const std::vector<ContactReference> contact = read_contact(options.reference_prefix + "_contact.csv");
    const std::vector<EnergyReference> energy = read_energy(options.reference_prefix + "_energy.csv");
    const std::size_t integration_points_per_element = options.reduced_integration ? 1 : 8;
    const cartesian::SpatialAssembly& spatial = cartesian::ProblemAccess::view(solved_problem);
    const bool has_contact = !definition.contacts.empty();
    const std::vector<std::size_t> contact_sources =
        has_contact ? cartesian::ProblemAccess::contact_secondary_source_nodes(solved_problem, 0)
                    : std::vector<std::size_t>{};
    if (snapshots.size() != options.expected_steps || nodes.size() != options.expected_steps * mesh.nodes().size()
        || integration.size() != options.expected_steps * mesh.elements().size() * integration_points_per_element
        || contact.size() != options.expected_steps * contact_sources.size()
        || energy.size() != options.expected_steps) {
        std::ostringstream message;
        message << options.case_name
                << " Abaqus full-field row counts do not match the Fuelsim case: steps=" << snapshots.size() << '/'
                << options.expected_steps << ", nodal_rows=" << nodes.size() << '/'
                << options.expected_steps * mesh.nodes().size() << ", integration_rows=" << integration.size() << '/'
                << options.expected_steps * mesh.elements().size() * integration_points_per_element
                << ", contact_rows=" << contact.size() << '/' << options.expected_steps * contact_sources.size()
                << ", energy_rows=" << energy.size() << '/' << options.expected_steps;
        throw std::invalid_argument(message.str());
    }

    std::map<std::size_t, std::size_t> source_to_global;
    for (std::size_t region = 0; region < spatial.region_count(); ++region)
        for (std::size_t local = 0; local < spatial.region_mesh(region).nodes().size(); ++local)
            source_to_global.emplace(spatial.region_mesh(region).source_node_ids().at(local),
                spatial.global_node(region, local));
    if (source_to_global.size() != mesh.nodes().size())
        throw std::invalid_argument(options.case_name + " Fuelsim source-node map is incomplete");

    const std::string prefix = options.case_name + "_";
    bool passed = true;
    std::array<FieldErrorMetrics, 8> nodal_metrics;
    GroupedFieldErrorMetrics displacement_vector, reaction_force_vector;
    TransientProblem reaction_problem(definition, mesh);
    const auto& reaction_dofs = cartesian::ProblemAccess::dof_map(reaction_problem);
    std::vector<unsigned char> constrained_reaction_dofs(reaction_problem.dof_count(), 0U);
    for (const DirichletCondition& condition : reaction_problem.dirichlet_conditions())
        constrained_reaction_dofs.at(condition.dof) = 1U;
    const std::array<Field, 4> fields = {Field::temperature,
        Field::displacement_x,
        Field::displacement_y,
        Field::displacement_z};
    std::size_t nodal_rows = 0;
    for (std::size_t increment = 1; increment <= options.expected_steps; ++increment) {
        const AbaqusHex8StepSnapshot& snapshot = snapshots.at(increment - 1);
        reaction_problem.begin_time_step({snapshot.time, snapshot.load_factor, true});
        const std::vector<double> reaction = raw_residual(reaction_problem, snapshot.state);
        for (const NodeReference& reference : nodes) {
            if (reference.increment != increment)
                continue;
            if (reference.node > mesh.nodes().size() || std::abs(reference.time - snapshot.time) > 1.0e-7)
                throw std::invalid_argument(options.case_name + " Abaqus nodal label or time is invalid");
            const std::size_t global = source_to_global.at(reference.node - 1);
            std::array<double, 3> actual_displacement{}, expected_displacement{}, actual_reaction{},
                expected_reaction{};
            for (std::size_t field = 0; field < fields.size(); ++field) {
                const std::size_t dof = reaction_dofs.dof(fields[field], global);
                nodal_metrics[field].add(snapshot.state[dof], reference.fields[field]);
                const double reaction_value = constrained_reaction_dofs[dof] != 0U ? reaction[dof] : 0.0;
                nodal_metrics[4 + field].add(reaction_value, reference.fields[4 + field]);
                if (field > 0) {
                    actual_displacement[field - 1] = snapshot.state[dof];
                    expected_displacement[field - 1] = reference.fields[field];
                    actual_reaction[field - 1] = reaction_value;
                    expected_reaction[field - 1] = reference.fields[4 + field];
                }
            }
            displacement_vector.add(actual_displacement.data(), expected_displacement.data(), 3);
            reaction_force_vector.add(actual_reaction.data(), expected_reaction.data(), 3);
            ++nodal_rows;
        }
        reaction_problem.commit_time_step(snapshot.state);
    }
    if (nodal_rows != nodes.size())
        throw std::invalid_argument(options.case_name + " did not map every Abaqus node");
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
                     field == 4 ? options.reaction_heat_flux_pointwise_absolute_tolerance : 0.0)
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
    std::vector<IsotropicThermoelasticMaterial> materials;
    for (const RegionDefinition& region : definition.regions)
        materials.emplace_back(region.material);
    double maximum_integration_coordinate_difference = 0.0;
    std::set<std::tuple<std::size_t, std::size_t, std::size_t>> mapped_integration_points;
    for (const IntegrationReference& reference : integration) {
        if (reference.increment > snapshots.size() || reference.element > mesh.elements().size()
            || reference.point > integration_points_per_element)
            throw std::invalid_argument(options.case_name + " Abaqus integration index is invalid");
        const AbaqusHex8StepSnapshot& snapshot = snapshots.at(reference.increment - 1);
        if (std::abs(reference.time - snapshot.time) > 1.0e-7)
            throw std::invalid_argument(options.case_name + " Abaqus integration time is invalid");
        const std::size_t source_element = reference.element - 1;
        const auto location = locate_source_element(spatial, source_element);
        const std::size_t region = location.first, local_element = location.second;
        const Hex8RegionMesh& region_mesh = spatial.region_mesh(region);
        const Hex8Element& element = region_mesh.elements().at(local_element);
        const Hex8Geometry& geometry = spatial.region_element_geometry(region, local_element);
        Hex8LocalValues local_state{};
        for (std::size_t local = 0; local < 8; ++local) {
            const std::size_t global = spatial.global_node(region, element.nodes[local]);
            local_state[local] = snapshot.state[spatial.dof(Field::temperature, global)];
            local_state[8 + local] = snapshot.state[spatial.dof(Field::displacement_x, global)];
            local_state[16 + local] = snapshot.state[spatial.dof(Field::displacement_y, global)];
            local_state[24 + local] = snapshot.state[spatial.dof(Field::displacement_z, global)];
        }
        std::size_t closest = 0;
        double closest_squared = std::numeric_limits<double>::max();
        std::array<double, 3> closest_position{};
        for (std::size_t q = 0; q < integration_points_per_element; ++q) {
            const Hex8QuadraturePoint& candidate =
                options.reduced_integration ? geometry.reduced_point : geometry.points[q];
            CartesianPoint3 current = candidate.position;
            for (std::size_t local = 0; local < 8; ++local) {
                current.x += candidate.shape[local] * local_state[8 + local];
                current.y += candidate.shape[local] * local_state[16 + local];
                current.z += candidate.shape[local] * local_state[24 + local];
            }
            const double distance_squared = std::pow(current.x - reference.position.x, 2)
                                            + std::pow(current.y - reference.position.y, 2)
                                            + std::pow(current.z - reference.position.z, 2);
            if (distance_squared < closest_squared) {
                closest = q;
                closest_squared = distance_squared;
                closest_position = {current.x, current.y, current.z};
            }
        }
        maximum_integration_coordinate_difference =
            std::max(maximum_integration_coordinate_difference, std::sqrt(closest_squared));
        if (!mapped_integration_points.emplace(reference.increment, source_element, closest).second)
            throw std::invalid_argument(options.case_name + " maps two Abaqus rows to one Fuelsim integration point");
        const Hex8QuadraturePoint& point =
            options.reduced_integration ? geometry.reduced_point : geometry.points[closest];
        const StrainFormulation formulation = definition.regions.at(region).strain_formulation;
        const auto diagnostics = options.reduced_integration
                                     ? test::recover_c3d8rt(geometry, local_state, {}, formulation)
                                     : test::recover_c3d8t(geometry, local_state, {}, formulation);
        const auto& point_diagnostics = diagnostics.points[closest];
        const double current_volume = diagnostics.current_volume;
        const double material_temperature = point_diagnostics.temperature;
        const auto& thermal_gradient = point_diagnostics.thermal_gradient;
        const double conductivity = materials.at(region)
                                        .conductivity(adlite::Scalar(material_temperature),
                                            {snapshot.time, point.position.x, point.position.y, point.position.z})
                                        .value();
        std::array<double, 3> actual_heat_flux{};
        for (std::size_t local = 0; local < 8; ++local)
            for (std::size_t component = 0; component < 3; ++component)
                actual_heat_flux[component] -= conductivity * thermal_gradient[local][component] * local_state[local];
        const CartesianMaterialPointState& material =
            snapshot.material_by_source_element.at(source_element).at(closest);
        const std::array<double, 6> actual_stress = components(material.stress),
                                    expected_stress = components(reference.stress);
        std::array<double, 6> actual_logarithmic = components(logarithmic_strain(point, local_state));
        if (formulation == StrainFormulation::finite && !options.reduced_integration) {
            const double local_trace = actual_logarithmic[0] + actual_logarithmic[1] + actual_logarithmic[2];
            const double selective_trace = std::log(current_volume / geometry.reference_volume);
            for (std::size_t component = 0; component < 3; ++component)
                actual_logarithmic[component] += (selective_trace - local_trace) / 3.0;
        }
        const std::array<double, 6> expected_logarithmic = components(reference.logarithmic_strain),
                                    expected_elastic = components(reference.elastic_strain),
                                    expected_plastic = components(reference.plastic_strain),
                                    expected_creep = components(reference.creep_strain);
        const std::array<double, 3> expected_position = {reference.position.x,
            reference.position.y,
            reference.position.z};
        integration_position.add(closest_position.data(), expected_position.data(), 3);
        for (std::size_t component = 0; component < 3; ++component)
            integration_metrics[component].add(actual_heat_flux[component], reference.heat_flux[component]);
        for (std::size_t component = 0; component < 6; ++component) {
            integration_metrics[3 + component].add(actual_stress[component], expected_stress[component]);
            integration_metrics[9 + component].add(actual_logarithmic[component], expected_logarithmic[component]);
            integration_metrics[15 + component].add(material.elastic_strain[component], expected_elastic[component]);
            integration_metrics[21 + component].add(material.plastic_strain[component], expected_plastic[component]);
            integration_metrics[28 + component].add(material.creep_strain[component], expected_creep[component]);
        }
        heat_flux_vector.add(actual_heat_flux.data(), reference.heat_flux.data(), 3);
        stress_tensor.add(actual_stress.data(), expected_stress.data(), 6);
        logarithmic_strain_tensor.add(actual_logarithmic.data(), expected_logarithmic.data(), 6);
        elastic_strain_tensor.add(material.elastic_strain.data(), expected_elastic.data(), 6);
        plastic_strain_tensor.add(material.plastic_strain.data(), expected_plastic.data(), 6);
        creep_strain_tensor.add(material.creep_strain.data(), expected_creep.data(), 6);
        integration_metrics[27].add(material.equivalent_plastic_strain, reference.equivalent_plastic_strain);
        integration_metrics[34].add(material.equivalent_creep_strain, reference.equivalent_creep_strain);
        integration_metrics[35].add(material_temperature, reference.temperature);
        const double volume = options.reduced_integration
                                  ? current_volume
                                  : (formulation == StrainFormulation::finite
                                            ? point.weighted_measure / geometry.reference_volume * current_volume
                                            : point.weighted_measure);
        integration_metrics[36].add(volume, reference.integration_volume);
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
                 1.0,
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
                 1.0e-8,
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
    if (maximum_integration_coordinate_difference >= 1.0e-3) {
        std::cerr << "[FAIL] " << options.case_name << " integration-point coordinates do not map uniquely\n";
        passed = false;
    }

    std::array<FieldErrorMetrics, 6> contact_scalar_metrics;
    GroupedFieldErrorMetrics contact_position, contact_slip, contact_normal_force, contact_shear_force,
        contact_shear_traction, contact_complete_force, contact_resultant, contact_moment, contact_center;
    FieldErrorMetrics replayed_contact_heat_rate, total_contact_heat_rate;
    double maximum_tangent_basis_error = 0.0, maximum_contact_heat_conservation_error = 0.0,
           maximum_replayed_contact_heat_conservation_error = 0.0;
    double maximum_summary_total_slip = 0.0, maximum_history_total_slip = 0.0;
    std::size_t maximum_history_total_slip_increment = 0, maximum_history_total_slip_index = 0;
    bool contact_states_match = true;
    std::array<std::size_t, 9> contact_state_pairs{};
    bool first_contact_state_mismatch_reported = false;
    std::array<FieldErrorMetrics, 8> energy_metrics;
    double maximum_abaqus_artificial_energy = 0.0, maximum_abaqus_artificial_energy_fraction = 0.0;
    double cumulative_elastic = 0.0, cumulative_plastic = 0.0, cumulative_creep = 0.0, cumulative_friction = 0.0,
           cumulative_external_work = 0.0;
    TransientProblem contact_replay_problem(definition, mesh);
    const auto& contact_replay_dofs = cartesian::ProblemAccess::dof_map(contact_replay_problem);
    const cartesian::SpatialAssembly& contact_replay_spatial = cartesian::ProblemAccess::view(contact_replay_problem);
    const std::array<Field, 4> replay_fields = {Field::temperature,
        Field::displacement_x,
        Field::displacement_y,
        Field::displacement_z};
    double replay_previous_time = 0.0;
    for (std::size_t increment = 1; increment <= options.expected_steps; ++increment) {
        const AbaqusHex8StepSnapshot& snapshot = snapshots.at(increment - 1);
        if (has_contact) {
            for (std::size_t history = 0; history < snapshot.contact_history.size(); ++history) {
                const std::array<double, 3>& total = snapshot.contact_history[history].cartesian_total_tangential_slip;
                const double magnitude = std::hypot(total[0], total[1], total[2]);
                if (magnitude > maximum_history_total_slip) {
                    maximum_history_total_slip = magnitude;
                    maximum_history_total_slip_increment = increment;
                    maximum_history_total_slip_index = history;
                }
            }
            if (snapshot.contact.size() != contact_sources.size())
                throw std::invalid_argument(options.case_name + " Fuelsim contact snapshot has an invalid node count");
            spatial.validate_state(snapshot.state);
            const std::vector<double> thermal_contact =
                thermal_contact_residual(spatial, snapshot.state, maximum_contact_heat_conservation_error);
            std::vector<double> abaqus_state = contact_replay_problem.committed_solution();
            std::size_t replayed_nodes = 0;
            for (const NodeReference& reference : nodes) {
                if (reference.increment != increment)
                    continue;
                if (std::abs(reference.time - snapshot.time) > 1.0e-7)
                    throw std::invalid_argument(options.case_name + " Abaqus replay time is invalid");
                const std::size_t global = source_to_global.at(reference.node - 1);
                for (std::size_t field = 0; field < replay_fields.size(); ++field)
                    abaqus_state[contact_replay_dofs.dof(replay_fields[field], global)] = reference.fields[field];
                ++replayed_nodes;
            }
            if (replayed_nodes != mesh.nodes().size() || !(snapshot.time > replay_previous_time))
                throw std::invalid_argument(options.case_name + " Abaqus thermal-contact replay is incomplete");
            contact_replay_problem.begin_time_step({snapshot.time, snapshot.load_factor, true});
            contact_replay_problem.validate_state(abaqus_state);
            const std::vector<double> abaqus_state_thermal_contact = thermal_contact_residual(contact_replay_spatial,
                abaqus_state,
                maximum_replayed_contact_heat_conservation_error);
            const std::vector<double> contact_area = secondary_nodal_areas(mesh,
                definition.contacts.at(0).secondary,
                spatial,
                source_to_global,
                snapshot.state);
            std::array<double, 3> actual_resultant{}, expected_resultant{}, actual_moment{}, expected_moment{},
                actual_center_sum{}, expected_center_sum{};
            double actual_center_weight = 0.0, expected_center_weight = 0.0, actual_heat_rate = 0.0,
                   replayed_heat_rate = 0.0;
            for (std::size_t output = 0; output < contact_sources.size(); ++output) {
                const std::size_t source = contact_sources.at(output);
                const auto found = std::find_if(contact.begin(), contact.end(), [&](const ContactReference& value) {
                    return value.increment == increment && value.node == source + 1;
                });
                if (found == contact.end())
                    throw std::invalid_argument(options.case_name + " Abaqus contact-node mapping is incomplete");
                if (std::abs(found->time - snapshot.time) > 1.0e-7)
                    throw std::invalid_argument(options.case_name + " Abaqus contact time is invalid");
                const CartesianContactNodeSummary& actual = snapshot.contact.at(output);
                const std::size_t global = source_to_global.at(source);
                const std::size_t primary_source = options.use_contact_summary_total_slip
                                                       ? std::numeric_limits<std::size_t>::max()
                                                       : matching_primary_source(mesh, contact_sources, source);
                const std::size_t primary_global = options.use_contact_summary_total_slip
                                                       ? std::numeric_limits<std::size_t>::max()
                                                       : source_to_global.at(primary_source);
                const std::array<Field, 3> displacement_fields = {Field::displacement_x,
                    Field::displacement_y,
                    Field::displacement_z};
                std::array<double, 3> actual_position = {mesh.nodes().at(source).x,
                    mesh.nodes().at(source).y,
                    mesh.nodes().at(source).z};
                const std::array<double, 3> expected_position = {found->position.x,
                    found->position.y,
                    found->position.z};
                std::array<double, 3> expected_slip{}, actual_slip{}, relative_displacement{}, actual_force{},
                    expected_force{}, actual_normal{}, expected_normal{}, actual_shear{}, expected_shear{},
                    actual_shear_stress{}, expected_shear_stress{};
                for (std::size_t component = 0; component < 3; ++component) {
                    actual_position[component] += snapshot.state[spatial.dof(displacement_fields[component], global)];
                    if (!options.use_contact_summary_total_slip)
                        relative_displacement[component] =
                            snapshot.state[spatial.dof(displacement_fields[component], global)]
                            - snapshot.state[spatial.dof(displacement_fields[component], primary_global)];
                    expected_slip[component] = found->slip_first * found->tangent_first[component]
                                               + found->slip_second * found->tangent_second[component];
                    actual_normal[component] = -actual.normal_contact_force[component];
                    expected_normal[component] = found->normal_force[component];
                    actual_shear[component] = -actual.tangential_contact_force[component];
                    expected_shear[component] = found->shear_force[component];
                    actual_force[component] = actual_normal[component] + actual_shear[component];
                    expected_force[component] = expected_normal[component] + expected_shear[component];
                    expected_shear_stress[component] =
                        -found->shear_traction_first * found->tangent_first[component]
                        - found->shear_traction_second * found->tangent_second[component];
                    if (actual.tributary_area > 0.0)
                        actual_shear_stress[component] = actual_shear[component] / actual.tributary_area;
                    actual_resultant[component] += actual_force[component];
                    expected_resultant[component] += expected_force[component];
                }
                if (options.use_contact_summary_total_slip)
                    actual_slip = actual.tangential_slip;
                else {
                    const double actual_slip_first = relative_displacement[0] * found->tangent_first[0]
                                                     + relative_displacement[1] * found->tangent_first[1]
                                                     + relative_displacement[2] * found->tangent_first[2];
                    const double actual_slip_second = relative_displacement[0] * found->tangent_second[0]
                                                      + relative_displacement[1] * found->tangent_second[1]
                                                      + relative_displacement[2] * found->tangent_second[2];
                    for (std::size_t component = 0; component < 3; ++component)
                        actual_slip[component] = actual_slip_first * found->tangent_first[component]
                                                 + actual_slip_second * found->tangent_second[component];
                }
                maximum_summary_total_slip =
                    std::max(maximum_summary_total_slip, std::hypot(actual_slip[0], actual_slip[1], actual_slip[2]));
                contact_position.add(actual_position.data(), expected_position.data(), 3);
                if (found->state != 0)
                    contact_slip.add(actual_slip.data(), expected_slip.data(), 3);
                contact_normal_force.add(actual_normal.data(), expected_normal.data(), 3);
                contact_shear_force.add(actual_shear.data(), expected_shear.data(), 3);
                contact_shear_traction.add(actual_shear_stress.data(), expected_shear_stress.data(), 3);
                contact_complete_force.add(actual_force.data(), expected_force.data(), 3);
                contact_scalar_metrics[0].add(actual.gap, found->opening);
                contact_scalar_metrics[1].add(actual.pressure, found->pressure);
                if (found->state != 0)
                    contact_scalar_metrics[2].add(std::hypot(actual_slip[0], actual_slip[1], actual_slip[2]),
                        std::hypot(expected_slip[0], expected_slip[1], expected_slip[2]));
                const double actual_normal_force = std::hypot(actual_normal[0], actual_normal[1], actual_normal[2]),
                             expected_normal_force =
                                 std::hypot(expected_normal[0], expected_normal[1], expected_normal[2]),
                             actual_shear_force = std::hypot(actual_shear[0], actual_shear[1], actual_shear[2]),
                             expected_shear_force = std::hypot(expected_shear[0], expected_shear[1], expected_shear[2]);
                contact_scalar_metrics[3].add(actual_normal_force, expected_normal_force);
                contact_scalar_metrics[4].add(actual_shear_force, expected_shear_force);
                const std::size_t temperature_dof = spatial.dof(Field::temperature, global);
                const std::size_t replay_temperature_dof = contact_replay_dofs.dof(Field::temperature, global);
                const double expected_nodal_heat_rate = found->heat_flux * contact_area.at(source);
                contact_scalar_metrics[5].add(thermal_contact.at(temperature_dof), expected_nodal_heat_rate);
                replayed_contact_heat_rate.add(thermal_contact.at(temperature_dof),
                    abaqus_state_thermal_contact.at(replay_temperature_dof));
                actual_heat_rate += thermal_contact.at(temperature_dof);
                replayed_heat_rate += abaqus_state_thermal_contact.at(replay_temperature_dof);
                const double actual_weight = actual_normal_force, expected_weight = expected_normal_force;
                actual_center_weight += actual_weight;
                expected_center_weight += expected_weight;
                for (std::size_t component = 0; component < 3; ++component) {
                    actual_center_sum[component] += actual_weight * actual_position[component];
                    expected_center_sum[component] += expected_weight * expected_position[component];
                }
                const auto add_moment = [](std::array<double, 3>& moment,
                                            const std::array<double, 3>& position,
                                            const std::array<double, 3>& force) {
                    moment[0] += position[1] * force[2] - position[2] * force[1];
                    moment[1] += position[2] * force[0] - position[0] * force[2];
                    moment[2] += position[0] * force[1] - position[1] * force[0];
                };
                add_moment(actual_moment, actual_position, actual_force);
                add_moment(expected_moment, expected_position, expected_force);
                maximum_tangent_basis_error = std::max({maximum_tangent_basis_error,
                    std::abs(found->tangent_first[0] * found->tangent_first[0]
                             + found->tangent_first[1] * found->tangent_first[1]
                             + found->tangent_first[2] * found->tangent_first[2] - 1.0),
                    std::abs(found->tangent_second[0] * found->tangent_second[0]
                             + found->tangent_second[1] * found->tangent_second[1]
                             + found->tangent_second[2] * found->tangent_second[2] - 1.0),
                    std::abs(found->tangent_first[0] * found->tangent_second[0]
                             + found->tangent_first[1] * found->tangent_second[1]
                             + found->tangent_first[2] * found->tangent_second[2])});
                const std::size_t actual_state = actual.pressure <= 0.0 ? 0 : actual.sliding ? 2 : 1;
                ++contact_state_pairs.at(3 * actual_state + found->state);
                if (actual_state != found->state && !first_contact_state_mismatch_reported) {
                    const double actual_elastic_slip = std::hypot(actual.elastic_tangential_slip[0],
                                     actual.elastic_tangential_slip[1],
                                     actual.elastic_tangential_slip[2]),
                                 actual_total_slip = std::hypot(actual_slip[0], actual_slip[1], actual_slip[2]),
                                 expected_total_slip = std::hypot(expected_slip[0], expected_slip[1], expected_slip[2]),
                                 friction = definition.contacts.at(0).friction_coefficient,
                                 actual_coulomb_ratio = actual.pressure > 0.0 && friction > 0.0
                                                            ? actual_shear_force / (friction * actual_normal_force)
                                                            : 0.0,
                                 expected_coulomb_ratio =
                                     found->pressure > 0.0 && friction > 0.0
                                         ? expected_shear_force / (friction * expected_normal_force)
                                         : 0.0;
                    std::cout << prefix << "first_contact_state_mismatch_increment=" << increment << '\n'
                              << prefix << "first_contact_state_mismatch_node=" << source + 1 << '\n'
                              << prefix << "first_contact_state_mismatch_actual=" << actual_state << '\n'
                              << prefix << "first_contact_state_mismatch_reference=" << found->state << '\n'
                              << prefix << "first_contact_state_mismatch_actual_elastic_slip=" << actual_elastic_slip
                              << '\n'
                              << prefix << "first_contact_state_mismatch_actual_total_slip=" << actual_total_slip
                              << '\n'
                              << prefix << "first_contact_state_mismatch_reference_total_slip=" << expected_total_slip
                              << '\n'
                              << prefix << "first_contact_state_mismatch_actual_coulomb_ratio=" << actual_coulomb_ratio
                              << '\n'
                              << prefix
                              << "first_contact_state_mismatch_reference_coulomb_ratio=" << expected_coulomb_ratio
                              << '\n';
                    first_contact_state_mismatch_reported = true;
                }
                contact_states_match = contact_states_match && actual_state == found->state;
            }
            total_contact_heat_rate.add(actual_heat_rate, replayed_heat_rate);
            contact_resultant.add(actual_resultant.data(), expected_resultant.data(), 3);
            contact_moment.add(actual_moment.data(), expected_moment.data(), 3);
            if (actual_center_weight > 0.0 && expected_center_weight > 0.0) {
                std::array<double, 3> actual_center{}, expected_center{};
                for (std::size_t component = 0; component < 3; ++component) {
                    actual_center[component] = actual_center_sum[component] / actual_center_weight;
                    expected_center[component] = expected_center_sum[component] / expected_center_weight;
                }
                contact_center.add(actual_center.data(), expected_center.data(), 3);
            }
            contact_replay_problem.commit_time_step(abaqus_state);
            replay_previous_time = snapshot.time;
        }

        const EnergyReference& expected = energy.at(increment - 1);
        if (expected.increment != increment || std::abs(expected.time - snapshot.time) > 1.0e-7)
            throw std::invalid_argument(options.case_name + " Abaqus energy time is invalid");
        cumulative_elastic += snapshot.conservation.elastic_energy_change;
        cumulative_plastic += snapshot.conservation.plastic_dissipation_increment;
        cumulative_creep += snapshot.conservation.creep_dissipation_increment;
        cumulative_friction += snapshot.conservation.friction_dissipation_increment;
        cumulative_external_work += snapshot.conservation.trapezoidal_pressure_traction_work_increment
                                    + snapshot.conservation.trapezoidal_dirichlet_reaction_work_increment
                                    + snapshot.conservation.body_force_work_increment;
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
        energy_metrics[6].add(snapshot.conservation.dirichlet_heat_input_rate, expected.boundary_heat_rate);
        energy_metrics[7].add(snapshot.conservation.mechanical_hourglass_energy, expected.artificial);
    }
    if (has_contact) {
        const std::array<std::string, 6> contact_scalar_names = {"contact_opening",
            "contact_pressure",
            "contact_slip_magnitude",
            "contact_normal_force_magnitude",
            "contact_shear_force_magnitude",
            "contact_recovered_heat_rate"};
        for (std::size_t field = 0; field < contact_scalar_metrics.size(); ++field)
            passed = report_metric(prefix + contact_scalar_names[field],
                         contact_scalar_metrics[field],
                         options.contact_relative_tolerance,
                         contact_pointwise_tolerance,
                         field == 0 || field == 2 ? 1.0e-10
                         : field == 1             ? 1.0
                                                  : 1.0e-2,
                         (field == 1 && options.gate_contact_pressure) || field == 3)
                     && passed;
        passed = report_metric(prefix + "contact_abaqus_state_replayed_heat_rate",
                     replayed_contact_heat_rate,
                     contact_replayed_heat_rate_relative_tolerance,
                     contact_replayed_heat_rate_pointwise_tolerance,
                     1.0e-2,
                     true,
                     options.contact_replayed_heat_rate_pointwise_absolute_tolerance)
                 && passed;
        passed = report_metric(prefix + "contact_total_heat_rate",
                     total_contact_heat_rate,
                     contact_total_heat_rate_relative_tolerance,
                     contact_total_heat_rate_pointwise_tolerance,
                     1.0e-2,
                     true)
                 && passed;
        passed = report_grouped(prefix + "contact_position",
                     contact_position,
                     options.contact_relative_tolerance,
                     contact_pointwise_tolerance,
                     options.coordinate_tolerance,
                     true)
                 && passed;
        passed = report_grouped(prefix + "contact_slip_vector",
                     contact_slip,
                     options.contact_relative_tolerance,
                     contact_slip_pointwise_tolerance,
                     1.0e-10,
                     options.gate_contact_slip,
                     options.contact_slip_pointwise_absolute_tolerance)
                 && passed;
        passed = report_grouped(prefix + "contact_normal_force_vector",
                     contact_normal_force,
                     options.contact_relative_tolerance,
                     contact_pointwise_tolerance,
                     1.0e-2,
                     true)
                 && passed;
        passed = report_grouped(prefix + "contact_complete_force_vector",
                     contact_complete_force,
                     options.contact_relative_tolerance,
                     contact_pointwise_tolerance,
                     1.0e-2,
                     true)
                 && passed;
        passed = report_grouped(prefix + "contact_resultant",
                     contact_resultant,
                     options.contact_relative_tolerance,
                     contact_pointwise_tolerance,
                     1.0e-2,
                     true)
                 && passed;
        passed = report_grouped(prefix + "contact_moment",
                     contact_moment,
                     options.contact_relative_tolerance,
                     contact_pointwise_tolerance,
                     1.0e-2,
                     true)
                 && passed;
        passed = report_grouped(prefix + "contact_normal_force_center",
                     contact_center,
                     options.contact_relative_tolerance,
                     contact_pointwise_tolerance,
                     options.coordinate_tolerance,
                     true)
                 && passed;
        report_grouped(prefix + "contact_shear_force_vector",
            contact_shear_force,
            options.contact_relative_tolerance,
            contact_pointwise_tolerance,
            1.0e-2,
            false);
        report_grouped(prefix + "contact_shear_traction_vector",
            contact_shear_traction,
            options.contact_relative_tolerance,
            contact_pointwise_tolerance,
            1.0,
            false);
        const std::size_t contact_state_count =
                              std::accumulate(contact_state_pairs.begin(), contact_state_pairs.end(), std::size_t{0}),
                          matching_contact_state_count =
                              contact_state_pairs[0] + contact_state_pairs[4] + contact_state_pairs[8];
        const double contact_state_match_fraction =
            contact_state_count > 0
                ? static_cast<double>(matching_contact_state_count) / static_cast<double>(contact_state_count)
                : 0.0;
        std::cout << prefix << "contact_states_match=" << contact_states_match << '\n'
                  << prefix << "contact_state_pairs_00_01_02_10_11_12_20_21_22=";
        for (std::size_t pair = 0; pair < contact_state_pairs.size(); ++pair)
            std::cout << (pair == 0 ? "" : ",") << contact_state_pairs[pair];
        std::cout << '\n'
                  << prefix << "contact_state_match_fraction=" << contact_state_match_fraction << '\n'
                  << prefix << "maximum_tangent_basis_error=" << maximum_tangent_basis_error << '\n'
                  << prefix << "maximum_summary_total_slip=" << maximum_summary_total_slip << '\n'
                  << prefix << "maximum_history_total_slip=" << maximum_history_total_slip << '\n'
                  << prefix << "maximum_history_total_slip_increment=" << maximum_history_total_slip_increment << '\n'
                  << prefix << "maximum_history_total_slip_index=" << maximum_history_total_slip_index << '\n'
                  << prefix << "maximum_contact_heat_conservation_error=" << maximum_contact_heat_conservation_error
                  << '\n'
                  << prefix << "maximum_replayed_contact_heat_conservation_error="
                  << maximum_replayed_contact_heat_conservation_error << '\n';
        if ((options.gate_contact_state && contact_state_match_fraction < options.minimum_contact_state_match_fraction)
            || maximum_tangent_basis_error >= options.tangent_basis_tolerance
            || maximum_contact_heat_conservation_error >= 1.0e-8
            || maximum_replayed_contact_heat_conservation_error >= 1.0e-8)
            passed = false;
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
                     field == 6 ? 1.0e-2 : 1.0e-8,
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
              << prefix << "compared_contact_rows=" << contact.size() << '\n'
              << prefix << "compared_energy_rows=" << energy.size() << '\n';
    if (passed)
        std::cout << "[PASS] " << options.case_name << " Abaqus full-field comparison\n";
    return passed;
}
} // namespace fuelsim::test
