#include "support/steady_single_region_problem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fuelsim {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

Quad4Coordinates element_coordinates(const StructuredRzMesh& mesh,
                                     const Quad4Element& element) {
    Quad4Coordinates coordinates{};
    for (std::size_t node = 0; node < element.nodes.size(); ++node)
        coordinates[node] = mesh.nodes().at(element.nodes[node]);
    return coordinates;
}

void append_boundary_conditions(const std::vector<std::size_t>& nodes,
                                Field field, double value,
                                const DofMap& dof_map,
                                std::vector<DirichletCondition>& output) {
    for (std::size_t node : nodes)
        output.push_back({dof_map.dof(field, node), value});
}

} // namespace

SteadySingleRegionProblem::SteadySingleRegionProblem(
    SteadySingleRegionParameters parameters)
    : SteadySingleRegionProblem(
          parameters, StructuredRzMesh::make_annulus(
                          parameters.inner_radius, parameters.outer_radius,
                          parameters.length, parameters.radial_elements,
                          parameters.axial_elements)) {}

SteadySingleRegionProblem::SteadySingleRegionProblem(
    SteadySingleRegionParameters parameters, StructuredRzMesh mesh)
    : _parameters(parameters), _mesh(std::move(mesh)),
      _dof_map(_mesh.nodes().size()),
      _kernel(IsotropicThermoelasticMaterial(parameters.fuel),
              parameters.volumetric_heat_source) {
    const auto same_geometry = [](double actual, double expected) {
        const double scale =
            std::max({1.0, std::abs(actual), std::abs(expected)});
        return std::abs(actual - expected) <= 1.0e-12 * scale;
    };
    if (_mesh.radial_elements() != _parameters.radial_elements ||
        _mesh.axial_elements() != _parameters.axial_elements ||
        !same_geometry(_mesh.inner_radius(), _parameters.inner_radius) ||
        !same_geometry(_mesh.outer_radius(), _parameters.outer_radius) ||
        !same_geometry(_mesh.length(), _parameters.length))
        throw std::invalid_argument(
            "SteadySingleRegionProblem imported mesh does not match "
            "SteadySingleRegionParameters");
    if (!(_parameters.inner_pressure >= 0.0) ||
        !(_parameters.outer_pressure >= 0.0))
        throw std::invalid_argument(
            "SteadySingleRegionProblem pressures must be nonnegative");

    _geometries.reserve(_mesh.elements().size());
    for (const Quad4Element& element : _mesh.elements())
        _geometries.push_back(
            make_quad4_rz_geometry(element_coordinates(_mesh, element)));

    build_dirichlet_conditions();
}

const SteadySingleRegionParameters&
SteadySingleRegionProblem::parameters() const noexcept {
    return _parameters;
}

const StructuredRzMesh& SteadySingleRegionProblem::mesh() const noexcept {
    return _mesh;
}

const DofMap& SteadySingleRegionProblem::dof_map() const noexcept {
    return _dof_map;
}

const Quad4RzThermoelasticKernel&
SteadySingleRegionProblem::kernel() const noexcept {
    return _kernel;
}

std::size_t SteadySingleRegionProblem::dof_count() const noexcept {
    return _dof_map.dof_count();
}

std::size_t SteadySingleRegionProblem::element_count() const noexcept {
    return _mesh.elements().size();
}

std::size_t SteadySingleRegionProblem::contribution_count() const noexcept {
    return element_count();
}

const std::vector<DirichletCondition>&
SteadySingleRegionProblem::dirichlet_conditions() const noexcept {
    return _dirichlet_conditions;
}

std::vector<double> SteadySingleRegionProblem::initial_state() const {
    std::vector<double> state(dof_count(), 0.0);
    for (std::size_t node = 0; node < _mesh.nodes().size(); ++node)
        state[_dof_map.temperature(node)] = _parameters.initial_temperature;

    for (const DirichletCondition& condition : _dirichlet_conditions)
        state[condition.dof] = condition.value;
    return state;
}

LocalDofs
SteadySingleRegionProblem::element_dofs(std::size_t element_index) const {
    return _dof_map.element_dofs(_mesh.elements().at(element_index));
}

LocalValues SteadySingleRegionProblem::element_state(
    std::size_t element_index, const std::vector<double>& global_state) const {
    return contribution_state(element_index, global_state);
}

LocalResidual
SteadySingleRegionProblem::element_residual(std::size_t element_index,
                                            const LocalValues& state) const {
    return _kernel.residual(_geometries.at(element_index), state);
}

LocalSystem
SteadySingleRegionProblem::linearize_element(std::size_t element_index,
                                             const LocalValues& state) const {
    return _kernel.linearize(_geometries.at(element_index), state);
}

LocalDofs SteadySingleRegionProblem::contribution_dofs(
    std::size_t contribution_index) const {
    return element_dofs(contribution_index);
}

LocalResidual SteadySingleRegionProblem::contribution_residual(
    std::size_t contribution_index, const LocalValues& state) const {
    return element_residual(contribution_index, state);
}

LocalSystem SteadySingleRegionProblem::linearize_contribution(
    std::size_t contribution_index, const LocalValues& state) const {
    return linearize_element(contribution_index, state);
}

void SteadySingleRegionProblem::add_state_independent_residual(
    std::vector<double>& residual) const {
    add_pressure_residual(residual);
}

const Quad4RzGeometry&
SteadySingleRegionProblem::element_geometry(std::size_t element_index) const {
    return _geometries.at(element_index);
}

void SteadySingleRegionProblem::build_dirichlet_conditions() {
    append_boundary_conditions(
        _mesh.boundary_nodes(BoundaryId::radial_outer), Field::temperature,
        _parameters.outer_temperature, _dof_map, _dirichlet_conditions);

    if (_parameters.inner_radius == 0.0) {
        append_boundary_conditions(
            _mesh.boundary_nodes(BoundaryId::radial_inner),
            Field::radial_displacement, 0.0, _dof_map, _dirichlet_conditions);
    }

    append_boundary_conditions(_mesh.boundary_nodes(BoundaryId::bottom),
                               Field::axial_displacement, 0.0, _dof_map,
                               _dirichlet_conditions);

    std::sort(_dirichlet_conditions.begin(), _dirichlet_conditions.end(),
              [](const DirichletCondition& lhs, const DirichletCondition& rhs) {
                  return lhs.dof < rhs.dof;
              });

    for (std::size_t index = 1; index < _dirichlet_conditions.size(); ++index) {
        const DirichletCondition& previous = _dirichlet_conditions[index - 1];
        const DirichletCondition& current = _dirichlet_conditions[index];
        if (previous.dof != current.dof)
            continue;
        if (previous.value != current.value)
            throw std::invalid_argument("SteadySingleRegionProblem has "
                                        "conflicting Dirichlet conditions");
        throw std::invalid_argument(
            "SteadySingleRegionProblem has duplicate Dirichlet conditions");
    }
}

void SteadySingleRegionProblem::add_pressure_residual(
    std::vector<double>& residual) const {
    const std::array<double, 2> gauss = {
        -0.577350269189625764509148780501957456,
        0.577350269189625764509148780501957456,
    };

    const auto add_boundary = [&](BoundaryId boundary, double pressure,
                                  double normal_r) {
        if (pressure == 0.0)
            return;

        for (const Line2BoundaryElement& edge :
             _mesh.boundary_elements(boundary)) {
            const RzPoint& first = _mesh.nodes().at(edge.nodes[0]);
            const RzPoint& second = _mesh.nodes().at(edge.nodes[1]);
            const double dr = second.r - first.r;
            const double dz = second.z - first.z;
            const double line_jacobian = 0.5 * std::sqrt(dr * dr + dz * dz);

            for (double xi : gauss) {
                const std::array<double, 2> shape = {
                    0.5 * (1.0 - xi),
                    0.5 * (1.0 + xi),
                };
                const double radius = shape[0] * first.r + shape[1] * second.r;
                const double measure = 2.0 * pi * radius * line_jacobian;

                for (std::size_t node = 0; node < 2; ++node) {
                    const std::size_t dof =
                        _dof_map.radial_displacement(edge.nodes[node]);
                    residual[dof] +=
                        measure * pressure * normal_r * shape[node];
                }
            }
        }
    };

    add_boundary(BoundaryId::radial_inner, _parameters.inner_pressure, -1.0);
    add_boundary(BoundaryId::radial_outer, _parameters.outer_pressure, 1.0);
}

} // namespace fuelsim
